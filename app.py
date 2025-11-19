import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import numpy as np
from io import StringIO
from converter import Converter

def zmat_string_to_cartesian_string(zmat_string):
    c = Converter()
    zmat = []
    lines = zmat_string.strip().split('\n')
    for i, line in enumerate(lines):
        items = line.split()
        name = items[0]
        mass = c.masses[name]
        if i == 0:
            zmat.append([name, [], mass])
        elif i == 1:
            name, atom1, distance = items
            zmat.append([name, [[int(atom1)-1, float(distance)],[],[]], mass])
        elif i == 2:
            name, atom1, distance, atom2, angle = items
            zmat.append([name, [[int(atom1)-1, float(distance)],[int(atom2)-1, np.radians(float(angle))],[]], mass])
        else:
            name, atom1, distance, atom2, angle, atom3, dihedral = items
            zmat.append([name, [[int(atom1)-1, float(distance)],
                                [int(atom2)-1, np.radians(float(angle))],
                                [int(atom3)-1, np.radians(float(dihedral))]], mass])
    c.zmatrix = zmat
    c.zmatrix_to_cartesian()
    return c.str_cartesian()


def parse_xyz(xyz_text):
    """
    Parse XYZ format or Z-matrix text and return element list and coordinates array.
    
    XYZ format:
    <number of atoms>
    <comment line>
    <element> <x> <y> <z>
    ...
    
    Z-matrix format (Fenske-Hall):
    atom1
    atom2  ref1  distance
    atom3  ref1  distance  ref2  angle
    atom4  ref1  distance  ref2  angle  ref3  dihedral
    ...
    
    Returns:
        tuple: (elements list, coords numpy array) or (None, None) if parsing fails
    """
    if not xyz_text.strip():
        return None, None
    
    lines = [line.strip() for line in xyz_text.strip().split('\n') if line.strip()]
    
    # Detect format: Z-matrix has 1, 3, 5, or 7 fields per line
    # XYZ has 4 fields (element x y z)
    first_data_line = lines[0]
    try:
        # Skip potential header in XYZ
        if first_data_line.split()[0].isdigit():
            first_data_line = lines[2] if len(lines) > 2 else lines[0]
    except:
        pass
    
    parts = first_data_line.split()
    
    # Check if it looks like Z-matrix (1, 3, 5, or 7 fields)
    if len(parts) in [1, 3, 5, 7]:
        # Likely Z-matrix format - try to convert it
        try:
            cartesian_string = zmat_string_to_cartesian_string(xyz_text)
            # Now parse the cartesian output (element x y z format)
            cart_lines = [line.strip() for line in cartesian_string.strip().split('\n') if line.strip()]
            elements = []
            coords = []
            for line in cart_lines:
                parts = line.split()
                if len(parts) >= 4:
                    elements.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            
            if elements:
                return elements, np.array(coords)
        except Exception as e:
            # If Z-matrix parsing fails, fall through to XYZ parser
            pass
    
    # Try XYZ format
    try:
        # Try standard XYZ format with number of atoms and comment
        num_atoms = int(lines[0].strip())
        # Skip comment line
        atom_lines = lines[2:2+num_atoms]
    except (ValueError, IndexError):
        # Try parsing without header (just element x y z lines)
        atom_lines = lines
    
    elements = []
    coords = []
    
    for line in atom_lines:
        parts = line.split()
        if len(parts) >= 4:
            elements.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    
    if not elements:
        return None, None
    
    return elements, np.array(coords)


def create_mol_from_xyz(elements, coords):
    """
    Create RDKit molecule from elements and coordinates.
    Infers bonds based on distance heuristics for visualization.
    
    Args:
        elements: list of element symbols
        coords: numpy array of shape (n_atoms, 3)
    
    Returns:
        RDKit Mol object or None
    """
    if elements is None or coords is None:
        return None
    
    # Create editable molecule
    mol = Chem.RWMol()
    
    # Add atoms
    for element in elements:
        atom = Chem.Atom(element)
        mol.AddAtom(atom)
    
    # Create conformer with coordinates
    conf = Chem.Conformer(len(elements))
    for i, coord in enumerate(coords):
        conf.SetAtomPosition(i, (float(coord[0]), float(coord[1]), float(coord[2])))
    
    mol.AddConformer(conf)
    
    # Infer bonds based on distances (for visualization)
    pt = Chem.GetPeriodicTable()
    for i in range(len(elements)):
        for j in range(i + 1, len(elements)):
            # Get covalent radii
            r1 = pt.GetRcovalent(elements[i])
            r2 = pt.GetRcovalent(elements[j])
            
            # Calculate distance
            pos1 = conf.GetAtomPosition(i)
            pos2 = conf.GetAtomPosition(j)
            dist = np.sqrt((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2)
            
            # Add bond if distance is within reasonable range (sum of covalent radii + tolerance)
            if dist < (r1 + r2) * 1.3:  # 30% tolerance
                mol.AddBond(i, j, Chem.BondType.SINGLE)
    
    # Convert to regular mol
    mol = mol.GetMol()
    
    return mol


def mol_to_xyz_string(mol):
    """Convert RDKit mol to XYZ format string."""
    if mol is None:
        return ""
    
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    
    xyz_lines = [str(num_atoms), ""]
    
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        element = atom.GetSymbol()
        xyz_lines.append(f"{element} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    
    return "\n".join(xyz_lines)


def calculate_distance(coord1, coord2):
    """Calculate distance between two points."""
    return np.linalg.norm(np.array(coord1) - np.array(coord2))


def calculate_angle(coord1, coord2, coord3):
    """Calculate angle (in degrees) between three points (1-2-3)."""
    v1 = np.array(coord1) - np.array(coord2)
    v2 = np.array(coord3) - np.array(coord2)
    
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
    angle_rad = np.arccos(cos_angle)
    return np.degrees(angle_rad)


def calculate_dihedral(coord1, coord2, coord3, coord4):
    """Calculate dihedral angle (in degrees) between four points (1-2-3-4)."""
    b1 = np.array(coord2) - np.array(coord1)
    b2 = np.array(coord3) - np.array(coord2)
    b3 = np.array(coord4) - np.array(coord3)
    
    # Normalize b2
    b2_norm = b2 / np.linalg.norm(b2)
    
    # Project b1 and b3 onto plane perpendicular to b2
    v1 = b1 - np.dot(b1, b2_norm) * b2_norm
    v2 = b3 - np.dot(b3, b2_norm) * b2_norm
    
    # Calculate angle
    x = np.dot(v1, v2)
    y = np.dot(np.cross(b2_norm, v1), v2)
    
    return np.degrees(np.arctan2(y, x))


def get_dark_cpk_color(element):
    """Get a darkened version of CPK colors for an element."""
    dark_colors = {
        'H': '#999999',
        'C': '#333333',
        'N': '#000080',
        'O': '#8B0000',
        'F': '#004D00',
        'P': '#CC6600',
        'S': '#999900',
        'Cl': '#006600',
        'Br': '#660000',
        'I': '#330066',
    }
    return dark_colors.get(element, '#404040')


def visualize_molecules(mol1, mol2, show_mol1, show_mol2, show_labels=False):
    """
    Create py3Dmol visualization of molecules.
    
    Returns:
        tuple: (py3Dmol view object, atom coordinate mapping)
    """
    view = py3Dmol.view(width=800, height=600)
    view.setBackgroundColor('#1a1a1a')  # Dark background
    
    model_count = 0
    atom_coords = {}
    global_atom_idx = 1
    
    if mol1 and show_mol1:
        xyz1 = mol_to_xyz_string(mol1)
        view.addModel(xyz1, 'xyz')
        
        # Molecule 1 with darkened CPK colors
        elements_in_mol = set([atom.GetSymbol() for atom in mol1.GetAtoms()])
        for element in elements_in_mol:
            color = get_dark_cpk_color(element)
            view.setStyle(
                {'model': model_count, 'elem': element},
                {
                    'stick': {'color': color, 'radius': 0.15},
                    'sphere': {'scale': 0.25, 'color': color}
                }
            )
        
        # Add labels and store coordinates
        conf = mol1.GetConformer()
        for i, atom in enumerate(mol1.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_coords[global_atom_idx] = (pos.x, pos.y, pos.z, atom.GetSymbol(), 1)
            if show_labels:
                view.addLabel(str(global_atom_idx), {
                    'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
                    'backgroundColor': 'black',
                    'fontColor': 'white',
                    'fontSize': 10,
                    'showBackground': True,
                    'backgroundOpacity': 0.7
                })
            global_atom_idx += 1
        
        model_count += 1
    
    if mol2 and show_mol2:
        xyz2 = mol_to_xyz_string(mol2)
        view.addModel(xyz2, 'xyz')
        
        # Molecule 2 with bright CPK colors
        view.setStyle({'model': model_count}, {
            'stick': {'colorscheme': 'default', 'radius': 0.15},
            'sphere': {'scale': 0.25, 'colorscheme': 'default'}
        })
        
        # Add labels and store coordinates
        conf = mol2.GetConformer()
        for i, atom in enumerate(mol2.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_coords[global_atom_idx] = (pos.x, pos.y, pos.z, atom.GetSymbol(), 2)
            if show_labels:
                view.addLabel(str(global_atom_idx), {
                    'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
                    'backgroundColor': 'navy',
                    'fontColor': 'white',
                    'fontSize': 10,
                    'showBackground': True,
                    'backgroundOpacity': 0.7
                })
            global_atom_idx += 1
        
        model_count += 1
    
    if model_count > 0:
        view.zoomTo()
    
    return view, atom_coords


def main():
    st.set_page_config(page_title="Double Molecule Viewer", layout="wide")
    
    # Reduce top margin
    st.markdown("""
        <style>
        .block-container {
            padding-top: 1rem;
            padding-bottom: 0rem;
        }
        </style>
        """, unsafe_allow_html=True)
    
    st.markdown("# Double Molecule Viewer")
    
    # Initialize session state
    if 'mol1' not in st.session_state:
        st.session_state.mol1 = None
        # Auto-load default molecule 1
        default_xyz1 = "O  0.0000  0.0000  0.0000\nH  0.7572  0.5865  0.0000\nH -0.7572  0.5865  0.0000"
        elements, coords = parse_xyz(default_xyz1)
        if elements is not None:
            st.session_state.mol1 = create_mol_from_xyz(elements, coords)
            st.session_state.mol1_original = create_mol_from_xyz(elements, coords)
    
    if 'mol2' not in st.session_state:
        st.session_state.mol2 = None
        # Auto-load default molecule 2
        default_xyz2 = "O  1.0000  1.0000  0.0000\nH  1.7572  1.5865  0.0000\nH  0.2428  1.5865  0.0000"
        elements, coords = parse_xyz(default_xyz2)
        if elements is not None:
            st.session_state.mol2 = create_mol_from_xyz(elements, coords)
            st.session_state.mol2_original = create_mol_from_xyz(elements, coords)
    if 'mol1_original' not in st.session_state:
        # Set to same as mol1 if it was auto-loaded
        st.session_state.mol1_original = st.session_state.mol1
    if 'mol2_original' not in st.session_state:
        # Set to same as mol2 if it was auto-loaded
        st.session_state.mol2_original = st.session_state.mol2
    if 'show_mol1' not in st.session_state:
        st.session_state.show_mol1 = True
    if 'show_mol2' not in st.session_state:
        st.session_state.show_mol2 = True
    if 'use_aligned' not in st.session_state:
        st.session_state.use_aligned = False
    if 'rmsd' not in st.session_state:
        st.session_state.rmsd = None
    if 'mol2_aligned' not in st.session_state:
        st.session_state.mol2_aligned = None
    if 'selected_atoms' not in st.session_state:
        st.session_state.selected_atoms = []
    if 'atom_coords_cache' not in st.session_state:
        st.session_state.atom_coords_cache = {}
    
    # Create left sidebar and right visualization area
    left_col, right_col = st.columns([1, 2])
    
    with left_col:
        # Molecule 1 inputs
        st.markdown("**Molecule 1** (Dark)")
        xyz_input1 = st.text_area(
            "XYZ or Z-matrix:",
            height=100,
            key="xyz1",
            value="O  0.0000  0.0000  0.0000\nH  0.7572  0.5865  0.0000\nH -0.7572  0.5865  0.0000",
            label_visibility="collapsed"
        )
        
        col1a, col1b = st.columns(2)
        with col1a:
            if st.button("Load", key="refresh1", use_container_width=True):
                elements, coords = parse_xyz(xyz_input1)
                if elements is not None:
                    st.session_state.mol1 = create_mol_from_xyz(elements, coords)
                    st.session_state.mol1_original = create_mol_from_xyz(elements, coords)
                    st.session_state.use_aligned = False
                else:
                    st.error("Parse failed")
        
        with col1b:
            st.session_state.show_mol1 = st.checkbox("Show", value=st.session_state.show_mol1, key="show1")
        
        # Show atom count for mol1
        if st.session_state.mol1:
            st.caption(f"✓ {st.session_state.mol1.GetNumAtoms()} atoms loaded")
        
        st.markdown("---")
        
        # Molecule 2 inputs
        st.markdown("**Molecule 2** (Bright)")
        xyz_input2 = st.text_area(
            "XYZ or Z-matrix:",
            height=100,
            key="xyz2",
            value="O  1.0000  1.0000  0.0000\nH  1.7572  1.5865  0.0000\nH  0.2428  1.5865  0.0000",
            label_visibility="collapsed"
        )
        
        col2a, col2b = st.columns(2)
        with col2a:
            if st.button("Load", key="refresh2", use_container_width=True):
                elements, coords = parse_xyz(xyz_input2)
                if elements is not None:
                    st.session_state.mol2 = create_mol_from_xyz(elements, coords)
                    st.session_state.mol2_original = create_mol_from_xyz(elements, coords)
                    st.session_state.use_aligned = False
                else:
                    st.error("Parse failed")
        
        with col2b:
            st.session_state.show_mol2 = st.checkbox("Show", value=st.session_state.show_mol2, key="show2")
        
        # Show atom count for mol2
        if st.session_state.mol2:
            st.caption(f"✓ {st.session_state.mol2.GetNumAtoms()} atoms loaded")
        
        st.markdown("---")
        
        # Alignment controls
        st.markdown("**Alignment & RMSD**")
        
        if st.button("Calculate RMSD & Align", type="primary", use_container_width=True):
            if st.session_state.mol1 and st.session_state.mol2:
                mol1 = st.session_state.mol1_original
                mol2 = st.session_state.mol2_original
                
                # Check if molecules have same number of atoms
                if mol1.GetNumAtoms() != mol2.GetNumAtoms():
                    st.error(f"Atom count mismatch! ({mol1.GetNumAtoms()} vs {mol2.GetNumAtoms()})")
                else:
                    # Create a copy of mol2 for alignment
                    mol2_copy = Chem.Mol(mol2)
                    
                    # Align mol2 to mol1 and get RMSD
                    rmsd = rdMolAlign.AlignMol(mol2_copy, mol1)
                    
                    st.session_state.rmsd = rmsd
                    st.session_state.mol2_aligned = mol2_copy
            else:
                st.error("Load both molecules first!")
        
        if st.session_state.mol2_aligned is not None:
            st.session_state.use_aligned = st.checkbox(
                "Use Aligned Coordinates",
                value=st.session_state.use_aligned,
                key="use_aligned_checkbox"
            )
        
        if st.session_state.rmsd is not None:
            st.metric("RMSD", f"{st.session_state.rmsd:.4f} Å")
    
    with right_col:
        # Visualization controls
        viz_col1, viz_col2 = st.columns([3, 1])
        with viz_col1:
            st.markdown("**3D Visualization**")
        with viz_col2:
            show_labels = st.checkbox("Show Labels", value=False, key="show_labels")
        
        # Determine which molecules to show
        mol1_display = st.session_state.mol1 if not st.session_state.use_aligned else st.session_state.mol1_original
        mol2_display = st.session_state.mol2 if not st.session_state.use_aligned else st.session_state.mol2_aligned
        
        if mol1_display or mol2_display:
            # Render viewer with labels option
            view, atom_coords = visualize_molecules(
                mol1_display,
                mol2_display,
                st.session_state.show_mol1,
                st.session_state.show_mol2,
                show_labels
            )
            st.session_state.atom_coords_cache = atom_coords
            st.components.v1.html(view._make_html(), height=540, scrolling=False)
            
            # Measurement tool at bottom (tight spacing)
            st.markdown("")
            meas_col1, meas_col2 = st.columns([1, 2])
            
            with meas_col1:
                st.markdown("**Measurement Tool**")
                atom_input = st.text_input(
                    "Atom numbers (comma-separated):",
                    value="",
                    key="atom_input",
                    placeholder="e.g., 1, 2, 3, 4",
                    help="Enable 'Show Labels' above to see atom numbers, then type them here"
                )
                
                # Parse atom input
                try:
                    if atom_input.strip():
                        atoms = [int(x.strip()) for x in atom_input.split(',') if x.strip()]
                        st.session_state.selected_atoms = atoms
                    else:
                        st.session_state.selected_atoms = []
                except ValueError:
                    st.error("Invalid input")
                
                if st.button("Clear", key="clear_atoms", use_container_width=True):
                    st.session_state.selected_atoms = []
                    st.rerun()
            
            with meas_col2:
                st.markdown("**Measurements**")
                selected = st.session_state.selected_atoms
                coords_cache = st.session_state.atom_coords_cache
                
                if len(selected) >= 2 and all(a in coords_cache for a in selected):
                    measurements = []
                    
                    # Distances
                    if len(selected) >= 2:
                        coord1 = coords_cache[selected[0]][:3]
                        coord2 = coords_cache[selected[1]][:3]
                        dist = calculate_distance(coord1, coord2)
                        measurements.append(f"**{selected[0]}-{selected[1]}:** {dist:.3f} Å")
                    
                    if len(selected) >= 3:
                        coord2 = coords_cache[selected[1]][:3]
                        coord3 = coords_cache[selected[2]][:3]
                        dist = calculate_distance(coord2, coord3)
                        angle = calculate_angle(coords_cache[selected[0]][:3], coord2, coord3)
                        measurements.append(f"**{selected[1]}-{selected[2]}:** {dist:.3f} Å")
                        measurements.append(f"**∠{selected[0]}-{selected[1]}-{selected[2]}:** {angle:.2f}°")
                    
                    if len(selected) >= 4:
                        coord3 = coords_cache[selected[2]][:3]
                        coord4 = coords_cache[selected[3]][:3]
                        dist = calculate_distance(coord3, coord4)
                        angle2 = calculate_angle(coords_cache[selected[1]][:3], coord3, coord4)
                        dihedral = calculate_dihedral(
                            coords_cache[selected[0]][:3],
                            coords_cache[selected[1]][:3],
                            coord3,
                            coord4
                        )
                        measurements.append(f"**{selected[2]}-{selected[3]}:** {dist:.3f} Å")
                        measurements.append(f"**∠{selected[1]}-{selected[2]}-{selected[3]}:** {angle2:.2f}°")
                        measurements.append(f"**⦜{selected[0]}-{selected[1]}-{selected[2]}-{selected[3]}:** {dihedral:.2f}°")
                    
                    for m in measurements:
                        st.markdown(m)
                elif len(selected) > 0:
                    st.caption("Need at least 2 atoms")
                else:
                    st.caption("📊 Supports XYZ & Z-matrix formats")
        else:
            st.info("👈 Load molecules to see visualization")


if __name__ == "__main__":
    main()

