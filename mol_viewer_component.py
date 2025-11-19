"""
Custom Streamlit component for clickable molecular viewer
"""
import streamlit.components.v1 as components
from rdkit import Chem
import json


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


def clickable_mol_viewer(mol1, mol2, show_mol1, show_mol2, show_labels, selected_atoms, height=600):
    """
    Create a clickable molecular viewer component.
    
    Returns:
        Clicked atom index (or None if no new click)
    """
    
    # Prepare molecule data
    # Clean selected_atoms - ensure it's a list of integers
    clean_selected = []
    if selected_atoms:
        for a in selected_atoms:
            try:
                clean_selected.append(int(a))
            except (ValueError, TypeError):
                pass  # Skip non-integer values
    
    mol_data = {
        'mol1': None,
        'mol2': None,
        'show_mol1': bool(show_mol1),
        'show_mol2': bool(show_mol2),
        'show_labels': bool(show_labels),
        'selected_atoms': clean_selected
    }
    
    atom_info = {}  # Global atom index -> (mol_id, local_idx, element, x, y, z)
    global_idx = 1
    
    if mol1 and show_mol1:
        xyz1 = mol_to_xyz_string(mol1)
        mol_data['mol1'] = xyz1
        
        # Get dark colors for each element
        dark_colors = {}
        for atom in mol1.GetAtoms():
            element = atom.GetSymbol()
            if element not in dark_colors:
                dark_colors[element] = get_dark_cpk_color(element)
        mol_data['mol1_colors'] = dark_colors
        
        # Store atom info
        conf = mol1.GetConformer()
        for i, atom in enumerate(mol1.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_info[global_idx] = {
                'mol': 1,
                'local_idx': i,
                'element': atom.GetSymbol(),
                'x': float(pos.x),
                'y': float(pos.y),
                'z': float(pos.z)
            }
            global_idx += 1
    
    if mol2 and show_mol2:
        xyz2 = mol_to_xyz_string(mol2)
        mol_data['mol2'] = xyz2
        
        # Store atom info
        conf = mol2.GetConformer()
        for i, atom in enumerate(mol2.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_info[global_idx] = {
                'mol': 2,
                'local_idx': i,
                'element': atom.GetSymbol(),
                'x': float(pos.x),
                'y': float(pos.y),
                'z': float(pos.z)
            }
            global_idx += 1
    
    # Convert atom_info keys to strings for JSON serialization
    mol_data['atom_info'] = {str(k): v for k, v in atom_info.items()}
    
    # Create HTML with embedded JavaScript
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                overflow: hidden;
            }}
            #viewer {{
                width: 100%;
                height: {height}px;
                position: relative;
            }}
        </style>
    </head>
    <body>
        <div id="viewer"></div>
        
        <script>
            const molData = {json.dumps(mol_data)};
            let selectedAtoms = new Set(molData.selected_atoms);
            let viewer;
            
            function initViewer() {{
                const element = document.getElementById('viewer');
                const config = {{ backgroundColor: 'white' }};
                viewer = $3Dmol.createViewer(element, config);
                
                let modelCount = 0;
                const atomMapping = {{}}; // Maps (model, serial) to global atom index
                let globalAtomIdx = 1;
                
                // Add molecule 1
                if (molData.mol1 && molData.show_mol1) {{
                    viewer.addModel(molData.mol1, 'xyz');
                    
                    // Count atoms in mol1
                    const lines = molData.mol1.trim().split('\\n');
                    const numAtoms = parseInt(lines[0]);
                    
                    // Apply dark colors by element
                    for (const [element, color] of Object.entries(molData.mol1_colors)) {{
                        viewer.setStyle(
                            {{model: modelCount, elem: element}},
                            {{
                                stick: {{color: color, radius: 0.15}},
                                sphere: {{scale: 0.25, color: color}}
                            }}
                        );
                    }}
                    
                    // Map atoms
                    for (let i = 0; i < numAtoms; i++) {{
                        const key = modelCount + '_' + i;
                        atomMapping[key] = globalAtomIdx;
                        globalAtomIdx++;
                    }}
                    
                    modelCount++;
                }}
                
                // Add molecule 2
                if (molData.mol2 && molData.show_mol2) {{
                    viewer.addModel(molData.mol2, 'xyz');
                    
                    // Count atoms in mol2
                    const lines = molData.mol2.trim().split('\\n');
                    const numAtoms = parseInt(lines[0]);
                    
                    // Apply bright CPK colors
                    viewer.setStyle(
                        {{model: modelCount}},
                        {{
                            stick: {{colorscheme: 'default', radius: 0.15}},
                            sphere: {{scale: 0.25, colorscheme: 'default'}}
                        }}
                    );
                    
                    // Map atoms
                    for (let i = 0; i < numAtoms; i++) {{
                        const key = modelCount + '_' + i;
                        atomMapping[key] = globalAtomIdx;
                        globalAtomIdx++;
                    }}
                    
                    modelCount++;
                }}
                
                // Add labels if requested
                if (molData.show_labels) {{
                    for (const [atomIdx, info] of Object.entries(molData.atom_info)) {{
                        const bgColor = info.mol === 1 ? 'black' : 'navy';
                        viewer.addLabel(atomIdx, {{
                            position: {{x: info.x, y: info.y, z: info.z}},
                            backgroundColor: bgColor,
                            fontColor: 'white',
                            fontSize: 10,
                            showBackground: true,
                            backgroundOpacity: 0.7
                        }});
                    }}
                }}
                
                // Highlight selected atoms
                for (const atomIdx of selectedAtoms) {{
                    if (molData.atom_info[atomIdx]) {{
                        const info = molData.atom_info[atomIdx];
                        viewer.addSphere({{
                            center: {{x: info.x, y: info.y, z: info.z}},
                            radius: 0.5,
                            color: 'yellow',
                            alpha: 0.5
                        }});
                    }}
                }}
                
                // Set up click handler
                viewer.setClickable({{model: -1}}, true, function(atom, viewer, event, container) {{
                    if (atom) {{
                        const modelIndex = atom.model;
                        const serial = atom.serial;
                        const key = modelIndex + '_' + serial;
                        const globalIdx = atomMapping[key];
                        
                        if (globalIdx) {{
                            // Send click event to Streamlit
                            window.parent.postMessage({{
                                type: 'streamlit:setComponentValue',
                                value: globalIdx
                            }}, '*');
                        }}
                    }}
                }});
                
                viewer.zoomTo();
                viewer.render();
            }}
            
            // Initialize when page loads
            window.addEventListener('load', initViewer);
        </script>
    </body>
    </html>
    """
    
    # Return the component with bidirectional communication
    clicked_atom = components.html(html_content, height=height, scrolling=False)
    
    return clicked_atom

