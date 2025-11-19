# Double Molecule Viewer with RMSD

A Streamlit web application for visualizing two molecules together and calculating their RMSD (Root Mean Square Deviation) alignment. Perfect for comparing molecular structures in computational chemistry!

## 🎯 Features

- **Dual Input Formats**: Load molecules from XYZ coordinates OR Z-matrix (Fenske-Hall format)
- **Interactive 3D Visualization**: Powered by py3Dmol with dark background
- **Show/Hide Controls**: Toggle visibility of each molecule independently
- **RMSD Alignment**: Calculate RMSD and align molecule 2 to molecule 1 using RDKit
- **Toggle Coordinates**: Switch between original and aligned positions
- **Measurement Tool**: Measure bond lengths, angles, and dihedral angles
- **Color-Coded Display**: Dark CPK colors (molecule 1) vs bright CPK colors (molecule 2)

## 🚀 Quick Start

### Local Installation

```bash
# Create virtual environment
python -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

Then open your browser to http://localhost:8501

### Deploy to Streamlit Cloud

1. Push this repo to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Click "New app" and select your repository
4. Done! Your app will be live at `your-app.streamlit.app`

## 📝 Input Formats

### XYZ Format

Standard format with header:
```
3
Water molecule
O  0.0000  0.0000  0.0000
H  0.7572  0.5865  0.0000
H -0.7572  0.5865  0.0000
```

Or simplified (without header):
```
O  0.0000  0.0000  0.0000
H  0.7572  0.5865  0.0000
H -0.7572  0.5865  0.0000
```

### Z-Matrix Format (Fenske-Hall)

```
atom1
atom2  ref1  distance
atom3  ref1  distance  ref2  angle
atom4  ref1  distance  ref2  angle  ref3  dihedral
```

Example (water):
```
O
H  1  0.96
H  1  0.96  2  104.5
```

Example (H2O2):
```
H
O 1 0.9
O 2 1.4 1 105.0
H 3 0.9 2 105.0 1 120.0
```

## 📊 Measurement Tool

1. Check **"Show Labels"** to display atom numbers
2. Type atom numbers in the measurement box (e.g., `1, 2, 3, 4`)
3. View calculated:
   - **Bond lengths** (2+ atoms)
   - **Bond angles** (3+ atoms)
   - **Dihedral angles** (4+ atoms)

## 🛠️ Technologies

- **Streamlit**: Web framework
- **RDKit**: RMSD calculation and molecular alignment
- **py3Dmol**: 3D molecular visualization
- **NumPy**: Numerical computations

## 📋 Requirements

- Python 3.8+
- See `requirements.txt` for package versions

## ⚠️ Notes

- Molecules must have the same number of atoms for RMSD calculation
- Atoms are matched by order in the input
- Coordinates should be in Ångströms (Å)
- Bonds are inferred automatically from distances

## 📄 License

Open source - feel free to use and modify!

## 🤝 Contributing

Contributions welcome! Please open an issue or PR.

