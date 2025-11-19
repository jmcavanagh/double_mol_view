# Double Molecule Viewer with RMSD

A Streamlit web application for visualizing two molecules together and calculating their RMSD (Root Mean Square Deviation) alignment.

## Features

- Load molecules from XYZ coordinates
- Interactive 3D visualization with py3Dmol
- Show/hide individual molecules
- Calculate RMSD and align molecule 2 to molecule 1
- Toggle between original and aligned coordinates
- Color-coded molecules (cyan and magenta)

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

Then open your browser to the URL shown (typically http://localhost:8501)

## XYZ Format

Standard XYZ format:
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

## Notes

- Molecules must have the same number of atoms for RMSD calculation
- Atoms are matched by order in the input
- Coordinates should be in Ångströms (Å)

