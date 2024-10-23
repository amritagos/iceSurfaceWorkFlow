from ase import Atoms
from ase.io import read, write
from pathlib import Path
import numpy as np


def parse_non_standard_xyz(file_content):
    lines = file_content.strip().split("\n")

    # Initialize lists to hold atomic data
    symbols = []
    positions = []

    # Variables to hold cell vectors
    cell = np.zeros((3, 3))

    # Iterate through lines
    for i, line in enumerate(lines):
        if line.startswith("#") or line.startswith("%PBC") or line.isdigit():
            # Skip comments, PBC, and the atom count line
            continue
        elif line.startswith("Vector1"):
            # Parse the first cell vector
            cell[0] = np.array([float(x) for x in line.split()[1:]])
        elif line.startswith("Vector2"):
            # Parse the second cell vector
            cell[1] = np.array([float(x) for x in line.split()[1:]])
        elif line.startswith("Vector3"):
            # Parse the third cell vector
            cell[2] = np.array([float(x) for x in line.split()[1:]])
        elif line.startswith("Offset"):
            continue
        elif len(line.split()) == 4:
            # Parse atomic data
            tokens = line.split()
            symbols.append(tokens[0])
            positions.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])

    # Convert positions to a numpy array
    positions = np.array(positions)
    # Create the ASE Atoms object
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    return atoms, cell


def main(input_path: Path, output_path: Path):
    with open(input_path) as f:
        file_content = f.read()
        atoms, cell = parse_non_standard_xyz(file_content)

    # Swap y and z
    for atom in atoms:
        new_y = atom.position[2]
        new_z = atom.position[1]
        atom.position[1] = new_y
        atom.position[2] = new_z
    old_cell_lengths = np.diagonal(cell)
    atoms.set_cell([old_cell_lengths[0], old_cell_lengths[2], old_cell_lengths[1]])
    atoms.set_pbc(True)
    write(output_path, atoms, format="extxyz")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_path", type=Path)
    parser.add_argument("output_path", type=Path)

    args = parser.parse_args()

    main(
        args.input_path,
        args.output_path,
    )
