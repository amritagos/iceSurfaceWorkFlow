from pathlib import Path
from typing import List, Tuple
from ase import Atoms
import numpy as np
import ase
from ase.io import read, write
from ase.build import molecule
import sys
import json
from ase.constraints import FixAtoms

sys.path.insert(1, "/home/amritagos/Git/Gitlab/scmecpp/")
from ase_interface import SCME_PS


def write_dipoles(xyz_file: Path, output_file: Path, add_vaccuum=0.0) -> None:
    """Writes out the dipole moments for each molecule in a .npy file

    Args:
        xyz_file (Path): Path to the input XYZ file
        output_file (Path): Path to the desired output file
        add_vaccuum (float): Vaccuum to be added to z dimension
    """
    atoms = read(xyz_file, format="extxyz")

    if add_vaccuum > 0:
        cell_lengths = atoms.cell.cellpar()[:3]
        cell_lengths[2] += add_vaccuum
        atoms.set_cell(cell_lengths)

    # Set up the calculator
    atoms.calc = SCME_PS(
        atoms,
        NC=np.array([0, 0, 0]),
        numerical=False,
        dms=True,
        irigidmolecules=False,
        system=np.ones(1) * 1,
    )
    """,
                     te=te,
                     td=td,
                     Ar=Ar,
                     Br=Br,
                     Cr=Cr,
                     rc_Disp=rc_Disp)"""

    energy = atoms.get_potential_energy()  # is this needed?
    dipoles = atoms._calc.dp * ase.units.Debye  # Calculate the dipole moment
    np.save(output_file, dipoles)


def main(
    initial_xyz: Path,
    final_xyz: Path,
    add_vaccuum_to_initial: float,
    initial_out: Path,
    final_out: Path,
):

    # Get the dipoles for the initial XYZ (before relaxation)
    if add_vaccuum_to_initial > 0.0:
        write_dipoles(initial_xyz, initial_out, add_vaccuum_to_initial)
    else:
        write_dipoles(initial_xyz, initial_out)
    # Get the dipoles for the final XYZ (after relaxation)
    write_dipoles(final_xyz, final_out)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--initial_xyz", type=Path)
    parser.add_argument("--final_xyz", type=Path)
    parser.add_argument("--add_vaccuum_to_initial", type=float)
    parser.add_argument("--initial_out", type=Path)
    parser.add_argument("--final_out", type=Path)

    args = parser.parse_args()

    main(
        args.initial_xyz,
        args.final_xyz,
        args.add_vaccuum_to_initial,
        args.initial_out,
        args.final_out,
    )
