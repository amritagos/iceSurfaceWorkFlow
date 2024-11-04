from pathlib import Path
from typing import List, Tuple
from ase import Atoms
import numpy as np
from ase.io import read, write
from ase.build import molecule
import sys
import json
from ase.constraints import FixAtoms

sys.path.insert(1, "/home/amritagos/Git/Gitlab/scmecpp/")
from ase_interface import SCME_PS


def shift_atom(
    point: np.ndarray, reference_point: np.ndarray, box_dimensions: np.ndarray
) -> np.ndarray:
    """Shift point coordinates, so that it is an unwrapped coordinate with respect to the reference point

    Args:
        point (np.ndarray): Position to be shifted
        reference_point (np.ndarray): Position of the reference point
        box_dimensions (np.ndarray): box lengths

    Returns:
        shifted_point: Shifted coordinates of the point
    """
    delta = point - reference_point
    for i in range(len(delta)):
        delta[i] -= box_dimensions[i] * round(delta[i] / box_dimensions[i])

    shifted_point = reference_point + delta
    return shifted_point


def sanitize_box(atoms: Atoms) -> None:
    """Make sure there are no broken H bonds. Assumes order is O H H where H atoms are bonded to O

    Args:
        atoms (Atoms): Atoms object
    """
    box_lengths = atoms.cell.cellpar()[:3]

    for atom in atoms:
        if atom.symbol == "O":
            o_idx = atom.index
            h1_idx = o_idx + 1
            h2_idx = o_idx + 2
            o_pos = atoms[o_idx].position
            atoms[h1_idx].position[:2] = shift_atom(
                atoms[h1_idx].position, o_pos, box_lengths
            )[:2]
            atoms[h2_idx].position[:2] = shift_atom(
                atoms[h2_idx].position, o_pos, box_lengths
            )[:2]
    all_pos = atoms.get_positions()
    # min_coord = np.min(all_pos, axis=0)
    # new_pos = all_pos - min_coord
    # atoms.set_positions(new_pos)


def shift_admolecule(admolecule_idx: int, atoms: Atoms) -> None:
    """Shift all coordinates so that the admolecule is in the center of the surface

    Args:
        admolecule_idx (int): O atom index in the water molecule to shift
        atoms (Atoms): Atoms object containing the admolecule
    """
    box_lengths = atoms.cell.cellpar()[:3]
    box_center = np.array(
        [0.5 * box_lengths[0], 0.5 * box_lengths[1], 0.5 * box_lengths[2]]
    )

    # Shift everything by moving the O atom into the center
    displacement = box_center - atoms[admolecule_idx].position
    displacement[2] = 0.0
    all_positions = atoms.get_positions()
    new_positions = all_positions + displacement
    new_positions[:, :2] = new_positions[:, :2] % box_lengths[:2]
    atoms.set_positions(new_positions)


def main(
    input_surface_path: Path,
    displacement_path: Path,
    input_metadata: Path,
    out_metadata: Path,
    out_xyz: Path,
):
    # Get the metadata
    with open(input_metadata, "r") as f:
        results = json.load(f)
    # Read the surface using ASE
    system = read(input_surface_path, format="extxyz")

    # Set PBCs to false and put admolecule in center
    system.pbc = [False, False, False]
    displacement = np.load(displacement_path)
    new_positions = system.get_positions() + displacement[:-3] # skip the last three, because we have no adatom
    system.set_positions(new_positions)
    # Also make sure there are no broken H bonds
    sanitize_box(system)

    # Set a constraint such that the frozen atoms are frozen (tag=2) to be safe
    freeze = FixAtoms(mask=[atom.tag == 2 for atom in system])
    system.set_constraint(freeze)

    system.calc = SCME_PS(
        system,
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

    energy = system.get_potential_energy()
    system._calc.get_dipole_moment()  # Calculate the dipole moment
    dipole_moment_magnitude = np.linalg.norm(system._calc.results["dipole"])

    # Write out the XYZ file after shifting etc
    out_xyz.parent.mkdir(exist_ok=True, parents=True)
    write(out_xyz, system)

    results["surface_energy"] = energy

    # Write results into metadata JSON file
    with open(out_metadata, "w") as f:
        f.write(json.dumps(results, indent=4))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_surface_path", type=Path)
    parser.add_argument("--displacement_path", type=Path)
    parser.add_argument("--input_metadata", type=Path)
    parser.add_argument("--out_metadata", type=Path)
    parser.add_argument("--out_xyz", type=Path)

    args = parser.parse_args()

    main(
        args.input_surface_path,
        args.displacement_path,
        args.input_metadata,
        args.out_metadata,
        args.out_xyz,
    )
