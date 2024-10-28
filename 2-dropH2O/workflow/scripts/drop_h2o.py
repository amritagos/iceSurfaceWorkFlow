from pathlib import Path
from typing import List
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


def get_dropped_water_coordinates(
    site_indices: List[int], surface_system: Atoms, shifting_distance=2.0
) -> np.ndarray:
    """Get the coordinates of the position where a water molecule should be dropped, onto a site defined by site_indices. The water molecule is placed shifting_distance Angstrom above this.
    WARNING: Assumes the box starts from 0,0,0

    Args:
        site_indices (List[int]): site indices that describe the site (a trio of O atoms)
        surface_system (Atoms): Atoms object for the surface

    Returns:
        np.ndarray: coordinates where the water molecule should go
    """
    surface_system_pos = surface_system.get_positions()
    pos0 = surface_system_pos[site_indices[0]]
    pos1 = surface_system_pos[site_indices[1]]
    pos2 = surface_system_pos[site_indices[2]]

    # Smallest box dimension
    box_lengths = surface_system.cell.cellpar()[:3]
    smallest_box_length = np.min(box_lengths)

    # Check distance between 0 and 1
    if (
        surface_system.get_distance(site_indices[0], site_indices[1], mic=False)
        > 0.5 * smallest_box_length
    ):
        # Shift point 1 around, keeping 0 constant
        pos1 = shift_atom(pos1, pos0, box_lengths)
    # Check the distance between 0 and 2
    if (
        surface_system.get_distance(site_indices[0], site_indices[2], mic=False)
        > 0.5 * smallest_box_length
    ):
        # Shift point 2 around, keeping 0 constant
        pos2 = shift_atom(pos2, pos0, box_lengths)
    # Get the average position :
    water_site_pos = np.mean(np.array([pos0, pos1, pos2]), axis=0)
    # Check that the new site position is inside the box
    # WARNING: assuming box starts from 0,0,0 (also don't do this for the z coordinate)
    for i in range(3):
        if surface_system.pbc[i] == True:
            water_site_pos[i] = water_site_pos[i] % box_lengths[i]
    # water_site_pos = water_site_pos % box_lengths
    # Shift the z coordinate by 2 Angstrom
    water_site_pos = water_site_pos + np.array([0, 0, shifting_distance])
    return water_site_pos


def main(
    input_surface_file: Path,
    site_file: Path,
    seed: int,
    site_type: str,
    site_idx: int,
    surface_energy: float,
    metadata_file: Path,
    out_xyz: Path,
):

    # Read the surface using ASE
    surface_system = read(input_surface_file, format="extxyz")

    # Do the computation
    all_sites = np.load(site_file)
    current_site_ind = all_sites[site_idx]
    # Drop the water molecule onto the site, 2 Angstrom above the averaged z coordinate of the trio (denoting the site)
    water_site = get_dropped_water_coordinates(
        current_site_ind, surface_system, shifting_distance=2.0
    )
    # Create the new water molecule
    single_water = molecule("H2O")
    displacement = (
        water_site - single_water[0].position
    )  # Assuming O is the center of mass
    single_water.translate(displacement)
    # Shift the water molecule to the desired position
    system = surface_system + single_water
    system.cell = surface_system.cell
    system.pbc = surface_system.pbc
    # Set a constraint such that the frozen atoms are frozen (tag=2)
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

    from ase.optimize import BFGS

    dyn = BFGS(system)
    dyn.run(fmax=0.01, steps=0)

    # Write out the final XYZ file
    write(out_xyz, system)

    # energy = 0
    energy = system.get_potential_energy()
    system._calc.get_dipole_moment()  # Calculate the dipole moment
    dipole_moment_magnitude = np.linalg.norm(system._calc.results["dipole"])

    # Write results into metadata JSON file
    with open(metadata_file, "w") as f:
        f.write(
            json.dumps(
                dict(
                    seed=seed,
                    site_type=site_type,
                    site_index=site_idx,
                    surface_energy=surface_energy,
                    system_energy=energy,
                    dipole_moment=dipole_moment_magnitude,
                ),
                indent=4,
            )
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_surface_path", type=Path)
    parser.add_argument("--site_path", type=Path)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--site_type", type=str)
    parser.add_argument("--site_index", type=int)
    parser.add_argument("--surface_energy", type=float)
    parser.add_argument("--out_metadata", type=Path)
    parser.add_argument("--out_xyz", type=Path)

    args = parser.parse_args()

    main(
        args.input_surface_path,
        args.site_path,
        args.seed,
        args.site_type,
        args.site_index,
        args.surface_energy,
        args.out_metadata,
        args.out_xyz,
    )
