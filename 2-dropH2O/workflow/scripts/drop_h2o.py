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


def get_dropped_water_coordinates(
    site_indices: List[int], surface_system: Atoms, shifting_distance
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


def check_final_site_conditions(
    dropped_o_index: int,
    initial_site_indices: np.ndarray,
    system: Atoms,
    o_o_cutoff: float = 3.5,
) -> Tuple[bool, str]:
    """Check whether, after the energy minimization, the water molecule has moved into another site or not

    Args:
        dropped_o_index (int): Index of the dropped water molecule
        initial_site_indices (np.ndarray): Site indices of the site where the water molecule was supposed to be dropped
        system (Atoms): Surface+dropped water molecule ASE Atoms object
        o_o_cutoff (float, optional): Cutoff for neighbouring water molecules. Defaults to 3.5.

    Raises:
        Exception: _description_

    Returns:
        Tuple[bool, str]: _description_
    """
    final_site_indices = []
    # find the oxygens in the top layer of the bilayer
    # oxygens in the top layer of the surface have a tag of 1.
    top_o_indices = [Atom.index for Atom in system if Atom.tag == 1]

    # Find the site for the water molecule (the indices)
    # The indices remain the same since the water molecule is added at the end of surface_system
    for o_ind in top_o_indices:
        if system.get_distance(o_ind, dropped_o_index, mic=True) <= o_o_cutoff:
            final_site_indices.append(o_ind)
    # Throw if the final_site_indices has more or less than three O atom indices
    if len(final_site_indices) != 3:
        raise Exception("The final site has more than three O atom indices in it.")
    # Sort in ascending order
    final_site_indices.sort()

    n_dangling_h = 0
    for o_ind in final_site_indices:
        h_indices = [o_ind + 1, o_ind + 2]
        for h_ind in h_indices:
            if system[h_ind].position[2] > system[o_ind].position[2]:
                n_dangling_h += 1
    # A sites: 1 dangling hydrogen
    if n_dangling_h == 1:
        final_site_type = "a"
    # B sites: 2 dangling hydrogens
    elif n_dangling_h == 2:
        final_site_type = "b"
    # C sites: no dangling hydrogens
    elif n_dangling_h == 0:
        final_site_type = "c"
    # D sites: all three are dangling hydrogens
    elif n_dangling_h == 3:
        final_site_type = "d"

    # Check that the final site indices and the initial site indices are the same
    site_moved = False
    if final_site_indices != initial_site_indices.tolist():
        site_moved = True

    return site_moved, final_site_type


def main(
    input_surface_file: Path,
    site_file: Path,
    seed: int,
    site_type: str,
    site_idx: int,
    surface_energy: float,
    distance_to_surface: float,
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
        current_site_ind, surface_system, distance_to_surface
    )
    # Create the new water molecule
    single_water = molecule("H2O")
    displacement = (
        water_site - single_water[0].position
    )  # Assuming O is the center of mass
    single_water.translate(displacement)
    # Shift the water molecule to the desired position
    dropped_o_index = len(
        surface_system
    )  # Index of the O in the dropped water molecule, in the system Atoms object
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
    dyn.run(fmax=0.01)

    # Write out the final XYZ file
    write(out_xyz, system)

    # energy = 0
    energy = system.get_potential_energy()
    system._calc.get_dipole_moment()  # Calculate the dipole moment
    dipole_moment_magnitude = np.linalg.norm(system._calc.results["dipole"])

    # Check the final site position, and site type
    site_moved, final_site_type = check_final_site_conditions(
        dropped_o_index, current_site_ind, system
    )

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
                    final_site_type=final_site_type,
                    site_moved=site_moved,
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
    parser.add_argument("--distance_to_surface", type=float)
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
        args.distance_to_surface,
        args.out_metadata,
        args.out_xyz,
    )
