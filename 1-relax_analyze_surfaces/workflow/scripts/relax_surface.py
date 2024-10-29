from pathlib import Path
import numpy as np
from ase.io import read, write
import sys
import json
from ase.constraints import FixAtoms

sys.path.insert(1, "/home/amritagos/Git/Gitlab/scmecpp/")
from ase_interface import SCME_PS


def sort_oxygen_by_z(o_symbol: str, atoms):
    """Get a list of sorted oxygen indices (with respect to the original atoms object)

    Args:
        o_symbol(str) : Oxygen atom symbol
        atoms (Atoms): ASE Atoms object
    """
    positions = atoms.get_positions()
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == o_symbol]
    oxygen_pos = positions[oxygen_indices]

    # Sort oxygen atoms based on their z coordinate
    ind = np.argsort(oxygen_pos[:, 2])

    sorted_oxygen_indices = [oxygen_indices[i] for i in ind]
    return sorted_oxygen_indices


def main(
    input_xyz: Path,
    seed: int,
    vaccuum_depth: float,
    n_bilayers: int,
    n_frozen_bilayers: int,
    o_symbol: str,
    output_xyz: Path,
    results_file: Path,
):

    # Read the surface using ASE
    surface = read(input_xyz, format="extxyz")
    # Add vaccuum in the z dimension
    cell_lengths = surface.cell.cellpar()[:3]
    cell_lengths[2] += vaccuum_depth
    surface.set_cell(cell_lengths)
    surface.pbc = [True, True, False]
    # Pre-processing
    n_water_molecules = len(surface) / 3
    n_waters_bilayer = (
        n_water_molecules / n_bilayers
    )  # Number of water molecules per bilayer
    n_waters_frozen = (
        n_frozen_bilayers * n_waters_bilayer
    )  # Number of water molecules to freeze
    # Frozen atoms will have tag 2 and moving atoms will have a tag of 0 (other moving atoms) or 1 (topmost top layer)
    sorted_o_ind = sort_oxygen_by_z(o_symbol, surface)
    # Top layer of the topmost bilayer
    top_indices = sorted_o_ind[-int(n_waters_bilayer / 2) :]
    # Set the tag of these atoms to 1
    for o_ind in top_indices:
        surface[o_ind].tag = 1
    # Find a distance cutoff for the frozen layer
    if n_frozen_bilayers > 0:
        frozen_bilayer_height = surface.get_positions()[
            sorted_o_ind[int(n_waters_frozen) - 1]
        ][
            2
        ]  # highest coordinate of the frozen bilayers
        next_bilayer_height = surface.get_positions()[
            sorted_o_ind[int(n_waters_frozen)]
        ][
            2
        ]  # lowest coordinate of the bilayer above it
        z_cutoff_frozen = 0.5 * (
            frozen_bilayer_height + next_bilayer_height
        )  # The average of this should be sufficient as a cutoff. All water molecules below this should be frozen
    else:
        z_cutoff_frozen = 0.0  # If no bilayers must be frozen
    # Set the tags of all frozen atoms to 2
    for atom in surface:
        if atom.tag == 0:
            if atom.position[2] < z_cutoff_frozen:
                atom.tag = 2
    # Set a constraint such that the frozen atoms are frozen (tag=2)
    freeze = FixAtoms(mask=[atom.tag == 2 for atom in surface])
    surface.set_constraint(freeze)

    surface.calc = SCME_PS(
        surface,
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

    dyn = BFGS(surface)  # Don't write out a trajectory
    dyn.run(fmax=0.01)

    energy = surface.get_potential_energy()

    surface._calc.get_dipole_moment()  # Calculate the dipole moment
    dipole_moment_magnitude = np.linalg.norm(surface._calc.results["dipole"])

    # Write out the final XYZ file
    write(output_xyz, surface)

    # Write results file as json or something easily parseable
    with open(results_file, "w") as f:
        f.write(
            json.dumps(
                dict(
                    seed=seed,
                    surface_energy=energy,
                    # dipole_moment_vector=system._calc.results["dipole"].tolist(),
                    dipole_moment=dipole_moment_magnitude,
                ),
                indent=4,
            )
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_xyz", type=Path)
    parser.add_argument("seed", type=int)
    parser.add_argument("vaccuum_depth", type=float)
    parser.add_argument("n_bilayers", type=int)
    parser.add_argument("n_frozen_bilayers", type=int)
    parser.add_argument("o_symbol", type=str)
    parser.add_argument("output_xyz", type=Path)
    parser.add_argument("results_file", type=Path)

    args = parser.parse_args()

    main(
        args.input_xyz,
        args.seed,
        args.vaccuum_depth,
        args.n_bilayers,
        args.n_frozen_bilayers,
        args.o_symbol,
        args.output_xyz,
        args.results_file,
    )
