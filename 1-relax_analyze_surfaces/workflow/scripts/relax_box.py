from pathlib import Path
import numpy as np
from ase.io import read, write
import sys
import json

sys.path.insert(1, "/home/amritagos/Git/Gitlab/scmecpp/")
from ase_interface import SCME_PS


def main(
    input_xyz: Path,
    seed: int,
    output_xyz: Path,
    results_file: Path,
):

    # Read the system using ASE
    system = read(input_xyz, format="extxyz")
    system.set_pbc(True)

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

    dyn = BFGS(system)  # Don't write out a trajectory
    dyn.run(fmax=0.01)

    energy = system.get_potential_energy()
    n_water_molecules = len(system) / 3
    cohesive_energy = energy / n_water_molecules

    system._calc.get_dipole_moment()  # Calculate the dipole moment
    dipole_moment_magnitude = np.linalg.norm(system._calc.results["dipole"])

    # Write out the final XYZ file
    write(output_xyz, system)

    # Write results file as json or something easily parseable
    with open(results_file, "w") as f:
        f.write(
            json.dumps(
                dict(
                    seed=seed,
                    cohesive_energy=cohesive_energy,
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
    parser.add_argument("output_xyz", type=Path)
    parser.add_argument("results_file", type=Path)

    args = parser.parse_args()

    main(
        args.input_xyz,
        args.seed,
        args.output_xyz,
        args.results_file,
    )
