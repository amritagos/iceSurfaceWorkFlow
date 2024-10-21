from pathlib import Path
import numpy as np
from ase.io import read
from ase.build import molecule
import sys

sys.path.insert(1, "/hpctheochem/amrita/scmecpp/")
from ase_interface import SCME_PS


def main(
    input_surface_file: Path,
    site_file: Path,
    results_file: Path,
    site_idx: int,
    surface_energy: float,
):

    # Read the surface using ASE
    surface_system = read(input_surface_file, format="extxyz")

    # Do the computation
    all_sites = np.load(site_file)
    current_site_ind = all_sites[site_idx]
    # Drop the water molecule onto the site, 2 Angstrom above the averaged z coordinate of the trio (denoting the site)
    surface_system_pos = surface_system.get_positions()
    site_avg_coord = np.mean(
        np.array(
            [
                surface_system_pos[current_site_ind[0]],
                surface_system_pos[current_site_ind[1]],
                surface_system_pos[current_site_ind[2]],
            ]
        ),
        axis=0,
    )
    # the desired z coordinate should be 2 Angstrom above this
    water_site = site_avg_coord + np.array([0.0, 0.0, 2.0])
    # Create the new water molecule
    single_water = molecule("H2O")
    displacement = (
        water_site - single_water[0].position
    )  # Assuming O is the center of mass
    single_water.translate(displacement)
    # Shift the water molecule to the desired position
    system = surface_system + single_water
    system.cell = surface_system.cell
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

    dyn = BFGS(system)
    dyn.run(fmax=0.01)

    energy = system.get_potential_energy()

    # TODO: write results file as json or something easily parseable
    with open(results_file, "w") as f:
        f.write(f"energy = {energy}\n")
        f.write(f"surface energy = {surface_energy}\n")
    pass


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_surface_path", type=Path)
    parser.add_argument("site_path", type=Path)
    parser.add_argument("results_file", type=Path)
    parser.add_argument("site_index", type=int)
    parser.add_argument("surface_energy", type=float)

    args = parser.parse_args()

    main(
        args.input_surface_path,
        args.site_path,
        args.results_file,
        args.site_index,
        args.surface_energy,
    )
