from pathlib import Path
import numpy as np
from ase.io import read

def main(input_surface_file : Path, site_file : Path, results_file : Path, site_idx : int, surface_energy: float):
    
    # Read the surface using ASE
    surface_system = read(input_surface_file, format="extxyz")
    
    # Do the computation

    energy = 0

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

    main(args.input_surface_path, args.site_path, args.results_file, args.site_index, args.surface_energy)