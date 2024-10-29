from pathlib import Path
import sys
import json
from typing import List
import pandas as pd


def main(
    result_paths: List[Path],
    isolated_water_energy: float,
    output_csv: Path,
):

    df = pd.DataFrame()

    results = []

    for p in result_paths:
        with open(p, "r") as f:
            res = json.load(f)
            site_type = res["site_type"]
            seed = res["seed"]
            site_index = res["site_index"]
            surface_energy = res["surface_energy"]
            system_energy = res["system_energy"]
            dipole_moment = res["dipole_moment"]
            final_site_type = (res["final_site_type"],)
            site_moved = res["site_moved"]
            binding_energy = system_energy - surface_energy - isolated_water_energy
            results.append(
                [
                    site_type,
                    seed,
                    site_index,
                    binding_energy,
                    dipole_moment,
                    system_energy,
                    surface_energy,
                    final_site_type,
                    site_moved,
                ]
            )

    df = pd.DataFrame(
        results,
        columns=[
            "site_type",
            "seed",
            "site_index",
            "binding_energy",
            "dipole_moment",
            "system_energy",
            "isolated_surface_energy",
            "final_site_type",
            "site_moved",
        ],
    )
    df.to_csv(output_csv)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--result_paths", type=Path, nargs="+")
    parser.add_argument("--isolated_water_energy", type=float)
    parser.add_argument("--output_csv", type=Path)

    args = parser.parse_args()

    main(args.result_paths, args.isolated_water_energy, args.output_csv)
