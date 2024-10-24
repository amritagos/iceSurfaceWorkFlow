from pathlib import Path
import sys
import json
from typing import List
import pandas as pd


def main(
    results_path: List[Path],
    output_path: Path,
):

    df = pd.DataFrame()

    results = []

    for p in results_path:
        with open(p, "r") as f:
            res = json.load(f)
            results.append([res["seed"], res["cohesive_energy"], res["dipole_moment"]])

    df = pd.DataFrame(results, columns=["seed", "cohesive_energy", "dipole_moment"])
    df.to_csv(output_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--results", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    main(args.results, args.output)
