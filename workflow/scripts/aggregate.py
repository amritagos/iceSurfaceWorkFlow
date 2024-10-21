from pathlib import Path
import sys
import json
from typing import List 
import pandas as pd

def main(
    results_path_a: List[Path],
    results_path_b: List[Path],
    output_path: Path,
):

    df = pd.DataFrame()
    water_energy = 3.6474808763708796e-06

    results = []

    print(results_path_b)

    for p in results_path_a:
        with open(p,"r") as f:
            res = json.load(f)
            relative_energy = res["energy"] - res["surface_energy"] - water_energy
            results.append( [ relative_energy, "A"] )
    
    for p in results_path_b:
        with open(p,"r") as f:
            res = json.load(f)
            relative_energy = res["energy"] - res["surface_energy"] - water_energy
            results.append( [ relative_energy, "B"] )

    df = pd.DataFrame(results, columns=["relative_energy", "site_type"])
    df.to_csv( output_path )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--results_a", type=Path, nargs="+")
    parser.add_argument("--results_b", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)

    args = parser.parse_args()

    main(
        args.results_a,
        args.results_b,
        args.output
    )
