from pathlib import Path
import sys
import json
from typing import List
import pandas as pd


def main(
    misc_metadata_path: List[Path],
    site_metadata_path: List[Path],
    output_metadata_path: Path,
    output_csv_path: Path,
):

    df = pd.DataFrame()

    results = []

    all_a_sites = 0
    all_b_sites = 0
    all_c_sites = 0
    all_d_sites = 0
    # Expected A, B, C, D site percentages
    out_metadata_dict = dict(
        {
            "exp_a_site_percent": 36.75,
            "exp_b_site_percent": 36.75,
            "exp_c_site_percent": 12.25,
            "exp_d_site_percent": 12.25,
        }
    )

    for p1, p2 in zip(misc_metadata_path, site_metadata_path):
        # read in the metadata from the surface relaxation
        with open(p1, "r") as f:
            misc = json.load(f)
        # Read in the metadata from the analysis
        with open(p2, "r") as f:
            sites = json.load(f)
        num_a_sites = sites["num_a_sites"]
        num_b_sites = sites["num_b_sites"]
        num_c_sites = sites["num_c_sites"]
        num_d_sites = sites["num_d_sites"]
        all_a_sites += num_a_sites
        all_b_sites += num_b_sites
        all_c_sites += num_c_sites
        all_d_sites += num_d_sites
        # Add to the results list of lists (will be written out to a CSV later)
        results.append(
            [
                misc["seed"],
                misc["surface_energy"],
                misc["dipole_moment"],
                num_a_sites,
                num_b_sites,
                num_c_sites,
                num_d_sites,
            ]
        )

    # Write out to a CSV file
    df = pd.DataFrame(
        results,
        columns=[
            "seed",
            "surface_energy",
            "dipole_moment",
            "num_a_sites",
            "num_b_sites",
            "num_c_sites",
            "num_d_sites",
        ],
    )
    df.to_csv(output_csv_path)

    # Calculate percentages
    all_sites = all_a_sites + all_b_sites + all_c_sites + all_d_sites
    out_metadata_dict["a_site_percent"] = all_a_sites * 100.0 / all_sites
    out_metadata_dict["b_site_percent"] = all_b_sites * 100.0 / all_sites
    out_metadata_dict["c_site_percent"] = all_c_sites * 100.0 / all_sites
    out_metadata_dict["d_site_percent"] = all_d_sites * 100.0 / all_sites
    out_metadata_dict["num_a_sites"] = all_a_sites
    out_metadata_dict["num_b_sites"] = all_b_sites
    out_metadata_dict["num_c_sites"] = all_c_sites
    out_metadata_dict["num_d_sites"] = all_d_sites
    # Write out to the metadata file
    with open(output_metadata_path, "w") as f:
        f.write(json.dumps(out_metadata_dict, indent=4))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--misc_metadata", type=Path, nargs="+")
    parser.add_argument("--site_metadata", type=Path, nargs="+")
    parser.add_argument("--out_metadata", type=Path)
    parser.add_argument("--out_csv", type=Path)

    args = parser.parse_args()

    main(args.misc_metadata, args.site_metadata, args.out_metadata, args.out_csv)
