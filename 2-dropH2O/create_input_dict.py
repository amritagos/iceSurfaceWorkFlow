from pathlib import Path
import re
import json
import pandas as pd
import numpy as np

relaxed_surface_results_dir = Path(
    "/home/amritagos/Git/Github/iceSurfaceWorkFlow/1-relax_analyze_surfaces/results"
)
dict_output_path = Path("./input_dict.json")


def extract_seeds_from_directory(base_dir):
    # Convert the base directory to a Path object
    base_path = Path(base_dir)
    seeds = []

    # Iterate over all files and subdirectories under the base directory
    for path in base_path.rglob("*"):
        # Match directories with the pattern 'surface_relax/1' and extract the integer at the end
        match = re.search(r"surface_relax/(\d+)$", str(path))

        if match:
            seeds.append(
                int(match.group(1))
            )  # Extract and store the seed as an integer

    return seeds


# Get the aggregated information for the surfaces
surface_aggregate_path = relaxed_surface_results_dir / "surfaces_aggregated.csv"
surface_data_df = pd.read_csv(surface_aggregate_path)

SEEDS = extract_seeds_from_directory(relaxed_surface_results_dir)

# Put everything you need (seed, site type, site index) etc into this dictionary
input_dict = dict()

# Loop through the seeds and update the dictionary
site_types = ["a", "b", "c"]

for seed in SEEDS:
    # Surface energy for each surface
    surface_energy = surface_data_df.loc[
        surface_data_df["seed"] == seed, "surface_energy"
    ].values[0]
    # Loop through every site type
    for site_type in site_types:
        # Get the number of sites for each site type
        path_to_site_file = (
            relaxed_surface_results_dir
            / f"surface_relax/{seed}/{site_type}_sites_seed_{seed}.npy"
        )
        n_sites = len(np.load(path_to_site_file))
        site_indices = [n for n in range(n_sites)]
        # Go through every site
        for site in site_indices:
            input_dict[f"{seed}_{site_type}_{site}"] = [
                seed,
                site_type,
                site,
                surface_energy,
            ]


with open(dict_output_path, "w") as f:
    f.write(json.dumps(input_dict, indent=4))
