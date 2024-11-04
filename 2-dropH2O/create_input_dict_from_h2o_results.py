from pathlib import Path
import re
import json
import pandas as pd
import numpy as np

drop_h2O_results_path = Path(
    "/home/amritagos/Git/Github/iceSurfaceWorkFlow/2-dropH2O/results/data"
)
dict_output_path = Path("./input_dict_from_h2o.json")

input_dict = dict()
for path in drop_h2O_results_path.glob("*/metadata.json"):
    with open(path, "r") as f:
        res = json.load(f)
        seed = res["seed"]
        site_type = res["site_type"]
        site_index = res["site_index"]
        surface_energy = res["surface_energy"]
        input_dict[f"{seed}_{site_type}_{site_index}"] = [
            seed,
            site_type,
            site_index,
            surface_energy,
        ]


with open(dict_output_path, "w") as f:
    f.write(json.dumps(input_dict, indent=4))
