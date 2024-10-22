from pathlib import Path
from spirit_extras.plotting import Paper_Plot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main(
    input_csv_path: Path,
    output_img_path: Path,
):  
    # Exclusion list 
    baddies = [[1024, 5], [500, 4], [12599, 4], [10000, 5], [1024, 6], [12599, 1], [10000, 2], [500, 9], [500, 6], [12599, 3], [500, 7], [12345, 5], [12599, 9], [10000, 7], [12345, 9], [10000, 11], [500, 9], [12599, 10], [10000, 8], [12345, 10], [10000, 12], [500, 6], [12345, 3], [12599, 7], [500, 10], [10000, 5], [12599, 11], [1024, 3], [500, 3], [1024, 7], [1024, 4], [10000, 10], [1024, 8]]
    # Input values
    cohesive_energy = 0.622 #  in eV 
    factor = -1.0 # to get positive values 
    n_bins = 75

    # Obtain the data
    df = pd.read_csv(input_csv_path)
    drop_indices = []

    for index, row in df.iterrows():
        seed = row["seed"]
        site_idx = row["site"]
        if [seed,site_idx] in baddies:
            drop_indices.append(index)
    
    df.drop(drop_indices)

    df_a_sites = df[df['site_type'] == 'A']
    df_b_sites = df[df['site_type'] == 'B']
    # All relative energies in eV
    energy_a_sites = factor*df_a_sites['relative_energy'].values
    energy_b_sites = factor*df_b_sites['relative_energy'].values
    
    # Lengths in inches in general
    CM = Paper_Plot.cm  # Use this to give lengths in cm...

    # Params
    params = {
        "font.size": 8,
        "font.family": ("Arial", "sans-serif"),
        "mathtext.fontset": "dejavuserif",
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.labelsize": 8,
    }

    # Either set the absolute height here or set the aspect ratio
    # inside apply_absolute_margins

    pplot = Paper_Plot(width=3.25, height=2.75, nrows=1, ncols=1, rcParams=params)

    # Vertical margin: bottom and then top
    # horizontal margin: left and then right
    # Golden ratio aspect ratio=1.618
    # hspace -> height space

    pplot.apply_absolute_margins(
        aspect_ratio=None,
        abs_horizontal_margins=[0.9 * CM, 0.05 * CM],
        abs_vertical_margins=[1.0 * CM, 0.05 * CM],
        abs_wspace=0.0 * CM,
        abs_hspace=0.2 * CM,
    )
    

    print(pplot.info_string())

    # Get the figure and gridpsec objects
    fig = pplot.fig()
    gs = pplot.gs()

    # Actually make the plot
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlabel(r"Binding Energy (eV)")
    ax1.set_ylabel(
        r"Frequency", labelpad=1
    )  # we already handled the x-label with ax1
    ax1.hist(energy_a_sites, bins=n_bins, ls='dashed', lw=3, fc=(0, 0, 1, 0.5), label="A sites") # blue, A sites
    ax1.hist(energy_b_sites, bins=n_bins, ls='dotted', lw=3, fc=(1, 0, 0, 0.5), label="B sites") # red, B sites 

    # ax1.set_xlim(0,1.1)

    # Vertical line through cohesive energy
    ax1.axvline(x=cohesive_energy, color="black", linestyle="--", linewidth=1, alpha=0.8, zorder=0)

    ax1.legend()

    # Save the image
    fig.savefig(output_img_path, dpi=300)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("input_csv_path", type=Path)
    parser.add_argument("output_img_path", type=Path)

    args = parser.parse_args()

    main(args.input_csv_path, args.output_img_path)
