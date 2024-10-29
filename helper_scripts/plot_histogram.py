from pathlib import Path
from spirit_extras.plotting import Paper_Plot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main(
    input_csv_path: Path,
    output_img_path: Path,
):  
    # Input values
    cohesive_energy = -0.618898383010492 #  in eV 
    cohesive_energy_stdev = 0.00230203355051586
    factor = -1.0 # to get positive values 
    n_bins = 75
    # Positive value for cohesive energy 
    cohesive_energy = factor*cohesive_energy

    # Obtain the data
    df = pd.read_csv(input_csv_path)

    df_a_sites = df[df['site_type'] == 'a']
    df_b_sites = df[df['site_type'] == 'b']
    df_c_sites = df[df['site_type'] == 'c']
    # All relative energies in eV
    energy_a_sites = factor*df_a_sites['binding_energy'].values
    energy_b_sites = factor*df_b_sites['binding_energy'].values
    energy_c_sites = factor*df_c_sites['binding_energy'].values
    
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
    ax1.hist(energy_c_sites, bins=n_bins, ls='dotted', lw=3, fc="tab:grey", label="B sites") # grey, C sites 

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
