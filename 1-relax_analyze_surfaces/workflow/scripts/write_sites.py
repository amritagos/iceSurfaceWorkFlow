from pathlib import Path
import sys
import json
from typing import Dict, List, Tuple
from ase import Atoms
import pandas as pd
import re
from ase.io import read, write
import numpy as np


def find_sites(
    o_o_neigh: Dict, o_h_neigh: Dict, top_o_indices: List, atoms: Atoms
) -> Tuple[List[List[int]], List[List[int]], List[List[int]], List[List[int]]]:
    """Find all the sites in the surface and classify them, returning A, B , C sites and the number of D sites

    Args:
        o_o_neigh (Dict): Neighbour list such that the Key -> oxygen index. Values -> List of neighbouring oxygen indices
        o_h_neigh (Dict): Neighbour list for intramolecular O-H bonds such that the Key -> oxygen index. Values -> List of hydrogen indices in the same molecule
        top_o_indices (List): Indices of topmost O atoms. Indices wrt original Atoms object for the surface
        atoms(Atoms): Atoms object for the entire surface

    Returns:
        Tuple[List[List[int]], List[List[int]], List[List[int]], List[List[int]]]: A sites, B sites, C sites, D sites. Each site is defined by a trio of O atom indices.
    """
    # Set of lists (in ascending order) to keep track of visited trios
    visited_trios = set()
    # List of lists for A sites, B sites and C sites
    a_sites = []
    b_sites = []
    c_sites = []
    d_sites = []

    for center_ind in top_o_indices:
        # Get neighbours of the center
        neigh_center = o_o_neigh[center_ind]
        # Now loop through neighbours
        for i_ind in neigh_center:
            # Find common neighbours between the center and i
            neigh_i = o_o_neigh[i_ind]
            common_ele = list(set(neigh_center).intersection(neigh_i))
            # Loop through the common elements to get the trios
            for j_ind in common_ele:
                trial_trio = [center_ind, i_ind, j_ind]
                trial_trio.sort()
                # If you've processed the trio already, skip it
                if frozenset(trial_trio) in visited_trios:
                    continue
                # Check the trio
                else:
                    visited_trios.add(frozenset(trial_trio))
                    check_site(
                        trial_trio, a_sites, b_sites, c_sites, d_sites, o_h_neigh, atoms
                    )

    # Return the A sites, B sites, C sites and the number of D sites
    return (a_sites, b_sites, c_sites, d_sites)


def add_neighbour(neigh_dict, center_ind, neigh_ind):
    """Adds O neighbour to O center in the oxygen-oxygen neighbour list (a dictionary)"""
    if neigh_dict.get(center_ind) is None:
        neigh_dict[center_ind] = [neigh_ind]
    else:
        neighbours = neigh_dict[center_ind]
        neighbours.append(neigh_ind)
        neigh_dict.update({center_ind: neighbours})


def check_site(
    trio: List,
    a_sites: List[List],
    b_sites: List[List],
    c_sites: List[List],
    d_sites: List[List],
    o_h_neigh: Dict,
    atoms: Atoms,
):
    """For a particular site (defined by a trio), check what kind of site it is

    Args:
        trio (List): Describes a particular site. Indices of a trio of O atoms in the topmost layer
        a_sites (List[List]): List of lists containing all trios (sites) that are A sites
        b_sites (List[List]): List of lists containing all trios (sites) that are B sites
        c_sites (List[List]): Contains all C sites
        d_sites (List[List]): Contains all D sites
        o_h_neigh (Dict): Neighbour list for hydrogen atoms bound to O atoms (intramolecular bonds)
        atoms (Atoms): Atoms object for the entire surface
    """
    # Number of dangling hydrogens (set to zero first)
    n_dangling_h = 0
    # Loop through the trio O sites
    for o_ind in trio:
        h_indices = o_h_neigh[o_ind]
        # Check the hydrogens
        # If the z coordinate of H is greater than the O, it be dangling
        for h_ind in h_indices:
            if atoms[h_ind].position[2] > atoms[o_ind].position[2]:
                n_dangling_h += 1
    # A sites: 1 dangling hydrogen
    if n_dangling_h == 1:
        a_sites.append(trio)
    # B sites: 2 dangling hydrogens
    elif n_dangling_h == 2:
        b_sites.append(trio)
    # C sites: no dangling hydrogens
    elif n_dangling_h == 0:
        c_sites.append(trio)
    # D sites: all three are dangling hydrogens
    elif n_dangling_h == 3:
        d_sites.append(trio)


def create_o_h_neighbour_list(
    top_o_indices: List, h_symbol: str, o_h_bond_cutoff: float, atoms: Atoms
) -> Dict:
    """Creates a dictionary containing the H atoms bonded to O atoms (only in the top most layer of the top bilayer).

    Args:
        top_o_indices (List): Indices of topmost O atoms. Indices wrt original Atoms object
        h_symbol (str): Symbol for Hydrogen
        o_h_bond_cutoff (float): Cutoff for neighbours
        atoms (Atoms): Original Atoms object for the entire surface

    Returns:
        Dict: Keys : oxygen indices, and values are a list of hydrogen indices
    """
    o_h_neigh = dict()  # Output neighbour list dictionary

    for i_ind in top_o_indices:
        h_indices = []
        for Atom in atoms:
            if Atom.symbol == h_symbol:
                if atoms.get_distance(i_ind, Atom.index, mic=True) < o_h_bond_cutoff:
                    h_indices.append(Atom.index)
            if len(h_indices) == 2:
                o_h_neigh[i_ind] = h_indices
                break

    return o_h_neigh


def create_o_o_neighbour_list(
    top_o_indices: List, o_o_top_cutoff: float, atoms: Atoms
) -> Dict:
    """Creates a neighbour list for the oxygens in the topmost layer (only O-O connections).

    Args:
        top_o_indices (List): Indices of the topmost bilayer
        o_o_top_cutoff (float): Distance cutoff for the topmost bilayer
        atoms (Atoms): Atoms object for the entire surface

    Returns:
        Dict: Keys -> oxygen index. Values -> List of neighbouring oxygen indices
    """
    o_o_neigh = dict()

    # Topmost O atoms connected to topmost O atoms
    for i in range(len(top_o_indices) - 1):
        i_ind = top_o_indices[i]
        for j in range(i + 1, len(top_o_indices)):
            j_ind = top_o_indices[j]
            dist_ij = atoms.get_distance(i_ind, j_ind, mic=True)
            if dist_ij < o_o_top_cutoff:
                # Add neighbours
                add_neighbour(o_o_neigh, i_ind, j_ind)
                add_neighbour(o_o_neigh, j_ind, i_ind)

    return o_o_neigh


def find_topmost_O_layer(atoms: Atoms) -> List[int]:
    """Find the topmost O layer atom indices. These have been processed already and have a tag of 1.

    Args:
        atoms (Atoms): Atoms object

    Returns:
        List[int]: list of oxygen atom indices
    """
    top_o_indices = [Atom.index for Atom in atoms if Atom.tag == 1]
    return top_o_indices


def main(
    input_xyz: Path,
    metadata: Path,
    h_symbol: str,
    o_o_cutoff: float,
    o_h_cutoff: float,
    a_site_path: Path,
    b_site_path: Path,
    c_site_path: Path,
):

    # Read in the XYZ file and analyze
    surface = read(input_xyz, format="extxyz")
    # The topmost O atoms have a tag of 1
    top_o_indices = find_topmost_O_layer(surface)

    # Get the neighbour lists for O-O connections in the top layer, and O-H connections
    o_h_neigh = create_o_h_neighbour_list(top_o_indices, h_symbol, o_h_cutoff, surface)
    o_o_neigh = create_o_o_neighbour_list(top_o_indices, o_o_cutoff, surface)

    # Find all the sites and get the A, B, C sites and the D sites
    a_sites, b_sites, c_sites, d_sites = find_sites(
        o_o_neigh, o_h_neigh, top_o_indices, surface
    )

    # Save the outputs
    np.save(a_site_path, a_sites)
    np.save(b_site_path, b_sites)
    np.save(c_site_path, c_sites)

    # Site info inside metadata
    res = dict()
    res["num_a_sites"] = len(a_sites)
    res["num_b_sites"] = len(b_sites)
    res["num_c_sites"] = len(c_sites)
    res["num_d_sites"] = len(d_sites)

    # Now overwrite the metadata
    with open(metadata, "w") as f:
        f.write(json.dumps(res, indent=4))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--input_xyz", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--h_symbol", type=str)
    parser.add_argument("--o_o_cutoff", type=float)
    parser.add_argument("--o_h_cutoff", type=float)
    parser.add_argument("--a_site_path", type=Path)
    parser.add_argument("--b_site_path", type=Path)
    parser.add_argument("--c_site_path", type=Path)

    args = parser.parse_args()

    main(
        args.input_xyz,
        args.metadata,
        args.h_symbol,
        args.o_o_cutoff,
        args.o_h_cutoff,
        args.a_site_path,
        args.b_site_path,
        args.c_site_path,
    )
