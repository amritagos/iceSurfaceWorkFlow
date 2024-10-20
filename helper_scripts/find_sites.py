
import ase 
from ase.io import read, Trajectory
import numpy as np
from pathlib import Path
import argparse

# CLI
parser = argparse.ArgumentParser()
parser.add_argument("seed")
args = parser.parse_args()

seed = args.seed

# Input conditions and file paths
z_layer_min = 19.29
# script_dir = Path(__file__).resolve().parent
# infilename = script_dir / "../res/ice-scme-1.000000-new.data"
# infilename = "/home/amritagos/Git/Github/iceSurfaceWorkFlow/res/seed_500_1.xyz"
infilename = '/hpctheochem/amrita/ice_surface_project/2-relax_surfaces/traj/relaxed_surface_seed_'+str(seed)+'.traj'
output_file_a_sites = '/hpctheochem/amrita/ice_surface_project/3-find_a_b_sites/outputs/a_sites_seed_'+str(seed)+'.npy'
output_file_b_sites = '/hpctheochem/amrita/ice_surface_project/3-find_a_b_sites/outputs/b_sites_seed_'+str(seed)+'.npy'

# Symbols for elements 
h_sym = 'H'
o_sym = 'O'
# Cutoffs 
# O-O cutoff in the topmost layer (Angstrom)
o_o_top_cutoff = 5.0
# O-H intramolecular cutoff 
o_h_bond_cutoff = 1.0

# Read in stuff 
# surface_system = read(infilename, format="extxyz") # lammps-data
surface_system = Trajectory(infilename)[-1]

# Find O atoms in the top part of the topmost bilayer
top_o_ind = []
for Atom in surface_system:
    if Atom.symbol == o_sym:
        if Atom.position[2] > z_layer_min:
            top_o_ind.append(Atom.index)
len(top_o_ind)

# Add neighbour to o_o_neigh dict
def add_neighbour(neigh_dict, center_ind, neigh_ind):
    if neigh_dict.get(center_ind) is None:
        neigh_dict[center_ind] = [neigh_ind]
    else:
        neighbours = neigh_dict[center_ind]
        neighbours.append(neigh_ind)
        neigh_dict.update({center_ind:neighbours})

# Get the Hydrogen atoms bonded to the topmost O atoms and other topmost O atoms 
o_h_neigh = dict()
o_o_neigh = dict()

# H atoms bonded to O atoms 
for i_ind in top_o_ind:
    h_indices = []
    for Atom in surface_system:
        if Atom.symbol == h_sym:
            if surface_system.get_distance(i_ind, Atom.index)<o_h_bond_cutoff:
                h_indices.append(Atom.index)
        if len(h_indices)==2:
            o_h_neigh[i_ind] = h_indices
            break

# Topmost O atoms connected to topmost O atoms 
for i in range(len(top_o_ind)-1):
    i_ind = top_o_ind[i]
    for j in range(i+1, len(top_o_ind)):
        j_ind = top_o_ind[j]
        dist_ij = surface_system.get_distance(i_ind, j_ind, mic=True)
        if dist_ij < o_o_top_cutoff:
            # Add neighbours 
            add_neighbour(o_o_neigh, i_ind, j_ind)
            add_neighbour(o_o_neigh, j_ind, i_ind)

# Set of lists (in ascending order) to keep track of visited trios
visited_trios = set()
# List of lists for A sites and B sites 
a_sites = []
b_sites = []

def check_site(trio, a_sites, b_sites, o_h_neigh, atoms):
    # Number of dangling hydrogens (set to zero first)
    n_dangling_h = 0 
    # Loop through the trio O sites 
    for o_ind in trio:
        h_indices = o_h_neigh[o_ind]
        # Check the hydrogens
        # If the z coordinate of H is greater than the O, it be dangling
        for h_ind in h_indices:
            if  atoms[h_ind].position[2] > atoms[o_ind].position[2]:
                n_dangling_h += 1
    # A sites: 1 dangling hydrogen 
    if n_dangling_h==1:
        a_sites.append(trio)
    # B sites: 2 dangling hydrogens
    elif n_dangling_h==2:
        b_sites.append(trio)

for center_ind in top_o_ind[:5]:
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
                check_site(trial_trio, a_sites, b_sites, o_h_neigh, surface_system)

print(len(a_sites))
print(len(b_sites))

# Save outputs
np.save(output_file_a_sites, a_sites)
np.save(output_file_b_sites, b_sites)

# a2_sites = np.load(output_file_a_sites)
# if np.array_equal(a2_sites, a_sites):
#     print("equal")


