from ase.io import read, write
from pathlib import Path

output_folder = Path('/home/amritagos/Git/Github/iceSurfaceWorkFlow/res/ice_surfaces_genIce')

infilepaths = list(Path('/home/amritagos/Git/Github/GenIce/').glob('*.xyz'))

for infilepath in infilepaths:
    outfilename = infilepath.name
    outfilepath = output_folder / outfilename
    
    atoms = read(infilepath, format='extxyz')
    # Switch y and z
    for atom in atoms:
        new_y = atom.position[2]
        new_z = atom.position[1]
        atom.position[1] = new_y
        atom.position[2] = new_z
    # from the CIF file
    atoms.set_cell([30.590228902373884, 26.502128397830692, 14.377611329389692])
    write(outfilepath, atoms, format='extxyz')

