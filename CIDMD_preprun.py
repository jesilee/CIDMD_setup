#!/usr/bin/env python

import numpy as np
import sys
import os 
import argparse
import glob




def to_xyz(xyz, fname):
    with open(fname, 'w') as f:
        f.write(f"{len(xyz)}\ncomment\n")
        for atom in xyz:
            sym, x, y, z = atom
            x, y, z = list(map(lambda a: float(a), [x,y,z]))
            f.write(f"{sym:4s} {x:8.5f} {y:8.5f} {z:8.5f}\n")

def main():

    mol_pos_xyz = sys.argv[1]
    col_pos_xyz = 'col_projection_pos.xyz'                                                             
    col_vel_xyz = 'col_projection_vel.xyz'                                                             
    
    mol_pos_array = np.loadtxt(mol_pos_xyz, skiprows=2, usecols=(1,2,3), dtype=float)
    mol_atom_array = np.loadtxt(mol_pos_xyz, skiprows=2, usecols=(0), dtype=str)                                        
    col_atom_array = np.loadtxt(col_pos_xyz, skiprows=2, usecols=(0), dtype=str)
    col_pos_array = np.loadtxt(col_pos_xyz,skiprows = 2, usecols=(1,2,3), dtype=float)
    col_vel_array = np.loadtxt(col_vel_xyz,skiprows = 2, usecols=(1,2,3), dtype=float)                                     
    for i, (col_pos, col_atom, col_vel) in enumerate(zip(col_pos_array, col_atom_array, col_vel_array)):
    
        start_xyz = np.vstack((mol_pos_array, col_pos))
        start_atom = np.hstack((mol_atom_array, [col_atom]))     
        start = np.hstack((np.atleast_2d(start_atom).T, start_xyz))
    
        dirname = "calcs/{:03d}/chunk_0000/".format(i)
        os.makedirs(dirname, exist_ok = True)
        xyz_fname = os.path.join(dirname, "start.xyz")
        to_xyz(start, xyz_fname)
    
        vel_xyz = np.vstack(( np.zeros_like(mol_pos_array), col_vel))
        vel = np.hstack((np.atleast_2d(start_atom).T, vel_xyz))
        vel_fname = os.path.join(dirname, "vel.xyz")
        to_xyz(vel, vel_fname)
        print("saved to", xyz_fname, vel_fname)
    



if __name__=="__main__":
    main()












    



