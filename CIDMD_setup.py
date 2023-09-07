#!/usr/bin/env python
##################################
#
#   Jesi Lee
#   CIDMD pilot study #4 
#   Jan. 2022
#
#################################
import numpy as np
import sys
import random
#import matplotlib.pyplot as plt
from itertools import product, combinations, repeat
from typing import Tuple
from lebedev import *
import argparse

def create_random_pvec(rng) -> np.ndarray:
    """
    Creates the random 3 numbers for the initial position.

    Parameters
    ----------
    rng : numpy.random._generator.Generator
        The random number generator.

    Returns
    -------
    random_pvec: numpy.ndarray
        The x, y, z position.
        The array is the size (3,). 
    """
    random_pvec = 2*(rng.random(3) - 0.5)
    return random_pvec

def create_random_dvec(rng) -> np.ndarray:
    """
    Creates the random 3 numbers for the initial direction (aka velocity).

    Parameters
    ----------
    rng : numpy.random._generator.Generator
        The random number generator.

    Returns
    -------
    random_dvec: numpy.ndarray
        The x, y, z direction.
        The array is the size (3,). 
    """
    while True:
        random_vec = 2*(rng.random(3) - 0.5)
        #random_vec[2] = np.sqrt(.5 - (random_vec[0])**2 - (random_vec[1])**2)
        mag = np.linalg.norm(random_vec)
        if mag <= 1.0:
            random_vec_norm = random_vec/mag
            return random_vec_norm

def print_stats(varname,array):
    min_val = np.min(array)
    max_val = np.max(array)
    mean_val = np.mean(array)
    print(f"{varname:16s} == Min: {min_val:9.5f} Mean: {mean_val:9.5f} Max: {max_val:9.5f}", end="\r")

def save_xyz(data, filename, symbols=None):
    if symbols is None:
        symbols = repeat("Ar")
    if type(filename) is str:
        f = open(filename, 'w')
    else:
        f = filename
    f.write(f"{len(data)}\n\n")

    for sym, (x,y,z) in zip(symbols, data):
        f.write(f"{sym:4s} {x:9.5f} {y:9.5f} {z:9.5f}\n")
    if type(filename) is str:
        f.close()

def sanity_check(mol_pos_array, new_pos, new_vel, set_distance):
    """
    Sanity check: new position should have one contact distance equal to set_distance, all others larger;
    further moving along v_hat should decrease the contact 

    Parameters
    ----------
    mol_pos_array : numpy.ndarray
        Array of atomic positions in input molecule.

    new_pos : np.ndarray
        New initial position to be added

    new_vel : np.ndarray
        New initial velocity to be added 

    set_distance : float
        Atomic radius of molecular envelope   
    """
    new_pos_dist = np.linalg.norm(mol_pos_array - new_pos[np.newaxis, :], axis=1)
    if np.abs(np.min(new_pos_dist) - set_distance) > 1e-6:
        raise RuntimeError('Minimum distance of collider position with atomic centers is not equal to set_distance')
    imin = np.argmin(new_pos_dist)
    new_pos_dist_displ = np.linalg.norm(mol_pos_array - (new_pos[np.newaxis, :] + 1e-3 * new_vel), axis=1)
    if new_pos_dist_displ[imin] > set_distance:
        raise RuntimeError('Moving collider along velocity direction appears to move away from molecular envelope')

def add_pos_vel(pos, vel, wts, new_pos, new_vel, new_wt, p_thre, v_thre, verbose, print_newline=True):
    """
    Add a pair of initial positions/velocities and associated weight to existing list,
    or if too similar to any existing position/velocity, increase the weight.

    Parameters
    ---------
    pos : list of np.ndarray
        Existing initial positions; may be modified by this function

    vel : list of np.ndarray
        Existing initial velocities; may be modified by this function

    wts : list of float
        Existing initial weights; may be modified by this function

    new_pos : np.ndarray
        New initial position to be added

    new_vel : np.ndarray
        New initial velocity to be added

    new_wt : float
        Weight associated with the new initial position/velocity

    p_thre : float
        Threshold for determining if two initial positions are the "same"

    v_thre : float
        Threshold for determining if two initial velocities are the "same"

    verbose : bool

    print_newline : bool

    Returns
    -------

    None (the variables pos, vel, wts are modified)
    """
    dupe = False
    for ipv, (p, v) in enumerate(zip(pos, vel)):
        if (np.linalg.norm(new_pos - p) < p_thre) and (np.linalg.norm(new_vel - v) < v_thre):
            if dupe: raise RuntimeError('Position/velocity vector is a duplicate of two existing ones?')
            dupe = True
            wts[ipv] += new_wt
            #print("Duplicate found, increasing weight of IC %i (%.3f + %.3f = %.3f)" % (ipv, wts[ipv]-new_wt, new_wt, wts[ipv]))
    if not dupe:
        pos.append(new_pos)
        vel.append(new_vel)
        wts.append(new_wt)
        if verbose:
            if print_newline:
                print()
            print("Adding IC %i with weight %.3f" % (len(pos), new_wt))

def gen_collider_pos_vel(mol_pos_array, p_vec, v_hat, set_distance, cutoff):
    """
    Given a molecular structure and a line representing the collider trajectory, 
    generate two pairs of initial positions and velocities, or None if it "misses".

    Parameters
    ---------
    mol_pos_array : numpy.ndarray
        Array of atomic positions in input molecule.

    p_vec : np.ndarray
        Cartesian coordinates of a point on the line representing the collider trajectory.

    v_hat : np.ndarray
        Unit vector representing the direction of the collider.

    set_distance : float
        Atomic radius for molecular envelope

    cutoff : float
        Keep initial condition only if collider passes within this distance of the molecule
    Returns
    -------
    (new_pos_fwd, new_vel_fwd, new_pos_bak, new_vel_bak) : Tuple or None
        each element of the tuple is a numpy.ndarray
    """

    # Array of displacement vectors pointing from Ar grid point to atoms in molecule
    r_1_arr = mol_pos_array - p_vec[np.newaxis, :]
    v_proj = np.dot(r_1_arr, v_hat)
    # Projection of displacement vectors along velocity direction 
    r_1v_arr = np.outer(v_proj, v_hat)
    # Displacement vector from nearest point on Ar trajectory to atomic center
    r_1p_arr = r_1_arr - r_1v_arr
    # Distances of closest approach of Ar trajectory to atomic center
    r_1p_norms = np.linalg.norm(r_1p_arr, axis=1)
    # Do not keep the initial condition if the minimum approach distance is greater than the cutoff
    if np.min(r_1p_norms) > cutoff: return None
    # The squared displacement along v_hat from the closest approach point to reach a
    # distance of set_distance from the atomic center. 
    a_squared = set_distance ** 2 - r_1p_norms**2
    # Shift p_vec by these amounts along v_hat to get distance of exactly set_distance from atomic center.
    shift_forward = []
    shift_backward = []
    for iatom, a2 in enumerate(a_squared):
        if a2 < 0: continue
        # Position vector of collider, translated such that it is now exactly set_distance from the atomic center.
        shift_forward.append(v_proj[iatom] - np.sqrt(a2))
        shift_backward.append(v_proj[iatom] + np.sqrt(a2))

    # This position and velocity corresponds to entering the molecular envelope in the "forward" direction.
    new_pos_fwd = p_vec + min(shift_forward)*v_hat
    new_vel_fwd = v_hat
    # This position and velocity corresponds to entering the molecular envelope in the "backward" direction.
    new_pos_bak = p_vec + max(shift_backward)*v_hat
    new_vel_bak = -v_hat
    return (new_pos_fwd, new_vel_fwd, new_pos_bak, new_vel_bak)

def generate_collider_ic_deterministic(mol_pos_array, lebedev_order, extent, set_distance, cutoff, cartgrid, verbose):
    """
    Generates bunch of collider positions and velocities using randomized algorithm.

    Parameters
    ---------
    mol_pos_array : numpy.ndarray
        Array of atomic positions in input molecule.

    lebedev_order : int
        Lebedev grid order. The higher the order, the more points. Only certain orders are valid (e.g. odd numbers under 31).
        The number of points for a given order are:
                   {3:6, 5:14, 7:26, 9:38, 11:50, 13:74, 15:86, 17:110, 19:146, 21:170, 23:194,
                   25:230, 27:266, 29:302, 31:350, 35:434, 41:590, 47:770, 53:974, 59:1202,
                   65:1454, 71:1730, 77:2030, 83:2354, 89:2702, 95:3074, 101:3470, 107:3890,
                   113:4334, 119:4802, 125:5294}

    extent : int
        Extent of Cartesian grid in each direction will be (-extent, +extent)

    set_distance : float
        Atomic radius for molecular envelope

    cutoff : float
        Keep initial condition only if collider passes within this distance of the molecule

    cartgrid : int
        Number of Cartesian grid points in each direction

    verbose : bool

    Returns
    -------
    (pos, vel) : Tuple
        pos : numpy.ndarray
        vel : numpy.ndarray
    """

    L = Lebedev(lebedev_order)
    # Array of distances between Lebedev grid points
    v_points = []
    for i1 in range(L.points.shape[0]):
        for i2 in range(0, i1):
            v_points.append(np.linalg.norm(L.points[i1] - L.points[i2]))

    # Threshold below which a pair of positions / directions are considered the same
    v_thre = 0.01*min(v_points)
    p_thre = 0.01*(2*extent)/(cartgrid-1)
    
    # The points on the Lebedev grid corresponding to inward-pointing velocities
    vecs = L.points
    vwts = L.weights / np.min(L.weights)
    for v_hat in L.points:
        if np.abs(np.linalg.norm(v_hat) - 1.0) > 1e-6:
            raise RuntimeError('Expected unit vectors from Lebedev grid')
    pos = []
    vel = []
    wts = []
    total_tries = len(np.linspace(-extent, extent, cartgrid))**3 * vecs.shape[0]
    print("Will try a total of %i collider position / velocity combinations" % total_tries)
    
    itrial = 0
    for grid_x in np.linspace(-extent, extent, cartgrid):
        for grid_y in np.linspace(-extent, extent, cartgrid):
            for grid_z in np.linspace(-extent, extent, cartgrid):
                p_vec = np.array([grid_x, grid_y, grid_z])
                for v_hat, v_wt in zip(vecs, vwts):
                    itrial += 1
                    print("\r Trial %i/%i" % (itrial, total_tries), end='')
                    # Generate two pairs of initial positions and velocities, or None if it "misses"
                    new_pos_vel = gen_collider_pos_vel(mol_pos_array, p_vec, v_hat, set_distance, cutoff)
                    if new_pos_vel is not None:
                        new_pos_fwd, new_vel_fwd, new_pos_bak, new_vel_bak = new_pos_vel
                        sanity_check(mol_pos_array, new_pos_fwd, new_vel_fwd, set_distance)
                        add_pos_vel(pos, vel, wts, new_pos_fwd, new_vel_fwd, v_wt, p_thre, v_thre, verbose)
                        sanity_check(mol_pos_array, new_pos_bak, new_vel_bak, set_distance)
                        add_pos_vel(pos, vel, wts, new_pos_bak, new_vel_bak, v_wt, p_thre, v_thre, verbose, print_newline=False)
    print("\nDone; generated %i position/velocity combinations" % len(pos))
    pos = np.array(pos)
    vel = np.array(vel)
    wts = np.array(wts)
    wts /= np.min(wts)
    print("After normalization, maximum weight is %.3f" % np.max(wts))
    return pos, vel, wts

def generate_collider_ic_random(mol_pos_array, n_points, extent, set_distance, cutoff, rng, verbose) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generates bunch of collider positions and velocities using randomized algorithm.

    Parameters
    ---------
    mol_pos_array : numpy.ndarray
        Array of atomic positions in input molecule.

    n_points : int
        The number of desired initial conditions for the collider atom

    extent : int
        Randomly generated collider positions will be in range (-extent, +extent)

    set_distance : float
        Atomic radius for molecular envelope

    cutoff : float
        Keep initial condition only if collider passes within this distance of the molecule

    rng : numpy.random._generator.Generator

    verbose : bool

    Returns
    -------
    (pos, vel) : Tuple
        pos : numpy.ndarray
        vel : numpy.ndarray
    """
    pos = []
    vel = []
    #wts is a dummy array of weights that will all be 1; we assume there will be no duplicates when generating random ICs.
    wts = []

    itrial = 0
    while len(pos) < n_points:
        # Generate randomized position and velocity direction
        col_pvec = create_random_pvec(rng)
        col_pvec *= extent
        col_dvec = create_random_dvec(rng)
        itrial += 1

        # Generate two pairs of initial positions and velocities, or None if it "misses"
        new_pos_vel = gen_collider_pos_vel(mol_pos_array, col_pvec, col_dvec, set_distance, cutoff)
        if new_pos_vel is not None:
            new_pos_fwd, new_vel_fwd, new_pos_bak, new_vel_bak = new_pos_vel
            sanity_check(mol_pos_array, new_pos_fwd, new_vel_fwd, set_distance)
            add_pos_vel(pos, vel, wts, new_pos_fwd, new_vel_fwd, 1.0, 0.0, 0.0, False)
            sanity_check(mol_pos_array, new_pos_bak, new_vel_bak, set_distance)
            add_pos_vel(pos, vel, wts, new_pos_bak, new_vel_bak, 1.0, 0.0, 0.0, False)
        
        # let's see the progress
        progress = len(pos)/n_points * 100. 
        print(f"{len(pos):4d}/{n_points}  {progress:5.2f}%  trials:{itrial:d}", end="\n" if (verbose and new_pos_vel is not None) else "\r")
        sys.stdout.flush()
        
    pos = np.array(pos)
    vel = np.array(vel)

    return pos, vel


def parse_user_input(*args):
    """ 
    Read user input from the command line interface. 
    """
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('xyz_in', type=str, help='REQUIRED positional argument: input .xyz structure of molecular ion\n ')
    parser.add_argument('--extent', type=float, default=50.0, 
                        help='Half the edge-length of cubic box in Angstrom in which initial positions are generated.'
                        'Should be larger than the spatial extent of the molecule.')
    parser.add_argument('--cutoff', type=float, default=2.0, 
                        help='Keep only the collider initial positions whose paths come to within this distance of any atom on the molecular ion (in Angstrom)')
    parser.add_argument('--set_distance', type=float, default=5.0, 
                        help='Atomic radius of molecular envelope surface on which the collider will be placed (in Angstrom).')
    parser.add_argument('--n_points', type=int, default=50, 
                        help='Number of initial conditions desired when using randomized algorithm.')
    parser.add_argument('--seed', type=int, default=2021, 
                        help='Number of initial conditions desired when using randomized algorithm.')
    parser.add_argument('--deterministic', action='store_true', 
                        help='If set, generate initial positions using deterministic algorithm instead of default randomized algorithm.')
    parser.add_argument('--cartgrid', type=int, default=51, 
                        help='Number of Cartesian grid points for deterministic algorithm. Default is set to generate 51 points as -50, -48, -46, ..., 48, 50.')
    parser.add_argument('--lebedev', type=int, default=5, 
                        help='Lebedev grid order for deterministic algorithm (use odd number in range [3,31]; higher orders not recommended)')
    
    
    parser.add_argument('--v_scale', type=float, default=1.) #like a vfac2 or 4

    parser.add_argument('-v', '--verbose', action='store_true', help='Turn on additional output for debugging.')
    
    return parser.parse_args(*args)

def main():
    ###=  Parse user input  =###
    args = parse_user_input(sys.argv[1:])
    assert args.cutoff <= args.set_distance, "args.cutoff should be <= args.set_distance"
    
    ###=  Reading the molecule position  =###
    mol_atom_array = np.loadtxt(args.xyz_in, skiprows = 2, usecols = (0), dtype = str)
    mol_pos_array = np.loadtxt(args.xyz_in, skiprows = 2, usecols = (1,2,3), dtype = float) 

    if args.deterministic:
        ###=  Deterministic algorithm; scan over Cartesian grid (position) and Lebedev grid (direction)  =###
        pos, vel, wts = generate_collider_ic_deterministic(mol_pos_array, args.lebedev, args.extent, args.set_distance, args.cutoff, args.cartgrid, args.verbose)
        np.savetxt("ar_points_weights.txt", wts)
    else:
        ###=  Generate ICs with random positions and directions  =###
        rng = np.random.default_rng(args.seed)
        pos, vel = generate_collider_ic_random(mol_pos_array, args.n_points, args.extent, args.set_distance, args.cutoff, rng, args.verbose)

    ###=  Save output to files  =###
    with open("colpts_mol_with_vel.xyz", 'w') as fid:
        for p,v in zip(pos, vel):
            pv = [p,p+v]
            save_xyz(pv, fid, repeat("Ar"))
    save_xyz(pos, "colpts_mol_projection.xyz", repeat("Ar"))
    save_xyz(pos, 'col_projection_pos.xyz', repeat("Ar"))
    save_xyz(vel*args.v_scale, "col_projection_vel.xyz", repeat("Ar"))

if __name__ == "__main__":
    main()
