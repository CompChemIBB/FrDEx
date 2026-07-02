#!/usr/bin/env python3
"""
Python script to replace all the geometries of fragments (Guanine/Cytosine) 
in a Quadruplex/Duplex by single ground state optimized structures 
in .com format for fragment diabatization.

Performs alignment using HEAVY ATOMS (C, N, O) ONLY.
Dynamically handles 16-atom (standard) and 17-atom (hemi-protonated) fragments
by matching them to their respective optimized reference structures provided via command line.
"""

import os
import sys
import argparse as arg
import numpy as np
import rmsd

# 1. Parse command-line arguments
parser = arg.ArgumentParser(description='Replace fragment geometries using heavy-atom alignment to 16-atom and 17-atom references.')
parser.add_argument('-i16', '--input16', required=True, help='Input 16-atom optimized reference geometry (standard cytosine)')
parser.add_argument('-i17', '--input17', required=True, help='Input 17-atom optimized reference geometry (hemi-protonated cytosine)')
parser.add_argument('-f', '--format', default='xyz', help='Format of the input reference files (default: xyz)')
parser.add_argument('-n', '--num_frags', type=int, default=6, help='Total number of fragments to process (default: 6)')
parser.add_argument('-s', '--system-name', default='test', help='System name for output files (default: test)')
args = parser.parse_args()

# Load 16-atom reference
ref16_atoms, ref16_coords = rmsd.get_coordinates(args.input16, args.format)
# Load 17-atom reference
ref17_atoms, ref17_coords = rmsd.get_coordinates(args.input17, args.format)

# Helper function to prepare reference objects for alignment
def prepare_reference(atoms, coords):
    heavy_idx = [i for i, atom in enumerate(atoms) if atom.strip().upper() in ['C', 'N', 'O']]
    heavy_atoms = atoms[heavy_idx]
    centroid = rmsd.centroid(coords[heavy_idx])
    centered_coords = coords - centroid
    return heavy_idx, heavy_atoms, centered_coords

# Precompute centered reference heavy-atom criteria
ref16_heavy_idx, ref16_heavy_atoms, ref16_coords_centered = prepare_reference(ref16_atoms, ref16_coords)
ref17_heavy_idx, ref17_heavy_atoms, ref17_coords_centered = prepare_reference(ref17_atoms, ref17_coords)

SYSTEM_NAME = args.system_name

# Open the job submission files and centroids file
single_frag_job_file = open(f'{SYSTEM_NAME}_single_frag_sub.sh', 'w')
single_frag_job_file.write("#!/bin/tcsh -f\n")

pair_frag_job_file = open(f'{SYSTEM_NAME}_pair_frag_sub.sh', 'w')
pair_frag_job_file.write("#!/bin/tcsh -f\n")

centroid_file = open(f'{SYSTEM_NAME}_mon_sub_centroids.out', 'w')
centroid_file.write("Frag Centroid(x) Centroid(y) Centroid(z)\n")

total_loops = args.num_frags
final_frag_atoms = [None] * total_loops
final_frag_coords = [None] * total_loops

# 2. Loop over fragments to align and replace geometries
for idx in range(total_loops):
    fragnum1 = idx + 1
    frag_file = f'frag{fragnum1}.xyz'
    
    if not os.path.exists(frag_file):
        print(f"Warning: Expected fragment file '{frag_file}' not found. Skipping.")
        continue

    print(f'Processing Fragment Number: {fragnum1}')
    f_atoms, f_coords = rmsd.get_coordinates(frag_file, 'xyz')
    num_atoms_target = len(f_atoms)

    # Automatically select the correct reference base based on fragment atom count
    if num_atoms_target == 16:
        ref_atoms = ref16_atoms
        ref_coords_centered = ref16_coords_centered
        ref_heavy_idx = ref16_heavy_idx
        ref_heavy_atoms = ref16_heavy_atoms
        print(" -> Detected: Standard Cytosine framework (16 atoms). Aligning to 16-atom reference.")
    elif num_atoms_target == 17:
        ref_atoms = ref17_atoms
        ref_coords_centered = ref17_coords_centered
        ref_heavy_idx = ref17_heavy_idx
        ref_heavy_atoms = ref17_heavy_atoms
        print(" -> Detected: Hemi-protonated Cytosine framework (17 atoms). Aligning to 17-atom reference.")
    else:
        print(f"Error: Fragment {fragnum1} has an unsupported number of atoms ({num_atoms_target}). Expected 16 or 17.")
        sys.exit(1)

    # Extract target fragment heavy atom indices
    frag_heavy_indices = [i for i, atom in enumerate(f_atoms) if atom.strip().upper() in ['C', 'N', 'O']]
    frag_heavy_atoms = f_atoms[frag_heavy_indices]
    
    # Isolate core heavy atom frameworks for rotation mapping
    ref_heavy_coords_centered = ref_coords_centered[ref_heavy_idx]
    frag_heavy_coords = f_coords[frag_heavy_indices]
    frag_heavy_centroid = rmsd.centroid(frag_heavy_coords)
    frag_heavy_coords_centered = frag_heavy_coords - frag_heavy_centroid

    # Align the core frameworks via brute force permutation matching
    reordered_heavy_idx = rmsd.reorder_brute(
        ref_heavy_atoms, frag_heavy_atoms, 
        ref_heavy_coords_centered, frag_heavy_coords_centered
    )
    frag_heavy_coords_centered = frag_heavy_coords_centered[reordered_heavy_idx]
    frag_heavy_atoms = frag_heavy_atoms[reordered_heavy_idx]

    if not all(ref_heavy_atoms == frag_heavy_atoms):
        print(f"Error: Heavy atom mapping could not be aligned for Fragment {fragnum1}.")
        sys.exit(1)
        
    # Generate the rotation matrix
    rotmat = rmsd.quaternion_rotate(ref_heavy_coords_centered, frag_heavy_coords_centered)
    print(f"...heavy atoms matched. Quaternion RMSD: {rmsd.quaternion_rmsd(ref_heavy_coords_centered, frag_heavy_coords_centered):.6f}")

    # Rotate the matching full reference geometry and reposition it to the snapshot space
    final_frag_atoms[idx] = ref_atoms
    final_frag_coords[idx] = np.dot(ref_coords_centered, rotmat) + frag_heavy_centroid


# 3. Generating Single Monomer Gaussian inputs
for idx in range(total_loops):
    fragnum1 = idx + 1
    if final_frag_coords[idx] is None:
        continue
        
    output_array = np.hstack((np.reshape(final_frag_atoms[idx], (-1, 1)), final_frag_coords[idx]))
    com_filename = f'{SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab.com'
    
    with open(com_filename, 'w') as output_file:
        output_file.write("%nproc=16\n")
        output_file.write(f"%chk={SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab.chk \n")
        output_file.write("#p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(3/33=4) td(nstates=10, conver=6) iop(9/40=5) nosymm \n\n")
        output_file.write(f"Fragment {fragnum1} diabatization snapshot\n\n")
        output_file.write("0 1\n")
        np.savetxt(output_file, output_array, fmt='%s  ', newline="\n")
        output_file.write("\n")

    single_frag_job_file.write(f"/home/roberto/IMOT/EFALBO/FrDEX/go-g16b-SSE4 {SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab\n")
    centroid_file.write(f"{fragnum1} {final_frag_coords[idx][0][0]:.6f} {final_frag_coords[idx][0][1]:.6f} {final_frag_coords[idx][0][2]:.6f}\n")

centroid_file.close()
single_frag_job_file.close()


# 4. Generating Multi-Fragment Block Pair Combinations
print("Building fragment pairing matrices...")
for idx1 in range(total_loops - 1):
    fragnum1 = idx1 + 1
    if final_frag_coords[idx1] is None: 
        continue
        
    for idx2 in range(idx1 + 1, total_loops):
        fragnum2 = idx2 + 1
        if final_frag_coords[idx2] is None: 
            continue
            
        output_array1 = np.hstack((np.reshape(final_frag_atoms[idx1], (-1, 1)), final_frag_coords[idx1]))
        output_array2 = np.hstack((np.reshape(final_frag_atoms[idx2], (-1, 1)), final_frag_coords[idx2]))
        
        pair_filename = f'{SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab.com'
        
        with open(pair_filename, 'w') as output_file:
            output_file.write("%nproc=16\n")
            output_file.write(f"%chk={SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab.chk \n")
            output_file.write("#p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(9/9=3) iop(3/33=4) td(nstates=40) iop(9/40=5) nosymm\n\n")
            output_file.write(f"Fragments {fragnum1} and {fragnum2} tracking interactions\n\n")
            output_file.write("0 1\n")
            np.savetxt(output_file, np.vstack((output_array1, output_array2)), fmt='%s  ', newline="\n")
            output_file.write("\n")

        pair_frag_job_file.write(f"/home/roberto/IMOT/EFALBO/FrDEX/go-g16b-SSE4 {SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab\n")

pair_frag_job_file.close()
print("All operational scripts configured successfully.")