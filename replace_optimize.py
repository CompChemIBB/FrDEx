#!/usr/bin/env python3
"""
Python script to replace all the geometries of fragments (Guanine/Cytosine) 
in a Quadruplex/Duplex by single ground state optimized structures 
in .com format for fragment diabatization.

Minimizes the RMSD via quaternion alignment to achieve this.
Requires files frag1.xyz ... fragN.xyz in the folder for each fragment.
Also needs the reference geometry optimized structure in .xyz format.
"""

import os
import sys
import argparse as arg
import numpy as np
import rmsd

# 1. Parse command-line arguments
parser = arg.ArgumentParser(description='Replace Guanine/Cytosine fragment geometries using an optimized reference structure.')
parser.add_argument('-i', '--input', default='cy.xyz', help='Input optimized reference geometry filename (default: cy.xyz)')
parser.add_argument('-f', '--format', default='xyz', help='Format of the input reference file (default: xyz)')
parser.add_argument('-n', '--num_frags', type=int, default=9, help='Total number of fragments to process (default: 9)')
parser.add_argument('-s', '--system-name', default='G4-1qdf', help='System name for output files (default: G4-1qdf)')
args = parser.parse_args()

# Get the coordinates of the atoms from the reference optimized file
meG_eq_atoms, meG_eq = rmsd.get_coordinates(args.input, args.format)

# Find its centroid and move reference structure to the origin
meG_eq_centroid = rmsd.centroid(meG_eq)
meG_eq -= meG_eq_centroid

# Get system name from command 
SYSTEM_NAME = args.system_name

# Open the job submission files and centroids file
single_frag_job_file = open(f'{SYSTEM_NAME}_single_frag_sub.sh', 'w')
single_frag_job_file.write("#!/bin/tcsh -f\n")

pair_frag_job_file = open(f'{SYSTEM_NAME}_pair_frag_sub.sh', 'w')
pair_frag_job_file.write("#!/bin/tcsh -f\n")

centroid_file = open(f'{SYSTEM_NAME}_mon_sub_centroids.out', 'w')
centroid_file.write("Frag Centroid(x) Centroid(y) Centroid(z)\n")

# Initialize list structures sized safely dynamically based on user parameter
total_loops = args.num_frags
frag_atoms = [None] * total_loops
frag = [None] * total_loops
frag_centroid = [None] * total_loops
meG_eq_rot = [None] * total_loops

# 2. Loop over fragments to align and replace geometries
for idx in range(total_loops):
    fragnum1 = idx + 1
    frag_file = f'frag{fragnum1}.xyz'
    
    if not os.path.exists(frag_file):
        print(f"Warning: Expected fragment file '{frag_file}' not found. Skipping.")
        continue

    print(f'Processing Fragment Number: {fragnum1}')
    frag_atoms[idx], frag[idx] = rmsd.get_coordinates(frag_file, 'xyz')

    # Find centroids for the extracted fragment snapshot
    frag_centroid[idx] = rmsd.centroid(frag[idx])

    # Move fragment to origin for optimal rotational matching
    frag[idx] -= frag_centroid[idx]

    # Reorder coordinates using brute force method if mapping doesn't match perfectly
    reordered_atoms = rmsd.reorder_brute(meG_eq_atoms, frag_atoms[idx], meG_eq, frag[idx])
    frag[idx] = frag[idx][reordered_atoms]
    frag_atoms[idx] = frag_atoms[idx][reordered_atoms]

    # Check alignment integrity
    if not all(meG_eq_atoms == frag_atoms[idx]):
        print(f"Error: Structure mapping not aligned for Fragment {fragnum1}. Use the -reorder option.")
        sys.exit(1)
    print("...structures are aligned matching elements perfectly.")  
    
    # Reshaping safely dynamically with -1 instead of hardcoded 19
    test_array = np.hstack((np.reshape(frag_atoms[idx], (-1, 1)), frag[idx]))
    
    # Find rotation matrix for reference optimized geometry to fragment framework via quaternion mapping
    rotmat = rmsd.quaternion_rotate(meG_eq, frag[idx])
    print(f"Rotation completed, RMSD: {rmsd.quaternion_rmsd(meG_eq, frag[idx]):.6f}")

    # Rotate reference optimized geometry at origin and move to physical coordinate space center
    meG_eq_rot[idx] = np.dot(meG_eq, rotmat)
    meG_eq_rot[idx] += frag_centroid[idx]


# 3. Generating Single Monomer Gaussian inputs
for idx in range(total_loops):
    fragnum1 = idx + 1
    if meG_eq_rot[idx] is None:
        continue
        
    # Formatting output geometry with dynamic shape scaling
    output_array = np.hstack((np.reshape(frag_atoms[idx], (-1, 1)), meG_eq_rot[idx]))
    com_filename = f'{SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab.com'
    
    with open(com_filename, 'w') as output_file:
        output_file.write("%nproc=12\n")
        output_file.write(f"%%chk={SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab.chk \n")
        output_file.write("#p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(3/33=4) td(nstates=10, conver=6) iop(9/40=5) nosymm \n\n")
        output_file.write(f"Fragment {fragnum1} diabatization snapshot\n\n")
        output_file.write("0 1\n")
        np.savetxt(output_file, output_array, fmt='%s  ', newline="\n")
        output_file.write("\n")

    # Append processing paths to cluster run script
    single_frag_job_file.write(f"/home/roberto/IMOT/EFALBO/FrDEX/go-g16b-SSE4 {SYSTEM_NAME}_frag{fragnum1}-mon_sub-diab\n")

    # Write centroid coordinate traces
    centroid_file.write(f"{fragnum1} {frag_centroid[idx][0]:.6f} {frag_centroid[idx][1]:.6f} {frag_centroid[idx][2]:.6f}\n")

# Close files before pair looping handles output arrays
centroid_file.close()
single_frag_job_file.close()


# 4. Generating Multi-Fragment Block Pair Combinations
print("Building fragment pairing matrices...")
for idx1 in range(total_loops - 1):
    fragnum1 = idx1 + 1
    if meG_eq_rot[idx1] is None: 
        continue
        
    for idx2 in range(idx1 + 1, total_loops):
        fragnum2 = idx2 + 1
        if meG_eq_rot[idx2] is None: 
            continue
            
        output_array1 = np.hstack((np.reshape(frag_atoms[idx1], (-1, 1)), meG_eq_rot[idx1]))
        output_array2 = np.hstack((np.reshape(frag_atoms[idx2], (-1, 1)), meG_eq_rot[idx2]))
        
        pair_filename = f'{SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab.com'
        
        with open(pair_filename, 'w') as output_file:
            output_file.write("%nproc=12\n")
            output_file.write(f"%%chk={SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab.chk \n")
            output_file.write("#p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(9/9=3) iop(3/33=4) td(nstates=40) iop(9/40=5) nosymm\n\n")
            output_file.write(f"Fragments {fragnum1} and {fragnum2} tracking interactions\n\n")
            output_file.write("0 1\n")
            np.savetxt(output_file, np.vstack((output_array1, output_array2)), fmt='%s  ', newline="\n")
            output_file.write("\n")

        # Record generation string directly to pair executor block
        pair_frag_job_file.write(f"/home/roberto/IMOT/EFALBO/FrDEX/go-g16b-SSE4 {SYSTEM_NAME}_frag{fragnum1}_frag{fragnum2}-mon_sub-diab\n")

pair_frag_job_file.close()
print("All operational scripts configured successfully.")