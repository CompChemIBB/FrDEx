####Python script to replace the geometry of a molecule by a ground state optimised structure
###Minimises the rmsd to achieve this using the python rmsd routine
###Takes as input the xyz file of the geometry you want to replace, and the geometry you want to replace it with

#Usage: python replace_molfrag_molopt.py molfrag.xyz molopt.xyz
#Both files must have the same number of atoms
#If the atoms are in a different order in each file, use the option -reorder to reorder the atoms via the brute force approach, as hungarian algorithm doesnt seem to work
#i.e.  python replace_molfrag_molopt.py -reorder molfrag.xyz molopt.xyz
#As output it will produce molsub.xyz

import sys
import numpy as np
import os
import argparse as arg

#add path to rmsd script
sys.path.append('/home/roberto/JAMES/Programs/rmsd-master/rmsd/')
#then import rmsd
#Doesnt work - need additional scripts
import calculate_rmsd


#Parse input arguments
parser=arg.ArgumentParser(description="Tool to replace local geometry of a molecule with geometry optimised structure. Both files must be provided in xyz format.")

parser.add_argument('-reorder',action='store_true',help='If the atoms in the xyz files are in different orders, use this option to reorder according to the brute force algorithm')

parser.add_argument('molfrag',help='The xyz file of the geometry you want to replace with the geometry optimised structure')

parser.add_argument('molopt',help='The xyz file of the geometry optimised structure')

args = parser.parse_args()

#Get the coordinates and atom labels from the files given as input arguments
molfrag_atoms,molfrag_xyz = rmsd.get_coordinates(args.molfrag,'xyz')
molopt_atoms, molopt_xyz = rmsd.get_coordinates(args.molopt,'xyz')

#check the files have the same number of atoms
if not (len(molfrag_atoms) == len(molopt_atoms)):
    print("error: different numbers of atoms in the two xyz files")
    quit()


#save number of atoms to a variable
tot_atoms = len(molopt_atoms)

#find centroids for molfrag and molopt
molfrag_centroid = rmsd.centroid(molfrag_xyz)
molopt_centroid = rmsd.centroid(molopt_xyz)

#move both to origin
molfrag_xyz -= molfrag_centroid
molopt_xyz -= molopt_centroid

#Reorder if necessary
if args.reorder == True :
    reordered_atoms = rmsd.reorder_brute(molfrag_atoms, molopt_atoms, molfrag_xyz, molopt_xyz)
    print(molopt_xyz)
    molopt_xyz = molopt_xyz[reordered_atoms]
    print(molopt_xyz)
    molopt_atoms = molopt_atoms[reordered_atoms]

#check that reorder has worked, and that it has been called if necessary
if not all(molopt_atoms == molfrag_atoms):
    print("error: Structure not aligned, use the -reorder option")
    quit()

#find rotation matrix for molopt to molfrag using the quaternion method and print the RMSD
rotmat = rmsd.quaternion_rotate(molopt_xyz,molfrag_xyz)
print("Rotation completed, RMSD: %s" % rmsd.quaternion_rmsd(molopt_xyz,molfrag_xyz))

#rotate molopt at the origin and add molfrag centroid to it
molopt_rot = np.dot(molopt_xyz,rotmat)
molopt_rot += molfrag_centroid

#write to file
output_array = np.hstack((np.reshape(molopt_atoms,(tot_atoms,1)),np.around(molopt_rot,decimals=8)))
output_file = open('molsub.xyz','w')
output_file.write('%s\n\n' % tot_atoms)
np.savetxt(output_file,output_array,fmt='%s     ', newline="\n")

output_file.close()
