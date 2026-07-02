#!/bin/bash
### Script to convert the PDB file to a gaussian COM file for G Quadruplex excitonic calculation
### The script converts the PDB file into a set of coordinates, then strips all nucleobases other than Guanine from the coordinates, removes the sugar and phosphate from Guanine, then adds a methyl group into the N9 position. The coordinates of the Me H's are based on an isolated Guanine geometry optimization in PCM water.

#Check whether the filename (without extension) is provided in the input to the script
if [ $# -ne 1 ]
then
  echo "Please supply filename of the pdb file"
  exit
fi

filename="$1"

#Make header for COM file
echo "%nproc=12" > ${filename}_wl.com
echo "%chk=${filename}" >> ${filename}_wl.com
echo "#p td=nstates=5 m052x 6-31g(d) eet=(fragment=22,fragmentcavity) nosymm scrf=(pcm,solvent=water)" >> ${filename}_wl.com
echo "" >> ${filename}_wl.com
echo "Excitonic calculation for G Quadruplex based on 143d.pdb" >> ${filename}_wl.com
echo "" >> ${filename}_wl.com
echo "0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1" >> ${filename}_wl.com

#Convert PDB file into coordinates
/shared/gaussian/g16.b01/SSE4/g16/newzmat -ipdb -oxyz ${filename}.pdb ${filename}_wl_coords.xyz

#sed command to strip sugars and phosphates off Nucleobases, leaving the C attached to N9 (this may need to be modified depending on ordering of atomic coords)

sed -i '1,7d;17,23d;28,38d;48,54d;59;161d;171,177d;183,291d;301,307d;312;318d;328,334d;340,448d;458,465d' ${filename}_wl_coords.xyz
#Dihedral angle C8-N9-CMe-H (phiA) and bond angle N9-CMe-H (thetaA) calculated on isolted Adenine in PCM for the 3 H's in the methyl

phiC[1]=2.0944	#methyl
phiC[2]=-2.0944	#methyl
phiC[3]=0.0 	#methyl
phiC[4]=0.012	#N-H for adding extra proton
thetaC[1]=1.91602245
thetaC[2]=1.91602245
thetaC[3]=1.89577663
thetaC[4]=1.9825

#C-H bond length in N9-Me group, calculated on isolated Adenine in PCM

CH_blC[1]=1.09
CH_blC[2]=1.09
CH_blC[3]=1.09
CH_blC[4]=1.01


#Dihedral angle C8-N9-CMe-H (phi) and bond angle N9-CMe-H (theta) calculated on isolted guanine in PCM for the 3 H's in the methyl

phiG[1]=-1.0581931
phiG[2]=1.0581931
phiG[3]=3.14159
thetaG[1]=1.92192167
thetaG[2]=1.92195657
thetaG[3]=1.87745068

#C-H bond length in N9-Me group, calculated on isalated guanine in PCM

CH_blG[1]=1.09227
CH_blG[2]=1.09231
CH_blG[3]=1.09230

#Line number counter for the coordinate xyz file
line_num=1

#Cartesian coordinates of C8-N9-Me atoms read from coordinates file
#First atom is CMe, second is N9 then third is C8, then 17 18 19, then 33 34 35 etc

for l in {1..5} #loop over number of Guanines/fragments
do
    for k in {1..13} #loop over number of atoms in 1 Thymine in the coordinate file (15 due to capping Me not having H's on C yet)
     do
      if [ $k -eq 1 ]   #coordinates of Methyl C
      then
       cme_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       cme_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       cme_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi
      if [ $k -eq 2 ]   #coordinates of N1
      then
       n1_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       n1_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       n1_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi
      if [ $k -eq 9 ]   #coordinates of C6
      then
       c6_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       c6_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       c6_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi

      if [ $k -eq 5 ]   #coordinates of N3
      then
       n5_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       n5_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       n5_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi
      if [ $k -eq 3 ]   #coordinates of C2
      then
       c3_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       c3_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       c3_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi
      if [ $k -eq 4 ]   #coordinates of O
      then
       o4_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_wl_coords.xyz)
       o4_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_wl_coords.xyz)
       o4_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_wl_coords.xyz)
      fi
  #append (fragment) onto the atom label in the COM file for the EET gaussian calc
      awk -v r="$line_num" -v frag="(fragment=$l)" 'NR==r {print $1 frag " " $2 " " $3 " " $4}' ${filename}_wl_coords.xyz >> ${filename}_wl.com
      printf "$l,$k \n" 
  #increment line number
      line_num=$((line_num+1))
#Excecute python script to calculate the coordinates of the H atoms in the Me group with variables passed as arguments
#Append coordinates with (fragment) to the COM file
     done
      for i in {1..4}
       do
        if [ $i -eq 4 ]
        then 
          hme[$i]=$(python calc_methyl_coords.py $o4_x $o4_y $o4_z $c3_x $c3_y $c3_z $n5_x $n5_y $n5_z ${phiC[$i]} ${thetaC[$i]} ${CH_blC[$i]})
        else
         hme[$i]=$(python calc_methyl_coords.py $c6_x $c6_y $c6_z $n1_x $n1_y $n1_z $cme_x $cme_y $cme_z ${phiC[$i]} ${thetaC[$i]} ${CH_blC[$i]})
        fi
        echo "H(fragment=$l) "${hme[$i]} >> ${filename}_wl.com
       done

done
#Insert blank line at the end of geometry specification
echo "" >> ${filename}_wl.com


