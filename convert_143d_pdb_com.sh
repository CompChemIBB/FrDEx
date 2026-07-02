#!/bin/bash
### Script to convert the 143d PDB file to a gaussian COM file for G Quadruplex excitonic calculation
### The script converts the PDB file into a set of coordinates, then strips all nucleobases other than Guanine from the coordinates, removes the sugar and phosphate from Guanine, then adds a methyl group into the N9 position. The coordinates of the Me H's are based on an isolated Guanine geometry optimization in PCM water.

#Check whether the filename (without extension) is provided in the input to the script
if [ $# -ne 1 ]
then
  echo "Please supply filename of the pdb file"
  exit
fi

filename="$1"

#Make header for COM file
echo "%nproc=12" > ${filename}.com
echo "%chk=${filename}" >> ${filename}.com
echo "#p td=nstates=10 m052x/6-31g(d) nosymm eet=(fragment=12,fragmentcavity) scrf=(pcm,solvent=water)" >> ${filename}.com
echo "" >> ${filename}.com
echo "Excitonic calculation for G Quadruplex based on 143d.pdb" >> ${filename}.com
echo "" >> ${filename}.com
echo "0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1" >> ${filename}.com

#Convert PDB file into coordinates
/shared/gaussian/g16.b01/SSE4/g16/newzmat -ipdb -oxyz ${filename}.pdb ${filename}_coords.xyz

#sed command to strip other nucleobases and sugars and phosphates off Guanine, leaving the C attached to N9 (this may need to be modified depending on ordering of atomic coords)
sed -i '1,7d;17,23d;28,38d;48,54d;59,161d;171,177d;183,291d;301,307d;312,318d;328,334d;340,448d;458,465d' ${filename}_coords.xyz


#Dihedral angle C8-N9-CMe-H (phi) and bond angle N9-CMe-H (theta) calculated on isolted guanine in PCM for the 3 H's in the methyl

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


#Line number counter for the coordinate xyz file
line_num=1

#Cartesian coordinates of C8-N9-Me atoms read from coordinates file
#First atom is CMe, second is N9 then third is C8, then 17 18 19, then 33 34 35 etc

for l in {1..6} #loop over number of Guanines/fragments
do
  if [[ l -eq 3 || l -eq 5 ]]
  then
   for k in {1..14} #loop over number of atoms in 1 Guanine in the coordinate file (16 due to capping Me not having H's on C yet)
   do
     if [ $k -eq 1 ]   #coordinates of Methyl C
     then
       cme_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       cme_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       cme_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi
     if [ $k -eq 2 ]   #coordinates of N9
     then
       n9_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       n9_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       n9_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi
     if [ $k -eq 3 ]   #coordinates of C8
     then
       c8_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       c8_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       c8_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi

  #append (fragment) onto the atom label in the COM file for the EET gaussian calc
    awk -v r="$line_num" -v frag="(fragment=$l)" 'NR==r {print $1 frag " " $2 " " $3 " " $4}' ${filename}_coords.xyz >> ${filename}.com

  #increment line number
    line_num=$((line_num+1))
   done
  else 
   for k in {1..13} #loop over number of atoms in 1 Guanine in the coordinate file (16 due to capping Me not having H's on C yet)
   do
     if [ $k -eq 1 ]   #coordinates of Methyl C
     then
       cme_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       cme_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       cme_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi
     if [ $k -eq 2 ]   #coordinates of N9
     then
       n9_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       n9_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       n9_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi
     if [ $k -eq 3 ]   #coordinates of C8
     then
       c8_x=$(awk -v r="$line_num" 'NR==r {print $2}' ${filename}_coords.xyz)
       c8_y=$(awk -v r="$line_num" 'NR==r {print $3}' ${filename}_coords.xyz)
       c8_z=$(awk -v r="$line_num" 'NR==r {print $4}' ${filename}_coords.xyz)
     fi

  #append (fragment) onto the atom label in the COM file for the EET gaussian calc
    awk -v r="$line_num" -v frag="(fragment=$l)" 'NR==r {print $1 frag " " $2 " " $3 " " $4}' ${filename}_coords.xyz >> ${filename}.com

  #increment line number
    line_num=$((line_num+1))
   done
fi
#Excecute python script to calculate the coordinates of the H atoms in the Me group with variables passed as arguments
#Append coordinates with (fragment) to the COM file

  for i in {1..3}
  do
    hme[$i]=$(python calc_methyl_coords.py $c8_x $c8_y $c8_z $n9_x $n9_y $n9_z $cme_x $cme_y $cme_z ${phiC[$i]} ${thetaC[$i]} ${CH_blC[$i]})
    echo "H(fragment=$l) "${hme[$i]} >> ${filename}.com
  done

done

#Insert blank line at the end of geometry specification
echo "" >> ${filename}.com

