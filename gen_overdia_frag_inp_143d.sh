# Generates diabatisation input files for all unique fragment pairs among 12 fragments.
# Writes one .inp file per pair and builds a run_diab.sh script to execute OVERDIA.
# Uses a fixed OVERDIA executable path and hard-coded diabatisation options/state definitions.
# Creates input filenames with fragment indices and redirects OVERDIA output to .out files.
# Iterates over frag1=1..12 and frag2=frag1+1..12 to avoid duplicate pairs and self-pairs.
#!/bin/bash
#Script to generate input files for the fragment diabatisation scheme for GQuadruplex, and a shell script to run the calculations

#Location of diabaisation excecutable
OVERDIA_EXE="/home/roberto/JAMES/Programs/OVERDIA-REF/FRAG-ORB/overdia-par.e"

#Make header for file to run diabatisation code
echo "#!/bin/bash" > run_diab.sh

for frag1 in {1..12}
do
  for frag2 in $(seq $((frag1+1)) 12)
  do 
    echo "'CONCISE'                        printing option" > diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "'143d_1_frag${frag1}_frag${frag2}-mon_sub-diab.log'  log file of pair" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo  "40                                 number of TD states" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "'FRAMO'     possible ops KSORB/MXORB/FRAMO/FRATD/" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "8                 number of diabatic states from input" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "'143d_1_frag${frag1}-mon_sub-diab.fchk' '143d_1_frag${frag2}-mon_sub-diab.fchk'     fchk files" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "2   4             number of occ and virtual orb from input"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "0 'FR1' 0 'FR2'             occ orbitals expressed on these KS occ orbitals 0 is HOMO" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1 'FR1' 2 'FR1' 1 'FR2' 2 'FR2'      virt orbitals expressed on these KS virtual orb   1 is LUMO" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1.0 0.0     coef  Occ    mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "0.0 1.0     coef  Occ    mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1.0 0.0 0.0 0.0    coef  Virt   mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "0.0 1.0 0.0 0.0    coef  Virt   mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "0.0 0.0 1.0 0.0    coef  Virt   mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "0.0 0.0 0.0 1.0     coef  Virt  mixing the orbitals" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components first state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components second state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components third state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components fourth state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components fifth state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components sixth state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components seventh state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1                no components eighth state"  >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1    1   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1    2   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1    3   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "1    4   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "2    1   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "2    2   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "2    3   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp
    echo "2    4   1.d0" >> diab_143d_1_frag${frag1}_frag${frag2}.inp

    echo "${OVERDIA_EXE} < diab_143d_1_frag${frag1}_frag${frag2}.inp > diab_143d_1_frag${frag1}_frag${frag2}.out" >> run_diab.sh
  done
done

