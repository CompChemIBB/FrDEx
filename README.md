# FrDEx → EXAT Workflow

A Python toolkit for converting **FrDEx** pairwise diabatic calculations into a complete excitonic Hamiltonian that can be directly used as input for **EXAT**.

The package reconstructs a global excitonic model starting from pairwise FrDEx calculations and automatically generates all input files required by EXAT for the calculation of absorption, Circular Dichroism (CD), and Linear Dichroism (LD) spectra.

---

# Overview

The complete workflow is

```
PDB Structure
      │
      ▼
fragment_build_xyz.py
      │
      ▼
Fragment geometries (.xyz, .com)
      │
      ▼
replace_optimize_heavy.py
      │
      ▼
Optimized fragment geometries
      │
      ▼
Gaussian TDDFT calculations
      │
      ▼
prep-overdia_mod.py
      │
      ▼
FrDEx input generation
      │
      ▼
FrDEx calculations
      │
      ▼
build_exat.py
      │
      ▼
Global excitonic Hamiltonian
      │
      ▼
EXAT input files
      │
      ▼
EXAT
      │
      ▼
Absorption / CD / LD spectra
```

---

# Features

The package automatically

- extracts fragment excited states from FrDEx calculations
- identifies identical diabatic states appearing in multiple fragment pairs
- reconstructs the global excitonic Hamiltonian
- averages transition properties over all available calculations
- computes transition centers from fragment geometries
- generates EXAT external input files
- interfaces directly with EXAT

---

# Requirements

## Python

Python **3.10** or newer is recommended.

Required packages

```bash
numpy
rmsd
```

Install them with

```bash
pip install numpy rmsd
```

---

# External Software

This workflow requires three external programs.

## Gaussian

Used for

- geometry optimization
- TDDFT excited-state calculations

---

## OVERDIA

Used to compute

- diabatic Hamiltonians
- transition dipoles
- magnetic transition dipoles
- velocity transition dipoles

for every fragment pair.

---

## EXAT

Used to compute

- excitonic Hamiltonian
- absorption spectra
- Circular Dichroism (CD)
- Linear Dichroism (LD)

---

# Installing EXAT

The current public version of EXAT requires **Python 2.7** together with older scientific libraries.

Since modern Linux distributions (Ubuntu ≥ 22.04) no longer provide Python 2 support, the recommended solution is to install EXAT inside an isolated Conda environment.

## Install Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

source ~/.bashrc
```

---

## Create a Python 2.7 environment

```bash
conda create -n py27 python=2.7

conda activate py27
```

---

## Install the required packages

```bash
conda install \
-c free \
-c defaults \
numpy=1.9.1 \
scipy=0.15.1 \
matplotlib=1.4.2 \
pygtk
```

---

## Fix the libgfortran issue

On modern Linux distributions you may encounter

```
ImportError: libgfortran.so.1
```

Locate your installed library

```bash
find /usr/lib -name "libgfortran.so.*"
```

Then create the symbolic link

```bash
ln -s \
/usr/lib/x86_64-linux-gnu/libgfortran.so.5 \
~/miniconda3/envs/py27/lib/libgfortran.so.1
```

(adjust the path according to your installation)

---

## Running EXAT

```bash
conda activate py27

python exat.py
```

---

# Workflow

---

# Step 1 — Build fragment geometries

The first step is to remove the sugar-phosphate backbone and unnecessary bases. This is is done automatically for any PDB which leave only cytosine and guanine bases removing all the rest. Currently, only these bases have been tested. The hydrogen atoms **MUST** be added beforehand according to your preferred protonation state. The program uses **calc_methyl_coords.py** module to remove all sugar base and replace them with methyl group.  


```
fragment_build_xyz.py
```

## Purpose

Reads a PDB structure and extracts every nucleobase into an independent fragment.

The script generates

- Gaussian input files
- optional XYZ files

---

## Usage

```bash
./fragment_build_xyz.py molecule.pdb
```

or

```bash
./fragment_build_xyz.py molecule.pdb --xyz
```

---

## Output

```
frag1.com
frag2.com
...

frag1.xyz
frag2.xyz
...
```

The generated XYZ files are later used to compute fragment centroids.

---

# Step 2 — Replace fragment geometries

```
replace_optimize_heavy.py
```

## Purpose

Replaces fragment geometries using heavy-atom RMSD alignment against optimized reference structures. It rotates optimised geometry to match the input geometry. The rotated optimized structure is saved in an array. The RMSD of this rotational transformation will be displayed. We expect the RMSD below 0.2. The saved rotated optimized coordinates will be used to generate TD-DFT input files. Monomer and pairs files can be generated. Currently, only pair base calculations are supported. Trimer or other moieties geometries will be explored in future. 

Currently supported references

- standard cytosine
- hemi-protonated cytosine

The procedure preserves the fragment orientation while replacing the coordinates with optimized geometries.

---

## Usage

```bash
./replace_optimize_heavy.py \
-i16 cyt16.xyz \
-i17 cyt17.xyz \
-n 8 \
-s quartet
```

---

## Output

Optimized Gaussian input files for excited-state calculations with the following header. 

```
1. Monomer calculation: #p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(3/33=4)
td(nstates=10, conver=6) iop(9/40=5) charge nosymm
```

```
2. Dimer: #p m052x/6-31G(d) scrf=(pcm,solvent=water) scf=tight iop(9/9=3) iop(3/33=4)
td(nstates=40, conver=6) iop(9/40=5) charge nosymm
Dimer calculations are generally used for Hinter Hamiltonian which involves CT states. There-
fore we need a minimum of 40 States to start with.
```

---

# Step 3 — Generate OVERDIA inputs

```
prep-overdia_mod.py
```

## Purpose
Once the TD-DFT results are ready, we can start our FrDEx Calculation. Generally, we will use Dimer calculation for calculating Excitonic coupling and CT couplings. Here our H-intra terms also comes from Pair calculation. Automatically generates FrDEx input files for every unique fragment pair. For **N** fragments the total number of calculations is N(N−1)/2. Example, 8 fragments → 28 FrDEx calculations.

---

## Usage

```bash
./prep-overdia_mod.py 8
```

or

```bash
./prep-overdia_mod.py 8 \
-s quartet
```
Here is a sample Input file for overdia 1.1 input file with comments for a dimer calculation with 2 Local Excitation in each Fragment and 4CT states.

```
{pair_log} log file of pair
0                       number of modes
0.1d0                   number of displacement size
40                      number of TD states
custom_states           possible ops KSORB/MXORB/FRAMO/FRATD/
8                       number of diabatic states from input
2                       number of fragments
{frag1_fchk}   fchk files without extension 
{frag2_fchk}   fchk files  
start_read_states       !Block to read state vectors from input files
d1 1 s1                 !read diabatic state 1 from fragment 1 (C1), s1 state (pipi1)
d2 1 s3                 !read diabatic state 2 from fragment 1 (C1), s3 state (pipi2)
d3 2 s1                 !read diabatic state 3 from fragment 2 (C2), s1 state (pipi1)
d4 2 s3                 !read diabatic state 4 from fragment 2 (C2), s5 state (pipi2)
end_read_states
start_define_states     !Block to define states based on orbital transitions
d5_start                !start the definition of diabatic state 5
1                       !number of components of the state
1H0 2L0 1.d0            !transition from HOMO fragment 1 to LUMO fragment 2
d5_end                  !end the definition of diabatic state 5
d6_start                !start the definition of diabatic state 6
1                       !number of components of the state
1H0 2L1 1.d0            !transition from HOMO fragment 1 to LUMO+1 fragment 2
d6_end                  !end the definition of diabatic state 6
d7_start                !start the definition of diabatic state 7
1                       !number of components of the state
2H0 1L0 1.d0            !transition from HOMO fragment 2 to LUMO fragment 1
d7_end                  !end the definition of diabatic state 7
d8_start                !start the definition of diabatic state 8
1                       !number of components of the state
2H0 1L1 1.d0            !transition from HOMO fragment 2 to LUMO+1 fragment 1
d8_end                  !end the definition of diabatic state 8
end_define_states
```


---

## Output

```
diab_fratd_resp_frag1_frag2.inp
diab_fratd_resp_frag1_frag3.inp
...
```


---

# Step 4 — Run OVERDIA

Run every generated input file with a submission script. This **MUST** be adapted to your own OS. The program is quite quick and it can be run locally. Only the current path is working **/home/roberto/JAMES/Programs/Overdia/1.1/source/overdia-par.e** 

```
#!/bin/tcsh -f
cd /home/roberto/IMOT/EFALBO/FrDEX/NEW
setenv OMP_NUM_THREADS 1

foreach f1 (1 2 3 4 5)
  @ start = $f1 + 1
  foreach f2 (`seq $start 6`)
/home/roberto/JAMES/Programs/Overdia/1.1/source/overdia-par.e  < diab_test_frag${f1}_frag${f2}.inp > diab_test_frag${f1}_frag${f2}.
out
  end
end
```

or submitting to a HPC-Scheduler (not tested). 

```
./run_diab.sh
```

Each calculation produces

```
diab_fratd_resp_fragX_fragY.out
```

Each output contains. Below there is an example 

- diabatic Hamiltonian
- transition dipoles
- velocity dipoles
- magnetic dipoles

```
Hamiltonian in the diabatic (overlap) basis (8x8 Matrix, Energy in eV):

          State 1    State 2    State 3    State 4    State 5    State 6    State 7    State 8
State 1:  5.09178   -0.00015    0.02015    0.01224    0.00115    0.00064   -0.00115   -0.00011
State 2: -0.00015    6.77737   -0.00079    0.01002    0.00035    0.00121   -0.00230   -0.00090
State 3:  0.02015   -0.00079    5.08167   -0.00016   -0.00114   -0.00006    0.00087    0.00018
State 4:  0.01224    0.01002   -0.00016    6.79181   -0.00010   -0.00096    0.00098    0.00129
State 5:  0.00115    0.00035   -0.00114   -0.00010    6.58998   -0.00225   -0.00000   -0.00000
State 6:  0.00064    0.00121   -0.00006   -0.00096   -0.00225    8.20203    0.00001   -0.00000
State 7: -0.00115   -0.00230    0.00087    0.00098   -0.00000    0.00001    6.77193    0.01585
State 8: -0.00011   -0.00090    0.00018    0.00129   -0.00000   -0.00000    0.01585    8.17386

TRANSITION DIPOLES DIABATIC STATES (Electric Dipoles, x y z format per state)
State 1:  -0.4286812   1.5114117  -0.0654172
State 2:   0.7309090   0.6993476   0.0350359
State 3:  -0.9073163   1.2679686   0.0675020
State 4:   0.4162707   0.9085905   0.0704559
State 5:   0.0009530  -0.0005924   0.0021263
State 6:   0.0007624  -0.0005075   0.0013506
State 7:  -0.0002055   0.0001566  -0.0036612
State 8:   0.0013849  -0.0001245  -0.0026196

TRANSITION VELOCITY DIPOLES DIABATIC STATES (Velocity Dipoles, x y z format per state)
State 1:   0.0911998  -0.2897161   0.0126361
State 2:  -0.1837120  -0.1699774  -0.0091568
State 3:   0.1817079  -0.2403075  -0.0126019
State 4:  -0.1059482  -0.2241269  -0.0172307
State 5:  -0.0004258   0.0002497  -0.0007042
State 6:  -0.0000392  -0.0001089  -0.0003190
State 7:   0.0010670  -0.0000941   0.0000189
State 8:  -0.0001890   0.0004077   0.0006046

TRANSITION MAGNETIC DIPOLES DIABATIC STATES (Magnetic Dipoles, x y z format per state)
State 1: -11.0611415  -3.3718956   2.5555534
State 2:  -6.3674059   6.9285226  -0.3884603
State 3:  -6.2079205  -4.8278801   2.6422359
State 4:  -5.7267382   2.7296817  -0.3795756
State 5:   0.0175498   0.0164066  -0.0035308
State 6:   0.0052206   0.0036829  -0.0007626
State 7:  -0.0058325   0.0032128  -0.0000573
State 8:  -0.0063028   0.0051785  -0.0019562
```


---

# Step 5 — Build the global excitonic Hamiltonian

```
build_exat.py
```

This is the main reconstruction program. It reads every FrDEx calculation and automatically reconstructs the global excitonic model. Internally it uses

```
parse_input.py
parse_output.py
pair_parser.py
state_registry.py
system_builder.py
centroid.py
exat_writer.py
models.py
```

---

## What it does

The program

- parses every FrDEx input/output pair
- identifies duplicate diabatic states
- builds a global state registry
- reconstructs the complete Hamiltonian
- averages repeated Hamiltonian elements
- averages transition dipoles
- averages magnetic dipoles
- averages velocity dipoles
- computes transition centers
- writes EXAT input files

---

## Global state registry

Every diabatic state receives a unique global index.

Example

```
0   LE_F1_S1

1   LE_F1_S3

2   LE_F2_S1

3   LE_F2_S3

4   CT_F1_F2

5   CT_F2_F1

...
```

---

## Transition centers

Fragment centers are computed automatically from

```
frag*.xyz
```

Rules

### Local Excited (LE) states

```
Center = fragment centroid
```

### Charge Transfer (CT) states

```
Center = midpoint(donor, acceptor)
```

---

# Step 6 — Generate EXAT input files

The program writes

```
chromlist.in

site.in

coup.in

dipo.in

dipovel.in

dipomag.in

cent.in
```

These correspond to the **external input** option of EXAT.

---

## chromlist.in

Chromophore list

Example

```
1 1
2 1
3 1
...
```

---

## site.in

Site energies

```
1 5.2034
2 6.7093
...
```

---

## coup.in

Upper-triangular excitonic couplings

Units

```
cm⁻¹
```

---

## dipo.in

Electric transition dipoles

```
μx  μy  μz
```

---

## dipovel.in

Velocity transition dipoles

---

## dipomag.in

Magnetic transition dipoles

---

## cent.in

Transition centers

LE

```
Fragment centroid
```

CT

```
Midpoint between donor and acceptor
```

---

# Step 7 — Run EXAT

Example in py27 environment. 

```
python exat.py \
-e \
-omag \
-o results \
-nloc 24
```

This computes

- excitonic states
- oscillator strengths
- rotatory strengths
- localized states

---

# Step 8 — Compute spectra

Example in py27 environment. 

```
python spectrum.py \
results.results.out \
-o spectrum \
--unitbroad eV \
--sigma 0.18 \
-s -0.65
```

The program generates

- Absorption spectrum
- Circular Dichroism
- Linear Dichroism

---

# Repository Structure

```
.
├── fragment_build_xyz.py
├── replace_optimize_heavy.py
├── prep-overdia_mod.py
│
├── build_exat.py
├── parse_input.py
├── parse_output.py
├── pair_parser.py
├── state_registry.py
├── system_builder.py
├── centroid.py
├── exat_writer.py
├── models.py
│
├── test_parser.py
├── test_registry.py
├── test_system_builder.py
│
└── README.md
```

---

# Scientific Workflow

```
PDB
 │
 ▼
Fragment extraction
 │
 ▼
Geometry replacement
 │
 ▼
Gaussian TDDFT
 │
 ▼
FrDEx
 │
 ▼
Pair diabatic Hamiltonians
 │
 ▼
Automatic reconstruction
 │
 ▼
Global excitonic Hamiltonian
 │
 ▼
EXAT
 │
 ▼
Absorption
Circular Dichroism
Linear Dichroism
```

---

# Citation

If this software contributes to published work, please cite the original publications describing

- OVERDIA
- EXAT
- FrDEX

and acknowledge the use of this workflow in the computational methods section.
