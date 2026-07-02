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

Replaces fragment geometries using heavy-atom RMSD alignment against optimized reference structures.

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

Optimized Gaussian input files for excited-state calculations.

---

# Step 3 — Generate OVERDIA inputs

```
prep-overdia_mod.py
```

## Purpose

Automatically generates FrDEx input files for every unique fragment pair.

For **N** fragments the total number of calculations is

```
N(N−1)/2
```

Example

```
8 fragments

↓

28 FrDEx calculations
```

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

---

## Output

```
diab_fratd_resp_frag1_frag2.inp

diab_fratd_resp_frag1_frag3.inp

...
```

---

# Step 4 — Run OVERDIA

Run every generated input file with a submission script. This **MUST** be adapted to your own OS.  


```
./run_diab.sh
```

Each calculation produces

```
diab_fratd_resp_fragX_fragY.out
```

Each output contains

- diabatic Hamiltonian
- transition dipoles
- velocity dipoles
- magnetic dipoles

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
