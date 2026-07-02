# FrDEx → EXAT Workflow

Python utilities for converting **FrDEx** pairwise diabatic calculations into a complete excitonic Hamiltonian suitable for **EXAT** calculations.

The workflow automates the complete procedure:

PDB → fragment geometries → optimized fragments → FrDEx input generation → FrDEx calculations → global excitonic Hamiltonian → EXAT input files → absorption/CD spectra.

---

# Workflow

```
PDB structure
      │
      ▼
fragment_build_xyz.py
      │
      ▼
fragment xyz files
      │
      ▼
replace_optimize_heavy.py
      │
      ▼
Gaussian excited-state inputs
      │
      ▼
prep-overdia_mod.py
      │
      ▼
FrDEx pair calculations
      │
      ▼
FrDEx parser
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

# Requirements

## Python

Python ≥ 3.10

Required packages

```bash
numpy
```

Install with

```bash
pip install numpy
```

---

## External software

The workflow assumes the following external programs are available.

### Gaussian

Used for TDDFT calculations and fragment optimization.

---

### FrDEx

Used to compute diabatic Hamiltonians for every fragment pair.

---

### EXAT

Used to compute

- Excitonic Hamiltonian
- Absorption spectra
- Circular Dichroism
- Linear Dichroism

---

# Step 1
## Build fragment geometries

```
fragment_build_xyz.py
```

### Purpose

Reads a PDB structure and extracts each nucleobase into an independent fragment.

Produces

- Gaussian input files
- optional XYZ files

### Usage

```bash
./fragment_build_xyz.py molecule.pdb
```

or

```bash
./fragment_build_xyz.py molecule.pdb --xyz
```

### Output

```
frag1.com
frag2.com
...

frag1.xyz
frag2.xyz
...
```

The XYZ files are later used to compute fragment centroids.

---

# Step 2
## Replace fragment geometries

```
replace_optimize_heavy.py
```

### Purpose

Replaces each fragment geometry using heavy-atom RMSD alignment against optimized reference structures.

Reference geometries

- standard cytosine
- hemi-protonated cytosine

The script preserves fragment orientation while replacing coordinates with optimized geometries.

### Usage

```bash
./replace_optimize_heavy.py \
-i16 cyt16.xyz \
-i17 cyt17.xyz \
-n 8 \
-s quartet
```

### Output

Optimized fragment geometries suitable for excited-state calculations.

---

# Step 3
## Generate FrDEx inputs

```
prep-overdia_mod.py
```

### Purpose

Creates FrDEx input files for every unique fragment pair.

For N fragments the script produces

```
N(N−1)/2
```

pair calculations.

Example

8 fragments

```
28 pair calculations
```

Usage

```bash
./prep-overdia_mod.py 8
```

or

```bash
./prep-overdia_mod.py 8 -s quartet
```

Example outputs

```
diab_fratd_resp_frag1_frag2.inp

diab_fratd_resp_frag1_frag3.inp

...
```

---

# Step 4
## Run FrDEx

Run every generated input.

Outputs

```
diab_fratd_resp_frag1_frag2.out

...

```

Each output contains

- diabatic Hamiltonian
- transition dipoles
- velocity dipoles
- magnetic dipoles

---

# Step 5
## Assemble the global excitonic Hamiltonian

The Python modules included in this repository automatically merge all pair calculations.

Modules

```
parse_input.py
parse_output.py
pair_parser.py
state_registry.py
system_builder.py
centroid.py
exat_writer.py
```

The code

- identifies identical diabatic states
- removes duplicates
- builds one global state registry
- averages repeated quantities
- reconstructs the complete Hamiltonian

---

## Global state registry

Each diabatic state receives a unique global index.

Example

```
0  LE_F1_S1
1  LE_F1_S3
2  LE_F2_S1
3  LE_F2_S3
4  CT_F1_F2
5  CT_F2_F1
...
```

---

## Averaged quantities

The program averages

- Hamiltonian elements
- transition dipoles
- velocity dipoles
- magnetic dipoles

over every occurrence of the same state in different pair calculations.

---

## Transition centers

Fragment centers are computed automatically from

```
frag*.xyz
```

Rules

LE states

```
center = fragment centroid
```

CT states

```
center = midpoint(donor,acceptor)
```

---

# Step 6
## Generate EXAT input

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

These files correspond to the **external input** mode of EXAT.

---

## chromlist.in

Maps chromophore IDs.

Example

```
1 1
2 1
3 1
...
```

---

## site.in

Diagonal Hamiltonian elements

```
1 5.20
2 6.71
...
```

---

## coup.in

Upper triangular couplings

Units

```
cm−1
```

---

## dipo.in

Electric transition dipoles

```
μx μy μz
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

For LE states

```
fragment centroid
```

For CT states

```
midpoint between donor and acceptor
```

---

# Step 7
## Run EXAT

Example

```bash
exat.py \
-e \
-omag \
-o test \
-nloc 24
```

This produces

- excitonic Hamiltonian
- oscillator strengths
- rotatory strengths
- localized states
- spectra input

---

# Step 8
## Compute spectra

Example

```bash
spectrum.py \
test.results.out \
-o spectrum \
--unitbroad eV \
--sigma 0.18 \
-s -0.65
```

Outputs

- Absorption
- Circular Dichroism
- Linear Dichroism

---

# Repository structure

```
.
├── fragment_build_xyz.py
├── replace_optimize_heavy.py
├── prep-overdia_mod.py
│
├── parse_input.py
├── parse_output.py
├── pair_parser.py
├── state_registry.py
├── system_builder.py
├── centroid.py
├── exat_writer.py
│
├── models.py
│
├── build_exat.py
│
├── test_parser.py
├── test_registry.py
├── test_system_builder.py
│
└── README.md
```

---

# Scientific workflow

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
Spectra
```

---

# Citation

If this workflow contributes to published work, please cite

- FrDEx
- EXAT
- the original theoretical methodology used for the diabatic Hamiltonian construction.
# FrDEx
