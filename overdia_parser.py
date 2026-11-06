#!/usr/bin/env python3
"""
overdia_parser.py

Generic parser for OVERDIA output.

Extracts all quantities required by EXAT.

Works for:
    - pair calculations
    - full-system calculations
"""

import argparse
from pathlib import Path
from dataclasses import dataclass, field
import numpy as np
import re
import argparse


# ==========================================================
# Data container
# ==========================================================

@dataclass
class OverdiaCalculation:

    hamiltonian: np.ndarray = None

    site_energies: np.ndarray = None

    electric_dipoles: np.ndarray = None

    velocity_dipoles: np.ndarray = None

    magnetic_dipoles: np.ndarray = None

    state_labels: list = field(default_factory=list)

    chromophores: list = field(default_factory=list)

    centers: list = field(default_factory=list)


# ==========================================================
# Utilities
# ==========================================================

_float = r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?"


def _read_matrix(lines, start):

    matrix = []

    i = start

    while i < len(lines):

        line = lines[i].strip()

        if line == "":
            break

        try:
            row = [float(x) for x in line.split()]
            matrix.append(row)
        except:
            break

        i += 1

    return np.array(matrix), i


def _read_vectors(lines, start):

    vec = []

    i = start

    while i < len(lines):

        line = lines[i].strip()

        if line == "":
            break

        toks = line.split()

        if len(toks) < 3:
            break

        try:
            vec.append([
                float(toks[-3]),
                float(toks[-2]),
                float(toks[-1]),
            ])
        except:
            break

        i += 1

    return np.array(vec), i


# ==========================================================
# Main parser
# ==========================================================

def parse_overdia(outfile):

    outfile = Path(outfile)

    with outfile.open() as f:
        lines = f.readlines()

    result = OverdiaCalculation()

    i = 0

    while i < len(lines):

        line = lines[i]

        #
        # Hamiltonian
        #

        if "DIABATIC HAMILTONIAN" in line.upper():

            H, i = _read_matrix(lines, i + 2)

            result.hamiltonian = H

            continue

        #
        # Site energies
        #

        if "DIABATIC ENERGIES" in line.upper():

            energies = []

            j = i + 1

            while j < len(lines):

                m = re.search(_float, lines[j])

                if m is None:
                    break

                energies.append(float(m.group()))

                j += 1

            result.site_energies = np.array(energies)

            i = j

            continue

        #
        # Electric dipoles
        #

        if "TRANSITION DIPOLES DIABATIC STATES" in line.upper():

            mu, i = _read_vectors(lines, i + 2)

            result.electric_dipoles = mu

            continue

        #
        # Velocity dipoles
        #

        if "TRANSITION VELOCITY DIPOLES" in line.upper():

            vel, i = _read_vectors(lines, i + 2)

            result.velocity_dipoles = vel

            continue

        #
        # Magnetic dipoles
        #

        if "TRANSITION MAGNETIC DIPOLES" in line.upper():

            mag, i = _read_vectors(lines, i + 2)

            result.magnetic_dipoles = mag

            continue

        i += 1

    #
    # Build default labels
    #

    if result.site_energies is not None:

        n = len(result.site_energies)

        result.state_labels = [f"d{i+1}" for i in range(n)]

    return result


# ==========================================================
# Test
# ==========================================================

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
    description="Generate OVERDIA output files."
    )

    parser.add_argument(
        "filename",
        help="OVERDIA output file"
    )

    args = parser.parse_args()  

    calc = parse_overdia(args.filename)

    print()

    print("Hamiltonian")
    print(calc.hamiltonian)

    print()

    print("Site energies")
    print(calc.site_energies)

    print()

    print("Electric dipoles")
    print(calc.electric_dipoles)

    print()

    print("Velocity dipoles")
    print(calc.velocity_dipoles)

    print()

    print("Magnetic dipoles")
    print(calc.magnetic_dipoles)