#!/usr/bin/env python3
"""
export_excitonic_system.py

Export the excitonic Hamiltonian and all associated quantities
from a complete OVERDIA calculation.

Usage
-----

python export_excitonic_system.py \
    diab_full.inp \
    diab_full.out \
    -o exported_system
"""

from pathlib import Path
import argparse
import numpy as np

from full_parser import parse_full_system


# ------------------------------------------------------------
# Writers
# ------------------------------------------------------------

def write_matrix(filename, matrix):

    np.savetxt(
        filename,
        matrix,
        fmt="%14.6f"
    )


def write_vectors(filename, vectors):

    np.savetxt(
        filename,
        vectors,
        fmt="%14.6f"
    )


def write_site_energies(filename, H):

    with open(filename, "w") as f:

        for i, e in enumerate(np.diag(H), start=1):

            f.write(f"{i:4d} {e:14.6f}\n")


def write_labels(filename, states):

    with open(filename, "w") as f:

        for i, state in enumerate(states, start=1):

            f.write(f"{i:4d} {state}\n")


def write_summary(filename, result):

    H = result.hamiltonian

    with open(filename, "w") as f:

        f.write("Excitonic system summary\n")
        f.write("========================\n\n")

        f.write(f"Number of diabatic states : {H.shape[0]}\n")

        f.write(f"Min energy : {np.min(np.diag(H)):.4f}\n")
        f.write(f"Max energy : {np.max(np.diag(H)):.4f}\n")

        off = H.copy()

        np.fill_diagonal(off, 0)

        f.write(
            f"Maximum |coupling| : {np.max(np.abs(off)):.4f}\n"
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("input")

    parser.add_argument("output")

    parser.add_argument(
        "-o",
        "--prefix",
        default="system"
    )

    args = parser.parse_args()

    result = parse_full_system(
        args.input,
        args.output
    )

    prefix = Path(args.prefix)

    write_matrix(
        prefix.with_suffix(".hamiltonian.dat"),
        result.hamiltonian
    )

    write_site_energies(
        prefix.with_suffix(".site.dat"),
        result.hamiltonian
    )

    write_vectors(
        prefix.with_suffix(".dipoles.dat"),
        result.dipoles
    )

    write_vectors(
        prefix.with_suffix(".velocity.dat"),
        result.velocity_dipoles
    )

    write_vectors(
        prefix.with_suffix(".magnetic.dat"),
        result.magnetic_dipoles
    )

    write_labels(
        prefix.with_suffix(".states.dat"),
        result.states
    )

    write_summary(
        prefix.with_suffix(".summary"),
        result
    )

    print()
    print("Export completed")
    print("----------------")
    print(prefix.with_suffix(".hamiltonian.dat"))
    print(prefix.with_suffix(".site.dat"))
    print(prefix.with_suffix(".dipoles.dat"))
    print(prefix.with_suffix(".velocity.dat"))
    print(prefix.with_suffix(".magnetic.dat"))
    print(prefix.with_suffix(".states.dat"))
    print(prefix.with_suffix(".summary"))


if __name__ == "__main__":
    main()