#!/usr/bin/env python3

"""
build_exat.py

Collect all FrDEx pair calculations and generate the EXAT input files.

Example
-------

Current directory

    build_exat.py

Specify input/output directories

    build_exat.py \
        -i pair_calculations \
        -o exat_input

"""

from pathlib import Path
import argparse

from pair_parser import parse_pair
from state_registry import build_registry
from system_builder import (
    build_system,
    print_system_summary,
)
from exat_writer import write_all


# ---------------------------------------------------------
# Command line
# ---------------------------------------------------------

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Build EXAT input files from FrDEx pair calculations."
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        default=".",
        help="Directory containing FrDEx *.inp/*.out files (default: current directory)"
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory where EXAT input files will be written (default: current directory)"
    )

    # parser.add_argument(
    #     "--xyz-pattern",
    #     default="frag*.xyz",
    #     help="Pattern used to locate fragment xyz files (default: frag*.xyz)"
    # )

    return parser.parse_args()


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    args = parse_arguments()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)

    print()
    print("FrDEx → EXAT")
    print("-------------")
    print(f"Input directory : {input_dir}")
    print(f"Output directory: {output_dir}")

    #
    # Locate FrDEx calculations
    #

    inp_files = sorted(input_dir.glob("*.inp"))

    print()
    print(f"Found {len(inp_files)} FrDEx input files.")

    pair_results = []

    #
    # Parse every pair
    #

    for inp in inp_files:

        out = inp.with_suffix(".out")

        if not out.exists():

            print(f"WARNING : missing {out.name}")

            continue

        pair_results.append(
            parse_pair(inp, out)
        )

    if len(pair_results) == 0:

        raise RuntimeError(
            "No complete FrDEx calculations (*.inp/*.out) were found."
        )

    #
    # Build registry
    #

    registry = build_registry(pair_results)

    print()
    print(f"Unique diabatic states : {len(registry)}")

    #
    # Assemble excitonic system
    #

    system = build_system(
        pair_results,
        registry,
        input_dir=input_dir,
    )
    #
    # Print summary
    #

    print_system_summary(system)

    #
    # Write EXAT files
    #

    write_all(
        system,
        output_dir,
    )

    print()
    print(f"EXAT input files written to")
    print(f"    {output_dir}")

    print()
    print("Done.")


if __name__ == "__main__":

    main()