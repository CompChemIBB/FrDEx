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
from full_parser import parse_full_system
from state_registry import build_registry
from system_builder import (
    build_system,
    print_system_summary,
)
from exat_writer import write_all

# ---------------------------------------------------------
# Detect input format
# ---------------------------------------------------------

def detect_input_type(inpfile):
    """
    Detect whether the input file is

      - a FrDEx pair calculation
      - a full OVERDIA calculation

    FrDEx inputs contain the keyword
        start_read_states

    whereas original OVERDIA inputs do not.
    """

    with open(inpfile) as f:

        text = f.read().lower()

    if "start_read_states" in text:
        return "pair"

    return "full"

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

    for inp in inp_files:
    
        out = inp.with_suffix(".out")
    
        if not out.exists():
    
            print(f"WARNING : missing {out.name}")
    
            continue
    
        input_type = detect_input_type(inp)
    
        print(f"Parsing {inp.name} ({input_type})")
    
        if input_type == "pair":
    
            result = parse_pair(inp, out)
    
        else:
    
            result = parse_full_system(inp, out)
    
        pair_results.append(result)

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