#!/usr/bin/env python3

from pathlib import Path

from pair_parser import parse_pair

from state_registry import build_registry

from system_builder import (
    build_system,
    print_system_summary
)

from exat_writer import write_all


def main():

    #
    # Locate all FrDEx calculations
    #

    pair_results = []

    inp_files = sorted(Path(".").glob("*.inp"))

    print()

    print(f"Found {len(inp_files)} FrDEx calculations.")

    #
    # Parse every pair
    #

    for inp in inp_files:

        out = inp.with_suffix(".out")

        if not out.exists():

            print(f"Missing {out.name}")

            continue

        pair_results.append(
            parse_pair(inp, out)
        )

    #
    # Build global registry
    #

    registry = build_registry(pair_results)

    print()

    print(f"Unique diabatic states : {len(registry)}")

    #
    # Assemble global excitonic system
    #

    system = build_system(
        pair_results,
        registry
    )

    #
    # Print summary
    #

    print_system_summary(system)

    #
    # Write EXAT input files
    #

    write_all(system)

    print()

    print("Done.")


if __name__ == "__main__":

    main()