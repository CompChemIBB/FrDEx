#!/usr/bin/env python3
"""
Generate OVERDIA input files for all unique fragment pairs.
"""

from pathlib import Path
import argparse

from state_classifier import (
    parse_gaussian_log,
    suggest_le_states,
    write_report,
)

# ------------------------------------------------------------
# Command line
# ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Generate OVERDIA input files."
)

parser.add_argument(
    "NFRAGS",
    type=int,
    help="Number of fragments"
)

parser.add_argument(
    "-s",
    "--system-name",
    default="test",
    help="System name"
)

parser.add_argument(
    "--auto-le",
    action="store_true",
    help="Automatically determine LE states from fragment TDDFT logs."
)

parser.add_argument(
    "--le-states",
    nargs="+",
    default=["s1", "s3"],
    help="Manual LE states."
)

parser.add_argument(
    "--ct-levels",
    type=int,
    default=2,
    help="Number of CT acceptor orbitals."
)

args = parser.parse_args()

NFRAGS = args.NFRAGS
SYSTEM_NAME = args.system_name
N_CT = args.ct_levels

OVERDIA_EXE = (
    "/home/roberto/JAMES/Programs/Overdia/1.1/source/overdia-par.e"
)

# ============================================================
# FIRST PASS
# Determine LE states for EVERY fragment
# ============================================================

fragment_le_states = {}

for frag in range(1, NFRAGS + 1):

    if args.auto_le:

        logfile = f"{SYSTEM_NAME}_frag{frag}-mon_sub-diab.log"

        print(f"Analyzing {logfile}")

        states = parse_gaussian_log(logfile)

        write_report(
            states,
            f"fragment{frag}.states.dat"
        )

        fragment_le_states[frag] = suggest_le_states(
            states,
            len(args.le_states)
        )

    else:

        fragment_le_states[frag] = list(args.le_states)

print()
print("Selected LE states")
print("------------------")

for frag in sorted(fragment_le_states):

    print(
        f"Fragment {frag}: {fragment_le_states[frag]}"
    )

# ------------------------------------------------------------
# Number of diabatic states
# ------------------------------------------------------------

N_LE = len(next(iter(fragment_le_states.values())))

TOTAL_STATES = 2 * N_LE + 2 * N_CT

# ============================================================
# SECOND PASS
# Generate all pair calculations
# ============================================================

run_script = Path("run_diab.sh")

with run_script.open("w") as runfile:

    runfile.write("#!/bin/bash\n\n")

    for frag1 in range(1, NFRAGS + 1):

        for frag2 in range(frag1 + 1, NFRAGS + 1):

            inp_name = (
                f"diab_{SYSTEM_NAME}_frag{frag1}_frag{frag2}.inp"
            )

            out_name = (
                f"diab_{SYSTEM_NAME}_frag{frag1}_frag{frag2}.out"
            )

            pair_log = (
                f"{SYSTEM_NAME}_frag{frag1}_frag{frag2}-mon_sub-diab.log"
            )

            frag1_fchk = (
                f"{SYSTEM_NAME}_frag{frag1}-mon_sub-diab"
            )

            frag2_fchk = (
                f"{SYSTEM_NAME}_frag{frag2}-mon_sub-diab"
            )

            lines = []

            # -------------------------------------------------
            # Header
            # -------------------------------------------------

            lines.append(f"'{pair_log}' log file of pair")
            lines.append("0                      number of modes")
            lines.append("0.1d0                  number of displacement size")
            lines.append("40                     number of TD states")
            lines.append("custom_states          possible ops KSORB/MXORB/FRAMO/FRATD/")
            lines.append(f"{TOTAL_STATES}                      number of diabatic states from input")
            lines.append("2                      number of fragments")
            lines.append(f"{frag1_fchk}   fchk files")
            lines.append(f"{frag2_fchk}   fchk files")

            # -------------------------------------------------
            # LE states
            # -------------------------------------------------

            lines.append("start_read_states")

            state_index = 1

            for state in fragment_le_states[frag1]:

                lines.append(
                    f"d{state_index} 1 {state}"
                )

                state_index += 1

            for state in fragment_le_states[frag2]:

                lines.append(
                    f"d{state_index} 2 {state}"
                )

                state_index += 1

            lines.append("end_read_states")

            # -------------------------------------------------
            # CT states
            # -------------------------------------------------

            lines.append("start_define_states")

            # frag1 -> frag2

            for level in range(N_CT):

                lines.append(f"d{state_index}_start")
                lines.append("1")
                lines.append(
                    f"1H0 2L{level} 1.d0"
                )
                lines.append(f"d{state_index}_end")

                state_index += 1

            # frag2 -> frag1

            for level in range(N_CT):

                lines.append(f"d{state_index}_start")
                lines.append("1")
                lines.append(
                    f"2H0 1L{level} 1.d0"
                )
                lines.append(f"d{state_index}_end")

                state_index += 1

            lines.append("end_define_states")

            # -------------------------------------------------

            with open(inp_name, "w") as inp:

                inp.write("\n".join(lines))
                inp.write("\n")

            runfile.write(
                f"{OVERDIA_EXE} < {inp_name} > {out_name}\n"
            )

print()

print(f"Fragments            : {NFRAGS}")
print(f"CT levels            : {N_CT}")
print(f"Diabatic states/pair : {TOTAL_STATES}")
print(f"Pair calculations    : {NFRAGS*(NFRAGS-1)//2}")

print()
print("Generated OVERDIA input files.")
print("Generated run_diab.sh")