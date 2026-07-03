#!/usr/bin/env python3

"""
Generate OVERDIA input files for all unique fragment pairs
and create a run_diab.sh script to execute them.

Equivalent to the original Bash script.
"""

from pathlib import Path
import argparse

# ------------------------------------------------------------------
# User settings
# ------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Generate OVERDIA input files for all unique fragment pairs."
)
parser.add_argument("NFRAGS", type=int, help="number of fragments")
parser.add_argument('-s', '--system-name', default='test', help='System name for output files (default: test)')
args = parser.parse_args()
NFRAGS = args.NFRAGS

OVERDIA_EXE = (
    # "/home/roberto/JAMES/Programs/OVERDIA-REF/FRAG-ORB/overdia-par.e"
    "/home/roberto/JAMES/Programs/Overdia/1.1/source/overdia-par.e"
)

SYSTEM_NAME = args.system_name

# ------------------------------------------------------------------
# Create run script
# ------------------------------------------------------------------

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
                # f"{SYSTEM_NAME}_frag{frag1}-mon_sub-diab.fchk"
            )

            frag2_fchk = (
                f"{SYSTEM_NAME}_frag{frag2}-mon_sub-diab"
                # f"{SYSTEM_NAME}_frag{frag2}-mon_sub-diab.fchk"
            )

            # Rebuilt exactly to mirror the sample formatting provided in the reference image
            lines = [
                f"'{pair_log}' log file of pair",
                "0                      number of modes",
                "0.1d0                  number of displacement size",
                "40                     number of TD states",
                "custom_states          possible ops KSORB/MXORB/FRAMO/FRATD/",
                "8                      number of diabatic states from input",
                "2                      number of fragments",
                f"{frag1_fchk} fchk files",
                f"{frag2_fchk} fchk files",
                "start_read_states      !Block to read state vectors from input files",
                "d1 1 s1                !read diabatic state 1 from fragment 1 (C1), s1 state (pipil)",
                "d2 1 s3                !read diabatic state 2 from fragment 1 (C1), s3 state (pipi2)",
                "d3 2 s1                !read diabatic state 3 from fragment 2 (C2), s1 state (pipil)",
                "d4 2 s3                !read diabatic state 4 from fragment 2 (C2), s5 state (pipi2)",
                "end_read_states",
                "start_define_states    !Block to define states based on orbital transitions",
                "d5_start               !start the definition of diabatic state 5",
                "1                      !number of components of the state",
                "1H0 2L0 1.d0           !transition from HOMO fragment 1 to LUMO fragment 2",
                "d5_end                 !end the definition of diabatic state 5",
                "d6_start               !start the definition of diabatic state 6",
                "1                      !number of components of the state",
                "1H0 2L1 1.d0           !transition from HOMO fragment 1 to LUMO+1 fragment 2",
                "d6_end                 !end the definition of diabatic state 6",
                "d7_start               !start the definition of diabatic state 7",
                "1                      !number of components of the state",
                "2H0 1L0 1.d0           !transition from HOMO fragment 2 to LUMO fragment 1",
                "d7_end                 !end the definition of diabatic state 7",
                "d8_start               !start the definition of diabatic state 8",
                "1                      !number of components of the state",
                "2H0 1L1 1.d0           !transition from HOMO fragment 2 to LUMO+1 fragment 1",
                "d8_end                 !end the definition of diabatic state 8",
                "end_define_states"
            ]

            with open(inp_name, "w") as inp:
                inp.write("\n".join(lines))
                inp.write("\n")

            runfile.write(
                f"{OVERDIA_EXE} < {inp_name} > {out_name}\n"
                # f"{OVERDIA_EXE}  {inp_name} > {out_name}\n"
            )

print(
    f"Generated {NFRAGS*(NFRAGS-1)//2} input files "
    f"and run_diab.sh"
)