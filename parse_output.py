#!/usr/bin/env python3
#  parse_output.py

import re
import numpy as np

def parse_hamiltonian(lines):

    blocks = []

    i = 0

    while i < len(lines):

        line = lines[i].strip()

        if line.startswith("ROOT NO."):

            cols = list(map(int, re.findall(r"\d+", line)))

            block = []

            i += 1

            while i < len(lines):

                line = lines[i].strip()

                # next ROOT block
                if line.startswith("ROOT NO."):
                    break

                # blank line
                if not line:
                    i += 1
                    continue

                toks = line.split()

                # skip eigenvalue line
                if not toks[0].isdigit():
                    i += 1
                    continue

                row = int(toks[0])

                vals = list(map(float, toks[1:]))

                block.append((row, vals))

                i += 1

            # <<< VERY IMPORTANT
            blocks.append((cols, block))

        else:

            i += 1

    n = max(max(cols) for cols, _ in blocks)

    H = np.zeros((n, n))

    for cols, block in blocks:

        for row, vals in block:

            for col, val in zip(cols, vals):

                H[row-1, col-1] = val

    return H


def extract_block(lines,start_marker,end_marker=None):

    start=None

    for i,line in enumerate(lines):

        if start_marker in line:
            start=i
            break

    if start is None:
        raise ValueError(start_marker)

    if end_marker:

        end=None

        for j in range(start+1,len(lines)):

            if end_marker in lines[j]:
                end=j
                break

        return lines[start:end]

    return lines[start:]

# parse_vector_block.py
def parse_vector_block(lines, header, nstates):

    start = None

    for i, line in enumerate(lines):

        if header in line:

            start = i + 1

            break

    if start is None:

        raise ValueError(
            f"Cannot find section: {header}"
        )

    data = []

    for line in lines[start:]:

        toks = line.split()

        if len(toks) != 3:
            continue

        try:

            data.append(
                list(map(float, toks))
            )

        except ValueError:
            continue

        if len(data) == nstates:
            break

    return np.array(data)

def parse_frdex_output(outfile):

    with open(outfile) as f:

        lines = f.readlines()

    ham_block = extract_block(
    lines,
    "Hamiltonian in the diabatic",
    "EIGENVALUES Hamiltonian matrix"
    )

    # This is a debug print to check the size of the Hamiltonian block
    # print("HAM_BLOCK SIZE =", len(ham_block))
    
    # for line in ham_block[:30]:
        # print(line.rstrip())

    # ham_block = extract_block(
    #     lines,
    #     "Hamiltonian in the diabatic"
    # )

    H = parse_hamiltonian(
        ham_block
    )

    n = H.shape[0]

    mu = parse_vector_block(
        lines,
        "TRANSITION DIPOLES DIABATIC STATES",
        n
    )

    vel = parse_vector_block(
        lines,
        "TRANSITION VELOCITY DIPOLES DIABATIC STATES",
        n
    )

    mag = parse_vector_block(
        lines,
        "TRANSITION MAGNETIC DIPOLES DIABATIC STATES",
        n
    )

    return H, mu, vel, mag

if __name__ == "__main__":

    H, mu, vel, mag = parse_frdex_output(
        "diab_fratd_resp_4ct_imot_frag1_frag2.out"
    )

    print("\nHamiltonian")
    print(H)

    print("\nShape:", H.shape)

    print("\nDiagonal:")
    print(np.diag(H))

    print("\nTransition dipoles:")
    print(mu)

    print("\nVelocity dipoles:")
    print(vel)

    print("\nMagnetic dipoles:")
    print(mag)