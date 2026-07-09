#!/usr/bin/env python3

"""
compare_spectra.py

Compare TDDFT and EXAT spectra.

Example
-------

python compare_spectra.py \
    MeGCGC-opt-wb97xd-gas-td70.log \
    results.results.out \
    --sigma 0.18
"""

import argparse
import re

import numpy as np
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Gaussian log parser
# ------------------------------------------------------------

def read_gaussian(logfile):

    energies = []
    osc = []
    rot = []

    with open(logfile) as f:
        lines = f.readlines()

    #
    # Excitation energies + oscillator strengths
    #

    for line in lines:

        if "Excited State" in line:

            m = re.search(
                r"Excited State\s+\d+.*?([0-9.]+)\s+eV.*?f=([0-9.]+)",
                line,
            )

            if m:

                energies.append(float(m.group(1)))
                osc.append(float(m.group(2)))

    #
    # Rotatory strengths
    #

    inside = False

    for line in lines:

        if "Rotatory Strengths" in line:

            inside = True
            continue

        if inside:

            if line.strip() == "":
                continue

            if "1/2[<0|r|b>" in line:
                continue

            toks = line.split()

            if len(toks) < 5:
                continue

            try:
                rot.append(float(toks[4]))
            except:
                pass

            if len(rot) == len(energies):
                break

    return (
        np.array(energies),
        np.array(osc),
        np.array(rot),
    )


# ------------------------------------------------------------
# EXAT parser
# ------------------------------------------------------------

def read_exat(filename):

    E = []
    osc = []
    rot = []

    with open(filename) as f:

        for line in f:

            if line.startswith("#"):
                continue

            toks = line.split()

            if len(toks) != 5:
                continue

            E.append(float(toks[1]))
            osc.append(float(toks[2]))
            rot.append(float(toks[4]))

    return (
        np.array(E),
        np.array(osc),
        np.array(rot),
    )


# ------------------------------------------------------------
# Gaussian broadening
# ------------------------------------------------------------

def broaden(E, I, sigma=0.18):

    xmin = min(E) - 1
    xmax = max(E) + 1

    x = np.linspace(xmin, xmax, 3000)

    y = np.zeros_like(x)

    pref = 1.0 / (sigma * np.sqrt(2*np.pi))

    for Ei, Ii in zip(E, I):

        y += Ii * pref * np.exp(
            -(x-Ei)**2/(2*sigma*sigma)
        )

    return x, y


# ------------------------------------------------------------

parser = argparse.ArgumentParser()

parser.add_argument("gaussian")

parser.add_argument("exat")

parser.add_argument(
    "--sigma",
    type=float,
    default=0.18,
)

args = parser.parse_args()

Eg, fg, rg = read_gaussian(args.gaussian)

Ee, fe, re = read_exat(args.exat)

#
# broaden
#

xg_uv, yg_uv = broaden(Eg, fg, args.sigma)
xe_uv, ye_uv = broaden(Ee, fe, args.sigma)

xg_cd, yg_cd = broaden(Eg, rg, args.sigma)
xe_cd, ye_cd = broaden(Ee, re, args.sigma)

#
# normalize
#

yg_uv /= yg_uv.max()
ye_uv /= ye_uv.max()

m = max(np.max(np.abs(yg_cd)),
        np.max(np.abs(ye_cd)))

yg_cd /= m
ye_cd /= m


# lam = 1239.84193 / energy

# order = np.argsort(lam)

# lam = lam[order]
# intensity = intensity[order]

#
# plotting
#

fig, ax = plt.subplots(2,2, figsize=(12,8))

#
# UV sticks
#

ax[0,0].vlines(Eg,0,fg,color="black",label="TDDFT")

ax[0,0].vlines(Ee,0,fe,color="red",label="EXAT")

ax[0,0].set_ylabel("Oscillator strength")

ax[0,0].set_title("Stick spectrum")

ax[0,0].legend()

#
# UV broadened
#


ax[0,1].plot(xg_uv,yg_uv,label="TDDFT")

ax[0,1].plot(xe_uv,ye_uv,label="EXAT")

ax[0,1].set_ylabel("Normalized intensity")

ax[0,1].set_title("Broadened UV")

ax[0,1].legend()

#
# CD sticks
#

ax[1,0].vlines(Eg,0,rg,color="black")

ax[1,0].vlines(Ee,0,re,color="red")

ax[1,0].axhline(0,color="k")

ax[1,0].set_ylabel("Rotatory strength")

ax[1,0].set_title("ECD sticks")

#
# CD broadened
#

ax[1,1].plot(xg_cd,yg_cd,label="TDDFT")

ax[1,1].plot(xe_cd,ye_cd,label="EXAT")

ax[1,1].axhline(0,color="k")

ax[1,1].set_ylabel("Normalized")

ax[1,1].set_title("Broadened ECD")

ax[1,1].legend()

for a in ax.flat:

    a.set_xlabel("Energy (eV)")

plt.tight_layout()

plt.show()