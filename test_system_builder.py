#!/usr/bin/env python3

from pathlib import Path
import numpy as np
from pair_parser import parse_pair
from state_registry import build_registry
from system_builder import (
    build_system,
    print_system_summary
)

#
# Read all pair calculations
#

pair_results = []

for inp in sorted(Path(".").glob("*.inp")):

    out = inp.with_suffix(".out")

    if out.exists():

        pair_results.append(

            parse_pair(inp, out)

        )

print()

print(f"Read {len(pair_results)} pair calculations.")

#
# Build registry
#

registry = build_registry(pair_results)

print(f"Unique diabatic states : {len(registry)}")

#
# Build Hamiltonian
#

system = build_system(

    pair_results,

    registry

)

print_system_summary(system)

print()

print("First five averaged transition dipoles")

for state, mu in zip(registry.states[:5], system.mu[:5]):

    print(f"{state.label:15s}", np.round(mu, 5))

print()

print("First five magnetic dipoles")

for state, mu in zip(registry.states[:5], system.mag[:5]):

    print(f"{state.label:15s}", np.round(mu, 5))

print("=======")
for state, n in zip(registry.states, system.mu_count):

    print(state.label, n)