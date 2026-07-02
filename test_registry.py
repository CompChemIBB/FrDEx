#!/usr/bin/env python3

from pathlib import Path

from pair_parser import parse_pair
from state_registry import (
    build_registry,
    print_registry
)

pair_results = []

for inp in sorted(Path(".").glob("*.inp")):

    out = inp.with_suffix(".out")

    if not out.exists():

        continue

    pair_results.append(

        parse_pair(inp, out)

    )

registry = build_registry(pair_results)

print_registry(registry)