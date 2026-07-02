#!/usr/bin/env python3

from models import StateRegistry


def build_registry(pair_results):
    """
    Build a unique registry of diabatic states from all pair calculations.
    """

    registry = StateRegistry()

    for pair in pair_results:

        for state in pair.states:

            key = state.key

            if key not in registry.lookup:

                state.global_index = len(registry.states)

                registry.lookup[key] = state.global_index

                registry.states.append(state)

            else:

                state.global_index = registry.lookup[key]

    return registry

def print_registry(registry):

    print("\nGLOBAL STATE REGISTRY\n")

    print(
        f"{'Idx':>4s} {'Type':>4s} {'Label':<15s}"
    )

    print("-"*35)

    for s in registry.states:

        print(
            f"{s.global_index:4d} "
            f"{s.state_type:>4s} "
            f"{s.label}"
        )

    print()

    print(f"Total states : {len(registry.states)}")