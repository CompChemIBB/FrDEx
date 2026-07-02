from dataclasses import dataclass
import numpy as np



@dataclass
class State:

    label: str

    state_type: str

    fragment: int | None = None

    donor: int | None = None

    acceptor: int | None = None

    root: int | None = None

    local_index: int | None = None

    global_index: int | None = None

    @property
    def key(self):
        """
        Unique identifier for a diabatic state.
        """

        if self.state_type == "LE":
            return ("LE", self.fragment, self.root)

        elif self.state_type == "CT":
            return ("CT", self.donor, self.acceptor)

        else:
            raise ValueError(
                f"Unknown state type: {self.state_type}"
            )
    

@dataclass
class FrDExPairResult:
    """
    Result of one FrDEx pair calculation.
    """

    pair: tuple[int, int]

    states: list[State]

    hamiltonian: np.ndarray

    dipoles: np.ndarray

    velocity_dipoles: np.ndarray

    magnetic_dipoles: np.ndarray

class StateRegistry:

    def __init__(self):

        self.states = []

        self.lookup = {}

    def __len__(self):
        return len(self.states)

    def __iter__(self):
        return iter(self.states)