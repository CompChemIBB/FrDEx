#!/usr/bin/env python3

# from os import system
import numpy as np
# from models import StateRegistry
from centroid import compute_fragment_centroids
from pathlib import Path



class ExcitonicSystem:
    """
    Global excitonic system assembled from all FrDEx pair calculations.
    """

    def __init__(self, registry):

        self.registry = registry

        n = len(registry)

        #
        # Hamiltonian
        #

        self.Hsum = np.zeros((n, n))
        self.Hcount = np.zeros((n, n), dtype=int)
        self.H = None

        #
        # Electronic Transition dipoles
        #

        self.mu_sum = np.zeros((n, 3))
        self.mu_count = np.zeros(n, dtype=int)
        self.mu = None

        #
        # Velocity dipoles
        #

        self.vel_sum = np.zeros((n, 3))
        self.vel_count = np.zeros(n, dtype=int)
        self.vel = None

        #
        # Magnetic dipoles
        #

        self.mag_sum = np.zeros((n, 3))
        self.mag_count = np.zeros(n, dtype=int)
        self.mag = None

        #
        # Geometry (filled later)
        #

        self.centroids = {}
        

        self.site_energies = None
        self.state_centers = np.zeros((n, 3))


# def build_system(pair_results, registry):

def build_system(
    pair_results,
    registry,
    input_dir=".",
    centroid_file="fragment_centroids.dat",
):

    system = ExcitonicSystem(registry)
    system.input_dir = input_dir

    input_dir = Path(input_dir)

    fragment_centroids = compute_fragment_centroids(
        directory=input_dir,
        outfile=centroid_file,
    )

    system.fragment_centroids = fragment_centroids

    #
    # Loop over all pair calculations
    #

    for pair in pair_results:

        Hpair = pair.hamiltonian

        mu = pair.dipoles
        vel = pair.velocity_dipoles
        mag = pair.magnetic_dipoles

        #
        # State vectors
        #

        for i, state in enumerate(pair.states):

            gi = state.global_index

            #
            # Transition dipoles
            #

            system.mu_sum[gi] += mu[i]
            system.mu_count[gi] += 1

            #
            # Velocity dipoles
            #

            system.vel_sum[gi] += vel[i]
            system.vel_count[gi] += 1

            #
            # Magnetic dipoles
            #

            system.mag_sum[gi] += mag[i]
            system.mag_count[gi] += 1

        #
        # Hamiltonian
        #

        for i, state_i in enumerate(pair.states):

            gi = state_i.global_index

            for j, state_j in enumerate(pair.states):

                gj = state_j.global_index

                system.Hsum[gi, gj] += Hpair[i, j]
                system.Hcount[gi, gj] += 1

    #
    # Average Hamiltonian
    #


    system.H = np.divide(
        system.Hsum,
        system.Hcount,
        out=np.zeros_like(system.Hsum),
        where=system.Hcount > 0
    )

     # no average
    # system.H = np.where(system.Hcount > 0, system.Hsum, 0)

    system.H = 0.5 * (system.H + system.H.T)

    #
    # Average transition dipoles
    #

    system.mu = np.divide(
        system.mu_sum,
        system.mu_count[:, None],
        out=np.zeros_like(system.mu_sum),
        where=system.mu_count[:, None] > 0
    )

    # system.state_centers = np.zeros((len(registry), 3))
    
    for state in registry.states:
    
        if state.state_type == "LE":
    
            center = fragment_centroids[state.fragment]
    
        else:
    
            c1 = fragment_centroids[state.donor]
            c2 = fragment_centroids[state.acceptor]
    
            center = 0.5 * (c1 + c2)
    
        system.state_centers[state.global_index] = center


    #
    # Average velocity dipoles
    #

    system.vel = np.divide(
        system.vel_sum,
        system.vel_count[:, None],
        out=np.zeros_like(system.vel_sum),
        where=system.vel_count[:, None] > 0
    )

    #
    # Average magnetic dipoles
    #

    system.mag = np.divide(
        system.mag_sum,
        system.mag_count[:, None],
        out=np.zeros_like(system.mag_sum),
        where=system.mag_count[:, None] > 0
    )
    

    # Site energies
    system.site_energies = np.diag(system.H).copy()

    return system

def print_system_summary(system):

    print()

    print("GLOBAL EXCITONIC SYSTEM")
    print("-----------------------")

    print(f"Number of states : {len(system.registry)}")

    print(f"Hamiltonian shape : {system.H.shape}")

    print(f"Symmetric : {np.allclose(system.H, system.H.T)}")

    print()

    print("Hamiltonian coverage")

    print(np.count_nonzero(system.Hcount),
          "/",
          system.Hcount.size)

    print()

    print("Average observations per state")

    print("Transition dipoles :",
          np.mean(system.mu_count))

    print("Velocity dipoles :",
          np.mean(system.vel_count))

    print("Magnetic dipoles :",
          np.mean(system.mag_count))

    print()

    print("Diagonal energies")

    print(np.round(np.diag(system.H), 5))


# def attach_centroids(system, centroid_dict):
#     """
#     Attach fragment centroids to the excitonic system.

#     Parameters
#     ----------
#     centroid_dict : dict

#         {
#             1 : np.array([x,y,z]),
#             2 : ...
#         }
#     """

#     system.centroids = centroid_dict

#     return system

# def write_cent_file(system, filename="cent.in"):
#     """
#     Write EXAT transition centers.
#     """

#     with open(filename, "w") as f:

#         for center in system.state_centers:

#             f.write(
#                 f"{center[0]:18.10f} "
#                 f"{center[1]:18.10f} "
#                 f"{center[2]:18.10f}\n"
#             )