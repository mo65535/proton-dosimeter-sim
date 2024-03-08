#!/usr/bin/python


from argparse import ArgumentParser
from math import exp, log, pi, sqrt
from random import random, seed
import sys
from time import time


# Define constants
########################################################################

# Wall thickness i.e. R_max(T_p_max = 1 MeV)
R_max_in_um = 21.22
R_max_in_m = R_max_in_um * (10**-6)

# Chamber radius
R_chamber_in_um = 0.5
R_chamber_in_m  = R_chamber_in_um * (10**-6)

# Media densities
rho_A_in_g_per_cm3 = 0.94
rho_B_in_g_per_cm3 = 6.6 * (10**-4)

rho_A_in_g_per_m3 = rho_A_in_g_per_cm3 * (10**6)
rho_B_in_g_per_m3 = rho_B_in_g_per_cm3 * (10**6)

# Media masses
m_A_in_g = rho_A_in_g_per_m3 * pi*(R_chamber_in_m**2)*R_max_in_m
m_B_in_g = rho_B_in_g_per_m3 * (4.0/3.0)*pi*(R_chamber_in_m**3)

# Probabilities of ionization, excitation in the event of a proton collision
P_ionization = 0.75
P_excitation = 0.25


"""
Describe problem geometry
########################################################################

We'll take Z to be the direction of propagation, i.e. along the long axis
of the pill-shaped geometry.

Because all particles generated travel only in the z direction, a particle
travels through the same length of plastic in region A regardless of its
s coordinate. That length is R_max(T_p = 1 MeV).

This behavior also means that the chord length, L, for a given value of
s is the only relevant aspect of region B; its spherical shape can be
disregarded and it can instead be treated as half of an ellipsoid with
minor radii of R_chamber and a major radius of 2*R_chamber.

Thus, to simplify the logic for determining which region a particle is in,
the simulation geometry can be considered to have the shape of a cartridge.

We'll also set z = 0 at the intersection between regions A and B for
convenience and take points with z = 0, s < R_chamber to be in region B.

                                        z = 0
                                          |
                                          v
*******************************************************
*                                         *           ********
*                                         *                  ****
*                                         *                     ****
*                                         *                        **
*                                         *                          *
*              Region A                   *        Region B           *
*     -R_max(T_p = 1 MeV) < z < 0         *      0 <= z < L(s)        *
*                                         *                          *
*                                         *                        **
*                                         *                     ****
*                                         *                  ****
*                                         *           ********
*******************************************************

"""


def parse_input():
    """
    Parses input supplied by the user on the command line.

    Returns an object containing the parsed input. This object is the
    output of ArgumentParser.parse_args(), so it has member variables
    with names that are the same as the argument names specified below.
    """

    parser = ArgumentParser(description="""Simulate micro Monte Carlo for
                                      neutron beam on Rossi chamber.""")

    parser.add_argument('-v', '--verbose', action='store_true',
                        dest='verbose', default=True,
                        help='Print additional output')

    parser.add_argument('-q', '--quiet', action='store_false',
                        dest='verbose',
                        help='Do not print additional output')

    parser.add_argument('--fixed-seed', type=int,
                        dest='fixed_seed', default=None,
                        help='Seed the PRNG with the given value when run')

    parser.add_argument('--Tn-mev', type=float,
                        dest='n_energy', default=1.0, choices=[1.0, 10.0],
                        help='Neutron energy in MeV')

    parser.add_argument('--Tp-ev', type=float,
                        dest='p_energy_thresh', default=1.0, #0.025,
                        help='Proton energy threshold for stopping in eV')

    parser.add_argument('--track-deltas', action='store_true',
                        dest='track_deltas', default=False,
                        help='Track delta rays arising from protons')

    parser.add_argument('--force-b', action='store_true',
                        dest='force_B', default=False,
                        help='Force all neutron collisions to occur in B')

    parser.add_argument('-o', '--out-file',
                        dest='out_file', default='./output.txt',
                        help='The file where output will be saved')

    parser.add_argument('--n-collide', type=int,
                        dest='n_collide', default=0,
                        help='The number of neutron collisions desired')


    input_args = parser.parse_args()
    return input_args


def run_sim(ia):
    """
    Simulates the emission of neutrons at the chamber until the user's
    threshold is reached. Tracks the energy deposited in region B by
    particles arising from neutron reactions.

    A short name is selected for the input arguments because we'll
    reference this object often, but what we really want is its member
    variables, e.g. ia.n_emit, ia.n_energy, etc.

    Prints results to output file.

    Returns 0 on success and 1 on failure
    """

    if ia.fixed_seed:
        seed(ia.fixed_seed)
    else:
        seed(int(time()))

    # Print parameters to output file. If verbose, also print to console.
    with open(ia.out_file,'w') as fp:
        print_sim_input(fp, ia)
    if ia.verbose:
        print_sim_input(sys.stdout, ia)

    # These lists will record energies deposited in region B by protons
    # (and deltas, if desired) with the respective types of trajectory
    insiders, starters, stoppers, crossers, no_shows = [],[],[],[],[]

    # Make a new list where energy lists are decorated with trajectory names
    traj_and_energies = [('insiders', insiders), ('starters', starters),
                         ('stoppers', stoppers), ('crossers', crossers),
                         ('no-shows', no_shows)]

    # Build a list of (particles simulated, percentage complete) tuples
    completion = [(ia.n_collide*pcc/100, pcc) for pcc in xrange(5,101,5)]

    start_time = time()
    neutrons_collided = 0

    while True:
        if ia.verbose:
            for n, pcc in completion:
                if neutrons_collided == n:
                    print " {}% complete".format(pcc)

        # If there have been enough neutron collisions, end sim.
        if neutrons_collided >= ia.n_collide:
            break

        (s_p, z_p, T_p) = emit_interacting_neutron(ia.n_energy, ia.force_B)
        neutrons_collided += 1

        (trajectory, E_dep_B) = track_proton(s_p, z_p, T_p,
                                             ia.p_energy_thresh,
                                             ia.track_deltas)

        if E_dep_B > 0: # if energy was deposited in region B, record it
            for traj, energies in traj_and_energies:
                if trajectory in traj:
                    energies.append(E_dep_B)

    run_time = time() - start_time

    # Append output info to output file. If verbose, also print to console.
    with open(ia.out_file,'a') as fp:
        print_sim_results(fp, run_time, traj_and_energies)
    if ia.verbose:
        print_sim_results(sys.stdout, run_time, traj_and_energies)

    # Append lists of energies to output file
    with open(ia.out_file,'a') as fp:
        for traj, energies in traj_and_energies:
            fp.write("E_dep_B by {} [eV]:\n".format(traj))
            for energy in energies:
                fp.write("{}\n".format(energy*(10**6))) # convert MeV to eV

    print 'Done!'
    return 0


def print_sim_input(fp, ia):
    """
    Prints the input args ia (as well as some globals defined higher up
    in the file) to the file pointer specified by fp (e.g. the output file
    or sys.stdout)
    """
    print >> fp, "Simulation input"
    print >> fp, "----------------------------------------"
    #print >> fp, "P(A)      = {:1.8e}".format(m_A_in_g/(m_A_in_g + m_B_in_g))
    #print >> fp, "P(B)      = {:1.8e}".format(m_B_in_g/(m_A_in_g + m_B_in_g))
    #print >> fp, "P(A)/P(B) = {:1.8e}".format(m_A_in_g/m_B_in_g)
    print >> fp, "N collisions: {}".format(ia.n_collide)
    print >> fp, "Track deltas: {}".format(ia.track_deltas)
    print >> fp, "Force in B:   {}".format(ia.force_B)
    print >> fp, "Fixed seed:   {}".format(ia.fixed_seed)
    print >> fp, "T_neutron:    {} [MeV]".format(ia.n_energy)
    print >> fp, "T_p_stop:     {} [eV]".format(ia.p_energy_thresh)
    print >> fp, "Out file:     {}".format(ia.out_file)
    print >> fp, ""


def print_sim_results(fp, run_time, traj_and_energies):
    """
    Prints results of simulation to the file specified by file pointer.
    """
    print >> fp, "\nSimulation results"
    print >> fp, "----------------------------------------"
    print >> fp, "Time taken: {} seconds".format(int(run_time))
    print >> fp, "Of the protons that deposited energy in B:"
    for traj, energies in traj_and_energies:
        print >> fp, " Number of {}: {}".format(traj, len(energies))
    print >> fp, ""


def emit_interacting_neutron(T_n_in_MeV, force_B):
    """
    Simulates a neutron emitted from the spatially uniform primary beam
    at energy T_n in MeV, which then goes on to interact in either region
    A or region B.

    Returns a tuple containing
     - The s coordinate
     - The z coordinate
     - The energy of the resultant proton
    """

    s, z, T_p = None, None, None

    # Determine whether interaction is in region A or B
    in_B = False
    if random() < m_B_in_g / (m_A_in_g + m_B_in_g) or force_B:
        in_B = True

    # Find valid coordinates for the site of interaction
    while True:
        # Map randoms from [0,1) to [-R_c,R_c) before assigning to x and y
        x = ((random()*2)-1)*R_chamber_in_m
        y = ((random()*2)-1)*R_chamber_in_m

        if in_B:
            z = random() * 2*R_chamber_in_m
        else:
            z = random() * -R_max_in_m

        s = sqrt((x**2) + (y**2))

        if s < R_chamber_in_m:
            if in_B and z >= L_in_m(s):
                continue  # z was invalid for region B, try again
            break # s and z were valid for for the given region, leave loop
        # Otherwise, resume loop and try a new possible interaction site

    # Get a random energy between 0 and T_n_in_MeV
    T_p = random() * T_n_in_MeV

    return (s, z, T_p)


def track_proton(s, z, T, T_threshold, track_deltas):
    """
    Tracks the given proton (described by coordinates s and z, and kinetic
    energy T [MeV]) until it comes to a stop or leaves the simulated regions.
    "Coming to a stop" occurs when the energy drops below the user's
    specified threshold, T_threshold [eV]

    If the deltas flag is True, then also track all delta rays that arise
    from this proton.

    Returns a tuple containing
     - A string that describes the proton track with respect to region B
        i.e. "insider", "starter", "stopper", "crosser" or "no-show"
        where "no-show" refers to a proton that stopped in region A but
        produced a delta ray that went on to deposit energy in region B.
     - Energy deposited in region B by this particle (and child deltas)
    """

    # If we're tracking deltas, when a collision creates a delta ray, store
    # tuples with the initial values for z_e and T_e in this list for later
    # processing by the track_electron() function.
    deltas = []

    z_i = z # Record initial z so we can later determine trajectory type.
    E_dep_B = 0

    while T*(10**6) > T_threshold and z < L_in_m(s):
        # Determine which region we're in, set density accordingly
        rho = None
        if z < 0:
            rho = rho_A_in_g_per_cm3
        else:
            rho = rho_B_in_g_per_cm3

        # Calculate distance, d, to next interaction (eq. 3)
        d = (-1/mu_in_per_m(T,rho))*log(1-random())

        if z+d >= L_in_m(s): # If the move would leave region B,
            z += d  # take the step so trajectory can be determined later
            break   # we're done tracking this proton.

        if z < 0 and z+d >= 0: # If the move would go between A and B,
            z = 0     # put the proton just inside B
            continue  # and return to the top of the loop

        # If flow reaches this point, there was a proton collision to
        # account for, with equation 5(6)
        T_lost = 57.14 * (T**0.045) * (10**-6) # converting eV to MeV
        z += d
        T -= T_lost

        if track_deltas:
            if random() < P_ionization:   # If ionization occurred,
                deltas.append((z,T_lost)) # store delta for later.
                continue

        # If flow reaches this point, energy was deposited locally.
        if z >= 0:   # So, if the event was in B, record it.
            E_dep_B += T_lost

    # Now that the proton has stopped or left region B, handle any deltas
    while deltas:
        (z_e, T_e) = deltas.pop()
        E_dep_B += track_electron(s, z_e, T_e)

    trajectory = None
    if z_i >= 0 and z < L_in_m(s):
        trajectory = "insider"
    if z_i >= 0 and z >= L_in_m(s):
        trajectory = "starter"
    if z_i < 0 and z < L_in_m(s):
        trajectory = "stopper"
    if z_i < 0 and z >= L_in_m(s):
        trajectory = "crosser"
    if z_i < 0 and z < 0:
        trajectory = "no-show"

    return (trajectory, E_dep_B)


def track_electron(s, z, T):
    """
    Tracks the given electron (described by coordinates s and z, and kinetic
    energy T) until it comes to a stop or leaves the simulated regions.

    Returns the amount of energy deposited in region B by this particle
    """

    if z < 0:
        # Calculate range of electron in region A to see if it will reach B
        R_e_A = R_e_in_m(T, rho_A_in_g_per_cm3)

        if z + R_e_A >= 0:
            # It will cross into B
            T_lost_A = T * (0 - z) / R_e_A

            # Set new energy and move to region B
            T -= T_lost_A
            z = 0
        else:
            return 0  # No energy deposited in B if it can't reach B

    # Calculate range of electron in region B
    R_e_B = R_e_in_m(T, rho_B_in_g_per_cm3)

    # Electron will traverse a distance in B that is the shorter of range
    # in B or remaining distance before B is left (chord length - z)
    E_dep_B = T * min(L_in_m(s)-z, R_e_B) / R_e_B
    return E_dep_B


def L_in_m(s_in_m):
    """
    Returns the chord length through a sphere for the given value of s
    """
    return 2*sqrt((R_chamber_in_m**2)-(s_in_m**2))


def mu_in_per_m(T_p_in_MeV, rho_in_g_per_cm3):
    """
    Returns the linear attenuation coefficient for a proton
    given proton energy and the density of the medium.

    Add missing term for density and correct weird exponent in eq. 4:
    (mu / rho) * 10^-6 [cm^2/g] = 4.94 * T^-0.821 [cm^2/g]
    ==>
    mu [1/m] = 4.94 * T^-0.821 [cm^2/g] * rho [g/cm^3] * 10^6 * 100 [cm/m]
    """
    return 4.94 * (T_p_in_MeV**-0.821) * rho_in_g_per_cm3  * (10**8)


def R_e_in_m(T_e_in_MeV, rho_in_g_per_cm3):
    """
    Returns the electron range in a medium defined by density rho, for the
    given electron energy. Based on equation 7
    """
    T = T_e_in_MeV * (10**3)   # Convert to keV
    R_e = 1.264*(10**-6)*(T**2) + 1.225*(10**-5)*T - 1.6045*(10**-5)

    # Divide out density to get range in cm and convert to m
    return (R_e / rho_in_g_per_cm3) * 10**-2   # 10^-2 [m/cm]




# Here is the entry point when the script is run from the command line
if __name__ == '__main__':
    input_args = parse_input()
    sys.exit(run_sim(input_args))

