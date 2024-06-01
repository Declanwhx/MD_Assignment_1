import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.constants as co

######################### CONSTANTS #########################
# Boltzmann Constant
k_B = 1.38065 * 10 ** -23  # [J/K]
# Lennard-Jones variables -> Depth of well and collisional diameter
epslj = 148 * co.k  # [J]
sigma = 3.73  # [Angstroms]
sigma_sq = sigma ** 2  # [Angstroms^2]
sigma_cu = sigma ** 3  # [Angstroms^3]
CH4_molar_mass = 16.04 * (10 ** -3)
avogadros_constant = 6.022 * 10 ** 23
# Degrees of freedom
dof = 3  # [-]
# System temperature in Kelvins
T = 150  # [K]
beta = 1 / (k_B * T)  # [1/J]


def convertMassDensity(mass_density):
    molecule_density = ((mass_density / CH4_molar_mass) * avogadros_constant) / (10 ** 30)
    return molecule_density


def readFile(file, l_domain):
    input_file = open(file, "r")
    # Reading and writing number of atoms/molecules in simulation box
    no_of_entities = int(input_file.readline())  # [-]
    # Simulation box size, note, we shall work with angstroms [Angstroms]
    domain_volume = l_domain ** 3
    # Calculating average number density of molecules in simulation
    molecule_density = no_of_entities / domain_volume  # [1/Angstroms^3]

    # Reading and writing coordinates into an array for Lennard-Jones Potential calculations
    molecules_coordinates = [[], [], []]
    input_coordinates = input_file.readlines()[1:]
    for line in input_coordinates:
        line = line.split()
        for i in range(0, 3):
            molecules_coordinates[i].append(float(line[i + 1]))
    return no_of_entities, molecule_density, molecules_coordinates


# 1.1
def totalEnergy(molecules_coordinates, l_domain):
    # molecules_coordinates = np.array(molecules_coordinates)
    no_of_entities = molecules_coordinates.shape[0]
    # Simulation box size, note, we shall work with angstroms [Angstroms]
    domain_volume = l_domain ** 3
    # Molecule density
    molecule_density = no_of_entities / domain_volume # [1/Angstroms^3]
    # Cutoff radius to prevent duplicate interactions [Angstroms]
    r_cut = 14

    U_lj = 0
    dU_dr = 0

    sr3 = sigma_cu / (r_cut ** 3)

    for (i, coordinates_i) in enumerate(molecules_coordinates):
        # list containing Δx, Δy, Δz for each pair
        d = (molecules_coordinates[i + 1:] - coordinates_i + l_domain / 2) % l_domain + l_domain / 2
        r_ij_sq = np.sum(d * d, axis=1)
        # print(type(r_ij_sq))

        sr2 = sigma_sq / r_ij_sq[r_ij_sq <= r_cut ** 2]  # filter out all r2 larger than rcut squared and get sigma^2/r^2 for all particles j>i
        sr6 = sr2 ** 3
        sr12 = sr6 ** 2
        print(sr2)

        U_lj += np.sum(sr12 - sr6)
        dU_dr += np.sum(sr6 - 2 * sr12)
    """
    for i in range(0, no_of_entities):
        # Starting from i+1 to prevent duplicate interactions
        for j in range(i + 1, no_of_entities):
            # The following 3 lines are a mathematical means of deriving the closest atom image to current atom per axis
            # Based on these 3 dimensions, a radial distance can be obtained
            x = (molecules_coordinates[j][0] - molecules_coordinates[i]
            [0] + l_domain / 2) % l_domain - l_domain / 2
            y = (molecules_coordinates[j][1] - molecules_coordinates[i]
            [1] + l_domain / 2) % l_domain - l_domain / 2
            z = (molecules_coordinates[j][2] - molecules_coordinates[i]
            [2] + l_domain / 2) % l_domain - l_domain / 2
            r_ij_sq = x ** 2 + y ** 2 + z ** 2
            # Implementation of the cut-off distance for duplicate interactions
            if r_ij_sq <= r_cut ** 2:
                sr2 = sigma_sq / r_ij_sq
                sr6 = sr2 * sr2 * sr2
                sr12 = sr6 * sr6
                U_lj += (sr12 - sr6)
                dU_dr += (sr6 - 2 * sr12)
                # Note: This should have been divided by the cutoff radius but because the virial equation
                # requires that we multiply by the cutoff radius again, that cancels out.
    """

    U_lj = 4 * epslj * U_lj  # [J]
    # print(U_lj)
    U_ljtail = no_of_entities * (8 / 3) * np.pi * epslj * molecule_density * sigma_cu * ((1 / 3) * (sr3 * sr3 * sr3) - sr3) # [J]
    # print(U_ljtail)
    # U_kin = 0.5 * co.k * T * dof * no_of_entities  # [J]
    # print(U_kin)
    U_tot = U_lj + U_ljtail # [J]
    # print(U_tot)

    # 1.3
    dU_dr = 24 * epslj * dU_dr
    P_tot = (molecule_density * k_B * T - (1 / (3 * domain_volume)) * dU_dr) * (10 ** 30)  # [Pa]
    # print(P)
    return U_tot, P_tot


# 1.2
# Note: particle no. specified here should be zero-indexed
def singleParticleEnergy(particle_no, molecules_coordinates, l_domain):
    no_of_entities = molecules_coordinates.shape[0]
    # print(no_of_entities)
    # Simulation box size, note, we shall work with angstroms [Angstroms]
    domain_volume = l_domain ** 3
    # Molecule density
    molecule_density = no_of_entities/domain_volume # [1/Angstroms^3]
    # Cutoff radius to prevent duplicate interactions [Angstroms]
    r_cut = 14

    sr3 = sigma_cu / (r_cut ** 3)

    # molecules_coordinates is a numpy array containing the coordinates of all the particles
    # this operation deducts the coordinate values from the respective coordinates of the other particles
    # and also reduces the original coordinates of the particle to a 0 list
    # to work around the 0 list, we set the rij value for that particle to be r_cut so that when we implement
    # the cut-off distance, this self-interaction is removed. In NumPy this is called broadcasting
    d = (molecules_coordinates - molecules_coordinates[particle_no] + l_domain / 2) % l_domain - l_domain / 2
    # print(d)
    r_ij_sq = np.sum(d * d, axis=1)
    r_ij_sq[particle_no] = r_cut ** 2
    # print(r_ij_sq)
    print(type(r_ij_sq))
    if np.any(r_ij_sq[:,] == 0):
        U_tot = 1e300
        return U_tot

    # returns list of sr2 values where r_ij_sq < r_cut^2
    sr2 = sigma_sq / r_ij_sq[r_ij_sq <= r_cut ** 2]
    # print(sr2)
    sr6 = sr2 ** 3
    sr12 = sr6 ** 2

    U_ljtail = (8 / 3) * np.pi * epslj * molecule_density * sigma_cu * ((1 / 3) * (sr3 * sr3 * sr3) - sr3)  # [J]
    # print(U_ljtail)
    U_lj = 4 * epslj * np.sum(sr12 - sr6) + U_ljtail  # [J]
    # print(U_lj)
    # U_kin = 0.5 * k_B * T * dof * no_of_entities  # [J]
    # print(U_kin)
    U_tot = U_lj + U_ljtail  # [J]
    # print(type(U_tot))

    # 1.3
    # Will not return this value since it's not needed for anything other than Q1.3
    # P = molecule_density * k_B * T - (1 / (3 * domain_volume)) * (24 * epslj * np.sum(sr6 - 2 * sr12))  # [Pa]
    # print(P)

    return U_tot


# 1.4
def translate(disp, molecules_coordinates, l_domain):
    # print(molecules_coordinates)
    no_of_entities = molecules_coordinates.shape[0]
    # Generates number between 0 and N, excludes N but has 0 -> Number here is the index of the particle we displace
    random_particle = np.random.randint(0, no_of_entities)

    accepted_move = 0

    # Calculate the single particle energy prior to displacement
    U_o = singleParticleEnergy(random_particle, molecules_coordinates, l_domain)
    U_tot, P_tot = totalEnergy(molecules_coordinates, l_domain)
    # print("Utotold: ", U_tot)
    # Generates a displacement array of 3
    # Delta setting -> 50% should be achieved with this value for LJ potential
    d = np.random.uniform(-disp, disp, size=3)
    # print(d)
    molecules_coordinates[random_particle,:] = (molecules_coordinates[random_particle] + d) % l_domain
    # print(molecules_coordinates)
    # Calculate the single particle energy post displacement
    U_n = singleParticleEnergy(random_particle, molecules_coordinates, l_domain)
    # print(U_o)
    # print(U_n)
    U_tot, P_tot = totalEnergy(molecules_coordinates, l_domain)
    # print("Utotnew: ", U_tot)

    # Energy change caused by displacement
    delta_U = U_n - U_o
    # print("ΔU = " + str(delta_U))
    if (delta_U < 0) or (np.exp(-delta_U / (co.k * T))) > np.random.uniform():
        # Counter for accepted move
        accepted_move += 1
    else:
        # Reverts coordinates back to original if move is rejected
        molecules_coordinates[random_particle,:] = (molecules_coordinates[random_particle] - d) % l_domain

    return molecules_coordinates, accepted_move


# 1.5
def startConf(disp, mass_density, l_domain):
    molecule_density = convertMassDensity(mass_density)
    # print(molecule_density)
    no_of_entities = int(np.ceil(molecule_density * (l_domain ** 3)))
    # print(no_of_entities)

    # Generating coordinates for molecules
    molecules_coordinates = np.random.uniform(0, l_domain, size=(no_of_entities, 3))
    # print(molecules_coordinates)
    # Recommended number of moves
    # Nint = no_of_entities * 50
    # Set Nint = 500000 when testing for disp value influence on acceptance rate
    Nint = 1
    accepted_moves = 0

    # Invoking Nint number of displacements to initialize the simulation box (Prevent overlaps)
    for i in range(0, Nint):
        print("Configuration: " + str(i + 1))
        # print(molecules_coordinates)
        molecules_coordinates = translate(disp, molecules_coordinates, l_domain)[0]
        accepted_moves += translate(disp, molecules_coordinates, l_domain)[1]

    rate_of_acceptance = round((accepted_moves / Nint) * 100,2)
    # print("The rate of acceptance is: ", rate_of_acceptance, "%")

    # Code can be modified to take delta values in an array and then generate plots
    """
    plt.scatter(disp, rate_of_acceptance)
    plt.title("Rate of acceptance vs Δ (N = 500,000)")
    plt.xlabel("Δ [Å]")
    plt.ylabel("Rate of acceptance [%]")
    plt.ylim(0,100)
    plt.savefig('Acceptance_vs_Delta_1.png')
    """

    """
    # Writing new .xyz file -> Might be wrong? But not important
    final_molecules_coordinates = open('generated_box.xyz', "w")
    final_molecules_coordinates.write(str(no_of_entities) + "\n")
    final_molecules_coordinates.write("box.pdb \n")
    for i in range(0, no_of_entities):
        final_molecules_coordinates.write(
            "{0:<6}{1:>12}{2:>15}{3:>16}".format("C", molecules_coordinates[0][i], molecules_coordinates[1][i],
                                                 molecules_coordinates[2][i]) + "\n")
    """

    return molecule_density, no_of_entities, molecules_coordinates


def averages(no_of_cycles, total_energy, total_pressure):
    ensemble_average_total_energy = total_energy / no_of_cycles
    ensemble_average_total_pressure = total_pressure / no_of_cycles
    return ensemble_average_total_energy, ensemble_average_total_pressure


def MC_NVT(disp, mass_density, l_domain):
    start = time.time()

    # Create box with specified mass_density and l_domain params using the startConf function
    molecule_density, no_of_entities, molecules_coordinates = startConf(disp, mass_density, l_domain)

    # No. of MC cycles
    Nint = 500000
    # Intervals for sampling
    intervals = 500
    # Total sample no. if we sample for every 500 cycles
    sample_no = int(Nint / intervals)
    # Array of cycles for graph plots
    cycles = np.zeros(sample_no+1)

    # Arrays to contain sample variables from the respective configuration
    U_tot = np.zeros(sample_no+1)
    P_tot = np.zeros(sample_no+1)

    for i in range(0, Nint):
        print("Configuration: " + str(i + 1))
        molecules_coordinates, accepted_move = translate(disp, molecules_coordinates, l_domain)
        # Conditional sampling -> Take samples at only every 500 cycles
        if i == 0:
            U_tot[0], P_tot[0] = totalEnergy(molecules_coordinates, l_domain)
        elif (i+1) % intervals == 0:
            print(i)
            cycle_index = (i+1)/intervals
            cycles[int(cycle_index)] = i+1
            U_tot[int(cycle_index)], P_tot[int(cycle_index)] = totalEnergy(molecules_coordinates, l_domain)

    print(U_tot)
    average_of_system_variables = averages(Nint, sum(U_tot), sum(P_tot))
    U_tot_average = average_of_system_variables[0]
    P_tot_average = average_of_system_variables[1]

    end = time.time()

    plt.scatter(cycles, U_tot)
    plt.title("Total Energy vs No. of cycles")
    plt.ylabel("U [J]")
    plt.xlabel("Number of cycles [-]")
    plt.savefig('Energy_vs_cycles.png')
    plt.show()

    print("The average total energy is: ", U_tot_average, "J")
    print("The average total pressure is: ", P_tot_average, "Pa")
    print("The time of execution of above program is :", (end - start) / 60, "mins")

    return U_tot_average, P_tot_average


# N = readFile('box.xyz', 30)[0]
# rho = readFile('box.xyz', 30)[1]
# coordinates = readFile('box.xyz', 30)[2]

# totalEnergy(N, rho, coordinates)
# singleParticleEnergy(361, N, rho, coordinates)
# translate(N, rho, coordinates)
startConf(0.45, 358.4, 30)

 #MC_NVT(0.45,358.4,30)
