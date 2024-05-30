import numpy as np
import time
import matplotlib.pyplot as plt

######################### CONSTANTS #########################
# Boltzmann Constant
k_B = 1.38065 * 10 ** -23  # [J/K]
# Lennard-Jones variables -> Depth of well and collisional diameter
epslj = 148 * k_B  # [J]
sigma = 3.73  # Note: this was 0.373 nm but is converted to Angstroms to maintain unit consistency
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


def readFile(file,l_domain):
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
def totalEnergy(no_of_entities, molecule_density, molecules_coordinates, l_domain):
    # Simulation box size, note, we shall work with angstroms [Angstroms]
    domain_volume = l_domain ** 3
    # Cutoff radius to prevent duplicate interactions [Angstroms]
    # r_cut = l_domain / 2
    r_cut = 14

    U_lj = 0
    U_ljtail = 0
    dU_dr = 0

    sr3 = sigma_cu / (r_cut ** 3)

    for i in range(0, no_of_entities):
        for j in range(i + 1, no_of_entities): # Starting from i+1 to prevent duplicate interactions
            # The following 3 lines are a mathematical means of deriving the closest atom image to current atom per axis
            # Based on these 3 dimensions, a radial distance can be obtained
            x = (molecules_coordinates[0][j]-molecules_coordinates[0][i]+l_domain/2) % l_domain-l_domain/2
            y = (molecules_coordinates[1][j]-molecules_coordinates[1][i]+l_domain/2) % l_domain-l_domain/2
            z = (molecules_coordinates[2][j]-molecules_coordinates[2][i]+l_domain/2) % l_domain-l_domain/2
            r_ij_sq = round((x ** 2 + y ** 2 + z ** 2), 4)
            r_ij = np.sqrt(r_ij_sq)
            # print(r_ij)
            # Implementation of the cut-off distance for duplicate interactions
            if r_ij <= r_cut:
                sr2 = sigma_sq / r_ij_sq
                sr6 = sr2 * sr2 * sr2
                sr12 = sr6 * sr6
                U_lj += (sr12 - sr6)
                dU_dr += (sr6 - 2 * sr12)
                # Note: This should have been divided by the cutoff radius but because the virial equation
                # requires that we multiply by the cutoff radius again, that cancels out.
        U_ljtail += ((1 / 3) * (sr3 * sr3 * sr3) - sr3)

    U_lj = 4 * epslj * U_lj  # [J]
    # print(U_lj)
    U_ljtail = (8 / 3) * np.pi * epslj * molecule_density * sigma_cu * U_ljtail  # [J]
    # print(U_ljtail)
    U_kin = 0.5 * k_B * T * dof * no_of_entities  # [J]
    # print(U_kin)
    U_tot = U_lj + U_ljtail + U_kin  # [J]
    # print(U_tot)

    # 1.3
    dU_dr = 24 * epslj * dU_dr
    P = (molecule_density * k_B * T - (1 / (3 * domain_volume)) * dU_dr) * (10 ** (30))  # [Pa]
    # print(P)
    return U_tot, P


# 1.2
# Note: particle no. specified here should be zero-indexed
def singleParticleEnergy(particle_no, no_of_entities, molecule_density, molecules_coordinates, l_domain):
    # Simulation box size, note, we shall work with angstroms [Angstroms]
    domain_volume = l_domain ** 3
    # Cutoff radius to prevent duplicate interactions [Angstroms]
    # r_cut = l_domain / 2
    r_cut = 14

    U_lj = 0
    dU_dr = 0

    for i in range(0, no_of_entities):
        # If condition implemented to prevent calculation of self-interaction
        if particle_no == i:
            continue
        # The following 3 lines are a mathematical means of deriving the closest molecule/image to current molecule per axis
        # Based on these 3 dimensions, a radial distance can be obtained
        x = (molecules_coordinates[0][i] - molecules_coordinates[0][particle_no] + l_domain / 2) % l_domain - l_domain / 2
        y = (molecules_coordinates[1][i] - molecules_coordinates[1][particle_no] + l_domain / 2) % l_domain - l_domain / 2
        z = (molecules_coordinates[2][i] - molecules_coordinates[2][particle_no] + l_domain / 2) % l_domain - l_domain / 2
        r_ij_sq = round((x ** 2 + y ** 2 + z ** 2), 4)
        r_ij = np.sqrt(r_ij_sq)
        # Implementation of the cut-off distance for duplicate interactions
        if r_cut >= r_ij > 0:
            sr2 = sigma_sq / r_ij_sq
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6
            U_lj += (sr12 - sr6)
            # dU_dr += (sr6 - 2 * sr12)
            # Note: This should have been divided by the cutoff radius but because the virial equation
            # requires that we multiply by the cutoff radius again, that cancels out.
        # Need to check for overlaps because when we generate random coordinates, overlaps can occur and if not
        # accounted for, the script will not run, so instead of letting U go to infinity in such instances, we impose
        # a large value such that the displacement function will displace it and clear the overlap.
        # elif r_ij < sigma:
        elif r_ij == 0:
            U_lj = 1 * 10**300

    U_lj = 4 * epslj * U_lj  # [J]
    # print(U_lj)
    sr3 = sigma_cu / (r_cut ** 3)
    U_ljtail = (8 / 3) * np.pi * epslj * molecule_density * sigma_cu * ((1 / 3) * (sr3 * sr3 * sr3) - sr3)  # [J]
    # print(U_ljtail)
    U_kin = 0.5 * k_B * T * dof * no_of_entities  # [J]
    # print(U_kin)
    U_tot = U_lj + U_ljtail + U_kin  # [J]
    # print(U_tot)

    # 1.3
    # dU_dr = 24 * epslj * dU_dr
    # P = molecule_density * k_B * T - (1 / (3 * domain_volume)) * dU_dr  # [Pa] CHECK UNITS AGAIN
    # print(P)
    return U_tot


# 1.4
def translate(disp, no_of_entities, molecule_density, molecules_coordinates, l_domain):
    # Generates number between 0 and N, excludes N but has 0 -> Number here is the index of the particle we displace
    random_particle = np.random.randint(0, no_of_entities)
    # Delta setting -> 50% should be achieved with this value for LJ potential
    disp_range = disp - (-disp)

    original_coordinates = []
    new_coordinates = []

    accepted_move = 0

    # Calculate the single particle energy prior to displacement
    U_o = singleParticleEnergy(random_particle, no_of_entities, molecule_density, molecules_coordinates, l_domain)

    for i in range(0, 3):
        # adding original coordinates of specific particle to an array
        original_coordinates.append(molecules_coordinates[i][random_particle])
        # Casting random displacement for each axis coordinate
        random_disp = (np.random.uniform(0.000,1.000) - 0.5) * disp_range
        new_coordinates.append(molecules_coordinates[i][random_particle] + random_disp)
        # PBC to bring back to original box
        if new_coordinates[i] > l_domain:
            new_coordinates[i] -= l_domain
        elif new_coordinates[i] < 0:
            new_coordinates[i] += l_domain
        molecules_coordinates[i][random_particle] = round(new_coordinates[i],5)
    # print(original_coordinates)
    # print(random_Displacements)
    # print(new_coordinates)

    # Calculate the single particle energy post displacement
    U_n = singleParticleEnergy(random_particle, no_of_entities, molecule_density, molecules_coordinates, l_domain)
    # print(U_o)
    # print(U_n)

    # Energy change caused by displacement
    delta_U = U_n - U_o
    # print("ΔU = " + str(delta_U))

    # Acceptance criteria
    acc_o_n = min(1, np.exp(-beta * delta_U))
    # print("The acceptance rate is: " + str(acc_o_n))

    # Random no. drawn to compare with acceptance criteria
    draw_random_No = np.random.uniform(0.000,1.000)
    # print("Random number drawn: " + str(draw_random_No))

    # Reverts coordinates back to original if move is rejected
    if acc_o_n < draw_random_No:
        for i in range(0,3):
            molecules_coordinates[i][random_particle] = original_coordinates[i]
    # Counter for accepted move
    else:
        accepted_move += 1

    return molecules_coordinates, accepted_move


# 1.5
def startConf(disp, mass_density, l_domain):
    molecule_density = convertMassDensity(mass_density)
    print(molecule_density)
    no_of_entities = int(np.ceil(molecule_density * (l_domain ** 3)))
    print(no_of_entities)

    # Generating coordinates for molecules
    molecules_coordinates = [[], [], []]
    for i in range(0, no_of_entities):
        for j in range(0,3):
            molecules_coordinates[j].append(round((np.random.random_sample() * l_domain), 5))

    # Recommended number of moves
    Nint = no_of_entities * 50
    # accepted_moves = 0

    # Invoking Nint number of displacements to initialize the simulation box (Prevent overlaps)
    for i in range(0,Nint):
        print("Configuration: " + str(i + 1))
        MC_cycle = translate(disp, no_of_entities, molecule_density, molecules_coordinates, l_domain)
        molecules_coordinates = MC_cycle[0]
        # accepted_moves += MC_cycle[1]

    # rate_of_acceptance = round((accepted_moves / Nint) * 100,2)
    # print("The rate of acceptance is: ", rate_of_acceptance, "%")

    # Code can be modified to take delta values in an array and then generate plots
    """
    plt.scatter(disp, rate_of_acceptance)
    plt.title("Rate of acceptance vs Δ (N = 500,000)")
    plt.xlabel("Δ [Å]")
    plt.ylabel("Rate of acceptance [%]")
    plt.ylim(0,100)
    plt.savefig('Acceptance_vs_Delta_2.png')
    plt.show()
    """

    # Writing new .xyz file -> Might be wrong? But not important
    final_molecules_coordinates = open('generated_box.xyz', "w")
    final_molecules_coordinates.write(str(no_of_entities) + "\n")
    final_molecules_coordinates.write("box.pdb \n")
    for i in range(0, no_of_entities):
        final_molecules_coordinates.write("{0:<6}{1:>12}{2:>15}{3:>16}".format("C", molecules_coordinates[0][i], molecules_coordinates[1][i], molecules_coordinates[2][i]) + "\n")

    return molecule_density, no_of_entities, molecules_coordinates

def averages(no_of_cycles, total_energy, total_pressure):
    ensemble_average_total_energy = total_energy/no_of_cycles
    ensemble_average_total_pressure = total_pressure/no_of_cycles
    return ensemble_average_total_energy, ensemble_average_total_pressure


def MC_NVT(disp, mass_density, l_domain):
    start = time.time()

    # Create box with specified mass_density and l_domain params using the startConf function
    initialize_box = startConf(disp, mass_density, l_domain)
    molecule_density = initialize_box[0]
    no_of_entities = initialize_box[1]
    molecules_coordinates = initialize_box[2]

    # No. of MC cycles
    Nint = 500000
    # Array of cycles for graph plots
    cycles = []

    # Arrays to contain sample variables from the respective configuration
    U_tot = []
    P_tot = []

    # Initial configuration variable samples
    U_tot.append(totalEnergy(no_of_entities, molecule_density, molecules_coordinates, l_domain)[0])
    P_tot.append(totalEnergy(no_of_entities, molecule_density, molecules_coordinates, l_domain)[1])

    for i in range(0,Nint):
        cycles.append(i+1)
        print("Configuration: " + str(i+1))
        MC_cycle = translate(disp, no_of_entities, molecule_density, molecules_coordinates, l_domain)
        molecules_coordinates = MC_cycle[0]
        # If condition to determine if particle was displaced, if not, take the U and P values from
        # the configuration before since nothing has changed
        if MC_cycle[1] == 0 and i != 0:
            U_tot.append(U_tot[i-1])
            P_tot.append(P_tot[i-1])
        # Else, recalculate total energy and pressure
        else:
            system_variables = totalEnergy(no_of_entities, molecule_density, molecules_coordinates, l_domain)
            U_tot.append(system_variables[0])
            P_tot.append(system_variables[1])
            
    average_of_system_variables = averages(i + 1, sum(U_tot), sum(P_tot))
    U_tot_average.append(average_of_system_variables[0])
    P_tot_average.append(average_of_system_variables[1])

    end = time.time()

    plt.scatter(U_tot, cycles)
    plt.title("Total Energy vs No. of cycles")
    plt.ylabel("U [J]")
    plt.xlabel("Number of cycles [-]")
    plt.savefig('Energy_vs_cycles.png')
    plt.show()

    print("The average total energy is: " + str(U_tot_average) + "J")
    print("The average total pressure is: " + str(P_tot_average) + "Pa")
    print("The time of execution of above program is :", (end - start) / 60, "mins")

    return U_tot_average, P_tot_average


# N = readFile('box.xyz', 30)[0]
# rho = readFile('box.xyz', 30)[1]
# coordinates = readFile('box.xyz', 30)[2]

# totalEnergy(N, rho, coordinates)
# singleParticleEnergy(361, N, rho, coordinates)
# translate(N, rho, coordinates)

# startConf(0.5, 358.4, 30)

MC_NVT(0.5,358.4,30)






