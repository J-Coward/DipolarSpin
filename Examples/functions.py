import quspin as qs
import numpy as np
import matplotlib.pyplot as plt
from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spin_basis_1d # Hilbert space spin basis
from quspin.operators import quantum_operator
import math
from quspin.basis import spin_basis_general
from matplotlib import colors
import random
import time


pulse_time = 0.0478e-3 # in s
free_time = 15000 * pulse_time # 0.6e-3 # in s
rabi_frequency = 2*np.pi*(1/pulse_time)/4 # From RbCs Ramsey Paper https://doi.org/10.1088/2058-9565/aaee35 Figure 7
# print((pulse_time*rabi_frequency)/(np.pi/2)) # Check pulse time and rabi frequency give a pi/2 pulse area 


####################################################################################


# Generates a matrix with the pairwise interaction strengths between a particle at the center and each lattice site
# Lx and Ly are the systems side lengths - Defult Ly = Lx
# theta and phi describe the orientation of the quantization axis - Defult pointing out of the plane of the system
# J is the interaction strength bewtween nearest neighbor particles
# xs and ys define how 'stretched' the system is in the x and y directions (could be simplified into 1 ratio variable in future)

def intmatrix(Lx, J, xs=1, ys=1, theta=0, phi=0, Ly='_', plot=False):
    # Define side length of simulated system as L
    # Handle odd and even side lengths such that every interaction in the system can be described with the interaction matrix
    if Lx % 2 == 1:
        pass
    else:
        Lx += 1

    if Ly == '_': # Defult equal length sides
        Ly = Lx
    else:
        if Ly % 2 == 1:
            pass
        else:
            Ly += 1
    
    mat = np.zeros((Lx, Ly))
    
    L_xrange, L_yrange = range(Lx), range(Ly)
    center_xindex = int(L_xrange[-1]/2)
    center_yindex = int(L_yrange[-1]/2)
    
    
    center = np.array([xs*center_xindex, ys*center_yindex, 0]) # Position of matrix center (in 3D)
    
    quan_axis = [np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)] # Quantization axis as an unit vector
    
    # Calculate each interaction strength
    for i in L_xrange:
        for j in L_yrange:
            displacement = np.array([xs*i, ys*j, 0]) - center
            r = np.linalg.norm(displacement)
            angle = np.arccos(displacement @ quan_axis/r)
            mat[i, j] = J*(1-3*(np.cos(angle)**2))/(r**3)
    mat[center_xindex, center_yindex] = 0.0 # No interaction with itself


    # Does not run by defult
    if plot == True:
        fig, ax = plt.subplots()

        mat, xside, yside = intmatrix(Lx, J, Ly=Ly, theta=theta, phi=phi)
        couplings = np.transpose(mat) # Transpose to show x and y axis correctly (this will also switch x and y corrdinate ordering when indexing)

        cmap = plt.cm.RdBu
        ax.matshow(couplings, cmap=cmap, norm=colors.CenteredNorm())

        for i in range(len(couplings)):
            for j in range(len(couplings[0])): 
                c = '{:g}'.format(float('{:.2g}'.format(couplings[i, j])))
                #c = '(' + str(j) + ', ' + str(i) + ')'
                    
                ax.text(j, i, str(c), va='center', ha='center')
    
    return mat, Lx, Ly #indexing interaction matrix is done with mat[x coordinate, y coordinate]


####################################################################################


# Uniformly distribute N particles into a lattice with x, y dimensions

def ran_occupied_sites(N, x, y): # Dimensions of simulated lattice (in units of lattice sites) and no of particles to distribute
    # Reqired x*y >= N
    total_sites = int(x*y)

    # Randomize which lattice sites are occupied
    occupied_sites = []
    free_sites = list(range(total_sites))
    
    for i in range(N):
        ran = random.randint(0, total_sites-i-1)
        #print('random int = ' + str(ran))
        occupied_sites.append(free_sites.pop(ran))
        
    return sorted(occupied_sites) # Returns ordered list of N unique numbers between 0 and the total sites - 1 inclusive


####################################################################################


# For each interaction between particles, defined in pairs, constructs the interaction in the format [interaction strength, 1st particle no., 2nd particle no.] in order to create the corresponding term in the Hamiltonian
# sites describes the occupied sites in the system
# interactions is the corresponding interaction matrix (see intmatrix() function) for the system
# the xy variable is a tuple/list of the dimensions of the interaction matrix which can be ignored as this is defult calculated from 'interactions'

def couples(sites, pairs, interactions, xy='_'):
    if xy == '_':    
        int_x, int_y = len(interactions)//2 + 1, len(interactions[0])//2 + 1 # dimensions of interaction matrix
    else:
        int_x, int_y = xy
    # use interaction matrix to get couplings between each occupied site and output this in the correct formatt with particle number
    couples = []
    for i in pairs:
        p1, p2 = sites[i[0]], sites[i[1]]
        
        couples.append([interactions[int_x - 1  + (p1 % int_x) - (p2 % int_x), int_y - 1 + (p1 // int_x) - (p2 // int_x)], i[0], i[1]]) 
    return couples


####################################################################################


# Periodic boundry conditions - interactions 'wrap around' the edges of the system to the opposite sides (thus we would expect stronger interactions)
# For each interaction between particles, defined in pairs, constructs the interaction in the format [interaction strength, 1st particle no., 2nd particle no.] in order to create the corresponding term in the Hamiltonian
# sites describes the occupied sites in the system
# interactions is the corresponding interaction matrix (see intmatrix() function) for the system

def pbc_couples(sites, pairs, interactions, xy='_'): # Constructs interaction strengths for periodic boundry conditions
    if xy == '_':    
        int_x, int_y = len(interactions)//2 + 1, len(interactions[0])//2 + 1 #centre of interaction matrix
    else:
        int_x, int_y = xy
    #use interaction matrix to get couplings between each occupied site and output this in the correct formatt with particle number
    couples = []
    for i in pairs:
        p1, p2 = sites[i[0]], sites[i[1]]
        dx, dy = (p1 % int_x) - (p2 % int_x), (p1 // int_x) - (p2 // int_x)
        
        if dx >= int_x:
            dx = dx - len(interactions)
        elif dx <= -int_x:
            dx = dx + len(interactions)
            
        if dy >= int_y:
            dy = dy - len(interactions[0])
        elif dy <= -int_y:
            dy = dy + len(interactions[0])
        
        couples.append([interactions[int_x - 1  + dx, int_y - 1 + dy], i[0], i[1]]) 
    return couples


####################################################################################


# Finds mirror opertion for some particle configuration
# *Function not used in the rest of the project
# Only relevent if lattice simulated has the symmetry which is tested on basis states
# Occupied sites labelled left to right up to down in basis
# For rectangular simulated systems (but sites occupied can be picked)

def symm_mat(x, y, occupied_sites, basis, axis=1): #defult finds mirror plane transformation which flips x coordinates
    symm = np.zeros((len(basis), len(basis)))
    count1 = 0
    asymmetrics = []
    flips = []
    for i in basis:
        basis_on_lattice = np.zeros((x, y))
        count2 = 0
        for n in i:
            if n == 1:
                basis_on_lattice[occupied_sites[count2] % x, occupied_sites[count2] // x] = 1

            count2 += 1

        flipped = np.flip(basis_on_lattice, axis)
        
        
        if flipped == basis_on_lattice:
            symm[count1, count1] = 1
        else:
            asymmetrics.append(basis_on_lattice)
            flips.append(flipped)
        count1 += 1

    for state in asymmetrics:
        pair = np.where(flips == state)[0]
        symm[pair[0], pair[1]], symm[pair[1], pair[0]] = 1, 1

    return symm


####################################################################################


# Tensor product for 2 strings A and B

def str_tens(A, B):
    new = []
    for i in A:
        for n in B:
            new.append(i + n)
    return new


####################################################################################


# Calculates the probability each particle is observed in the up state for each time evolution (descrbed in populations)
# L is the number of particles in the system

def sparse_up_prob(L, populations):
    basis = spin_basis_general(N=L, pauli=False)
    #col_len = int(math.ceil(L/row_len))

    # Contruct outer product of single particle up spin state operator to find the probability it is in up
    outer_products = []
    for i in range(L):
        factors = [[[0.5, i]], [[-1, i]]]
        opers = [['I', factors[0]], ['z', factors[1]]]
        input_dict = {'up': opers}
        outer_products.append(quantum_operator(input_dict, basis=basis, check_herm=False, check_symm=False, check_pcon=False))
        #__ = quantum_operator(input_dict, basis=basis).tohamiltonian()
        #print(__.toarray())


    # Calculate the probability each particle is spin up for each time evolved
    results = []
    for i in np.transpose(populations):
        probs = []
        for n in range(L):
            probs.append(outer_products[n].expt_value(i))
        results.append(sum(probs)/L)

    #print(results)
    return results


####################################################################################


# Time dependant pulse sequences
# All functions have the same input parameters for easy switches between pulses in further programs

# Define pulse sequences
def pulse_ramsey(t, pulse_time, free_time, rabi_frequency):
    if t <= pulse_time:
        return rabi_frequency
    elif t < pulse_time + free_time:
        return 0 
    elif t <= 2 * pulse_time + free_time:
        return rabi_frequency
    else:
        return 0

def pulse_single(t, pulse_time, free_time, rabi_frequency):
    if t <= pulse_time:
        return rabi_frequency
    else:
        return 0

def no_pulse(t, pulse_time, free_time, rabi_frequency):
    return 0


####################################################################################


# Constructing the Hamitonian - Returns QuSpin Hamitonain object
# N - Number of particles in the system
# J - Interaction strength between nearest neighbors
# det - Detuning term coefficient (Pauli z operator)
# ax, ay - How 'stretched' the lattice of particle sites is in the x and y directions respectively
# Lx - System x side length (Ly=Lx defult for y side)
# params - Parameters for the pulse function - formatted [pulse_time, free_time, rabi_frequency]
# pulse - Time dependant pulse sequence function
# theta, phi - Define the direction of the quantization axis
# occupied_sites - Defines the occupied sites in the lattice, defult RANDOM occupancy of lattice - formatted as a list of integers (size N) indicating the occupied sites labeled in the lattice from left to right top to bottom (e.g. 1 is top left, 2 is the site to its right or the first site on the next row)
# pbc - Calculate interaction terms with periodic boundry conditions, defult FALSE
# int_plot - Plot the interaction strength matrix, defult false

def construct_H(N, J, det, ax, ay, Lx, params, pulse, theta=0, phi=0, Ly='_', int_plot=False, occupied_sites='_', pbc=False):
    basis = spin_basis_general(N=N, pauli=False)
    if Ly == '_':
        Ly = 1

    # Construct interaction
    sites = range(N)
    pairs = [[i, k] for idx, i in enumerate(sites) for k in sites[idx + 1:]]
    
    #interaction strengths between each particle
    interaction_matrix, _, _ = intmatrix(Lx*2 - 1, 1, xs=ax, ys=ay, theta=theta, phi=phi, Ly=Ly*2 - 1, plot=int_plot) 
    #above outputs 3 variables (also set J=1 in this function as it xJ at the end anyway) 

    #determin sites that are occupied (random or not)
    if occupied_sites == '_':
        sites = ran_occupied_sites(N, Lx, Ly)
    else:
        sites = sorted(occupied_sites) #make sure occupied sites is the correct length with numbers under total sites

    #contruct infomation in to reqired format
    if pbc == True: #For periodic boundry conditions
        intr = pbc_couples(sites, pairs, interaction_matrix)
    else:    
        intr = couples(sites, pairs, interaction_matrix)
    
    # Construct resonant and detuning terms
    drive = [[1.0, i] for i in range(N)]
    detuning = [[1.0, i] for i in range(N)]

    operator_list_0 = [['z', detuning]]
    operator_list_1 = [['x', drive]]

    # Construct full Hamiltonian
    interaction = [['+-', intr], ['-+', intr]]
    operator_dict = {'detuning': operator_list_0, 'drive': operator_list_1, 'interaction': interaction}
    H = quantum_operator(operator_dict, basis=basis, check_herm=False, check_symm=False, check_pcon=False)
    params_dict = dict(detuning=det, drive=(pulse, params), interaction=J)
    H_lmbda1 = H.tohamiltonian(params_dict)
    return H_lmbda1


####################################################################################


# Run many time evolutions of a smaller system according to some particle distribution to approximate the dynamics in a large system with the same particle distribution
# 

def ensemble_Ramsey(samples, time, p_max, r_max, initial_state, N, J, det, ax, ay, params, pulse, theta=0, phi=0):
    results = []
    for i in range(samples):

        # Sample the system size to distribute particles to mimic the particle distribution of the larger system being studied
        r = random.randint(0, int(3000*r_max/4))/1000
        p = p_max*(1 - (r/r_max)**2)
        side_length = int(np.sqrt(N/p))+1
        #print(side_length)
        #side_length = 2

        # Build Hamiltonian
        Hamiltonian = construct_H(N, J, det, ax, ay, side_length, params, pulse, theta=theta, phi=phi, Ly=side_length)

        # Evolve system
        populations = Hamiltonian.evolve(initial_state, 0, time)

        # Extract the probability each particle is in the up state at each time evolution
        ups = sparse_up_prob(N, populations)

        # Record the final result
        results.append(ups[-1])
        #print('Sample ' + str(i) + ' completed')
    return sum(results)/samples