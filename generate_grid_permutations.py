import numpy as np
from itertools import product
import pickle

# This file GENERATES and SAVES its result.
# The output is already saved in "Orbits.pickle" so generally this file does not need to be run again.

####################################################################################


# In this file we aim to find the RELAVENT smaller particle configurations to simulate the quantum system and format the result such that our program can use it.
# This code ONLY works for finding all relevent configuratiuons in a 3x3 grid, however the method may be adapted for other limits (e.g. configurations in a 4x4 grid, 3x2 grid etc.) with considerations of their symmetry
# This means the output repersents all the dynamics which are possible within a SQUARE 3x3 grid.

# We remove configurations which are DYNAMICALLY EQUIVALENT to others from symmetry (e.g. rotations and reflections) then save the output into a numpy list of 3x3 matricies corresponding to each position in the grid.
# Rotational and mirror symmetries are found and handled by code, however translational symmetries and configurations which are smaller then a 3x3 grid are manually found/removed then added in the appropiate formatt.


####################################################################################


# Generate all particle configurations within a 3x3 grid, repersenting each as a 9 element string of 0's and 1's corresponding to each sites particle occupancy

open_sites = 9
x = np.array([i for i in product(range(2), repeat=open_sites)])


####################################################################################


# Build 3x3 symmetry matricies to identify the rotational and mirror symmetries between configurations

#4-fold symmetry, mirror plane ### format of below array is [symmetry matrix, number of opertations to give the identity matrix]
symmertry_matricies = [[np.zeros((9, 9)), 4], [np.kron(np.eye(3), [[0, 0, 1], [0, 1, 0], [1, 0, 0]]), 1]]
symmertry_matricies[0][0][0, 6], symmertry_matricies[0][0][1, 3], symmertry_matricies[0][0][2, 0], symmertry_matricies[0][0][3, 7], symmertry_matricies[0][0][4, 4], symmertry_matricies[0][0][5, 1], symmertry_matricies[0][0][6, 8], symmertry_matricies[0][0][7, 5], symmertry_matricies[0][0][8, 2] = np.ones((9))
#do this with tensor products?

# Create all symmetry matricies which may create unique configurations
unique_symm_mat = []
for i in symmertry_matricies:
    o = i[0]
    for n in range(i[1]):
        unique_symm_mat.append(o)
        o = np.matmul(o, i[0])
rot_mat = np.matmul(np.array(unique_symm_mat[0: 3]), unique_symm_mat[4])

for i in rot_mat:
    unique_symm_mat.append(i)
        
unique_symm_mat = np.array(unique_symm_mat)


####################################################################################


# Apply EVERY symmetry matrix to EVERY configuration (currently repersented by a 9D vector).
# Produces an array where every configuration has every other configuration connected by rotational and mirror symmetry grouped together (i.e. every configurations ORBIT)

g = np.inner(unique_symm_mat, x)
l = g
v = np.swapaxes(g, -1, -2)
v = x - v
v = np.swapaxes(v, 0, 1)

g = np.moveaxis(g, -1, 0)


####################################################################################


# Remove all duplicate orbits (configurations connected by dyanmical symmetry) and use a single configuration from each orbit to REPERSENT the dynamics of the entire orbit
# Remove all single particle configurations (repersentative added back later)

count = 0
orbits = []
for config in v:
    checked = 0
    for symm_map in config:
        non_zero = symm_map[next((i for i, x in enumerate(symm_map) if x), 0)]
        if non_zero < 0:
            break
        elif non_zero > 0:
            checked += 1
            continue
        checked += 1
    if checked == open_sites-1:
        orbits.append(np.unique(g[count], axis=0))
    count += 1
    
# Find all configurations with only a single particle then delete them
t = 0
onesandzeros = [0]
for i in orbits: 
    if np.sum(i[0]) == 1:
        onesandzeros.append(t)
    t += 1
for i in sorted(onesandzeros, reverse=True):
    del orbits[i]


####################################################################################


# Listing indicies for manually removing duplicate configurations due to translational symmetry or unoccupied rows/columns for a 3x3 grid

# Manually found indicies of relavent repersentatives
trans_symm = [0, 16, 1, 24, 2, 62, 3, 12, 6, 50, 37, 7, 68, 26, 28, 14, 31, 39, 74, 43, 80, 4, 20, 41, 52, 58, 86]
rem_row = [92, 56, 44, 30, 10]
trans_symm, rem_row = sorted(trans_symm), sorted(rem_row)
remove_indicies = trans_symm + rem_row
remove_indicies = sorted(remove_indicies)


####################################################################################


# Listing and formatting configurations which are smaller then 3x3 to manually add back into our resulting repersentative list

symm_2 = [ [[1, 0, 1]], [[1, 1]], [[1, 1, 1]], [[1, 0], [0, 1]], [[1, 1, 1], [1, 1, 1]], [[1, 1], [0, 0], [1, 1]] ]
symm_4 = [ [[0, 1, 0], [1, 0, 1]], [[1, 1], [1, 0]], [[0, 1, 0], [1, 1, 1]], [[1, 1, 1], [1, 0, 1]] ]
symm_4_m = [ [[1, 1], [0, 0], [1, 0]], [[1, 1, 0], [0, 0, 1]], [[1, 1], [0, 1], [0, 1]], [[1, 1, 0], [1, 0, 1]], [[1, 1, 1], [1, 1, 0]]  ]
symm_2_m = [ [[1, 1, 0], [0, 1, 1]], [[1, 0, 0], [0, 0, 1]] ]

for i in range(len(symm_2)):
    new = []
    new.append(np.array(symm_2[i]))
    new.append(np.rot90(symm_2[i]))
    symm_2[i] = new

for i in range(len(symm_4)):
    new = []
    transformed = np.array(symm_4[i])
    for u in range(4):
        new.append(transformed)
        transformed = np.rot90(transformed)
    symm_4[i] = new

for i in range(len(symm_4_m)):
    new = []
    transformed = np.array(symm_4_m[i])
    for s in range(2):
        for u in range(4):
            new.append(transformed)
            transformed = np.rot90(transformed)
        transformed = np.fliplr(symm_4_m[i])
    symm_4_m[i] = new

for i in range(len(symm_2_m)):
    new = []
    transformed = np.array(symm_2_m[i])
    for s in range(2):
        for u in range(2):
            new.append(transformed)
            transformed = np.rot90(transformed)
        transformed = np.fliplr(symm_2_m[i])
    symm_2_m[i] = new

corrected_orbits = []

new_orbits = [ symm_2, symm_4, symm_4_m, symm_2_m, [[np.array([[1, 1], [1, 1]])]], [[np.array([[1]])]] ]
for w in new_orbits:
    for p in w:
        corrected_orbits.append(p)


####################################################################################


# Manage smaller configurations into our repersentatives and format the results for further use

# Convert to 3x3 matricies
for u in range(len(orbits)):
    orbits[u] = orbits[u].reshape((-1, 3, 3))

# Remove orbits with at least one row or column full of 0's (as we simplify or combine these)
for i in remove_indicies[::-1]:
    orbits.pop(i)

# Re add removed orbits with the more optimized correct orbits calculated above
for j in corrected_orbits:
    orbits.append(j)

# Pad all configurations with a single layer of unoccupied sites (0's)

padded_orbits = []
x = 0
for r in orbits[:-1]: # Treat final single particle configuration seperately below
    y = 0
    orb = []
    for l in r:
        kernel = np.pad(l, 1, 'constant', constant_values=((0, 0), (0, 0)))
        orb.append(kernel)
        y += 1
    padded_orbits.append(orb)
    x += 1

# Pad the single particle configuration twice
# Isolates the particle such that there is no overlap with other configurations/orbits

isolated_particle = np.pad(orbits[-1][0], 2, 'constant', constant_values=((0, 0), (0, 0)))
padded_orbits.append([isolated_particle])

####################################################################################


# Save the orbit repersentatives

with open('Orbits.pickle', 'wb') as handle:
    pickle.dump(padded_orbits, handle, protocol=pickle.HIGHEST_PROTOCOL)


####################################################################################
####################################################################################
####################################################################################