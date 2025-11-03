import numpy as np
from itertools import product
import pickle

open_sites = 9
x = np.array([i for i in product(range(2), repeat=open_sites)])

#4-fold symmetry, mirror plane ### format of below array is [symmetry matrix, number of opertations to give the identity matrix]
symmertry_matricies = [[np.zeros((9, 9)), 4], [np.kron(np.eye(3), [[0, 0, 1], [0, 1, 0], [1, 0, 0]]), 1]]
symmertry_matricies[0][0][0, 6], symmertry_matricies[0][0][1, 3], symmertry_matricies[0][0][2, 0], symmertry_matricies[0][0][3, 7], symmertry_matricies[0][0][4, 4], symmertry_matricies[0][0][5, 1], symmertry_matricies[0][0][6, 8], symmertry_matricies[0][0][7, 5], symmertry_matricies[0][0][8, 2] = np.ones((9))
#do this with tensor products?

#create all symmetry matricies which may create unique configurations
unique_symm_mat = []
for i in symmertry_matricies:
    o = i[0]
    for n in range(i[1]):
        unique_symm_mat.append(o)
        o = np.matmul(o, i[0])
rot_mat = np.matmul(np.array(unique_symm_mat[0: 3]), unique_symm_mat[4])
# print(rot_mat)
# print(len(unique_symm_mat))

for i in rot_mat:
    unique_symm_mat.append(i)
        
unique_symm_mat = np.array(unique_symm_mat)
# print(unique_symm_mat.shape)

g = np.inner(unique_symm_mat, x)
l = g
v = np.swapaxes(g, -1, -2)
v = x - v
v = np.swapaxes(v, 0, 1)
# f = 288
# print(x[f])
# print(v[f])

# print(v.shape)

g = np.moveaxis(g, -1, 0)
# print(g.shape)
# print(np.unique(v[f], axis=0))
# print(np.unique(g[f], axis=0))

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
    
#orbits = np.unique(orbits, axis=0)


#print(len(orbits))
unique = 0
for i in orbits:
    unique += len(i)

print('Total permutations found = ' + str(unique))

t = 0
onesandzeros = [0]
for i in orbits: 
    if np.sum(i[0]) == 1:
        print(t)
        onesandzeros.append(t)
        print(i)
    t += 1

for i in sorted(onesandzeros, reverse=True):
    del orbits[i]
    
for i in orbits: 
    if np.sum(i[0]) == 1:
        print(i)

unique = 0
for i in orbits:
    unique += len(i)

print('Total permutations found = ' + str(unique))
print('Unique configurations found = ' + str(len(orbits)))

with open('Orbits.pickle', 'wb') as handle:
    pickle.dump(orbits, handle, protocol=pickle.HIGHEST_PROTOCOL)



#Manually removing duplicate configurations due to translational symmetry or unoccupied rows/columns for a 3x3 grid
trans_symm = [0, 16, 1, 24, 2, 62, 3, 12, 6, 50, 37, 7, 68, 26, 28, 14, 31, 39, 74, 43, 80, 4, 20, 41, 52, 58, 86]
rem_row = [92, 56, 44, 30, 10]
trans_symm, rem_row = sorted(trans_symm), sorted(rem_row)
remove_indicies = trans_symm + rem_row
remove_indicies = sorted(remove_indicies)

# print(trans_symm)
# print(len(trans_symm))
# print(rem_row)
# print(len(rem_row))
# print(remove_indicies)
# print(len(remove_indicies))


with open('Orbits.pickle', 'rb') as handle:
    b = pickle.load(handle)

# print(type(b))
# print(type(b[0]))
# print(type(b[0][0]))
# print(len(b))
# print(b)

g = [[1, 0, 0], [9, 1, 2]]
print(type(np.fliplr(np.rot90(g))))


#Manually add back in smaller configurations which dont use 
symm_2 = [ [[1, 0, 1]], [[1, 1]], [[1, 1, 1]], [[1, 0], [0, 1]], [[1, 1, 1], [1, 1, 1]], [[1, 1], [0, 0], [1, 1]] ]
symm_4 = [ [[0, 1, 0], [1, 0, 1]], [[1, 1], [1, 0]], [[0, 1, 0], [1, 1, 1]], [[1, 1, 1], [1, 0, 1]] ]
symm_4_m = [ [[1, 1], [0, 0], [1, 0]], [[1, 1, 0], [0, 0, 1]], [[1, 1], [0, 1], [0, 1]], [[1, 1, 0], [1, 0, 1]], [[1, 1, 1], [1, 1, 0]]  ]
symm_2_m = [ [[1, 1, 0], [0, 1, 1]], [[1, 0, 0], [0, 0, 1]] ]
print(len(symm_2) + len(symm_4) + len(symm_2_m) + len(symm_4_m))

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

# print(corrected_orbits[-2])

#load previous data
with open('Orbits.pickle', 'rb') as handle:
    b = pickle.load(handle)

#convert to 3x3 matricies
# for u in range(len(b)):
#     b[u] = b[u].reshape((-1, 3, 3))

#remove orbits with at least one row or column full of 0's (as we simplify or combine these)
for i in remove_indicies[::-1]:
    b.pop(i)

#re add removed orbits with the more optimized correct orbits calculated above
for j in corrected_orbits:
    b.append(j)

#pad all kernels
# optimized_orbits = []
# x = 0
# for r in b:
#     y = 0
#     orb = []
#     for l in r:
#         kernel = np.pad(l, 1, 'constant', constant_values=((0, 0), (0, 0)))
#         orb.append(kernel)
#         y += 1
#     optimized_orbits.append(orb)
#     x += 1

# print(b[0])


#save orbits
with open('Orbits.pickle', 'wb') as handle:
    pickle.dump(b, handle, protocol=pickle.HIGHEST_PROTOCOL)


