import numpy as np
import random

# Return a radial distribution function with some max density in the center which extends radially based on r_max

def generate_radial_distribution(lattice_filling_max, r_max):
    def generator_dist(r):
        return lattice_filling_max*(1 - (r/r_max)**2)
    return generator_dist


####################################################################################


# Return a uniform distribution function with some desity

def generate_uniform_distribution(lattice_filling):
    def uniform_dist(r):
        return lattice_filling
    return uniform_dist


####################################################################################


# 

def prob_mat(function, lattice_dimensions, samples=1, ax=1, ay=1, polar=True):
    
    dist = np.zeros(lattice_dimensions)
    cen_x, cen_y = lattice_dimensions[0]/2, lattice_dimensions[1]/2
    
    cy = -cen_y
    for i in range(lattice_dimensions[1]):
        cx = -cen_x
        for n in range(lattice_dimensions[0]):
            sam_xpoints, sam_ypoints = np.linspace(cx - ax/2, cx + ax/2, samples + 1), np.linspace(cy - ay/2, cy + ay/2, samples + 1)

            if polar == True:
                values = []
                for u in sam_xpoints:
                    sam_radii = np.sqrt(u**2 + sam_ypoints**2)
                    values.append(function(sam_radii))
            ###################### next 'else' code is untested #####################
            else:
                values = []
                for u in sam_xpoints:
                    values.append(function(u, sam_ypoints))

            segments = 0
            for l in range(samples**2):
                avg_func = sum([values[ l % samples][ l // samples], values[ l % samples + 1][ l // samples], values[ l % samples][ l // samples + 1], values[ l % samples + 1][ l // samples + 1]])/4 
                segments += avg_func * (ax/samples) * (ay/samples)
            
            dist[n, i] = segments

            cx += ax
        cy += ay

    return dist


####################################################################################


# Generates sample uniform lattice for N particles with a lattice filling/desity specified

def unif_random_config(filling, N):
    side_length = int(np.sqrt(N/filling)) + 1
    positions = random.sample(range(1, side_length**2), N)
    config = np.zeros((side_length, side_length))
    for i in positions:
        config[i // side_length, i % side_length] = 1 
    return config


####################################################################################


# Measures the average and standard deviation of the frequency the patterns/particle configurations in config appear in a N particle uniform distibution with some latticle filling
# Samples the distribution then searches and counts the number of specific patterns within, increasing samples increases accuray and reduces varience of result
# unif_random_config could be repalaced to apply this method to other more complex distribtuions
    
def config_freqs(configs, filling, N, samples):
    frequencies = []
    for _ in range(samples):
        sample_config = unif_random_config(filling, N)
        sample_config = np.pad(sample_config, ((1, 1), (1, 1)), 'constant', constant_values=((0, 0), (0, 0)))

        instance = []
        for m in configs:
            found = 0
            particle_no = m[0].sum()
            for c in m:
                found += count_geometry(sample_config, c, particle_no=particle_no)
            instance.append(particle_no*found/N)
            #print(instance[-1])
        frequencies.append(instance)
    # returns list of averages and list of standard deviations
    return np.mean(frequencies, axis=0), np.std(frequencies, axis=0)

def count_geometry(image, kernel, particle_no='_'):
    kernel = np.array(kernel)
    m, n = kernel.shape
    if particle_no == '_':
        particle_no = m*n
        #print(particle_no)
    y, x = image.shape
    y = y - m + 1
    x = x - n + 1
    
    found = 0
    for i in range(y):
        for j in range(x):
            if (image[i:i+m, j:j+n] == kernel).all():
                found += 1
    return found

def count_geometry(image, kernel, particle_no='_'):
    #if particle_no == '_' :
        #particle_no = np.sum(kernel)
    kernel = np.select([kernel == 0], [-1], kernel)
    kernel = np.pad(kernel, ((1, 1), (1, 1)), 'constant', constant_values=((-1, -1), (-1, -1)))
    #print(kernel)
    g = signal.convolve2d(image, kernel)
    return np.count_nonzero(g == particle_no)

def find_geometry(image, kernel, particle_no='_'):
    #if particle_no == '_' :
        #particle_no = np.sum(kernel)
    kernel = np.select([kernel == 0], [-1], kernel)
    kernel = np.pad(kernel, ((1, 1), (1, 1)), 'constant', constant_values=((-1, -1), (-1, -1)))
    #print(kernel)
    g = signal.convolve2d(image, kernel, mode='same')
    g = np.select([g != particle_no], [0], g)
    g = np.select([g == particle_no], [1], g)
    return g

def rotated(array_2d):
    list_of_tuples = zip(*array_2d[::-1])
    return [list(elem) for elem in list_of_tuples]

def pattern_search(image, kernel, particle_no = '_', symmetrical=False):
    rot = 4
    if symmetrical == True:
        rot = 2
    patterns_found = 0
    for i in range(rot):
        patterns_found += count_geometry(image, kernel, particle_no=particle_no)
        kernel = rotated(kernel)
    return patterns_found

def random_config(discrete_dist):
    y, x = discrete_dist.shape
    config = np.zeros((y, x))
    for i in range(y):
        for n in range(x):
            roll = random.randint(0, 1000)/1000
            if roll <= discrete_dist[i, n]:
                config[i, n] = 1
    return config