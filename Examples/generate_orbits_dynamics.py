ALL_RAMSEYS = []
progress = 0
leno = len(orbits)

for i in orbits:
    config = i[-1]
    occ_sites = np.nonzero(config.flatten())[0]
    up_probabilities = []

    Ly, Lx = config.shape
    particles = int(np.sum(config))
    v0 = np.zeros(2**particles)
    v0[-1] = 1
    
    for i in holds:
        #time = [0, pulse_time, pulse_time + i, 2*pulse_time + i]
        time = np.linspace(0, 2*pulse_time + i, 3)
        params = [pulse_time, i, rabi_frequency]
    
        Hamiltonian = construct_H(particles, J, det, ax, ay, Lx, params, pulse_ramsey, theta=theta, phi=phi, Ly=Ly, occupied_sites=occ_sites)
        populations = Hamiltonian.evolve(v0, 0, time)
        ups = sparse_up_prob(particles, particles, 1, populations)
        up_probabilities.append(ups[-1])

    progress += 1
    print(str(progress) + '/' + str(leno) + ' Completed')
    
    ALL_RAMSEYS.append(np.array(up_probabilities))

# Save specific dynamics
with open('Orbit_Ramseys.pickle', 'wb') as handle:
    pickle.dump(ALL_RAMSEYS, handle, protocol=pickle.HIGHEST_PROTOCOL)