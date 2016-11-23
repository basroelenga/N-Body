'''
Created on Nov 16, 2016

@author: roelenga
'''

import matplotlib.pyplot as plt

import scipy.constants as scp

import BarnesHut as bh
import numpy as np

import pickle as pk

particles = []

# Sun
x = 0
y = 2
z = 4

vx = 0
vy = 0
vz = 0

mass = 1.989 * 10**30

particles.append(bh.particle3D(0, x, y, z, vx, vy, vz, mass))

# Earth
x = scp.astronomical_unit
y = 1
z = 0.2 * scp.astronomical_unit

vx = 0
vy = 29700
vz = 3000

mass = 6 * 10**24

particles.append(bh.particle3D(1, x, y, z, vx, vy, vz, mass))

# Moon
x = scp.astronomical_unit + 380000000
y = 3
z = 0.2 * scp.astronomical_unit

vx = 0
vy = 30700
vz = 3000

mass = 7.3 * 10**22 

particles.append(bh.particle3D(2, x, y, z, vx, vy, vz, mass))

# Jupiter
x = 1.3 * scp.astronomical_unit
y = 2
z = 0.5 * scp.astronomical_unit

vx = 0
vy = 20000
vz = -2000

mass = 6 * 10**28

particles.append(bh.particle3D(3, x, y, z, vx, vy, vz, mass))

# Jupiter 2
x = -1.4 * scp.astronomical_unit
y = 4
z = 0.5 * scp.astronomical_unit

vx = 0
vy = -20000
vz = 2000

mass = 6 * 10**28

particles.append(bh.particle3D(4, x, y, z, vx, vy, vz, mass))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

steps = 4000
dt = 25000

particle_list = []

for i in range(0, steps):

    print(i)
    
    particle_list.append([i * dt,particles])
    
    tree = bh.octo_tree(particles)
    particles = bh.calculate_step3D(tree, dt)

    print(particles[1].getX(), particles[1].getY(), particles[1].getZ(), len(particles))

    tree.save_plots(i, False)
    
# Save the data
catalogdType = {'names': ('t', 'x', 'y', 'z'), 'formats': ('f4', 'f4', 'f4', 'f4')}
catalog = np.recarray((len(particle_list) * len(particles)), dtype = catalogdType)

t = []
x = []
y = []
z = []

for i in range(0, len(particle_list)):
    for j in range(len(particles)):
        
        t.append(particle_list[i][0])
        x.append(particle_list[i][1][j].getX())
        y.append(particle_list[i][1][j].getY())
        z.append(particle_list[i][1][j].getZ())
        
catalog['t'] = t
catalog['x'] = x
catalog['y'] = y
catalog['z'] = z

# Pickle it
obj = open('nbody.pydat', "wb")
pk.dump(catalog, obj)

obj.close()