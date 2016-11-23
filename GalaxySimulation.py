'''
Created on Oct 30, 2016
@author: Bas
'''

from __future__ import division

import scipy.constants as spc
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import BarnesHut as bh

import numpy as np
import time 

def kuzmin_disk(M, r, a):

    pdf = (M / (2 * np.pi * a**2)) * (1 + (r**2/a**2))**2.5
    return pdf

def plummer(M, r, a):
    
    rho = (3 * M * a**2) / (4 * np.pi * (a**2 + r**2)**2.5)
    return rho

def hernquist(M, r, a):
    
    rho = (M * a) / (2 * np.pi * r * (a + r)**3)
    return rho

def pdf_plummer(M, r, a):
    
    mass = (M * r**3) / ((a**2 + r**2)**1.5)
    return mass

def velocity(r, m):
    
    v = np.sqrt((2 * spc.G * m) / r)
    return v

sol_mass = 1.989 * 10**30
mass_gal = 5.8 * 10**9 * sol_mass

r_gal = 5 * 1000 * spc.parsec
a_gal = 3 * 1000 * spc.parsec

# Getting radius from the pdf (acception-rejection method)
points = 1000
complete = 0

r_list = []
m_list = []

while(complete < points):
    
    r = np.random.uniform(0, 1) * r_gal
    m = np.random.uniform(0, 1) * mass_gal
    
    mass = pdf_plummer(mass_gal, r, a_gal)
    
    if(m <= mass):
        
        r_list.append(r)
        m_list.append(m)
        
        complete += 1
        
# Assign positions
x_list = []
y_list = []
z_list = []

theta_list = []
phi_list = []

for i in range(0, len(r_list)):
    
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.random.uniform(0, 1 * np.pi)
    
    x = r_list[i] * np.cos(theta) * np.sin(phi)
    y = r_list[i] * np.sin(theta) * np.sin(phi)
    z = r_list[i] * np.cos(phi)
        
    x_list.append(x)
    y_list.append(y)
    z_list.append(z)

    theta_list.append(theta)
    phi_list.append(phi)

# Get the mass in each shell
shells = 50

radius = np.linspace(0, r_gal, shells)
r_step = radius[1] - radius[0]

mass_list = []

for i in range(0, len(radius) - 1):
    
    mass_per_shell = pdf_plummer(mass_gal, radius[i + 1], a_gal) - pdf_plummer(mass_gal, radius[i], a_gal)
    mass_list.append(mass_per_shell)

# Check how many galaxies there are in each shell
galaxies_per_shell = np.zeros(len(radius) - 1)

for i in range(0, len(radius) - 1):
    for j in range(0, len(r_list)):
    
        if(r_list[j] > radius[i] and r_list[j] < radius[i + 1]):
            galaxies_per_shell[i] += 1

# Assign a mass to each point based on radial distance
mass_points = []

for i in range(0, len(r_list)):
    for j in range(0, len(radius) - 1):
        
        if(r_list[i] > radius[j] and r_list[i] < radius[j + 1]):
            
            mass_points.append(mass_list[j] / galaxies_per_shell[j])
            continue

# Assign veloctities to each particle
vel_x_list = []
vel_y_list = []
vel_z_list = []

for i in range(0, len(r_list)):
    
    mass_in = pdf_plummer(mass_gal, r_list[i], a_gal)
    vel = velocity(r_list[i], mass_in)
    
    vel_x = vel * np.cos(theta_list[i]) * np.sin(2 * phi_list[i])
    vel_y = vel * np.sin(theta_list[i]) * np.sin(2 * phi_list[i])
    vel_z = vel * np.sin(2 * phi_list[i])
    
    vel_x_list.append(vel_x)
    vel_y_list.append(vel_y)
    vel_z_list.append(vel_z)

# Create all the 3D particles
particle_list = []

# Add a super massive black hole
black_hole = bh.particle3D(0, 0, 0, 0, 0, 0, 0, 0.5 * 10**9 * sol_mass, 0, 0)
particle_list.append(black_hole)

for i in range(0, len(r_list)):
    
    particle = bh.particle3D(i + 1, x_list[i], y_list[i], z_list[i], vel_x_list[i], vel_y_list[i], vel_z_list[i], mass_points[i], theta_list[i], phi_list[i])
    particle_list.append(particle)

# Start the simulation
steps = 1000
dt = np.pi * 10**13

for i in range(0, steps):
    
    print(i)
    
    tree = bh.octo_tree(particle_list)
    particle_list = bh.calculate_step3D(tree, dt)
    
    tree.save_plots(i, False)