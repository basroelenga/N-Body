'''
Created on Nov 1, 2016
Purpose: The barnes-hut algorithm
@author: Bas
'''

from __future__ import division
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

import scipy.constants as spc

import numpy as np
import time

# A cell used in creating the tree
class cell2D:
    
    def __init__(self, ID, parentNode, nodeRange, level, objects, xList, yList):
        
        self.ID = ID
        self.PN = parentNode
        
        self.NR = nodeRange
            
        self.level = level
        self.objects = objects
        
        xList.append(xList[0])
        yList.append(yList[0])
        
        self.xList = xList
        self.yList = yList
        
        self.center_of_mass = self.CoM()
        
        # This binds a particle to a cell
        if(len(objects) == 1):
            self.PID = objects[0].getID()
    
    # Function that calculates the center of mass for each cell    
    def CoM(self):
        
        center_of_mass = 0

        if(len(self.objects) > 0):
            
            x = 0
            y = 0

            total_mass = 0
            
            for i in range(0, len(self.objects)):
                
                x += self.objects[i].getX() * self.objects[i].getMass()
                y += self.objects[i].getY() * self.objects[i].getMass()

                total_mass += self.objects[i].getMass()

            x = x / total_mass
            y = y / total_mass

            center_of_mass = [x, y, total_mass]
            
        return center_of_mass
     
    def getID(self):
        return self.ID
        
    def getPID(self):
        return self.PID
        
    def getParentNode(self):
        return self.PN
        
    def getNodeRange(self):
        return self.NR
        
    def getLevel(self):
        return self.level
    
    def getObjects(self):
        return self.objects
    
    def getCoM(self):
        return self.center_of_mass
    
    def getXList(self):
        return self.xList
    
    def getYList(self):
        return self.yList
    
    def plot(self):
        plt.plot(self.xList, self.yList, 'r')

# The main difference between these cells is the in which the boxes are plotted.
# In the 2 dimensional case the cells can be plotted as a 2 dimensional box which requires just one plotting command
# In 3 dimensions there needs to be a specific order in plotting commands the ensure that the 3 dimensional boxes look right.
class cell3D:
    
    def __init__(self, ID, parentNode, nodeRange, level, objects, xList, yList, zList):
        
        self.ID = ID
        self.PN = parentNode
        
        self.NR = nodeRange
            
        self.level = level
        self.objects = objects
        
        self.xList = xList
        self.yList = yList
        self.zList = zList
        
        self.center_of_mass = self.CoM()
        
        # This binds a particle to a cell
        if(len(objects) == 1):
            self.PID = objects[0].getID()
    
    # Function that calculates the center of mass for each cell    
    def CoM(self):
        
        center_of_mass = 0

        if(len(self.objects) > 0):
            
            x = 0
            y = 0
            z = 0

            total_mass = 0
            
            for i in range(0, len(self.objects)):
                
                x += self.objects[i].getX() * self.objects[i].getMass()
                y += self.objects[i].getY() * self.objects[i].getMass()
                z += self.objects[i].getZ() * self.objects[i].getMass()

                total_mass += self.objects[i].getMass()

            x = x / total_mass
            y = y / total_mass
            z = z / total_mass

            center_of_mass = [x, y, z, total_mass]
            
        return center_of_mass
     
    def getID(self):
        return self.ID
        
    def getPID(self):
        return self.PID
        
    def getParentNode(self):
        return self.PN
        
    def getNodeRange(self):
        return self.NR
        
    def getLevel(self):
        return self.level
    
    def getObjects(self):
        return self.objects
    
    def getCoM(self):
        return self.center_of_mass
    
    def getXList(self):
        return self.xList
    
    def getYList(self):
        return self.yList
    
    def getZList(self):
        return self.zList
    
    # The plot function has a specified order in which to plot
    def plot(self):
        
        # Front
        xF = self.xList[0:4]
        xF.append(self.xList[0])
        
        yF = self.yList[0:4]
        yF.append(self.yList[0])
        
        zF = self.zList[0:4]
        zF.append(self.zList[0])
        
        plt.plot(xF, yF, zF, 'r')
        
        # Back
        xF = self.xList[4:8]
        xF.append(self.xList[0])
        
        yF = self.yList[4:8]
        yF.append(self.yList[0])
        
        zF = self.zList[4:8]
        zF.append(self.zList[4])

        plt.plot(xF, yF, zF, 'r')

        # Connecting points between front and back
        
        for i in range(0, 4):
            
            xList = [self.xList[i], self.xList[i + 4]]
            yList = [self.yList[i], self.yList[i + 4]]
            zList = [self.zList[i], self.zList[i + 4]]

            plt.plot(xList, yList, zList, 'r')
        
        
# Particle as an object in 2D
class particle2D:
    
    def __init__(self, ID, x, y, vx, vy, mass):
        
        self.x = x
        self.y = y
        
        self.vx = vx
        self.vy = vy
        
        self.ID = ID
        
        self.mass = mass
        
    def getX(self):
        return self.x
    
    def getY(self):
        return self.y
    
    def getR(self):
        r = np.sqrt(self.x**2 + self.y**2)
        return r
    
    def getVX(self):
        return self.vx
    
    def getVY(self):
        return self.vy

    def getMass(self):
        return self.mass
    
    def setMass(self, mass):
        self.mass = mass
    
    def getID(self):
        return self.ID

# Particle as an object in 3D
class particle3D:
    
    def __init__(self, ID, x, y, z, vx, vy, vz, mass, theta, phi):
        
        self.ID = ID
        
        self.x = x
        self.y = y
        self.z = z
        
        self.vx = vx
        self.vy = vy
        self.vz = vz
        
        self.mass = mass

        self.theta = theta
        self.phi = phi

    def getID(self):
        return self.ID

    def getX(self):
        return self.x

    def getY(self):
        return self.y
    
    def getZ(self):
        return self.z

    def getVX(self):
        return self.vx

    def getVY(self):
        return self.vy

    def getVZ(self):
        return self.vz

    def getMass(self):
        return self.mass

    def setMass(self, mass):
        self.mass = mass

    def getPhi(self):
        return self.phi
    
    def getTheta(self):
        return self.theta

class quad_tree:
    
    # The pointlist should contain the particles
    def __init__(self, pointList):
        
        time_start = time.time()
        
        self.pointList = pointList
        self.cellList = self.createTree(self.pointList)
    
        self.tree_calc_time = time.time() - time_start
        print('Computing time (tree): ' + str(self.tree_calc_time))
    
    # This function creates a list of cells created by a tree algorithm.
    # Each cell contains a center of mass
    def createTree(self):
    
        # Check for right dimensions
        if(isinstance(self.pointList[0], particle2D) != True):
            print('Incorrect type, should be of particle2D type.')
            return

        # Create the tree
        level = 0
        completeCells = 0
        
        cellList = []
        ID = 0
        
        # This loop will run until there is only 1 object per cell
        while(True):
 
            # Exit condition
            if(completeCells >= len(self.pointList)): break

            # Initial condition, clockwise rotation
            if(level == 0):
                
                self.xList = []
                self.yList = []
                
                for i in range(0, len(self.pointList)):
                    self.xList.append(self.pointList[i].getX())
                    self.yList.append(self.pointList[i].getY())
                
                # A correction factor is added to prevent particles from being outside of the box (on the box)
                corr_x = np.abs(np.max(self.xList) + np.min(self.xList)) / 100
                corr_y = np.abs(np.max(self.yList) + np.min(self.yList)) / 100
                
                x_max = np.max(self.xList) + corr_x
                y_max = np.max(self.yList) + corr_y

                x_min = np.min(self.xList) - corr_x
                y_min = np.min(self.yList) - corr_y
                
                # The root node is its own parent node
                cellList.append(cell2D(ID, ID, [0, 0, 0, 0], level, self.pointList, [x_max, x_max, x_min, x_min], [y_max, y_min, y_min, y_max]))
                ID += 1
                
            for i in range(0, len(cellList)):
                if(cellList[i].getLevel() == level - 1 and len(cellList[i].getObjects()) > 1):
                    
                    x_min = np.min(cellList[i].getXList())
                    y_min = np.min(cellList[i].getYList())
                    
                    x_max = np.max(cellList[i].getXList())
                    y_max = np.max(cellList[i].getYList())
                    
                    x_dist_step = (x_max - x_min) / 2
                    y_dist_step = (y_max - y_min) / 2
        
                    obj_in_cell = cellList[i].getObjects()
                    cell_ID_range = np.linspace(ID, ID + 3, 4)
        
                    for j in range(0, 2):
                        for k in range(0, 2):
                        
                            x = [x_min + (j + 1) * x_dist_step, x_min + (j + 1) * x_dist_step, x_min + j * x_dist_step, x_min + j * x_dist_step]
                            y = [y_min + (k + 1) * y_dist_step, y_min + k * y_dist_step, y_min + k * y_dist_step, y_min + (k + 1) * y_dist_step]
        
                            obj_in_new_cell = []
        
                            for m in range(0, len(obj_in_cell)):
                                
                                # Check if the object is in the newly defined cell
                                if(obj_in_cell[m].getX() > np.min(x) and obj_in_cell[m].getX() < np.max(x) and obj_in_cell[m].getY() > np.min(y) and obj_in_cell[m].getY() < np.max(y)):
                                    obj_in_new_cell.append(obj_in_cell[m])
                            
                            if(len(obj_in_new_cell) == 1):
                                completeCells += 1
                                
                            cellList.append(cell2D(ID, cellList[i].getID(), cell_ID_range, level, obj_in_new_cell, x, y))
                            ID += 1
                  
            level += 1
            
        print('Tree created')
        return cellList
    
    def plot_cells(self):
        
        for i in range(0, len(self.cellList)):
            self.cellList[i].plot()
    
    def plot_points(self):
        plt.plot(self.xList, self.yList, 'b.')
    
    def save_plots(self, i, boxes):
        
        fig = plt.figure()
        
        plt.plot(self.xList, self.yList, 'b.')
        
        if(boxes):
            for j in range(0, len(self.cellList)):
                self.cellList[j].plot()
        
        plt.xlim(-2 * 10**4 * spc.parsec, 2 * 10**4 * spc.parsec)
        plt.ylim(-2 * 10**4 * spc.parsec, 2 * 10**4 * spc.parsec)

        if(i < 10 and i>= 0): plt.savefig('plots/000' + str(i) + '.png')
        if(i < 100 and i >= 10): plt.savefig('plots/00' + str(i) + '.png')
        if(i < 1000 and i >= 100): plt.savefig('plots/0' + str(i) + '.png')
        if(i < 10000 and i >= 1000): plt.savefig('plots/' + str(i) + '.png')
    
        fig.clf()
    
    def getPoints(self):
        return self.pointList

    def getCellList(self):
        return self.cellList

    def getCalcTime(self):
        return self.tree_calc_time

# When looking at three dimension every box can be subdivided into 8 boxes
class octo_tree:
    
    def __init__(self, pointList):
        
        time_start = time.time()
        
        self.pointList = pointList
        self.cellList = self.createTree()
    
        self.tree_calc_time = time.time() - time_start
        print('Computing time (tree): ' + str(self.tree_calc_time))
        
    def createTree(self):
        
        # Check for right dimensions
        if(isinstance(self.pointList[0], particle3D) != True):
            print('Incorrect type, should be of particle3D type.')
            return
        
        # Create the tree
        level = 0
        completeCells = 0
        
        cellList = []
        ID = 0
        
        while(True):
       
            # Exit condition
            if(completeCells >= len(self.pointList)): break

            # Initial condition, clockwise rotation
            if(level == 0):
                
                self.xList = []
                self.yList = []
                self.zList = []
                
                self.phiList = []
                self.thetaList = []
                
                for i in range(0, len(self.pointList)):
                    
                    self.xList.append(self.pointList[i].getX())
                    self.yList.append(self.pointList[i].getY())
                    self.zList.append(self.pointList[i].getZ())
                    
                    self.phiList.append(self.pointList[i].getPhi())
                    self.thetaList.append(self.pointList[i].getTheta())
                
                # A correction factor is added to prevent particles from being outside of the box (on the box)
                corr_x = (np.abs(np.max(self.xList) + np.min(self.xList))) / 100
                corr_y = (np.abs(np.max(self.yList) + np.min(self.yList))) / 100
                corr_z = (np.abs(np.max(self.zList) + np.min(self.zList))) / 100
                
                x_max = np.max(self.xList) + corr_x
                y_max = np.max(self.yList) + corr_y
                z_max = np.max(self.zList) + corr_z 

                x_min = np.min(self.xList) - corr_x
                y_min = np.min(self.yList) - corr_y
                z_min = np.min(self.zList) - corr_z
        
                # It is a bit harder to create a 3 dimensional box than a 2 dimensional one (4 ---> 8 points)
                xBox = [x_max, x_max, x_min, x_min, x_max, x_max, x_min, x_min]
                yBox = [y_max, y_min, y_min, y_max, y_max, y_min, y_min, y_max]
                zBox = [z_min, z_min, z_min, z_min, z_max, z_max, z_max, z_max]
                
                # The root node is its own parent node
                cellList.append(cell3D(ID, ID, [0, 0, 0, 0, 0, 0, 0, 0], level, self.pointList, xBox, yBox, zBox))
                
                ID += 1
        
            for i in range(0, len(cellList)):
                if(cellList[i].getLevel() == level - 1 and len(cellList[i].getObjects()) > 1):
                    
                    x_min = np.min(cellList[i].getXList())
                    y_min = np.min(cellList[i].getYList())
                    z_min = np.min(cellList[i].getZList())
                    
                    x_max = np.max(cellList[i].getXList())
                    y_max = np.max(cellList[i].getYList())
                    z_max = np.max(cellList[i].getZList())
                    
                    x_dist_step = (x_max - x_min) / 2
                    y_dist_step = (y_max - y_min) / 2
                    z_dist_step = (z_max - z_min) / 2
        
                    # There are now 8 cells in each parent cell
                    obj_in_cell = cellList[i].getObjects()
                    cell_ID_range = np.linspace(ID, ID + 7, 8)
                    
                    # Divide the cells
                    for j in range(0, 2):
                        for k in range(0, 2):
                            for m in range(0, 2):
                                
                                # Generate the new 8 coordinates for each box
                                x = [x_min + (j + 1) * x_dist_step, x_min + (j + 1) * x_dist_step, x_min + j * x_dist_step, x_min + j * x_dist_step, x_min + (j + 1) * x_dist_step, x_min + (j + 1) * x_dist_step, x_min + j * x_dist_step, x_min + j * x_dist_step]
                                y = [y_min + (k + 1) * y_dist_step, y_min + k * y_dist_step, y_min + k * y_dist_step, y_min + (k + 1) * y_dist_step, y_min + (k + 1) * y_dist_step, y_min + k * y_dist_step, y_min + k * y_dist_step, y_min + (k + 1) * y_dist_step]
                                z = [z_min + m * z_dist_step, z_min + m * z_dist_step, z_min + m * z_dist_step, z_min + m * z_dist_step, z_min + (m + 1) * z_dist_step, z_min + (m + 1) * z_dist_step, z_min + (m + 1) * z_dist_step, z_min + (m + 1) * z_dist_step]
            
                                obj_in_new_cell = []
            
                                for l in range(0, len(obj_in_cell)):
                                    
                                    # Check if the object is in the newly defined cell
                                    if(obj_in_cell[l].getX() > np.min(x) and obj_in_cell[l].getX() < np.max(x) and obj_in_cell[l].getY() > np.min(y) and obj_in_cell[l].getY() < np.max(y) and obj_in_cell[l].getZ() > np.min(z) and obj_in_cell[l].getZ() < np.max(z)):
                                        obj_in_new_cell.append(obj_in_cell[l])
                                
                                if(len(obj_in_new_cell) == 1):
                                    completeCells += 1
                                    
                                cellList.append(cell3D(ID, cellList[i].getID(), cell_ID_range, level, obj_in_new_cell, x, y, z))
                                ID += 1
                                
            level += 1                    
            
        print('Tree created')   
        return cellList

    def plot_cells(self):
        
        for i in range(0, len(self.cellList)):
            self.cellList[i].plot()
    
    def plot_points(self):
        
        plt.plot(self.xList, self.yList, self.zList, 'b.')
    
    def save_plots(self, i, boxes):
        
        fig = plt.figure(figsize=(9, 4.5))
        ax = fig.add_subplot(121, projection='3d')
        
        bx = []
        by = []
        bz = []
        
        rx = []
        ry = []
        rz = []
        
        # Different colors check
        for j in range(0, len(self.phiList)):
            
            if(self.phiList[j] < 1/4 * np.pi or self.phiList[j] > 3/4 * np.pi):
                
                bx.append(self.xList[j])
                by.append(self.yList[j])
                bz.append(self.zList[j])
            
            else:
                
                rx.append(self.xList[j])
                ry.append(self.yList[j])
                rz.append(self.zList[j])
        
        ax.plot(bx, by, bz, 'b.')
        ax.plot(rx, ry, rz, 'r.')
        
        if(boxes):
            for j in range(0, len(self.cellList)):
                self.cellList[j].plot()
        
        ax2 = fig.add_subplot(122)
        plt.plot(self.xList, self.yList, 'b.')
        
        limit = 10 * 1000 * spc.parsec
        
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_zlim(-limit, limit)

        ax.view_init(elev=0, azim=0)

        ax.set_xlabel('x - axis')
        ax.set_ylabel('y - axis')
        ax.set_zlabel('z - axis')

        ax2.set_xlim(-limit, limit)
        ax2.set_ylim(-limit, limit)

        if(i < 10 and i>= 0): plt.savefig('plots/000' + str(i) + '.png', dpi=200)
        if(i < 100 and i >= 10): plt.savefig('plots/00' + str(i) + '.png', dpi=200)
        if(i < 1000 and i >= 100): plt.savefig('plots/0' + str(i) + '.png', dpi=200)
        if(i < 10000 and i >= 1000): plt.savefig('plots/' + str(i) + '.png', dpi=200)
    
        fig.clf()
        plt.close()
    
    def getPoints(self):
        return self.pointList

    def getCellList(self):
        return self.cellList

    def getCalcTime(self):
        return self.tree_calc_time


def getAcceleration2D(part1, cell_CoM):

    r = ((part1.getX() - cell_CoM[0])**2 + (part1.getY() - cell_CoM[1])**2)**0.5
    acc = cell_CoM[2] * spc.G / r**2

    accX = acc * (cell_CoM[0] - part1.getX()) / r
    accY = acc * (cell_CoM[1] - part1.getY()) / r
    
    return accX, accY

def getVelocity2D(accX, accY, velX, velY, dt):

    velX += accX * dt
    velY += accY * dt
    
    return velX, velY

def getPosition2D(velX, velY, posX, posY, dt):

    posX += velX * dt
    posY += velY * dt
    
    return posX, posY

# Calculate a time step in the simulation
def calculate_step2D(tree, dt):
    
    time_start = time.time()
    
    # Check if tree is a quad_tree instance
    if(isinstance(tree, quad_tree) == False):
        print('Incorrect type, should be of quad_tree type.')
        return
    
    cell_smallest = []
    
    # Keep only the cells which contain 1 particle
    for i in range(0, len(tree.getCellList())):
        if(len(tree.getCellList()[i].getObjects()) == 1):
            cell_smallest.append(tree.getCellList()[i])
    
    # List for storing all the new particles
    particle_list = []    

    for i in range(0, len(cell_smallest)):
        
        acc_x = 0
        acc_y = 0
        
        current_cell = cell_smallest[i]
        
        # This ups the level
        for j in range(0, current_cell.getLevel()):
            
            # Calculate the gravitational acceleration from the remaining nodes
            for k in range(0, 4):
            
                # Skip own node
                if(current_cell.getNodeRange()[k] != current_cell.getID()):
                    
                    CoM = tree.getCellList()[int(current_cell.getNodeRange()[k])].getCoM()
                    
                    # If the CoM is 0 it means that there are no particles in the cell, thus it can be skipped
                    if(CoM == 0): continue
                    
                    acc = getAcceleration2D(cell_smallest[i].getObjects()[0], CoM)

                    acc_x += acc[0]
                    acc_y += acc[1]
                    
            # Now the current node needs to be set to the parent node
            current_cell = tree.getCellList()[current_cell.getParentNode()]
        
        # Update velocity and position
        vel = getVelocity2D(acc_x, acc_y, cell_smallest[i].getObjects()[0].getVX(), cell_smallest[i].getObjects()[0].getVY(), dt)
        pos = getPosition2D(vel[0], vel[1], cell_smallest[i].getObjects()[0].getX(), cell_smallest[i].getObjects()[0].getY(), dt)
        
        # Update the particles
        particle_new = particle2D(cell_smallest[i].getObjects()[0].getID(), pos[0], pos[1], vel[0], vel[1], cell_smallest[i].getObjects()[0].getMass())
        particle_list.append(particle_new)
    
    calc_time = time.time() - time_start
    print('Computing time (step): ' + str(calc_time))
     
    # Return list with all the particles 
    return particle_list

def getAcceleration3D(part1, cell_CoM):

    r = np.sqrt((part1.getX() - cell_CoM[0])**2 + (part1.getY() - cell_CoM[1])**2 + (part1.getZ() - cell_CoM[2])**2)
    acc = cell_CoM[3] * spc.G / r**2

    accX = acc * (cell_CoM[0] - part1.getX()) / r
    accY = acc * (cell_CoM[1] - part1.getY()) / r
    accZ = acc * (cell_CoM[2] - part1.getZ()) / r
    
    return accX, accY, accZ

def getVelocity3D(accX, accY, accZ, velX, velY, velZ, dt):

    velX += accX * dt
    velY += accY * dt
    velZ += accZ * dt
    
    return velX, velY, velZ

def getPosition3D(velX, velY, velZ, posX, posY, posZ, dt):

    posX += velX * dt
    posY += velY * dt
    posZ += velZ * dt
    
    return posX, posY, posZ

# Calculate a time step in the simulation
def calculate_step3D(tree, dt):
    
    time_start = time.time()
    
    # Check if tree is a quad_tree instance
    if(isinstance(tree, octo_tree) == False):
        print('Incorrect type, should be of quad_tree type.')
        return
    
    cell_smallest = []
    
    # Keep only the cells which contain 1 particle
    for i in range(0, len(tree.getCellList())):
        if(len(tree.getCellList()[i].getObjects()) == 1):
            cell_smallest.append(tree.getCellList()[i])
    
    # List for storing all the new particles
    particle_list = []    

    for i in range(0, len(cell_smallest)):
        
        acc_x = 0
        acc_y = 0
        acc_z = 0
        
        current_cell = cell_smallest[i]
        
        # This ups the level
        for j in range(0, current_cell.getLevel()):
            
            # Calculate the gravitational acceleration from the remaining nodes
            for k in range(0, 8):
            
                # Skip own node
                if(current_cell.getNodeRange()[k] != current_cell.getID()):
                    
                    CoM = tree.getCellList()[int(current_cell.getNodeRange()[k])].getCoM()
                    
                    #if(cell_smallest[i].getObjects()[0].getID() == 1): print(CoM)
                    
                    # If the CoM is 0 it means that there are no particles in the cell, thus it can be skipped
                    if(CoM == 0): continue
                    
                    acc = getAcceleration3D(cell_smallest[i].getObjects()[0], CoM)

                    acc_x += acc[0]
                    acc_y += acc[1]
                    acc_z += acc[2]
                    
            # Now the current node needs to be set to the parent node
            current_cell = tree.getCellList()[current_cell.getParentNode()]
        
        # Update velocity and position
        vel = getVelocity3D(acc_x, acc_y, acc_z, cell_smallest[i].getObjects()[0].getVX(), cell_smallest[i].getObjects()[0].getVY(), cell_smallest[i].getObjects()[0].getVZ(), dt)
        pos = getPosition3D(vel[0], vel[1], vel[2], cell_smallest[i].getObjects()[0].getX(), cell_smallest[i].getObjects()[0].getY(), cell_smallest[i].getObjects()[0].getZ(), dt)

        phi = np.arccos(pos[2] / np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2))
        theta = np.arctan2(pos[1], pos[0])

        # Update the particles
        particle_new = particle3D(cell_smallest[i].getObjects()[0].getID(), pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], cell_smallest[i].getObjects()[0].getMass(), theta, phi)
        particle_list.append(particle_new)
    
    calc_time = time.time() - time_start
    print('Computing time (step): ' + str(calc_time))

    # Return list with all the particles 
    return particle_list