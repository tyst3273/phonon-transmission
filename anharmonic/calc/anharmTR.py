#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 4 13:34:46 2019

@author: Ty Sterling
ty.sterling@colorado.edu

This code calculates the anharmonic phonon transmission across interfaces
in crystals. It uses force constants and velocities from MD simulations 
following the method in PHYSICAL REVIEW B 95, 115313 (2017)
"""
import numpy as np
import anharmFX as fx

##### SIMULATION PARAMETERS #####
dT = 20 #NEMD bath temperature difference

dtMD = 0.5e-15 #MD timestep
dn = 2**4 #frequency velocities are printed
dt = dtMD*dn #effective timestep for velocity data
steps = 2**20 #number of simulation steps

split = 2 #number of chunks to split data into for averaging
tn = (steps/dn/split) #number of steps per block of data

win = 0.05 #gaussian smoothing width = THz
kb = 1.38063e-23 #Boltzmann's constant, J/K

forcefile = 'Fijk.dat'
velsfile = 'vels.compact.dat'
outfile = 'anh.tr.TEST.dat'

#####################################################################
#--------------------------------------------------------------------
#####################################################################

print('\n\tUsing MD timestep:\t\t'+str(dtMD*1e12)+'\t\tps')
print('\tEffective timestep:\t\t'+str(dt*1e12)+'\t\tps')
print('\tTemperature bias:\t\t'+str(dT)+'\t\tK')
print('\tTotal Number of Steps:\t\t'+str(steps))
print('\tBlocks for averaging:\t\t'+str(split))
print('\tTime per block:\t\t\t'+str(np.round(tn*dt*1e9,3))+'\t\tns')
print('\tMaximum Frequency:\t\t'+str(0.5/dt*1e-12)+'\t\tTHz')
print('\tFrequency Resolution:\t\t'+str(np.round(1/dt*1e-9/tn,3))+'\t\tMHz\n')
print('\t--------------------------------------------------------\n')

## Get force constants from file 
print('\tNow reading force constants from '+str(forcefile)+'\n')
fijk, du, idsl, idsr, nl, nr, n = fx.readFijk(forcefile)

## Read velocites 
for i in range(split): #loop over blocks for averaging
    print('\n\tNow reading velocities from block '+str(i)+'')
    vels = np.zeros((n,3)) #vx,vy,vz
    
    

