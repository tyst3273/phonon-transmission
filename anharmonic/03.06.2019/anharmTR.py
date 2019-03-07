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
import sys
import anharmFX as fx

##### SIMULATION PARAMETERS #####
dT = 20 #NEMD bath temperature difference

dtMD = 0.5e-15 #MD timestep
dn = 2**4 #frequency velocities are printed
dt = dtMD*dn #effective timestep for velocity data
steps = 2**19 #number of simulation steps

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
## Print parameters to screen
fx.printParams(dtMD,dt,dT,steps,split,tn)

## Get force constants from file 
#fx.tic()
#print('\tNow reading force constants from '+str(forcefile)+'\n')
#fijk, du, idsl, idsr, ids, nl, nr, n = fx.readFijk(forcefile)
#fx.toc()

## Generate time and freq. arrays
om, thz, dom = fx.makeTime(dt,tn)

## Read velocites 
with open(velsfile,'r') as fid:
    if int(fid.readline().strip().split()[1]) != n: #error checking
                sys.exit('LAMMPS ERROR: Numer of atoms in vels run doesn\'t'
                         'match forces run')
    if int(fid.readline().strip().split()[1]) != dn: #error checking
        sys.exit('LAMMPS ERROR: Different stride given in '
                 +str(velsfile))
    fid.readline() #skip comment
    for j in range(n): #read the ids for error checking
        if ids[j] != int(fid.readline().strip().split()[0]):
            sys.exit('LAMMPS ERROR: ids in '+str(velsfile)+' don\'t'
                 ' match'+str(forcefile))
    tmp = fid.readline() #skip comment
            
    for i in range(1):#split): #loop over blocks for averaging
        fx.tic()
        vels = np.zeros((tn,n*3)) # [time,(1vx,1vy,1vz,2vx,2vy,...)]
        print('\n\tNow reading velocities from block '+str(i+1)+'\n')                
        num = tn/10
        for j in range(tn): 
            if j != 0 and j%(num) == 0: #print progress updates
                print('\t\tNow '
                      +str(10*np.round(j/num,decimals=1))+
                      '% done reading block'+str(i+1))
            for k in range(n*3):
                vels[j,k] = float(fid.readline().strip().split()[0])
        
        lvf, rvf, vf = fx.computeAnh(vels,fijk,idsl,idsr,dt,tn,nl,nr)
        #fft of left atoms, right atoms, and ALL atoms respectively
        #might get memory heavy but don't want to index in the loop
        
        fx.toc()
            
        

            
                
            
            
    
    

