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
fx.tic()
print('\tNow reading force constants from '+str(forcefile)+'\n')
fijk, du, idsl, idsr, ids, nl, nr, n, = fx.readFijk(forcefile)
fx.toc()

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

     
    ## Get velocities for each block
    for s in range(1):#split): #loop over blocks for averaging
        fx.tic()
        vels = np.zeros((tn,n*3)) # [time,(1vx,1vy,1vz,2vx,2vy,...)]
        print('\n\tNow reading velocities from block '+str(s+1)+'\n')                
        num = tn/10
        for j in range(tn): 
            if j != 0 and j%(num) == 0: #print progress updates
                print('\t\tNow '
                      +str(10*np.round(j/num,decimals=1))+
                      '% done reading block '+str(s+1))
            for k in range(n*3):
                vels[j,k] = float(fid.readline().strip().split()[0])
        fx.toc()
        
        
        ## Compute the spectral conductance terms
        fx.tic()
        print('\n\tNow computing anharmonic contributions for block '
              +str(s+1)+'\n')
        vfft = np.fft.fft(vels,axis=0)*dt #time scaling, have a link somewhere
        for i in range(1): #loop over left atoms
            
            if i != 0 and i%(nl/10) == 0: #print progress updates
                print('\t\tNow '+str(10*np.round(i/(nl/10),1))+'% done '
                      'comptuting spectral conductance for block '+str(s+1))
                
            vi = vfft[:,idsl[i]*3:idsl[i]*3+3] #vx on left side
            for j in range(1):#nr): #loop over right atoms
                vj = vfft[:,idsr[j]*3:idsr[j]*3+3] #vx on left side
                for k in range(1):#n): #loop over all atoms contributing to two
                    # moderating the two phonon contribution
                    vk = vfft[:,k*3:k*3+3] #vx on left side
                    
                    phix = fijk[i*3,j*3:j*3+3,k*3:k*3+3] #phi on i in x
                    phiy = fijk[i*3+1,j*3:j*3+3,k*3:k*3+3] #phi on i in y
                    phiz = fijk[i*3+2,j*3:j*3+3,k*3:k*3+3] #phi on i in z
#                    
#                    for w in range(tn): #omega
#                        for wp in range(tn): #omega prime

        fx.toc()
            
#qom[j,0] = -(np.dot(lvelsfft[j,:],np.dot(kij,rvelsfft[j,:].transpose().conj()))/om[j])

            
                
            
            
    
    

