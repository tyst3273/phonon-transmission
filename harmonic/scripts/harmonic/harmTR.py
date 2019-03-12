#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:36:56 2019

@author: Ty Sterling
ty.sterling@colorado.edu

Another go at harmonic tr, bump hasn't gone away

This is the main function calling from the module
"""

import numpy as np
import sys
import harmFX as fx
import matplotlib.pyplot as plt
fx.tic()

##### SIMULATION PARAMETERS #####
dT = 6 #NEMD bath temperature difference

dtMD = 4.3e-15 #MD timestep
dn = 2**5 #frequency velocities are printed
dt = dtMD*dn #effective timestep for velocity data
steps = 2**22 #number of simulation steps

split = 4 #number of chunks to split data into for averaging
tn = (steps/dn/split) #number of steps per block of data

win = 0.9 #gaussian smoothing width = THz
kb = 1.38063e-23 #Boltzmann's constant, J/K

forcefile = 'Fij.dat'
velsfile = 'vels.compact.dat'
outfile = 'tr.TEST.dat'

conv = 1.602e-19/1e-24 #convert: v*dF/du*v = [A/ps]*[eV/A]/[A]*[A/ps]  
# eV/ps^2; 1 eV = 1.602e-19 J, 1 ps = 1e-12 s

#####################################################################
#--------------------------------------------------------------------
#####################################################################
## Print parameters to screen
fx.printParams(dtMD,dt,dT,steps,split,tn)
## Generate time and freq. arrays
om, thz, dom = fx.makeTime(dt,tn)
## Get force constants from file 
print('\tNow reading force constants from '+str(forcefile)+'')
kij, du, idsl, idsr, ids, nl, nr, n = fx.readFij(forcefile)
fx.toc()

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
    fid.readline() #skip comment
    
    ## Get velocities for each block
    qomRaw = np.zeros((tn,split)).astype(complex)
    for s in range(split): #loop over blocks for averaging
        vels = np.zeros((tn,n*3)) # [time,(1vx,1vy,1vz,2vx,2vy,...)]
        print('\n\tNow reading velocities from block '+str(s+1)+'\n')                
        num = tn/10
        for j in range(tn): 
            if j != 0 and j%(num) == 0: #print progress updates
                print('\t\tNow '+str(10*np.round(j/num,decimals=1))+
                      '% done reading block '+str(s+1))
            for k in range(n*3):
                vels[j,k] = float(fid.readline().strip().split()[0])
        fx.toc()
         
        qomRaw[:,s] = fx.computeQ(kij,vels,om,idsl,idsr,nl,nr,tn,dt,s) 
        fx.toc()
        #save each split into an array
 
qomRaw = -2*conv*np.imag(qomRaw)/(tn*dt) 
qom = qomRaw.mean(axis=1) #avg across splits
tr = qom/kb/dT #transmission
  
## Smoothen transmission
tr = fx.gsmooth(tr,win,dom) #gaussian convolution
## Write to output file
fx.writeTr(outfile,thz,tr,qom,dT,dtMD,tn,win)

print('\n\tALL DONE!\n')
fx.toc()

## Plot the results
fig,ax = plt.subplots()
ax.plot(thz,tr)
plt.xlabel('THz')
plt.ylabel('Transmission')
plt.title('Interfacial Phonon Transmission')
plt.show()




        
        
        
        
    
