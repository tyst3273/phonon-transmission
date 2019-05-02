#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 10:21:12 2019

@author: ty
"""

import numpy as np
import copy as cp
sig = 3.4
a = 1.5496*sig/2.0 #scaling factor for integer basis array

nx1 = int(raw_input('Enter number of unit cells of Ar in x:\n\t>'))
nx2 = int(raw_input('Enter number of unit cells of heavy-Argon in x:\n\t>'))
ny = int(raw_input('Enter width in y:\n\t>'))
nz = int(raw_input('Enter width in z:\n\t>'))
nx = nx1+nx2
lammps = str(raw_input('Enter the name of the lammps data file or "none" if '
                       'you dont want one\n\t>'))
xyz = str(raw_input('Enter the name of the xyz data file or "none" if '
                       'you dont want one\n\t>'))

basis = np.array([[0,0,0],
                  [0,1,1],
                  [1,0,1],
                  [1,1,0]])

pos = cp.deepcopy(basis)
tmp = cp.deepcopy(basis)
for i in range(nx-1):
    tmp[:,0] = tmp[:,0]+2
    pos = np.append(pos,tmp,axis=0)
    
tmp = cp.deepcopy(pos)
for i in range(ny-1):
    tmp[:,1] = tmp[:,1]+2
    pos = np.append(pos,tmp,axis=0)

tmp = cp.deepcopy(pos)
for i in range(nz-1):
    tmp[:,2] = tmp[:,2]+2
    pos = np.append(pos,tmp,axis=0)
    

    
size = nx*ny*nz*4
pos = pos[np.lexsort((pos[:,2],pos[:,1],pos[:,0]))]
pos = np.append(np.ones((size,1)),pos,axis=1) #types
pos = np.append(np.arange(1,size+1,1).reshape(size,1),pos,axis=1) #ids
if nx2 != 0:
    num = nx1*ny*nz*4 #number of Ar atoms
    pos[num:size,1] = 2 
    masses = np.array([39.948,4*39.948])
    types = 2
else:
    masses = np.array([39.948])
    types = 1

pos[:,2:5] = pos[:,2:5]*a

buff = a/2.0
xmin = pos[:,2].min()-buff
xmax = pos[:,2].max()+buff
ymin = pos[:,3].min()-buff
ymax = pos[:,3].max()+buff
zmin = pos[:,4].min()-buff
zmax = pos[:,4].max()+buff
    
if lammps != 'none':
    with open(lammps,'w') as fid:
        fid.write('LAMMPS DATA FILE\n\n'+str(size)+' atoms\n\n'+str(types)+
                  ' atom types\n\n')
        fid.write(str(xmin)+' '+str(xmax)+' xlo xhi\n')
        fid.write(str(ymin)+' '+str(ymax)+' ylo yhi\n')
        fid.write(str(zmin)+' '+str(zmax)+' zlo zhi\n')
        fid.write('\nMasses\n\n')
        for i in range(types):
            fid.write(str(i+1)+' '+str(masses[i])+'\n')
        fid.write('\nAtoms\n\n')
        for i in range(size-1):
            fid.write(str(int(pos[i,0]))+' '+str(int(pos[i,1]))+' '+str(pos[i,2])+
                      ' '+str(pos[i,3])+' '+str(pos[i,4])+'\n')
        fid.write(str(int(pos[-1,0]))+' '+str(int(pos[-1,1]))+' '+str(pos[-1,2])+
                      ' '+str(pos[-1,3])+' '+str(pos[-1,4]))
        
if xyz != 'none':
    with open(xyz,'w') as fid:
        fid.write(str(size)+'\n')
        fid.write(str(xmin)+' '+str(xmax)+' '+str(ymin)+' '+str(ymax)+' '+
                  str(zmin)+' '+str(zmax)+'\n')
        for i in range(size-1):
            fid.write(str(int(pos[i,1]))+' '+str(pos[i,2])+
                      ' '+str(pos[i,3])+' '+str(pos[i,4])+'\n')
        fid.write(str(int(pos[-1,1]))+' '+str(pos[-1,2])+
                      ' '+str(pos[-1,3])+' '+str(pos[-1,4]))

