#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 20:25:02 2019

@author: Ty Sterling
ty.sterling@colorado.edu

Functions for computing the anharmonic phonon transmission 
"""
import numpy as np
import sys

def tic():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    import numpy as np
    import time
    if 'startTime_for_tictoc' in globals():
        print("\n\tElapsed time is "+
              str(np.round(time.time()-
                       startTime_for_tictoc,decimals=3))+" seconds.")
    else:
        print("\n\t\tToc: start time not set") 

#def readFijk(infile):
#    """
#    This function reads in the [(d^2Fi)/(duj duk)] force constanst from 
#    LAMMPS compute commands. The structure of the infile should be blocks
#    of data that are NR atoms long. Here, NR means the number of atoms on
#    the right side of the interface. Each atom of the NR atoms in the block 
#    is the 'i' index in d^2Fi/dujduk. The LAMMPS code loops over all atoms j
#    on the otherside of the interface; j in the d^2Fi/dujduk expression. There
#    are NL = NR atoms on the left side of the inteface, i.e. the number on 
#    each side should (and MUST) be the same. Each atom on the left side, j, is 
#    looped over; in each loop iteration the atom is moved in the +x, then -x, 
#    then back to equilibrium. Then the atom is moved in the +y, then -y, then
#    back to equilibrium. Finally, the atom is moved in the +z, then -z, then 
#    back to equilibrium.
#    
#    At each of the NL loop steps, ANOTHER LOOP is intitiated that loop over all
#    NL + NR = N atoms in the whole interface region. Each of the atoms, k in
#    the d^2Fi/dujduk expression, is moved first the positive x direction then
#    the force change dFi of all the atoms on the right side due the movement of
#    j in one direction and k in another is written to the file. 
#    
#    Thus, each line in each block in the file corresponds to an atom i, on the 
#    right side of the interface. Each element on each side is the force felt
#    by the atom in the x, y, and z direction (the force vector elements). For
#    each block, the atom j is moved iteratively in the +x, then-x, then +y ...
#    etc. Each time the atom j is moved, k is moved to in all the same fashion.
#    Thus, if the interface has NL = NR = 64, then N = 128. Each block in the 
#    file will have NR = 64 lines, and there will be (NL*6)*(N*6) = 64*6*128*6 
#    = 294912 blocks. #     
#    """
tic()

filename = 'Fijk.dat' #Force constant data
print('\n\tNow reading force constants from '+str(filename)+'\n')
with open(filename) as fid:
    nl = int(fid.readline().strip().split()[1]) #number of atoms on left side
    nr = int(fid.readline().strip().split()[1]) #number of atoms on right side
    n = nl+nr #total number of atoms
    
    ids = np.zeros(n) #unique lammps ID
    side = np.zeros(n) #1 if on the left side, 2 if on the right
    for i in range(n): #read in the side data
        tmp = fid.readline().strip().split()
        ids[i] = int(tmp[0])
        side[i] = int(tmp[1])
    
    idsl = ids[np.argwhere(side[:] == 1)] #left ids
    idsr = ids[np.argwhere(side[:] == 2)] #left right
    
    if (nl != len(idsl) or nr != len(idsr) or nl != nr):
        sys.exit('\n\tLAMMPS SETUP ERROR: Number of atoms on left and '
                 'right side don\'t match!\n')
    
    du = float(fid.readline().strip().split()[1]) #dr step size
    fijk = np.zeros((3*nr,3*nl,3*n)) #see docstring
    
    num = 3*nl*3*n
    for j in range(3*nl):
        fjplus = np.zeros((nr*3,n*3)) 
        fjminus = np.zeros((nr*3,n*3))
        for k in range(3*n):
            fkplus = np.zeros((nr,3)) #df on atom i in x, y, and z
            fkminus = np.zeros((nr,3)) #df on atom i in x, y, and z       
            
            if j*k != 0 and (j*k)%int(num/10) == 0:
                print('\t\tNow '
                      +str(np.round(100*j*k/float(num/2),decimals=0))+
                      '% done reading force constants')
                
            ####################################
            for l in range(9): #skip comment lines
                fid.readline()
            for l in range(nl):
                tmp = fid.readline().strip().split() 
                fkplus[:] = tmp[1:4]
                
            for l in range(9): #skip comment lines
                fid.readline()
            for l in range(nl):
                tmp = fid.readline().strip().split() #df from minus movement
                fkminus[:] = tmp[1:4]
                
            fk = np.subtract(fkplus,fkminus)
            fjplus[:,k] = np.reshape(fk,3*nr)
            
            ###################################
            for l in range(9): #skip comment lines
                fid.readline()
            for l in range(nl):
                tmp = fid.readline().strip().split() 
                fkplus[:] = tmp[1:4]
        
            for l in range(9): #skip comment lines
                fid.readline()
            for l in range(nl):
                tmp = fid.readline().strip().split() #df from minus movement
                fkminus[:] = tmp[1:4]
                
            fk = np.subtract(fkplus,fkminus)
            fjminus[:,k] = np.reshape(fk,3*nr)
            
        fijk[:,j,:] = np.subtract(fjplus,fjminus)
        
toc()

