#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:36:56 2019

@author: Ty Sterling
ty.sterling@colorado.edu

Another go at harmonic tr, bump hasn't gone away

This script contains functions called by main
"""

import numpy as np
import sys
import scipy.linalg.blas as blas

######################################################################
def makeTime(dt,tn):
    """
    This function takes in the step size and number of steps in a block
    and returns the corresponding angular frequency and THz array
    """
    om = np.arange(0,tn)*2*np.pi/(tn*dt) #angular frequency, positive val
    thz = om/2/np.pi*1e-12 #frequency in THz
    dom = om[1]-om[0]
    
    return [om, thz, dom]

########################################################################
def printParams(dtMD,dt,dT,steps,split,tn):
    """
    Prints params to screen
    """
    print('\n\tUsing MD timestep:\t\t'+str(dtMD*1e12)+'\t\tps')
    print('\tEffective timestep:\t\t'+str(dt*1e12)+'\t\tps')
    print('\tTemperature bias:\t\t'+str(dT)+'\t\tK')
    print('\tTotal Number of Steps:\t\t'+str(steps)+'\t\t--')
    print('\tBlocks for averaging:\t\t'+str(split)+'\t\t--')
    print('\tTime per block:\t\t\t'+str(np.round(tn*dt*1e9,3))+'\t\tns')
    print('\tMaximum Frequency:\t\t'+str(np.round(0.5/dt*1e-12,3))+'\t\tTHz')
    print('\tFrequency Resolution:\t\t'+str(np.round(1/dt*1e-9/tn,3))+'\t\tMHz\n')
    print('\t--------------------------------------------------------\n')

############################################################################
def gsmooth(Raw, win, dom):
    """
    This function takes an numpy array defined over some domain and returns
    a gaussian smoothened version. "Raw" is the input numpy array. "win" is 
    the width of the gaussian smoothing window used to initialize the 
    smoothing window; e.g. if you have a signal from 0 to 20 THz and you want
    to smooth with a window of 1/3 THz, set "win = 1/3.0". dom is the constant
    spacing between values in the domain; e.g. if you have a numpy array,
    "freq", ranging from 0 to 20 THz and with length 10000, "dom = freq[1] - 
    freq[0]" = 0.002 THz. If none of this makes sense, do like I did and
    figure it our yourself ;)
    """
    
    gwin = round(win*1e12*2*np.pi/dom) #number of array elements in window
    if gwin % 2 == 0: #make sure its odd sized array
        gwin = gwin+1
    if gwin == 1:
        gauss = np.array([1]) #if array is size 1, convolve with self
    else:
        n = 2*np.arange(0,gwin)/(gwin-1)-(1) #centered at 0, sigma = 1
        n = (3*n) #set width of profile to 6*sigma i.e. end values ~= 0
        gauss = np.exp(-np.multiply(n,n)/2.0) #gaussian fx to convolve with
        gauss = gauss/np.sum(gauss) #normalized gaussian

    smooth = np.convolve(Raw,gauss,mode='same')
    return smooth


############################################################################
def readFij(infile):
    """
    This function reads in the [dFi/duj] force constanst from 
    LAMMPS compute commands. The structure of the infile should be blocks
    of data that are N atoms long. Here, N means the number of atoms in
    the interface region. Each atom of the N atoms in the block 
    is the 'i' index in d^2Fi/duj. The LAMMPS code loops over all atoms j
    on the left side of the interface; j in the dFi/duj expression. There
    are NL = NR atoms on the left side of the inteface, i.e. the number on 
    each side should (and MUST) be the same. Each atom on the left side, j, is 
    looped over; in each loop iteration the atom is moved in the +x, then -x, 
    then back to equilibrium. Then the atom is moved in the +y, then -y, then
    back to equilibrium. Finally, the atom is moved in the +z, then -z, then 
    back to equilibrium.
    
    Thus, each line in each block in the file corresponds to an atom i. 
    Each element in each line is the force felt by the atom in the x, y, and 
    z direction (the force vector elements). For each block, the atom j is 
    moved iteratively in the +x, then-x, then +y ... etc.
    
    The function returns kij, the matrix elements of the force constants 
    between atoms on either side of the interface.
    kij has the shape [3*nr,3*nl] -> [dFi-x,y,z,duj-x,y,z]. Each individual 
    force constant is defined as:
        
         dFi_a
        ----------
         drj_b
        
    Which is interpreted as the change in force on atom i in the direction
    a = (x,y,z) due to the movement of atom j in b=(x,y,z).
    
    Each element of first dimension of the kij matrix corresponds to the 
    force on atom i due to the movement of atom j in each direction 
        
         dF1_x            dF1_x             dF1_x          dF1_x        
        ----------       ----------       ----------     ----------  .....
         dr1_x            dr1_y            dr1_z           dr2_x
        
        
          dF1_y           dF1_y             dF1_y          dF1_y        
        ----------       ----------       ----------     ----------  .....
          dr1_x           dr1_y             dr1_z          dr2_x
        
        
           dF1_z           dF1_z             dF1_z          dF1_z        
        ----------       ----------       ----------     ----------  .....
           dr1_x           dr1_y             dr1_z          dr2_x
          
          
           dF2_x           dF2_x             dF2_x          dF2_x        
        ----------       ----------       ----------     ----------  .....
          dr1_x           dr1_y             dr1_z          dr2_x
        
            :               :                  :               :
            :               :                  :               :    
    """
    with open(infile,'r') as fid:
        nl = int(fid.readline().strip().split()[1]) #number of atoms on left side
        nr = int(fid.readline().strip().split()[1]) #number of atoms on right side
        n = nl+nr #total number of atoms
        
        ids = np.zeros(n) #unique lammps ID
        side = np.zeros(n) #1 if on the left side, 2 if on the right
        for i in range(n): #read in the side data
            tmp = fid.readline().strip().split()
            ids[i] = int(tmp[0])
            side[i] = int(tmp[1])
            
        idsl = np.argwhere(side[:] == 1) #left ids
        idsr = np.argwhere(side[:] == 2) #left right
        
        if (nl != len(idsl) or nr != len(idsr) or nl != nr):
            sys.exit('\n\tLAMMPS SETUP ERROR: Number of atoms on left and '
                     'right side don\'t match!\n')
            
        du = float(fid.readline().strip().split()[1]) #dr step size
        
        kii = np.zeros((nr*3,nl*3)) #force constants
        kij = np.zeros((nr*3,nl*3)) #force constants
        
        for i in range(nl*3):
            fplus = np.zeros((n,3)) # j in plus x,y,z
            fmin = np.zeros((n,3)) # j in minus x,y,z
        
            for j in range(9): #skip comments
                fid.readline()
            for j in range(n): #read force dF
                fplus[j,:] = fid.readline().strip().split()[1:4]
                
            for j in range(9): #skip comments
                fid.readline()
            for j in range(n): #read force dF
                fmin[j,:] = fid.readline().strip().split()[1:4]
                
            kii[:,i] = (fplus[idsl,:]-fmin[idsl,:]).reshape(nl*3)
            kij[:,i] = (fplus[idsr,:]-fmin[idsr,:]).reshape(nr*3)
            
    kii = -kii/(2*du)
    kij = -kij/(2*du)

    return [kij, du, idsl.reshape(nl), idsr.reshape(nr), ids, nl, nr, n]
##########################################################################
    


##########################################################################
def computeQ(kij,vels,om,idsl,idsr,nl,nr,tn,dt,s):
    """
    Compute the spectral heat current across the interface
    """
    vfft = np.fft.fft(vels,axis=0)*dt
    
    qom = np.zeros(tn).astype(complex)
    vl = np.zeros((tn,nl*3)).astype(complex)
    vr = np.zeros((tn,nr*3)).astype(complex)
    
    vl[:,0::3] = vfft[:,idsl*3]
    vl[:,1::3] = vfft[:,idsl*3+1]
    vl[:,2::3] = vfft[:,idsl*3+2]
    vr[:,0::3] = vfft[:,idsr*3]
    vr[:,1::3] = vfft[:,idsr*3+1]
    vr[:,2::3] = vfft[:,idsr*3+2]

    print('\n\tNow computing heat current for block '+str(s+1)+'\n')
    num = tn/10
    for j in range(tn):
        if j != 0 and j%(num) == 0: #print progress updates
            print('\t\tNow '+str(10*np.round(j/num,decimals=1))+
                  '% done computing heat current for block '+str(s+1))
        if j == 0:
            qom[j] = 0+0j
        else:
            qom[j] = -(blas.cgemm((1+0j),vr[j,None].conj(),
               blas.cgemm((1+0j),kij,vl[j,None].T)))/om[j]
            
    return qom

############################################################################
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
    import time
    if 'startTime_for_tictoc' in globals():
        print("\n\t\tElapsed time is "+
              str(np.round(time.time()-
                       startTime_for_tictoc,decimals=3))+" seconds.")
    else:
        print("\n\t\tToc: start time not set") 
        
############################################################################
def writeTr(outfile,thz,tr,qom,dT,dt,tn,win):
    ""
    ""
    with open(outfile, 'w') as fid:
        fid.write('delta T = '+str(dT)+' K\nMD timestep = '+str(dt)+' s\n'
                  'Gaussian window = '+str(win)+' THz\n-------\n')
        fid.write('THz tr qomRaw\n--------------------------------------\n')
        for i in range(tn-1):
            fid.write(str(thz[i])+'\t'+str(tr[i])+'\t'+str(qom[i])+'\n')
        fid.write(str(thz[-1])+'\t'+str(tr[-1])+'\t'+str(qom[-1]))
        
        
        
#########################################################################
def readTr(infile):
    """
    """
    with open(infile, 'r') as fid:
        num_lines = sum(1 for line in fid)-6 #get number of line in file
        fid.seek(0) #return cursor to start of file
        
        for i in range(4):
            fid.readline() #skip comments
        
        columns = fid.readline().strip().split()
        ncol = len(columns)
        data = np.zeros((num_lines,ncol))
        fid.readline() #skip comment
        
        for i in range(num_lines):
            tmp = fid.readline().strip().split()
            for j in range(ncol):
                data[i,j] = float(tmp[j])
            
    return [data, columns]