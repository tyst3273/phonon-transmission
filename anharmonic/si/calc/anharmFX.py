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

######################################################################
#def computeAnh(vels,fijk,idsl,idsr,dt,tn,nl,nr):
#    """
#    See PHYSICAL REVIEW B 95, 115313 (2017) ... :)
#    """
#    vf = np.fft.fft(vels,axis=0)*dt #time scaling, i have a link about it
#    #somewhere
#    
#    lvf = np.zeros((tn,3*nl)).astype(complex) #fft of left atoms
#    rvf = np.zeros((tn,3*nr)).astype(complex) #fft of right atoms
    
#    lvf[:,0::3] = vf[:,idsl*3] #vx on left side
#    lvf[:,1::3] = vf[:,idsl*3+1] #vy on left side 
#    lvf[:,2::3] = vf[:,idsl*3+2] #vz on left side
#    rvf[:,0::3] = vf[:,idsr*3] #vx on right side
#    rvf[:,1::3] = vf[:,idsr*3+1] #vy on right side 
#    rvf[:,2::3] = vf[:,idsr*3+2] #vz on right side
    
    #lvps = (abs(lvf)**2).mean(axis=1)
    #rvps = (abs(rvf)**2).mean(axis=1)
    #vps = (abs(vf)**2).mean(axis=1)

#    return [lvf, rvf, vf] #, lvps, rvps, vps]
######################################################################
def makeTime(dt,tn):
    """
    This function takes in the step size and number of steps in a block
    and returns the corresponding angular frequency and THz array
    """
    om = np.arange(0,tn)*2*np.pi/(tn*dt) #angular frequency
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
    print('\tMaximum Frequency:\t\t'+str(0.5/dt*1e-12)+'\t\tTHz')
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
def readFijk(infile):
    """
    This function reads in the [(d^2Fi)/(duj duk)] force constanst from 
    LAMMPS compute commands. The structure of the infile should be blocks
    of data that are NR atoms long. Here, NR means the number of atoms on
    the right side of the interface. Each atom of the NR atoms in the block 
    is the 'i' index in d^2Fi/dujduk. The LAMMPS code loops over all atoms j
    on the otherside of the interface; j in the d^2Fi/dujduk expression. There
    are NL = NR atoms on the left side of the inteface, i.e. the number on 
    each side should (and MUST) be the same. Each atom on the left side, j, is 
    looped over; in each loop iteration the atom is moved in the +x, then -x, 
    then back to equilibrium. Then the atom is moved in the +y, then -y, then
    back to equilibrium. Finally, the atom is moved in the +z, then -z, then 
    back to equilibrium.
    
    At each of the NL loop steps, ANOTHER LOOP is intitiated that loop over all
    NL + NR = N atoms in the whole interface region. Each of the atoms, k in
    the d^2Fi/dujduk expression, is moved first the positive x direction then
    the force change dFi of all the atoms on the right side due the movement of
    j in one direction and k in another is written to the file. 
    
    Thus, each line in each block in the file corresponds to an atom i, on the 
    right side of the interface. Each element in each line is the force felt
    by the atom in the x, y, and z direction (the force vector elements). For
    each block, the atom j is moved iteratively in the +x, then-x, then +y ...
    etc. Each time the atom j is moved, k is moved in all the same fashion.
    Thus, if the interface has NL = NR = 64, then N = 128. Each block in the 
    file will have NR = 64 lines, and there will be (NL*6)*(N*6) = 64*6*128*6 
    = 294912 blocks. #     
    
    The function returns fijk, the matrix elements of the anharmonic forces.
    fijk has the shape [3*nr,3*nl,3*n] -> [j,i,k]. Each individual force 
    constant is defined as:
        
         d^2 Fi_a
        ----------
        drj_b drk_c
        
    Which is interpreted as the change in force on atom i in the direction
    a = (x,y,z) due to the movement of atoms j in b=(x,y,z) and k in c=(x,y,z).
    
    Each element of first dimension of the fijk matrix corresponds to the 
    force on atom i due to the movement of atom j in each direction and k in 
    each direction. I.e. the first element is a matrix with the force 
    constant elements:
        
         d^2 F1_x         d^2 F1_x         d^2 F1_x      d^2 F1_x        
        ----------       ----------       ----------     ----------  .....
        dr1_x dr1_x      dr1_x dr1_y      dr1_x dr1_z    dr2_x dr1_z
        
        
          d^2 F1_x        d^2 F1_x         d^2 F1_x      d^2 F1_x        
        ----------       ----------       ----------     ----------  .....
        dr1_y dr1_x      dr1_y dr1_y      dr1_y dr1_z    dr2_x dr1_z
        
        
        d^2 F1_x         d^2 F1_x           d^2 F1_x      d^2 F1_x        
        ----------       ----------       ----------     ----------  .....
        dr1_z dr1_x      dr1_z dr1_y       dr1_z dr1_z    dr2_x dr1_z
        
        
         d^2 F1_x         d^2 F1_x          d^2 F1_x      d^2 F1_x     
        ----------       ----------        ----------     ---------  .....
        dr2_x dr1_x      dr2_x dr1_y       dr2_x dr1_z   dr2_x dr1_z 
        
            :               :                  :               :
            :               :                  :               :
        
    It seems counter intuitive to have the dimensions ordered this way but it
    is more straightforward to extract them from the file and to view them in
    Spyder. 
    
    """
    
    filename = 'Fijk.dat' #Force constant data
    with open(filename,'r') as fid:
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
        fijk = np.zeros((3*nl,3*nr,3*n)) #see docstring
        
        num = 3*nl/10
        for j in range(3*nl): #loop over x,y,z for all atoms in the left side
            
            if j != 0 and j%(num) == 0: #print progress updates
                print('\t\tNow '
                      +str(10*np.round(j/num,decimals=1))+
                      '% done reading force constants')
                
            fikplus = np.zeros((nr*3,n*3)) 
            fikminus = np.zeros((nr*3,n*3))
            
            for k in range(3*n): #loop over x,y,z all atoms at interface
                fkplus = np.zeros((nr,3)) #df on atom i in x, y, and z
                fkminus = np.zeros((nr,3)) #df on atom i in x, y, and z       
                ####################################
                for l in range(9): #skip comment lines
                    fid.readline()
                for l in range(nr): #+xk, +yk, +xk
                    tmp = fid.readline().strip().split() 
                    fkplus[l,:] = tmp[1:4]
                for l in range(9): #skip comment lines
                    fid.readline()
                for l in range(nr):
                    tmp = fid.readline().strip().split() #df from neg. dir
                    fkminus[l,:] = tmp[1:4]
                dfk = np.subtract(fkplus,fkminus) #dfi due to movement of k
                fikplus[:,k] = np.reshape(dfk,3*nr)/(2*du)
                #dFi/drk, each matrix corresponds to the '+j' in d^2Fi/drk/drj
                #rows are ix,iy,iz columns kx,ky,kz
                
            for k in range(3*n): #loop over x,y,z all atoms at interface
                fkplus = np.zeros((nr,3)) #df on atom i in x, y, and z
                fkminus = np.zeros((nr,3)) #df on atom i in x, y, and z       
                ####################################
                for l in range(9): #skip comment lines
                    fid.readline()
                for l in range(nr): #+xk, +yk, +xk
                    tmp = fid.readline().strip().split() 
                    fkplus[l,:] = tmp[1:4]
                for l in range(9): #skip comment lines
                    fid.readline()
                for l in range(nr):
                    tmp = fid.readline().strip().split() #df from pos. dir
                    fkminus[l,:] = tmp[1:4]
                dfk = np.subtract(fkplus,fkminus) #dfi due to movement of k
                fikminus[:,k] = np.reshape(dfk,3*nr)/(2*du) 
                #dFi/drk, each matrix corresponds to the '-j' in d^2Fi/drk/drj
                #rows are ix,iy,iz columns kx,ky,kz
                
            fijk[:,j,:] = -np.subtract(fikplus,fikminus)/(2*du) #dfi due to 
            #movement of j. matrix elements are third order terms, see docstring

    return [fijk, du, idsl.reshape(nl), idsr.reshape(nr), ids, nl, nr, n]


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