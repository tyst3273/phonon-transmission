#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Oct, 10 2018

This is my version of the tr scipt.

This script assumes the MD simulation ran for 10,000,000 steps with
dt = 0.5e-15 s -> total time = 5 ns

@author: Ty Sterling <ty.sterling@colorado.edu>
"""
import numpy as np
import sys
import matplotlib.pyplot as plt

## GLOBAL VARIABLES ##
dT = 60 #NEMD bath temperature difference
split = 20 #number of chunks to split data into for averaging
dt = 0.5e-15 #MD timestep
dn = 10 #frequency velocities are printed
steps = 10000000 #number of simulation steps
tn = (steps/dn/split) #total number of times data was printed divided by split
win = 2/3.0 #gaussian smoothing width = THz

kb = 1.38063e-23 #Boltzmann's constant, J/K

forcefile = 'Fij.dat'
velsfile = 'vels.compact.dat'

outfile = 'tr.dat'

print('\tTemp difference is '+str(dT)+' K\n\n\tVelocity data print interval '
      'is '+str(dn)+' steps\n\n\tMD timestep '+str(dt)+' ps\n\t-----------\n')
####

## READ IN FORCES ##
with open(forcefile, 'r') as fid:
    print('\tReading forces from file...\n')
    NL = int(fid.readline().strip().split()[1]) #number of atoms on left
    NR = int(fid.readline().strip().split()[1]) #number of atoms on right
    N = NL+NR #total number of atoms
    kii = np.zeros((NL*3,NL*3)) #force constant matrix of left atoms 
    kij = np.zeros((NR*3,NR*3)) #force constant matrix or right atoms
    idsL = np.zeros((1,1)) #ids of left atoms
    idsR = np.zeros((1,1)) #ids of right atoms
    idsfij = np.zeros((N,1)) #LAMMPS ids in force file
    
    for i in range(N):
        tmp = fid.readline().strip().split()
        idsfij[i] = int(tmp[0])
        if int(tmp[1]) == 1:
            idsL = np.append(idsL,i)
        else:
            idsR = np.append(idsR,i)
    idsL = np.delete(idsL,0,axis=0).astype(int)
    idsR = np.delete(idsR,0,axis=0).astype(int)
    
    if (NL != NR) or (len(idsL) != len(idsR)): #make sure NL == NR
        sys.exit('Number of atoms on left must equal number on right!')
    
    du = float(fid.readline().strip().split()[1]) #atomic displacement size
    
    for i in range(NL*3): #move left atoms (i) x, y, and z
        fplus = np.zeros((N,3)) #dFij from positive shift of atom i
        fminus = np.zeros((N,3)) #dFij from negative shift of atom i
        
        for j in range(9): #skip comments
            fid.readline()
        for k in range(N): #dF on all atoms, positive shift
            tmp = fid.readline().strip().split()
            fplus[k,:] = tmp[1:4]
            
        for j in range(9): #skip comments
            fid.readline()
        for k in range(N): #dF on all atoms, negative shift of atom i
            tmp = fid.readline().strip().split()
            fminus[k,:] = tmp[1:4]
            
        fii = np.subtract(fplus[idsL,:],fminus[idsL,:]) #total dFii on atoms i
        fij = np.subtract(fplus[idsR,:],fminus[idsR,:]) #total dFij on atoms j
        
        kii[i,:] = (np.reshape(fii,NL*3)).transpose()
        kij[i,:] = (np.reshape(fij,NL*3)).transpose()

#Force constant matricies Kij = [d^2U/duiduj] = [-dFj/dui]
kii = -kii/(2*du) #finite difference derivative: dFplus-dFminus/2*du
kij = -kij/(2*du) 
del fij, fii, fminus, fplus, forcefile #wont need them again
####

## READ IN VELOCITIES AND CALCULATE HEAT CURRENT ##
with open(velsfile, 'r') as fid:
    print('\tDone reading forces, now reading velocities \n\tand calculating ' 
          'heat current...\n'
          '\t\tNumber of chunks is '+str(split))
    
    ## READ COMMENTS AND IDS ##
    tmp = fid.readline().strip().split() #number of atoms in velocity run
    if int(tmp[1]) != N:
        sys.exit('Number of atoms in force calulcation must match number in'
                 'velocity calculation!')
        
    tmp = fid.readline().strip().split()
    if int(tmp[1]) != dn: #check is velocity print frequency matches run
        print('\tVelocity timestep in file doesn\'t match dn in script.\n'
              '\t\tResetting timestep to '+str(tmp[1]))
        dn = int(tmp[1])
    tmp = fid.readline().strip().split() #skip comment
        
    for i in range(N): #read in LAMMPS ids from velocity file
        tmp = fid.readline().strip().split()
        if idsfij[i] != int(tmp[0]): #if ids dont match, exit
            sys.exit('Atom ids in force file don\'t match ids in '
                     'velocity file!')     
    tmp = fid.readline().strip().split() #skip comment
    ####
    
    om = np.arange(0,tn)*2*np.pi/(tn*dt*dn) #angular frequency
    thz = om/2/np.pi*1e-12 #frequency in THz
    dom = om[1]-om[0]
    qom = np.zeros((tn,1)).astype(complex) #spectral heat current function
    qomraw = np.zeros((tn,split)).astype(complex)
    
    dos = np.zeros((tn,1)) #vdos of whole interface region
    dosx = np.zeros((tn,1)) #x-direction vdos of interface region
    dosy = np.zeros((tn,1)) #y-direction vdos of interface region
    dosz = np.zeros((tn,1)) #z-direction vdos of interface region
    ldos = np.zeros((tn,1)) #dos of left region 
    ldosx = np.zeros((tn,1)) # ...
    ldosy = np.zeros((tn,1))
    ldosz = np.zeros((tn,1))
    rdos = np.zeros((tn,1)) #dos of right region
    rdosx = np.zeros((tn,1)) # ...
    rdosy = np.zeros((tn,1))
    rdosz = np.zeros((tn,1))
    
    for i in range(split): #average across split number of chunks
        print('\n\t\tNow on chunk: '+str(i+1))        
        vels = np.zeros((N*3*tn)) #velocity data is single column,
        #v1x,v1y,v1z,v2x,v2y... Read all atoms tn times
               
        for j in range(N*3*tn): #velocities of all N atoms in all 3 directions
            if (j/(N*3.0))%(tn/10) == 0:
                print('\t\t'+str(((j/(N*3)/(tn/10))+1)*10)+'% done with chunk '+
                                 str(i+1))
            vels[j] = float(fid.readline().strip().split()[0]) #velocity data
            
        vels = vels.reshape(tn,N*3) #rows are timestep, 
        #columns are v1x,v1y,v1x,v2x,v2y,v2z,...
        
        ## SPECTRAL HEAT CURRENT ##
        velsfft = np.fft.fft(vels,axis=0)*dt*dn #compute Fourier transform 
        #of velocities #not sure why *dt*dn is needed but matches matlab code. 
        #Something to do with dicrete fourier transforms.
        
        lvelsfft = np.zeros((tn,NL*3)).astype(complex)
        lvelsfft[:,0::3] = velsfft[:,idsL*3] #vx of left atoms
        lvelsfft[:,1::3] = velsfft[:,idsL*3+1] #vy
        lvelsfft[:,2::3] = velsfft[:,idsL*3+2] #vz
        
        rvelsfft = np.zeros((tn,NL*3)).astype(complex)
        rvelsfft[:,0::3] = velsfft[:,idsR*3] #vx of right atoms
        rvelsfft[:,1::3] = velsfft[:,idsR*3+1] #vy
        rvelsfft[:,2::3] = velsfft[:,idsR*3+2] #vz
        
        for j in range(tn):
            if j == 0: #avoid divide by zero
                qom[j,0] = 0+0j
            else:
                qom[j,0] = -(np.dot(lvelsfft[j,:], 
                   np.dot(kij,rvelsfft[j,:].transpose().conj()))/om[j])
        qom = -2*1j*qom/(tn*dt*dn)
        qom = qom*1.602e-19/1e-24 #convert to J/s. Units are F*v^2 = 
        #eV/A^2*(A/ps)^2 = eV/ps^2 #1 eV = 1.602e-19 J, 1 ps = 1e-12 s 
        qomraw[:,i] = qom[:,0] #save each split
        ####
        
        ## VIBRATIONAL DENSITY OF STATES ##
        Dos = (np.multiply(abs(velsfft),abs(velsfft))/
               np.tile(np.multiply(vels,vels).mean(axis=0),(tn,1))/(tn*dt*dn))       
        #square of ffts is autocorrelation. Divide each atom by its 
        #average temperature and by the total simulation time
        
        lDosx = Dos[:,idsL*3] #left atoms in x
        lDosy = Dos[:,idsL*3+1] #y
        lDosz = Dos[:,idsL*3+2] #z
        rDosx = Dos[:,idsR*3] #left atoms in x
        rDosy = Dos[:,idsR*3+1] #y
        rDosz = Dos[:,idsR*3+2] #z
        
        Dosx = Dos[:,0::3] #joint in x direction
        Dosy = Dos[:,1::3] #joint in y direction
        Dosz = Dos[:,2::3] #joint in z direction
        
        Dos = Dos.mean(axis=1) #avg across all atoms
        Dosx = Dosx.mean(axis=1) # ...
        Dosy = Dosy.mean(axis=1)
        Dosz = Dosz.mean(axis=1)
        
        lDosx = lDosx.mean(axis=1)
        lDosy = lDosy.mean(axis=1)
        lDosz = lDosz.mean(axis=1)
        lDos = (lDosx+lDosy+lDosz)/3
        rDosx = rDosx.mean(axis=1)
        rDosy = rDosy.mean(axis=1)
        rDosz = rDosz.mean(axis=1)
        rDos = (rDosx+rDosy+rDosz)/3
        
        dos[:,0]=Dos[:]+dos[:,0] #avg across splits
        dosx[:,0]=Dosx[:]+dosx[:,0] # ...
        dosy[:,0]=Dosy[:]+dosy[:,0]
        dosz[:,0]=Dosz[:]+dosz[:,0]
        ldos[:,0]=lDos[:]+ldos[:,0]
        ldosx[:,0]=lDosx[:]+ldosx[:,0]
        ldosy[:,0]=lDosy[:]+ldosy[:,0]
        ldosz[:,0]=lDosz[:]+ldosz[:,0]
        rdos[:,0]=rDos[:]+rdos[:,0]
        rdosx[:,0]=rDosx[:]+rdosx[:,0]
        rdosy[:,0]=rDosy[:]+rdosy[:,0]
        rdosz[:,0]=rDosz[:]+rdosz[:,0]
        ####
            
## AVERAGE VDOS CHUNKS ##
dos=dos/split #avg across all splits
dosx=dosx/split # ...
dosy=dosy/split
dosz=dosz/split
ldos=ldos/split
ldosx=ldosx/split
ldosy=ldosy/split
ldosz=ldosz/split
rdos=rdos/split
rdosx=rdosx/split
rdosy=rdosy/split
rdosz=rdosz/split
####

## AVERAGE TRANSMISSION FUNCTION ##
qom = qomraw.mean(axis=1) #heat current function
trRaw = qom/kb/dT #transmission function    

### GAUSSIAN SMOOTHING ##
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

tr = np.convolve(trRaw,gauss,mode='same') #smoothened transmission
dos = np.convolve(dos[:,0],gauss,mode='same') #smoothened dos
dosx = np.convolve(dosx[:,0],gauss,mode='same') #smoothened dosx
dosy = np.convolve(dosy[:,0],gauss,mode='same') #smoothened dosy
dosz = np.convolve(dosz[:,0],gauss,mode='same') #smoothened dosz
ldos = np.convolve(ldos[:,0],gauss,mode='same') #smoothened ldos
ldosx = np.convolve(ldosx[:,0],gauss,mode='same') #smoothened ldosx
ldosy = np.convolve(ldosy[:,0],gauss,mode='same') #smoothened ldosy
ldosz = np.convolve(ldosz[:,0],gauss,mode='same') #smoothened ldosz
rdos = np.convolve(rdos[:,0],gauss,mode='same') #smoothened rdos
rdosx = np.convolve(rdosx[:,0],gauss,mode='same') #smoothened rdosx
rdosy = np.convolve(rdosy[:,0],gauss,mode='same') #smoothened rdosy
rdosz = np.convolve(rdosz[:,0],gauss,mode='same') #smoothened rdosz
####

## WRITE TO FILE ##
with open(outfile, 'w') as fid:
    fid.write('delta T = '+str(dT)+' K\nMD timestep = '+str(dt)+' s\n'
              'Gaussian window = '+str(win)+' THz\n-------\n')
    fid.write('THz trRaw tr dos dosx dosy dosz ldos ldosx ldosy ldosz'
              ' rdos rdosx rdosy rdosz\n--------\n')
    for i in range(tn-1):
        fid.write(str(thz[i])+'\t'+str(trRaw[i])+'\t'+str(tr[i])+'\t'+
                  str(dos[i])+'\t'+str(dosx[i])+'\t'+str(dosy[i])+'\t'+
                  str(dosz[i])+'\t'+str(ldos[i])+'\t'+str(ldosx[i])+'\t'+
                  str(ldosy[i])+'\t'+str(ldosz[i])+'\t'+str(rdos[i])+'\t'+
                  str(rdosx[i])+'\t'+str(rdosy[i])+'\t'+str(rdosz[i])+'\n')
    fid.write(str(thz[-1])+'\t'+str(trRaw[-1])+'\t'+str(tr[-1])+'\t'+
              str(dos[-1])+'\t'+str(dosx[-1])+'\t'+str(dosy[-1])+'\t'+
              str(dosz[-1])+'\t'+str(ldos[-1])+'\t'+str(ldosx[-1])+'\t'+
              str(ldosy[-1])+'\t'+str(ldosz[-1])+'\t'+str(rdos[-1])+'\t'+
              str(rdosx[-1])+'\t'+str(rdosy[-1])+'\t'+str(rdosz[-1]))
####
    
## PLOTS ##
#fig, ax = plt.subplots()
#ax.plot(thz[0:6521],-np.real(tr[0:6521]),'k')
#fig.legend(('avg','x','y','z'))
#fig.savefig('tr.jpg',dpi=300)
####
    
## CLEAN UP VARIABLES ##    
del Dos, Dosx, Dosy, Dosz, lDos, lDosx, lDosy #clean up variables
del lDosz, rDos, rDosx, rDosy, rDosz, idsfij, tmp, vels, velsfile
del idsL, idsR, ldosx, ldosy, ldosz, rdosx, rdosy, rdosz, 
del split, steps, i, j, k, du, dn, #gwin, n, gauss, dom
del kij, kii, lvelsfft, rvelsfft, velsfft, qomraw, tn, win, outfile
####

print('\n\tAll Done!')
    


