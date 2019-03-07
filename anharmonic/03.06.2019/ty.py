#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This file is Ty's custom module to import functions into scritps

DATE STAMP: 01.08.2019 MM.DD.YYYY
"""

############################################################
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
    import numpy as np
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



def smoothSED(Raw, win, dom):
    """
    Same as gsmooth but smooths 2d array along axis=0; designed for smoothing
    SED.
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
    import numpy as np
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
    
    smooth = np.zeros((len(Raw[:,0]),len(Raw[0,:])))
    for i in range(len(Raw[0,:])):
        smooth[:,i] = np.convolve(Raw[:,i],gauss,mode='same')
    return smooth


############################################################
def findNN(pos,index=0):
    """
    This function takes a position array with 5 columns = [id, type, x, y, z]
    as its argument. It returns a list containing [nl , dl, nn] where nl, 
    dl, and nn are the neighbor list, distance to neighbor lists, and number 
    of nearest neighbors. Each row in the neighbor list correspond to the atom
    at the same position in the position array. Each colum is that atoms 
    nearest neighbors sorted in ascending order. The first element is the atom. 

    Optional argument "index" can be either intdex=0 or index=1. 0 is default.
    index=0 returns neighbor lists whose elements correspond to the POSITION
    of the neighbor in the position array. index=1 returns the unique id of 
    the neighbor atom provided in the positions array.
    
    For each atom and atom in nl, dl provides the distance between the two.
    
    nn provides the number SELF=1, 1st nn, 2nd nn, 3rd nn, and so on
    for each atom in the position array
    """
    import numpy as np
    import sys
    if index != 0 and index != 1:
        sys.exit('Index must be 0 (default) of 1')
    num = len(pos[:,0])
    
    nl = np.zeros((num,num))
    dl = np.zeros((num,num))
    nn = np.zeros((num,num))
    
    for i in range(num): #loop over atoms
        tmp = np.zeros((num,2))
        for j in range(num): #loop over neighbors
            tmp[j,0] = j
            tmp[j,1] = np.round(np.sqrt((pos[i,2]-pos[j,2])**2+
               (pos[i,3]-pos[j,3])**2+(pos[i,4]-pos[j,4])**2),decimals=2) 
                #euclidean metric
        if index == 1: #unique ids or indicies
            tmp[:,0] = pos[:,0]
        tmp = tmp[np.argsort(tmp[:,1]),:] #sort
        nl[i,:] = tmp[:,0] #neighbor list
        dl[i,:] = tmp[:,1] #distance list
        
    for i in range(num): #find number of neighbors
        count = 0 #accumulate number of neighbors, dont recount
        n = 0 #number of neighbors
        for j in range(num):
            if count == num: 
                break #cant have more neighbors than there are atoms
            dist = dl[i,count] #neighbor distance
            n = np.sum(dl[i,:] == dist) #number of neighbors
            nn[i,j] = n #number of neighbors
            count = count+n
            
    return [nl,dl,nn]


############################################################
def readData(filename):
    """
    This function reads a lammps data file and returns:1: the number of atoms,
    "num", in the file. 2: the number of atom types, "types". 3: the masses of
    each atom type and 4: an array, "pos", that contains 5 columns that are 
    num elements long. The columns are [id, type, x, y, z].
    """
    import numpy as np
    with open(filename, 'r') as fid:
        for i in range(2): #skip comments
                fid.readline()
                
        num = int(fid.readline().strip().split()[0]) #number of atoms
        fid.readline() #skip comments
        
        types = int(fid.readline().strip().split()[0]) #number of atom types
        masses = np.zeros((types))
        for i in range(7): #skip comments
            fid.readline()
        for i in range(types):
            masses[i] = float(fid.readline().strip().split()[1]) #atomic masses
        for i in range(3):
            fid.readline() #skip comments
            
        pos = np.zeros((num,5))
        for i in range(num):
            tmp = fid.readline().strip().split()
            pos[i,0] = int(tmp[0]) #atom id
            pos[i,1] = int(tmp[1]) #atom type
            pos[i,2] = float(tmp[2]) #x-coordinate
            pos[i,3] = float(tmp[3]) #y-coordinate
            pos[i,4] = float(tmp[4]) #z-coordinate
            
    return [num, types, masses, pos]


############################################################
def readFij(filename,num,pos,basis,nl,nn,neighbor=2,tol=0):
    """
    This function reads a file, 'forcefile', of lammps force data calculated
    by finite displacements. Some assumptions are made about the structure of 
    the lammps output file. The first line should be 'du ###' where ### is the
    size of the displacement step. The lammps code loops over all atoms, num,
    by displacing them first in the positive x direction then the negative 
    x and then the positive y, negative y, etc. It does this for each atom.
    Each time an atom is displaced, the forces acting on all atoms is printed
    to the file. The difference in forces in each direction, e.g. (F_+x)-(F_-x)
    is calculated and interpreted as the 'finite difference derivative' dF. 
    The force constant d^2V_ij/du_i/du_j is the second order term of the
    taylor expansion of the potential energy calculated as d^2V_ij/du_i/du_j =
    -dF_i/du_j. dF_i is the difference in force on each atom du to displacement
    of atom j in the x,y,z directions. NOTE: its acutally -dF/(2*du) since 
    the derivative is taken about the equilibrium position in +/- minus x,y,z.
    
    The matrix fij this creates is formatted as follows:
        
               i1x      i1y      i1z       i2x ...
        ____|______________________________________
            |
        j1x |  -dFj1x   -dFj1x  -dFj1x   -dFj1x   
            |   ----     ----     ----     ----   ... 
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j1y |  -dFj1y   -dFj1y   -dFj1y   -dFj1y  
            |   ----     ----     ----     ----   ...
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j1z |  -dFj1z   -dFj1z   -dFj1z   -dFj1z   
            |   ----     ----     ----     ----   ...
            |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
            |
            |
        j2x |  -dFj2x  -dFj2x   -dFj2x    -dFj2x
          . |   ----    ----     ----      ----   ...
          . |  2*dui1x  2*dui1y  2*dui1z  2*dui2x
          . |     .       .         .        .
                  .       .         .        .
                  .       .         .        .
                  
    The optional argument 'tol' can be set to some small float value.
    if the magnitude of the value in the matrix is lower than tol, it is set
    to 0. This could be used to sparsify the matrix but I dont care.
 
    The rest of the arguments are used to find which force constant matricies
    correcpond to which atoms. This function returns phi and rvec. Each of 
    this objects contains a set of item for each basis atom and each item 
    contains 1 item per nearest neighbor (included the basis atom) up to the 
    2nd neighbor. E.g. for si, this is 17 items per basis atom, 2 basis atoms.
    
    The optional argument neighbor is how many neighbors to include .. e.g. 
    by default, neighbor=2 includes only up to the second nearest neighbor.
    
    phi contains the force constant matricies and rvec contains the direction
    vectors to the neighbor atoms ; this are needed to compute the dynamical
    matrix.
          
    """
    with open(filename, 'r') as fid:
        import numpy as np
#        import matplotlib.pyplot as plt
        
        fij = np.zeros((num*3,num*3))
        du = float(fid.readline().strip().split()[1]) #angstroms
        
        for i in range(num*3): #iteratively displace each atom, du_i, in x,y,z
            fplus = np.zeros((num,3)) #fx, fy, fz; positive displacement
            fminus = np.zeros((num,3)) #fx, fy, fz; positive displacement
            
            for j in range(9): 
                fid.readline() #skip comments
            for j in range(num): #forces on all atoms due to du_i plus x,y,z
                tmp = fid.readline().strip().split() #get data
                fplus[j,:] = tmp[1:4]
                
            for j in range(9): 
                fid.readline() #skip comments
            for j in range(num): #forces on all atoms due to du_i minus x,y,z
                tmp = fid.readline().strip().split() #get data
                fminus[j,:] = tmp[1:4]
                
            df = fplus-fminus #dF/(2*du)
            df = np.reshape(df,num*3) #row vector for each dui_x,y,z ; each 
            #element is the force constant between atom i in x,y,z (each row)
            #and each other atom j in x,y,z (each column)
    
            if tol != 0:
                df[np.argwhere(np.abs(df) < tol)] = 0 
                
            fij[i,:] = df #force constant matrix

            
        fij = -fij/(2*du) #force constant = d^2Uij/dui/dij = -dFj/dui
    
    nb = len(basis) #mumber of basis atoms
    phi = [0]*nb #force constant matricies
    rvec = [0]*nb #direction vectors
    for i in range(nb):
        atom = basis[i] #basis atom
        no = int(np.sum(nn[atom,0:neighbor+1])) #number of nn's
        phi[i] = [0]*no #fij matrix for each neighbor
        rvec[i] = [0]*no #direction vector to each neighbor
        tmp = nl[atom,0:no].astype(int)
        for j in range(no): #look up the matrix for each atom
            rvec[i][j] = pos[tmp[j],2:5]-pos[atom,2:5]
            phi[i][j] = fij[atom*3:(atom*3+3),tmp[j]*3:(tmp[j]*3+3)]
        
    #fig, ax = plt.subplots() #visualize the fij matrix
    #plot = ax.imshow(fij,cmap='plasma')
    #fig.colorbar(plot,ax=ax)
    #fig.tight_layout()
    #plt.show()
    
    return phi, rvec


############################################################
def makeKpoints(prim,specialk,dk):
    """
    This function takes primitive lattice vectors, prim, and constructs
    reciprocal lattice vectors from them. It also takes an arbitrarily long
    (n,3) array of normalized special k points and an integer, dk.
    Reciprocal lattice vectors are constructed from the primitive lattice 
    vector and then those reciprocal lattice vectors are projected onto the
    special k points. An array of k points containing the special k points
    is then created and dk points between each k point are populated linearly.

    This function returns 2 variables: kpoints and kdist. kpoints is the array
    containing ((n-1)*dk,3) k points. kdist is cumulative sum of the euclidean
    k-space distance between each special k point - this is used later for
    plotting.
    """
    import numpy as np
    kvec = np.zeros((3,3)) #k space vectors 
    vol = np.dot(prim[0,:],np.cross(prim[1,:],prim[2,:])) #real space cell volume
    kvec[0,:] = 2*np.pi*np.cross(prim[1,:],prim[2,:])/vol
    kvec[1,:] = 2*np.pi*np.cross(prim[2,:],prim[0,:])/vol
    kvec[2,:] = 2*np.pi*np.cross(prim[0,:],prim[1,:])/vol
    #reciprocal lattice vectors with cubic symmetry are parallel to real 
    #space vectors
    for i in range(len(specialk)): #project into reciprocal space 
        specialk[i,:] = np.dot(kvec,specialk[i,:])
        
    nk = (len(specialk)-1)*dk
    kpoints = np.zeros((nk,3)) 
    for i in range(len(specialk)-1): #populate space between special k points
        for j in range(3): #kx, ky, kz
            kpoints[i*dk:(i+1)*dk,j] = np.linspace(specialk[i,j],
                    specialk[(i+1),j],dk) 
            
    kdist = np.zeros((len(specialk)))
    for i in range(len(specialk)-1): #euclidean distance between k points
        kdist[(i+1)] = np.sqrt((specialk[(i+1),0]-specialk[i,0])**2+
             (specialk[(i+1),1]-specialk[i,1])**2+
             (specialk[(i+1),2]-specialk[i,2])**2) 
    kdist = np.cumsum(kdist) #cumulative distance between special k points
    
    return [kpoints, kdist]


############################################################
def dynamicalMat(phi,rvec,pos,basis,masses,kpoints,nl):
    """
    This function computes the phonon dynamical matrix by space fourier
    transforming the force constants matricies. The eigenvales of the 
    dynamical matrix are the frequencies of each branch; the eigenvectors
    are the polarization vectors. Input arguments are the force constant 
    matricies around each bases atom, phi, the displacement vectors to each
    neighbor, rvec. It returns eigen frequencies, om, and could be 
    adapted to return polarization vectors. The frequencies are in THZ
    """
    import numpy as np
    nb = len(basis) #number of basis atoms
    nk = len(kpoints) #number of k points to sum over
    
    om = np.zeros((nk,nb*3)) #eigenvalues, 3 solutions per basis atom
    
    for i in range(nk):
        kvec = kpoints[i,:] #k point
        dyn = np.zeros((nb*3,nb*3)).astype(complex) #dynamical matrix
        for j in range(nb): #sum over basis atoms
            atom1 = basis[j] #basis atom
            type1 = int(pos[atom1,1])-1
            mass1 = masses[type1]
            index1 = np.arange(0,3)+type1*3 #indicies in dyn matrix
            for k in range(len(phi[j][:])):
                atom2 = int(nl[atom1,k]) #neighbor atom
                type2 = int(pos[atom2,1])-1
                mass2 = masses[type2]
                index2 = np.arange(0,3)+type2*3 #indicies in dyn matrix
                dyn[index1[0]:index1[-1]+1,index2[0]:index2[-1]+1] = (
                        dyn[index1[0]:index1[-1]+1,index2[0]:index2[-1]+1]+
                        phi[j][k]*np.exp(1j*np.dot(kvec,rvec[j][k]))/
                        np.sqrt(mass1*mass2)) #space FFT of force constant
                        #matrix, dived by sqrt of the masses comes from
                        #the eigenvalue equation
                        
        eigVal, eigVec = np.linalg.eig(dyn) #eigen values are frequency, vectors
        #are displacement vectors
        om[i,:] = np.sqrt(np.real(eigVal))*1000/2.0/np.pi/10.18
        
    return om #thz


############################################################
def makeFCCdiamond(nx,ny,nz,lammps='no',element='si'):
    """
    This function takes arguments x,y,z: the length of the structure (in 
    number of unit cells) in each direction. Optional argument element
    is default to si but can be c or ge (carbon or germanium).
    If its c or ge, the mass and lattice constant of c and ge are used.
    This function returns num, pos, an array with nx*ny*nz*8 atoms, 5 columns,
    masses, a 2 element array with the mass of the element duplicated for each
    basis site. it also returns uc, a list of the unit cell index of each atom.
    The format of pos is [id, type, x, y, z]. type 1 is the atom at the fcc 
    basis site and 2 is the neighbor at the diamond site. each element in 
    uc correspond to the atom in the same position in pos. uc is needed to
    comput SED. num is the number of atoms. lammps can be set to a string
    the is the output filename.
    """
    import numpy as np
    import sys
    import copy as cp
    
    if element == 'si':
        a = 5.431 #Si lattice constant
        masses = np.array([28.0855,28.0855]) #atomic mass of Si
    if element == 'ge':
        a = 5.658 #Ge lattice constant
        masses = np.array([72.6400,72.6400]) #atomic mass of Ge
    if element == 'c':
        a = 3.57 #C lattice constant
        masses = np.array([12.0107,12.0107]) #atomic mass of C
    if (nx%2.0 != 0) and (ny%2.0 != 0) and (nz%2.0 != 0):
        sys.exit('Number of unit cells in each direction must be an even'
                 ' integer')
        
    nT = np.array([nx,ny,nz]).astype(int) #times to translate in each direction
    unit = np.array(([1,0,0,0], #FCC-diamond unit cell (basis-index,x,y,z,uc)
                     [1,0,2,2], #Index 0 are FCC sites, 1 are diamond sites,
                     [1,2,0,2], #uc is the particular unit cell
                     [1,2,2,0],
                     [2,1,1,1],
                     [2,1,3,3],
                     [2,3,1,3],
                     [2,3,3,1])).astype(float)
    
    uc = np.array([0,1,2,3,0,1,2,3])
    
    tmp = cp.deepcopy(unit)
    pos = cp.deepcopy(unit) #basis-index and coordinates of atoms
    tuc = cp.deepcopy(uc)
    for i in range(nT[0]-1): #translate in x-direction 
        tmp[:,1] = tmp[:,1]+4 #array positions in x
        pos = np.append(pos,tmp,axis=0)
        tuc = tuc+4
        uc = np.append(uc,tuc)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(nT[1]-1): #translate in y-direction
        tmp[:,2] = tmp[:,2]+4 
        pos = np.append(pos,tmp,axis=0)
        tuc2 = tuc+max(uc)+1
        uc = np.append(uc,tuc2)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(nT[2]-1): #translate in z-direction
        tmp[:,3] = tmp[:,3]+4
        pos = np.append(pos,tmp,axis=0)
        tuc2 = tuc+max(uc)+1
        uc = np.append(uc,tuc2)
    
    num = len(pos[:,0])
    ids = np.zeros((num,1))
    ids[:,0] = np.arange(0,num,1)+1 #will be used later    
    pos = np.append(ids,pos,axis=1)
    pos[:,2:5] = pos[:,2:5]*a/4.0 #rescale coordinates
    
    if lammps != 'no':
        buff = a/8.0
        with open(lammps, 'w') as fid:
            xmin = min(pos[:,2]); xmax = max(pos[:,2])
            ymin = min(pos[:,3]); ymax = max(pos[:,3])
            zmin = min(pos[:,4]); zmax = max(pos[:,4])

            fid.write(str('LAMMPS DATA FILE\n'))
        
            fid.write('\n' + str(len(pos)) + ' atoms\n')
            fid.write('\n' + str(len(masses)) + ' atom types\n')
            fid.write('\n' + str(xmin-buff)+' '+str(xmax+buff)+' xlo'+' xhi\n')
            fid.write(str(ymin-buff)+' '+str(ymax+buff)+' ylo'+' yhi\n')
            fid.write(str(zmin-buff)+' '+str(zmax+buff)+' zlo'+' zhi\n')
            fid.write('\nMasses\n')
            for i in range(len(masses)):
                fid.write('\n' + str(i+1) + ' ' + str(float(masses[i])))
            fid.write('\n\nAtoms\n\n')
            for i in range(len(pos)-1):
                fid.write(str(int(i+1)) + ' ' + str(int(pos[i,1])) + ' ' 
                          + str(pos[i,2]) + ' ' +
                        str(pos[i,3]) + ' ' + str(pos[i,4]) + '\n')
            fid.write(str(len(pos)) +  ' ' + str(int(pos[-1,1])) + ' ' 
                      + str(pos[-1,2]) + ' ' +
                    str(pos[-1,3]) + ' ' + str(pos[-1,4]))

    return [num, pos, masses, uc, a]     



def makeTriclinic(n1,n2,n3,lammps='no',element='si'):
    """
    See docstring for makeFCCdiamond
    """
    import numpy as np
    import sys
    import copy as cp
    
    if element == 'si':
        a = 5.431 #Si lattice constant
        masses = np.array([28.0855,28.0855]) #atomic mass of Si
    if element == 'ge':
        a = 5.658 #Ge lattice constant
        masses = np.array([72.6400,72.6400]) #atomic mass of Ge
    if element == 'c':
        a = 3.57 #C lattice constant
        masses = np.array([12.0107,12.0107]) #atomic mass of C
    if (n1%2.0 != 0) and (n2%2.0 != 0) and (n3%2.0 != 0):
        sys.exit('Number of unit cells in each direction must be an even'
                 ' integer')
        
    nT = np.array([n1,n2,n3]).astype(int) #times to translate in each direction
    unit = np.array(([1,0,0,0], #FCC-diamond unit cell (basis-index,x,y,z,uc)
                     [2,1,1,1])).astype(float)
    
    uc = np.array([0,0])
    
    tmp = cp.deepcopy(unit)
    pos = cp.deepcopy(unit) #basis-index and coordinates of atoms
    tuc = cp.deepcopy(uc)
    for i in range(nT[0]-1): #translate in a1-direction 
        tmp[:,1] = tmp[:,1]+2
        pos = np.append(pos,tmp,axis=0)
        tuc = tuc+1
        uc = np.append(uc,tuc)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(nT[1]-1): #translate in a2-direction
        tmp[:,1] = tmp[:,1]+1
        tmp[:,2] = tmp[:,2]+3
        pos = np.append(pos,tmp,axis=0)
        tuc2 = tuc+max(uc)+1
        uc = np.append(uc,tuc2)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(nT[2]-1): #translate in a3-direction
        tmp[:,1] = tmp[:,1]+1
        tmp[:,2] = tmp[:,2]+1
        tmp[:,3] = tmp[:,3]+4
        pos = np.append(pos,tmp,axis=0)
        tuc2 = tuc+max(uc)+1
        uc = np.append(uc,tuc2)
    
    num = len(pos[:,0])
    ids = np.zeros((num,1))
    ids[:,0] = np.arange(0,num,1)+1 #will be used later    
    pos = np.append(ids,pos,axis=1)
    pos[:,2] = pos[:,2]*a*np.sqrt(2.0)/4.0 #rescale coordinates
    pos[:,3] = pos[:,3]*a/np.sqrt(24.0) #rescale coordinates
    pos[:,4] = pos[:,4]*a/np.sqrt(3.0)/4.0 #rescale coordinates
    
    if lammps != 'no':
        with open(lammps, 'w') as fid:
            xbuff = (pos[np.argwhere(pos[:,3] == 0)[:,0],2][1]-
                     pos[np.argwhere(pos[:,3] == 0)[:,0],2][0])/2.0
            xmax = pos[np.argwhere(pos[:,3] == 0)[:,0],2].max()+xbuff
            xmin = 0-xbuff #by construction
            ybuff = (np.unique(pos[np.argwhere(pos[:,4] == 0)[:,0],3])[1]-
                     np.unique(pos[np.argwhere(pos[:,4] == 0)[:,0],3])[0])/2.0
            ymax = pos[np.argwhere(pos[:,4] == 0)[:,0],3].max()+ybuff
            ymin = 0-ybuff #by construction
            zbuff = (np.unique(pos[:,4])[2]-np.unique(pos[:,4])[1])/2.0
            zmax = pos[:,4].max()+zbuff
            zmin = 0-zbuff #by construction
            
            xy = pos[np.argwhere(pos[:,4] == pos[:,4].max())[:,0],:]
            xy = xy[np.argwhere(xy[:,3] == xy[:,3].min())[:,0],2].min()
            xz = xy
            yz = pos[np.argwhere(pos[:,4] == pos[:,4].max())[:,0],3].min()
    
            fid.write(str('LAMMPS DATA FILE\n'))
        
            fid.write('\n' + str(len(pos)) + ' atoms\n')
            fid.write('\n' + str(len(masses)) + ' atom types\n')
            fid.write('\n' + str(xmin)+' '+str(xmax)+' xlo'+' xhi\n')
            fid.write(str(ymin)+' '+str(ymax)+' ylo'+' yhi\n')
            fid.write(str(zmin)+' '+str(zmax)+' zlo'+' zhi\n')
            fid.write(str(xy)+' '+str(xz)+' '+
                          str(yz)+' xy xz yz\n')
            fid.write('\nMasses\n')
            for i in range(len(masses)):
                fid.write('\n' + str(i+1) + ' ' + str(float(masses[i])))
            fid.write('\n\nAtoms\n\n')
            for i in range(len(pos)-1):
                fid.write(str(int(i+1)) + ' ' + str(int(pos[i,1])) + ' ' 
                          + str(pos[i,2]) + ' ' +
                        str(pos[i,3]) + ' ' + str(pos[i,4]) + '\n')
            fid.write(str(len(pos)) +  ' ' + str(int(pos[-1,1])) + ' ' 
                      + str(pos[-1,2]) + ' ' +
                    str(pos[-1,3]) + ' ' + str(pos[-1,4]))
            
    return [num, pos, masses, uc, a] 



def makeGaN(nx,ny,nz,lammps='no'):
    """ 
    See docstring for makeFCCdiamond
    """
    import numpy as np
    import copy as cp
    masses = np.array([69.723,14.007,69.723,14.007,69.723,14.007,69.723,14.007])
    a = 3.189 #a = b = diagonal length
    c = np.round(np.sqrt(8/3.0)*a,decimals=3) #5.185/2.0 #ideal c = sqrt(8/3)*a
    
    basis = np.array(([1,0,0,0], #Ga
                      [2,0,0,3], #N
                      [3,2,0,4], #Ga
                      [4,2,0,7], #N
                      [5,3,1,0], #Ga
                      [6,3,1,3], #N
                      [7,5,1,4], #Ga
                      [8,5,1,7])).astype(float) #N
    uc = np.array([0,0,0,0,0,0,0,0])
    
    pos = cp.deepcopy(basis)
    tmp = cp.deepcopy(basis)
    tuc = cp.deepcopy(uc)
    for i in range(nx-1):
        tmp[:,1] = tmp[:,1]+(6)
        pos = np.append(pos,tmp,axis=0)
        tuc = tuc+1
        uc = np.append(uc,tuc)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(ny-1):
        tmp[:,2] = tmp[:,2]+2
        tuc = tuc+nx
        uc = np.append(uc,tuc)
        pos = np.append(pos,tmp,axis=0)
    tmp = cp.deepcopy(pos)
    tuc = cp.deepcopy(uc)
    for i in range(nz-1):
        tmp[:,3] = tmp[:,3]+8
        pos = np.append(pos,tmp,axis=0)
        tuc = tuc+nx*ny
        uc = np.append(uc,tuc)
    pos[:,1] = pos[:,1]*np.round(np.sqrt(3)*a/6.0,decimals=4)
    pos[:,2] = pos[:,2]*a/2.0
    pos[:,3] = pos[:,3]*c/8.0
    num = len(pos)
    ids = np.zeros((num,1))
    ids[:,0] = np.arange(1,num+1)
    pos = np.append(ids,pos,axis=1)
    
    if lammps != 'no':    
        with open(lammps,'w') as fid:
            xbuff = 0.4603
            ybuff = 0.79725
            zbuff = 0.3255
            fid.write('LAMMPS DATA FILE FOR SED\n')
            fid.write('\n'+str(num)+' atoms\n')
            fid.write('\n'+str(len(masses))+' atom types\n')
            fid.write('\n'+str(pos[:,2].min()-xbuff)+' '+
                      str(pos[:,2].max()+xbuff)+' xlo xhi\n')
            fid.write(str(pos[:,3].min()-ybuff)+' '+str(pos[:,3].max()+ybuff)+
                      ' ylo yhi\n')
            fid.write(str(pos[:,4].min()-zbuff)+' '+str(pos[:,4].max()+zbuff)+
                      ' zlo zhi\n')
            fid.write('\nMasses\n\n')
            for i in range(len(masses)):
                fid.write(str(i+1)+' '+str(masses[i])+'\n')
            fid.write('\nAtoms\n\n')
            for i in range(num-1):
                fid.write(str(int(pos[i,0]))+' '+str(int(pos[i,1]))+' '
                          +str(pos[i,2])+' '+str(pos[i,3])+' '+
                          str(pos[i,4])+'\n')
            fid.write(str(int(pos[-1,0]))+' '+str(int(pos[-1,1]))+' '
                      +str(pos[-1,2])+' '+str(pos[-1,3])+' '+str(pos[-1,4]))
    return [num, pos, masses, uc, a, c] 


##########################################################
def writeSED(outfile,thz,kpoints,sed,dos):
    """
    This function is simple. It writes the frequency data array, k points, 
    SED matrix, and DOS array to file 'outfile'. Read file with readSED
    """
    import numpy as np
    nf = len(thz)
    nk = len(kpoints[:,0])
    sed = np.reshape(sed,(nf*nk,1))

    with open(outfile, 'w') as fid:
        fid.write('nf = '+str(nf)+'\n')
        fid.write('nk = '+str(nk)+'\n')
        for i in range(nf):
            fid.write(str(thz[i])+'\n')
        for i in range(nk):
            fid.write(str(kpoints[i,0])+'\t'+str(kpoints[i,1])+'\t'+
                      str(kpoints[i,2])+'\n')
        for i in range(nk*nf):
            fid.write(str(sed[i,0])+'\n')
        for i in range(nf):
            fid.write(str(dos[i,0])+'\n')
         
def readSED(infile):
    """
    This function reads in the SED outout file and returns THz, kpoints, SED, 
    and DOS written to file using writeSED
    """
    with open(infile,'r') as fid:
        import numpy as np
        nf = int(fid.readline().strip().split()[2])
        nk = int(fid.readline().strip().split()[2])
        thz = np.zeros((nf,1))
        kpoints = np.zeros((nk,3))
        sed = np.zeros((nf*nk,1)).astype(complex)
        dos = np.zeros((nf,1))
        for i in range(nf):
            thz[i,0] = float(fid.readline())
        for i in range(nk):
            kpoints[i,:] = fid.readline().strip().split()[:]
        kpoints = kpoints.astype(float)
        for i in range(nf*nk):
            sed[i,0] = complex(fid.readline())
        sed = np.reshape(sed,(nf,nk))
        for i in range(nf):
            dos[i,0] = float(fid.readline())
    return thz, kpoints, sed, dos


##########################################################            
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

def toc(logFlag='yes'):
    """
    Same as MATLAB tic and toc functions. Use ty.tic() at the beginning of
    code you want to time and ty.toc() at the end. Once ty.toc() is reached,
    elapsted time will be printed to screen and optionally (by default) written
    to 'log.txt' file.
    """
    import numpy as np
    import time
    if 'startTime_for_tictoc' in globals():
        if logFlag == 'yes':
            log("\n\tElapsed time is "+str(np.round(time.time()-
                           startTime_for_tictoc,decimals=3))+" seconds.")
        else:
            print("\n\tElapsed time is "+
                  str(np.round(time.time()-
                           startTime_for_tictoc,decimals=3))+" seconds.")

    else:
        print("\n\t\tToc: start time not set") 
        
        
##########################################################     
def vdos(vels,tn,num,dt,dn,win,thz):
    """
    This function calculates vibrational density of states from EMD velocity
    data. Intended to be used adjunct to SED code. 
    """
    import numpy as np
    vels = vels.reshape(tn,num*3)
    dos = np.zeros((tn,num*3))
    velsfft = np.fft.fft(vels,axis=0)*dt*dn
    dos = (np.multiply(abs(velsfft),abs(velsfft))/
               np.tile(np.multiply(vels,vels).mean(axis=0),(tn,1))/(tn*dt*dn))
    dosx = gsmooth(dos[:,0::3].mean(axis=1),win,(thz[1]-thz[0])*2*np.pi*1e12)
    dosy = gsmooth(dos[:,1::3].mean(axis=1),win,(thz[1]-thz[0])*2*np.pi*1e12)
    dosz = gsmooth(dos[:,2::3].mean(axis=1),win,(thz[1]-thz[0])*2*np.pi*1e12)
    dos = gsmooth(dos.mean(axis=1),win,(thz[1]-thz[0])*2*np.pi*1e12)

    return [dos, dosx, dosy, dosz]


##########################################################
def log(string,outfile='log.txt',suppress='no',new='no'):
    """
    This function prints output to a file called log.txt and to the screen. 
    Useful for tracking whats happening when submitted using qsub or slurm etc.
    If you don't want to print to screen, enter supress='yes'
    """
    if new == 'no':
        with open(outfile,'a') as fid:
            fid.write(string)
    else:
        with open(outfile,'w') as fid:
            fid.write(string)
    if suppress == 'no':
        print(string)
    
        
