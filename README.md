# phonon-transmission
Harmonic phonon transmission calculations from molecular dynamics trajectories

# harmonic
Within harmonic are the codes 'compute_transmission.py' and 'subroutines.py'. 'subroutines.py' is a module that 'compute_transmission.py' calls. The module does most of the work, the compute code is where you provide simulations parameters that the code needs. The it mainly reads data in and plots it. The subdirectory scripts contains another directory, 'matlab', that has a matlab version of the code and some other functions that are useful. ORIGINAL_CODES contains the original code that my codes are based off. Also, compactify.cpp shoulde by compiled and used (USAGE: a.x vels.dat vels.dat.compact) to 'compactify' the velocity data into a more easy to read (for the code!) version. The subdirectory 'lammps_input_files' has lammps input scripts that produce the data for the transmission calculation.

