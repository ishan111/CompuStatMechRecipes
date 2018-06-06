# CompuStatMechRecipes
This will contain all the code I write for molecular simulation. Only getting started here, so I don't expect it to be a professional, efficient  and readable code just yet.
Hopefully, the final versions of our code will be professional, efficient, readable and capable of being run in parallel in various interfaces

28-08-2017 UPDATE: LenardJones_2D_MonteCarlo_v0.m A basic program to generate random samples according to boltzmann statistics, assuming Canonical ensemble, and 6-12-LJ potential interaction between particles in a 2D square box and Periodic Boundary conditions.

29-08-2017 UPDATE: modifiable_2D_MonteCarlo_v0.m, modifiableLJ_2D_MonteCarlo_v0.m, modifiableHD_2D_MonteCarlo_v0.m have been added to the toy program folder. They are more modifiable in that the particle interaction can be modified by modifying the "Potential_function" anonymous function. modifiableLJ and modifiableHD are already modified for 6-12 Lenard-Jones and Hard Disk potentials respectively

23-10-2017 UPDATE: GC_TMMCv0.m is added. Its a toy GC-TMMC code, verified by similar plots of density of states of macro variable number of particles, obtained by histogram and by TMMC.

20-11-2017 UPDATE: NEW REWRITTEN CODE ADDED FOR GC-TMMC. This code is more easily modified.

21-11-2017 UPDATE: SPECIAL CODE  ADDED FOR HARD DISK SIMULATION. This turned out to be very straightforward and small code due to prior experiences and special nature of code.

28-11-2017 UPDATE: SPECIAL CODE ADDED FOR LJ FLUID SIMULATION. I am now using the single decision statement to do everything. No compartmentalization here.

20-11-2017 UPDATE: MODULAR UNBIASED CODE AND ONE FILE CODE FOR UNBIASED SIMULATION ADDED too much fluctuation in the particle histogram still. Concluded that scripts are much slower to be called than functions on MATLAB

01-01-2018 UPDATE: FORTRAN INCOMPLETE CODE many functions still need to be defined, and too many variables are global as yet.

06-06-2017 UPDATE: MATLAB ISING MODEL TMMC CODE WAS INCLUDED
