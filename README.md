# FindMole
Find molecular structrues or fragments in the energetical matrial trajectories, mainly suport the lammps software


To use "findmole", you need cite my paper:
    "F. Guo et al., J. Phys. Chem. A 2012, 116, 3514−3520"



The running windows are like this:


 * Input in file name:
equ.lammpstrj
 * Coordinate Format:
 * 1. Lammmpstrj
 * 2. xyz
1
 * In put time interval between trajectory frames (ps):
5.0
**********************************************************************
      Cell parameters (Angstroms/Degrees):
**********************************************************************
      a=     61.2429  alpha=     90.0000
      b=     86.9832  beta =    106.6220
      c=     92.5457  gamma=     90.0000
**********************************************************************
Totally read in partical number: 31104
Totally read in  5184 C   atoms...
Totally read in  5184 H   atoms...
Totally read in 10368 O   atoms...
Totally read in 10368 N   atoms...
* The total mass in the simulation box:
*      378463.09989166  g/mole
* The total density of the system
*           0.80115071  g/mole/A^3
*           1.33025605      g/cm^3
**********************************************************************
 1. analysis MS among all the trajectorys
 2. analysis MS from one of the trajectorys
 3. analysis the Mean squared displacement of the trajectorys
 4. compute the Mean lattice parameters of the system
2
 * please input the step number of the trajectorys to analysis:
100000
 * The last frame!
 * exit now!
