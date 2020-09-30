# EQdyna  

A parallel Finite Element software to model earthquake spontaneous dynamic ruptures. The software is designed to use high-performance computing. It is written in FORTRAN 90 and MPI.

## Features

EQdyna is based on Finite Element Method of solid mechanics and features 
1) 3D under-integrated hexahedral elements stabilized by hourglass control; 
2) explicit central difference time integration for efficient and accurate modeling of wave propagation; 
3) the traction-at-split-node technique to model earthquake faults or fractures;
4) the perfectly matched layer absorbing boundary;
5) the coarse-grained Q model for attenuation of seismic waves;
6) various forms of friction laws including slip-weakening, time-weakening, rate- and state- friction with the aging law, slip law and dynamic weakening;
7) 3D MPI. 

## Initializing environments and compiling on ADA, Texas A&M University, HPC.

```
Module load intel/2018a
make # run makefile to complie EQdyna. An excutable file eqdyna-hyb will be generated. 
```

## Batch file preparation on ADA

Here is an example of the batch file, runeqdyna.txt to submit jobs on ADA. The key word is #BSUB and contexts after # is commented. 

```
#BSUB -J job.name

## send stderr and stdout to the same file 
#BSUB -o mpitest.%J

## login shell to avoid copying env from login session
## also helps the module function work in batch jobs
#BSUB -L /bin/bash

## 4 hours of walltime ([HH:]MM)
#BSUB -W 04:00

#BSUB -n 72 -M 2700    # 72: total number of cores allocated to the job
#BSUB -R rusage[mem=2700] span[ptile=12] 
# mem: the memory allocated to each core
# span[ptile]: the number of cores used on each computing node

# load intel toolchain
ml intel/2018a
ml
export OMP_NUM_THREADS=1

mpirun -np 72 ./eqdyna-hyb # excuting EQdyna
```

## Job submission on ADA
After submission, the job will be quened and excuated when resources are avalible.
```
Bsub <runeqdyna.txt
```

## Postporcessing

1) run rtp3DMPI.m to plot rupture time contour.
It coverts results frt.txt* to cplot.txt.  
2) run timeanalysis.m to analyze computational times used by various modules.
It processes timeinfo* files.

## History
3.1
3.1.1
3.1.2
3.2.1
4.0
4.1.1
4.1.2
4.1.3
4.2

## Ajustable parameters for Version 4.1.3
### in globalvar.f90
```
 C_elastic: 0: plastic model; stress tensors should be assigned to every element in the volume in meshgen.f90
            1: elastic model; normal and shear stresses should be assigned to on-fault split nodes;
 
 C_Nuclea:  0: disabled;
            1: allow artificial nucleation;
 
 friclaw:   1: slip weakening;
            2: time weakening;
            3: rate- and state- friction with the aging law;
            4: rate- and state- friction with strong rate-weakening;
 
 C_Q:       0: attenuation disabled;
            1: allow frequency-dependent/-independet attenuation of seismi waves; currently only works with uniform element                 size (rat == 1.0) and elastic models;
 
 C_hg:      0：visocus hourglass control;
            1: KF78 hourglass control; recommended;
 
 nPML:      6: default value; number of layers in the Perfetly Matched Layers absorbing boundaries;
 
 rat:       1.025: bufferring ratio
 
 dx:        element size;
 
 dis4uniF/dis4uniB: the number of elements with dx edge along y axis (-y/+y);
 
 critt0:    temporal scale in time weakening for the nucleation phase in slip-weakening law;
 
 srcrad0:   the radius of the nucleation patch for the slip-weakening law;
 
 vrupt0:    forced rupture velocity inside the nucleation patch for the slip-weakening law;
 rhow:      fluid density in the volume in kg/m3;
 bulk:      bulk modulus for Drucker-Prager Plasticity;
 coheplas:  cohesion for Drucker-Prager Plasticity in Pa;
 tv:        temporal scale for viscoplasticity helping smoothing the pulling back to the yielding curve in s;
 xsource/ysource/zsource:
            hypotencer locations along x, y, z directions, respectively, in km for atificial nucleation;
 vmaxPML:   maximum P wave velocity in the whole model used for PML layers.
 term:      termination time in s.
 dt:        time step in s.
 npx/npy/npz:
            MPI discretization along x/y/z directions, respectively.            
```

### in meshgen.f90
```
  fric:   1: mus
          2: mud
          3: D0 in m
          4: cohesion in Pa
          5: time for forced rupture in s; void;
          6: pore pressure in Pa for the elastic case;
          7: initial normal stress in Pa for the elastic case; negative for compressional;
          8: initial shear stress in Pa for the elastic case; positive for right-lateral slipp;
          9: a in rsf;
         10: b in rsf;
         11: L in rsf;
         12: V0 in rsf;
         13: f0 in rsf;
         14: fw in rsf with strong rate-weakening;
         15: Vw in rsf with strong rate-weakening;
         16: initial slip-rate along the x axis;
         17: initial slip-rate along the y axis;
         18: initial slip-rate along the z axis;
         19: initial slip-rate magnitude;
```

```
  mat:    1: P wave velocity in m/s;
          2: S wave velocity in m/s;
          3: material density in km/m3;
```
## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Benchun Duan**-*Initial work*-[Earthquake Numerical Simulation Lab](https://geogeo.tamu.edu/people/faculty/duanbenchun)
* **Dunyu Liu**
* **Bin Luo**

See also the list of [contributors]() who participated in this project.

## License

To be determined (An exmaple may look like: This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details)

## Reference
Duan, B. (2010). Role of initial stress rotations in rupture dynamics and ground
  motion: A case study with implications for the Wenchuan earthquake, J. Geophys.
  Res. 115.

Duan, B. (2012). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles
  of a possible subducting seamount, J. Geophys. Res. 117

Duan, B., D. Liu, and A. Yin (2017). Seismic shaking in the North China Basin
  expected from ruptures of a possible seismic gap, Geophys. Res. Lett. 44 4855-4862.
  
* Hat tip to anyone whose code was used
* Inspiration
* etc
