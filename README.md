JeDi
====

#### Cloning JeDi from GitHub

To clone JeDi to your machine from GitHub:
```bash
git clone https://github.com/rpavlick/JeDi.git
```

#### Compiling JeDi

To compile JeDi, use the ```makejedi``` script.
```bash

> ./makejedi --help

usage: ./makejedi [options]

  options:
    --build  <buildname>  - provide a name for this build (default is git commit SHA-1)
    --compiler <ifort,pgf90,gcc>  - choose compiler (default is ifort) 
    --fast   - compile with optimizations
    --debug  - compile with debug options
    --nompi  - compile without MPI libraries for single processor
    --clean  - removes old compile files
    --netcdf - compile with netcdf libraries

makejedi version:heads/master-0-g79705a1-dirty
Copyright © 2014-2015 Ryan Pavlick
This is free software and comes with ABSOLUTELY NO WARRANTY
For more information, see: https://github.com/rpavlick/JeDi
```

As an example, let's compile JeDi with gfortran (now part of gcc) and name the build ```testbuild```:
```bash

>./makejedi --clean --compiler gcc --build testbuild

************************************************************
/home/rpavlick/JeDi/makejedi --clean --compiler gcc --build testbuild
------------------------------------------------------------
git revision: heads/master-0-g79705a1-dirty
host: graphene.jpl.nasa.gov: Linux 2.6.32-696.6.3.el6.x86_64 x86_64 x86_64
date: Wed May  2 09:44:35 PDT 2018
************************************************************
 WARNING: Git repository is dirty! 
COMPILER=gcc
BUILD=testbuild
BUILD_DIR=/home/rpavlick/JeDi/build/testbuild
MODULES_FILE=/home/rpavlick/JeDi/build/testbuild/modules
Saving git diff to /home/rpavlick/JeDi/build/testbuild/testbuild.diff

Cleaning old build files

Copying source files to build directory
Saving modules files to /home/rpavlick/JeDi/build/testbuild/modules
Currently Loaded Modulefiles:
  1) rocks-openmpi

GNU Fortran compiler found at /usr/bin/gfortran
Using built-in specs.
Target: x86_64-redhat-linux
Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-languages=c,c++,objc,obj-c++,java,fortran,ada --enable-java-awt=gtk --disable-dssi --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-1.5.0.0/jre --enable-libgcj-multifile --enable-java-maintainer-mode --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --disable-libjava-multilib --with-ppl --with-cloog --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
Thread model: posix
gcc version 4.4.7 20120313 (Red Hat 4.4.7-18) (GCC) 

MPI Fortran compiler found at /opt/openmpi/bin/mpif90
Using built-in specs.
Target: x86_64-redhat-linux
Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-languages=c,c++,objc,obj-c++,java,fortran,ada --enable-java-awt=gtk --disable-dssi --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-1.5.0.0/jre --enable-libgcj-multifile --enable-java-maintainer-mode --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --disable-libjava-multilib --with-ppl --with-cloog --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
Thread model: posix
gcc version 4.4.7 20120313 (Red Hat 4.4.7-18) (GCC) 

creating version.f90
compiling ./globe/globe_mod_func.F90
compiling ./globe/globe_mod_flux.F90
compiling ./globe/globe_mod_stat.F90
compiling ./globe/globe_mod.F90
compiling ./globe/globe.F90
compiling ./globe/globe_sub.F90
compiling ./globe/globe_fio.F90
compiling ./globe/globe_surf.F90
compiling ./globe/globe_functions.F90
compiling ./globe/globe_mpimod.F90
compiling ./jam/jam_mod.F90
compiling ./jam/jam.F90
compiling ./jam/jam_fio.F90
compiling ./jam/jam_sub.F90
compiling ./jam/jam_globe.F90
compiling ./jam/jam_carbon.F90
compiling ./jedi/jedi_mod.F90
compiling ./jedi/jedi_mod_dyn.F90
compiling ./jedi/jedi.F90
compiling ./jedi/jedi_sub.F90
compiling ./jedi/jedi_fio.F90
compiling ./jedi/jedi_globe.F90
compiling ./jedi/jedi_dyn.F90
compiling ./jedi/jedi_opti.F90
compiling ./jedi/jedi_spec.F90
compiling ./jedi/jedi_kristin.f90
compiling ./nosoil/nosoil_globe.f90
linking jedi.x

Success! ☺
jedi.x is now in /home/rpavlick/JeDi/build/testbuild

saving the output of this script...
/home/rpavlick/JeDi/build/testbuild/makejedi_testbuild.log
```

As you can see in the output from ```makejedi``` above, the script has
1. logged the command line arguments used when calling the script
2. logged the current git revision of the JeDi repository you are working with
3. logged the hostname and operating system of the machine you are compiling on
4. logged the date and time
5. displays a warning if the Git repository is "dirty" (i.e. there are uncommitted changes anywhere in the repository)
6. created a directory ```build\testbuild```
7. created a git diff file ```build\testbuild\testbuild.diff``` with those uncommitted changes (if any) 
8.* created a modules file ```build/testbuild/modules``` that will used later by ```runjedi``` to load modules necessary to run JeDi
9. compiled JeDi with the gcc compiler, found at ```JeDi/build/testbuild```
10. saved a log file of the ```makejedi``` output to ```build/testbuild/makejedi_testbuild.log```.

Notes about compiling with ```makejedi```
* If your cluster doesn't use the `module` system or PBS queueing system, then you will need to modify `makejedi` and `runjedi` to suit those needs.
* You will likely need to make modifications to `makejedi` due the specifics of your compiler and the specific modules used on your cluster. 

#### Running JeDi
To run JeDi, use the ```runjedi``` script:

```bash
> ./runjedi --help

Usage: ./runjedi <configuration file> [ parameter1=value parameter2=value ... ]

runjedi version:heads/master-0-g79705a1-dirty
Copyright © 2014-2015 Ryan Pavlick
This is free software and comes with ABSOLUTELY NO WARRANTY
For more information, see: https://github.com/rpavlick/JeDi
```

```runjedi``` requires a configuration file. Below is an example configuration file ```testrun.cfg``` that is included in the JeDi git repository.

```bash
> cat testrun.cfg

    # experiment name
    EXP=testrun
    RUN=1

    ### restart setting
    RESTART=0
    ### directory containing restart files
    RESTART_DIR=

    ### location of model executable
    BUILD=""

    BASE_DIR=$PWD
    ### daily meteorological forcing directory
    FORCING_DIR="/home/rpavlick/nobackup0/FORCING/ISIMIP/0.5degree_amazon/HadGEM2-ES/spinup"

    ### namelist directory
    NAMELIST_DIR="${BASE_DIR}/data/testrun"

    ### directory containing surface description files
    SURFACE_DIR="${BASE_DIR}/data/0.5degree/amazon"

    ### directory containing species parameter file
    SPECPARM_FILE="${BASE_DIR}/data/testrun/jedi_specparms.txt"

    ### atmospheric CO2 file
    PCO2_FILE=""

    ### jedi parameter table file
    PARTAB_FILE="${BASE_DIR}/data/jedi_partab"

    ### pbs cluster options
    NODES=1
    NCPUS=12
    QUEUE=mediumq
    WALLTIME=12:00:00
```


Input files

paw.srv
- unit plant available water (mm water per mm soil)
- calculated following Kleidon and Heimann 2000, Batjes 1996

landsea.srv
- land mask (1=land/0=ocean/glacier)

tas.srv
- 2m air temperature
- Kelvin

rlns.srv
- net longwave radiation at the surface
- W/m^2
 
rsds.srv
- downward shortwave radiation at the surface
- W/m^2

pr.srv
- daily total precipitation
- kg/m^2/s

jedi_specparm.txt
- pseudorandom trait values generated by create_jedi_specparm


Postprocessing

jedi_partab
landsea.nc
landfraction.nc
