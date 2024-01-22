# PlutoBenchmarks
A few benchmarks for the PLUTO MHD code. The benchmarks have been tested and assembled by S. Orlando at INAF-OAPa.

PLUTO is a freely-distributed software for the numerical solution of mixed hyperbolic/parabolic systems of partial differential equations.

The software is developed at the [Dipartimento di Fisica, Torino University](http://www.unito.it/) in a joint collaboration with [INAF, Osservatorio Astronomico di Torino](http://www.oato.inaf.it/index.php?lang=en) and the [SCAI Department of CINECA](http://www.hpc.cineca.it/).

Please see the [HomePage of Pluto](http://plutocode.ph.unito.it/)

System Requirements
-------------------------
PLUTO can run on most platforms but some software prerequisites
must be met, depending on the specific configuration you intend to
use. The minimal set to get PLUTO running on a workstation with a
static grid (no AMR) requires Python, a C compiler and the make
utility. These are usually installed by default on most Linux/Unix
platforms. 

The minimum software requirements for PLUTO applications are the following 
- Python (> 2.0)
- C compiler
- GNU make
- MPI library

The Chombo library is required for computations making use of
Adaptive Mesh Refinement. The HDF5 library is required for
I/O with the Chombo library and may also be used with the static
grid version of the code. In case Adaptive Mesh Refinement is used,
the software requirements are:
- C, C++ and Fortran compilers;
- the MPI library (for parallel runs).
- GNU make
- the Chombo library1 (version 3.2 is recommended);
- the HDF5 library2 (version < 1.8.14 is recommended);
- the chombo-3.2-patch.tar.gz provided with the PLUTO distribution, which replaces some of the library source files.

For the installation of the Chombo library, refer to the userguide pages 123-126

=====================================================================================

Set environment variables
-------------------------
Set the environment variable PLUTO_DIR to point to the code directory.
Depending on your shell (e.g. tcsh or bash) use either one of

     ~/MyWorkDir> export PLUTO_DIR=/home/user/PLUTO # If you’re using the bash shell;
     ~/MyWorkDir> setenv PLUTO_DIR /home/user/PLUTO # If you’re using the tcsh shell;


=====================================================================================
BENCHMARK SHOCK_CLOUD (BENCH_TESTS/SHOCK_CLOUD)

1) Change directory to test problem SHOCK_CLOUD under BENCHMARK/BENCH_TESTS

     ~/MyWorkDir> cd BENCHMARK/BENCH_TESTS/SHOCK_CLOUD

2) Run the Python script using

     ~/MyWorkDir> python $PLUTO_DIR/setup.py

   select “Change makefile”, choose a suitable makefile (e.g.
   Linux.mpicc.defs) and press enter.

   (It is possible to edit and change the configuration file for
   mpicc PLUTO/Config/Linux.mpicc.defs to specify the compiler
   and/or attributes in the compilation)

   All the information relevant to the specific problem should now be
   stored in four files: init.c (assigns initial condition and
   user-supplied boundary conditions), pluto.ini (sets the number of
   grid zones, Riemann solver, output frequency, etc.), definitions.h
   (specifies the geometry, number of dimensions, reconstruction, time
   stepping scheme, and so forth) and the makefile.

3) Exit from the main menu (“Quit” or press ’q’) and type

       ~/MyWorkDir> make

   to compile the code. 

4) Run the code

       ~/MyWorkDir> mpirun [...] ./pluto [args]

   [...] are options given to MPI, such as number of processors,
   etc, while [args] are command line options specific to PLUTO
   (see Table 1.3 in the user guide). For example, 

       ~/MyWorkDir> ./pluto -restart 5 -maxsteps 840

   will restart from the 5-th double precision output file and stop
   computation after 840 steps.

   Useful arguments:

       -restart n          Restart computations from the n-th output file
       -maxsteps n         Stop computations after n steps
       -no-write           Do not write data to disk
       -no-x#par           Do not perform parallel domain decomposition
                           along the x# direction (where # = 1, 2, 3)

5) It is possible to change the size of the numerical grid by editing
   the file pluto.ini and changing the [Grid]. Examples:

   (for a grid 512^3):
      X1-grid    1    0.0    512    u    1.0
      X2-grid    1    0.0    512    u    1.0
      X3-grid    1    0.0    512    u    1.0

   (for a grid 1024^3):
      X1-grid    1    0.0    1024    u    1.0
      X2-grid    1    0.0    1024    u    1.0
      X3-grid    1    0.0    1024    u    1.0


=====================================================================================
BENCHMARK SN_EXPLOSION (BENCH_TESTS/SN_EXPLOSION)

In this setup, the initial conditions are stored in external files
under the subdirectories EXT and CSM. The files in these directories
are read, and the data are interpolated in the numerical grid.

1) Change directory to test problem SN_EXPLOSION under BENCHMARK/BENCH_TESTS

     ~/MyWorkDir> cd BENCHMARK/BENCH_TESTS/SN_EXPLOSION

2) Run the Python script using

     ~/MyWorkDir> python $PLUTO_DIR/setup.py

   select “Change makefile”, choose a suitable makefile (e.g.
   Linux.mpicc.defs) and press enter.

   (It is possible to edit and change the configuration file for
   mpicc PLUTO/Config/Linux.mpicc.defs to specify the compiler
   and/or attributes in the compilation)

   All the information relevant to the specific problem should now be
   stored in five files: init.c (assigns initial condition and
   user-supplied boundary conditions), mean_mol_weight.c (calculate
   the mean molecular weight), pluto.ini (sets the number of grid
   zones, Riemann solver, output frequency, etc.), definitions.h
   (specifies the geometry, number of dimensions, reconstruction,
   time stepping scheme, and so forth) and the makefile.

3) Exit from the main menu (“Quit” or press ’q’) and type

       ~/MyWorkDir> make

   to compile the code.

4) Run the code

       ~/MyWorkDir> mpirun [...] ./pluto [args]

   [...] are options given to MPI, such as number of processors,
   etc, while [args] are command line options specific to PLUTO
   (see Table 1.3 in the user guide). For example,

       ~/MyWorkDir> ./pluto -restart 5 -maxsteps 840

   will restart from the 5-th double precision output file and stop
   computation after 840 steps.

   Useful arguments:

       -restart n          Restart computations from the n-th output file
       -maxsteps n         Stop computations after n steps
       -no-write           Do not write data to disk
       -no-x#par           Do not perform parallel domain decomposition
                           along the x# direction (where # = 1, 2, 3)

5) It is possible to change the size of the numerical grid by editing
   the file pluto.ini and changing the [Grid]. Examples:

   (for a grid 256^3):
      X1-grid    1    0.0    256    u    1.0
      X2-grid    1    0.0    256    u    1.0
      X3-grid    1    0.0    256    u    1.0

   (for a grid 1024^3):
      X1-grid    1    0.0    1024    u    1.0
      X2-grid    1    0.0    1024    u    1.0
      X3-grid    1    0.0    1024    u    1.0

=====================================================================================
BENCHMARK PWN (BENCH_TESTS/PWN_AMR)

The steps to run the setup are basically the same illistrated before.
A summary is provided below

How to run a setup:

     ~/MyWorkDir> cd BENCHMARK/BENCH_TESTS/PWN_AMR

     ~/MyWorkDir> python $PLUTO/DIR/setup.py --with-chombo: MPI=TRUE

     ~/MyWorkDir> make

     ~/MyWorkDir> mpirun [...] ./pluto [args]


------ SETUP1 ---->
This setup is though to run on a small number of processors (tested
on 64). It has a limited grid with base level (128^3 cells) + 3
Adaptive Mesh Refinement levels (for an equivalent resolution of
1024^3 cells at maximum resolution).
The simulation runs up to a time Tfinal, and outputs are produced
in .hdf5 format every Tfinal/10.


------ SETUP2 ---->
This setup is though to run on a large number of processors (tested
on 1100-2000). It has a full extended grid with base level (136^3
cells) + 7 Adaptive Mesh Refinement levels (for an equivalent
resolution of 17408^3 cells at maximum resolution).
The simulation runs up to a time Tfinal, and outputs are produced
in .hdf5 format every Tfinal/10.
