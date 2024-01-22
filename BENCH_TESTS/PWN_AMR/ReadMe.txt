How to run a setup:

> python $PLUTO/DIR/setup.py --with-chombo: MPI=TRUE

> make

> mpirun -n NN ./pluto &


------ SETUP1 ---->
This setup is though to run on a small number of processors (tested on 64). It has a limited grid with base level (128^3 cells) + 3 Adaptive Mesh Refinement levels (for an equivalent resolution of 1024^3 cells at maximum resolution). 
The simulation runs up to a time Tfinal, and outputs are produced in .hdf5 format every Tfinal/10.


------ SETUP2 ---->
This setup is though to run on a large number of processors (tested on 1100-2000). It has a full extended grid with base level (136^3 cells) + 7 Adaptive Mesh Refinement levels (for an equivalent resolution of 17408^3 cells at maximum resolution). 
The simulation runs up to a time Tfinal, and outputs are produced in .hdf5 format every Tfinal/10.

