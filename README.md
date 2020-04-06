# Finite_element_HPC
Finite element code written in C++ in parallel using OpenMP
Used BLAS and SCALAPACK for offloading linear algebra operations

The project is a deliverable for High Performance Computing Course at Imperial College London.

To run the programme make;
task1: for static single processor
task2: for dynamic single processor
task3: for static multi-core
task4: for dynamic multi-core

the user will be prompted for the following inputs
length, Elements, Area, Moment of Inertia, Young Modulus, time, timesteps, density


Files are output in a text file and post-processed in python scripts not included here.
