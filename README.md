# Axonal-PathFinding
This repository contains the ABAQUS VUMAT subroutine code for modeling axonal growth and pathfinding during brain development. Included in the folder "Abaqus_inp_sample" is a sample geometrical input file. To launch the Abaqus job, use the following command:

Abaqus job=G_6_Case_0.inp user=vumat_axongrowth_G_6.f cpus=32 interactive

You may adjust the number of CPUs as needed. Before launching the job, ensure that Fortran is properly linked with ABAQUS.
