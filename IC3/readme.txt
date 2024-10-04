
In this directory you can find all the code necesary to reconstruct all the simulations perfomed from the IC2
To create the initial conditions, run de commands:

chmod 777 Make_ensemble_IC3.sh
./Make_ensemble_IC3.sh

The initial conditions for each simulation will appear in the directory /bin
They will be named with the following code:

IC3.dat00001_00043_00043.modi.dat.gz

where:
  - IC3.dat stands for the initial condition IC3
  - the three numbers correspond to the properties/behaviors activaded in each territory (see the codes in the file ranges_ENSEMBLE.dat)

Note that these files are compressed and must be extracted before using them to run the simulations with EMaker.


