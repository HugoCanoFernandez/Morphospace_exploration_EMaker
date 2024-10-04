
In this directory you can find all the code necesary to reconstruct all the simulations perfomed from the IC1
To create the initial conditions, run de commands:

chmod 777 Make_ensemble_IC1.sh
./Make_ensemble_IC1.sh

The initial conditions for each simulation will appear in the directory /bin
They will be named with the following code:

IC1.dat00006_00010_00001_00002_.modi.dat.gz

where:
  - IC1.dat stands for the initial condition IC1
  - the first number (00006 in this case) stands for the property or behavior activated in the animal pole (see the codes in the file ranges_ENSEMBLE.dat)
  - the second number (00010 in this case) stands for the property or behavior activated in the dorsal stripe (see the codes in the file ranges_ENSEMBLE.dat)
  - the third number (00001 in this case) stands for sign of the regulation of properties/behaviors in the animal pole (00001 is increasing and 00002 is decreasing)
  - the fourth number (00002 in this case) stands for sign of the regulation of properties/behaviors in the dorsal stripe (00001 is increasing and 00002 is decreasing)

Note that these files are compressed and must be extracted before using them to run the simulations with EMaker.


