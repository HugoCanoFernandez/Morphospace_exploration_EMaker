
In this directory you can find all the code necesary to reconstruct all the simulations perfomed from the IC2
To create the initial conditions, run de commands:

chmod 777 Make_ensemble_IC2.sh
./Make_ensemble_IC2.sh

The initial conditions for each simulation will appear in the directory /bin
They will be named with the following code:

IC2.dat00006_00001_00006_00001_00001_00043.modi.dat.gz

where:
  - IC2.dat stands for the initial condition IC2
  - the first number  stands for the property or behavior activated in the epithelial animal pole (see the codes in the file ranges_ENSEMBLE.dat)
  - the second number  stands for sign of the regulation of properties/behaviors in epithelial the animal pole (00001 is increasing and 00002 is decreasing)

  - the third number  stands for the property or behavior activated in the epithelial dorsal stripe (see the codes in the file ranges_ENSEMBLE.dat)
  - the fourth number  stands for sign of the regulation of properties/behaviors in epithelial the dorsal stripe (00001 is increasing and 00002 is decreasing)
   
  - the fifth number stands for the property or behavior activated in the mesenchymal animal pole (see the codes in the file ranges_ENSEMBLE.dat)
  - the sixth number  stands for the property or behavior activated in the mesenchymal dorsal stripe (see the codes in the file ranges_ENSEMBLE.dat)

Note that these files are compressed and must be extracted before using them to run the simulations with EMaker.


