#!/bin/bash
#SBATCH --mail-type=end          
#SBATCH --mail-user=hugo.cano@telefonica.net 
#SBATCH --job-name="compilEMaker"
#SBATCH --qos=debug              
#SBATCH --time 00:05:00          
#SBATCH --error=error.dat        
#SBATCH --output=output.dat                   
#SBATCH --ntasks=1   

cd src/core


echo 'descending'
pwd

aleas=aleas.mod.f90

#making object files
echo 'making object files'

gfortran -S -w -fallow-invalid-boz  OpenGL_gl.f90 OpenGL_glu.f90 OpenGL_glut.f90
gfortran -w -c $aleas
gfortran -w -c general.mod.f90 
gfortran -w -c geompack3.f90 gnuplotter.mod.f90
gfortran -w -c shell.mod.f90 neighboring.mod.f90
gfortran -w -c genetic.mod.f90
gfortran -w -c energy.mod.f90 io.mod.f90 
gfortran -w -c polarization.mod.f90
gfortran -w -c biomechanic_pola.mod.f90 
gfortran -w -c biomechanic.mod.f90 death.mod.f90 pola.mod.f90 ecm.mod.f90 growth.mod.f90
gfortran -w -c ic.mod.f90 mitosis.mod.f90
gfortran -w -c single_node.mod.f90
gfortran -w -c nexus.mod.f03
gfortran -w -c inicial.mod.f90
gfortran -w -c model.mod.f90
gfortran -w -c automaticon.mod.f90
gfortran -w -c elli_MN4.f90



#linking
echo 'linking EMaker'

gfortran -w -O2 -fexceptions -fno-underscoring -fcheck=all gnuplotter.mod.o polarization.mod.o biomechanic_pola.mod.o aleas.mod.o general.mod.o neighboring.mod.f90 genetic.mod.o energy.mod.o shell.mod.o io.mod.o pola.mod.o mitosis.mod.o growth.mod.o death.mod.o single_node.mod.o ic.mod.o ecm.mod.o nexus.mod.o biomechanic.mod.o model.mod.o inicial.mod.o automaticon.mod.f90 geompack3.f90 elli_MN4.o -o EMaker  

echo 'linking ens.e'

gfortran general.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 ens0_mes.f90 -o ens.e

#cleaning
echo 'cleaning'
rm *.s *.o *.mod

cd -
echo 'moving to bin'
mv src/core/EMaker bin
mv src/core/ens.e bin

echo 'executables installed in bin/'

cd bin
echo 'Making the ensemble'
./ens.e IC1.dat ranges_ENSEMBLE.dat

echo 'finished'


