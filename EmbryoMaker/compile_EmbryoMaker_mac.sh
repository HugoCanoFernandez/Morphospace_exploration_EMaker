#!/bin/bash


#echo 'pre'
#pwd


cd src/core



echo 'descending'
pwd

aleas=aleas.mod.f90

#making object files
echo 'making object files'

gfortran-mp-4.6 -S -w OpenGL_gl.f90 OpenGL_glu.f90 OpenGL_glut.f90
gfortran-mp-4.6 -w -c $aleas
gfortran-mp-4.6 -w -c general.mod.f90 
gfortran-mp-4.6 -w -c geompack3.f90 gnuplotter.mod.f90
gfortran-mp-4.6 -w -c shell.mod.f90 neighboring.mod.f90
gfortran-mp-4.6 -w -c genetic.mod.f90
gfortran-mp-4.6 -w -c energy.mod.f90 io.mod.f90 
gfortran-mp-4.6 -w -c biomechanic.mod.f90 death.mod.f90 pola.mod.f90 ecm.mod.f90 growth.mod.f90
gfortran-mp-4.6 -w -c ic.mod.f90 mitosis.mod.f90
gfortran-mp-4.6 -w -c single_node.mod.f90
gfortran-mp-4.6 -w -c nexus.mod.f03
gfortran-mp-4.6 -w -c inicial.mod.f90
gfortran-mp-4.6 -w -c editor.mod.f90 model.mod.f90
gfortran-mp-4.6 -w -c drawer.mod.f90 automaticon.mod.f90
gfortran-mp-4.6 -w -c elli.f90


#linking
echo 'linking'

gfortran-mp-4.6 -w -fexceptions -fno-underscoring -fbounds-check gnuplotter.mod.o aleas.mod.o general.mod.o neighboring.mod.f90 genetic.mod.o energy.mod.o shell.mod.o io.mod.o pola.mod.o mitosis.mod.o growth.mod.o death.mod.o single_node.mod.o ic.mod.o ecm.mod.o nexus.mod.o biomechanic.mod.o model.mod.o inicial.mod.o editor.mod.o drawer.mod.f90 automaticon.mod.f90 geompack3.f90 elli.o -o EMaker -framework OpenGL -framework glut

#cleaning
echo 'cleaning'
rm *.s *.o *.mod

cd .. && cd ..
mv src/core/EMaker bin/


echo 'executables installed in bin/'


