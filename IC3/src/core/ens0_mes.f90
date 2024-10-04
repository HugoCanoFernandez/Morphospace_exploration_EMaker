program rewrite


! gfortran general.mod.f90 genetic.mod.f90 io.mod.f90 geompack3.f90 ens0_mes.f90 -o ens.e

use io
use general
use genetic
!use neighboring

 implicit none
  character*400 :: input,output, rangfile, when, used
  character*500 :: inpu1
  character*152 :: name_dat2
  integer :: ihc, iihc, jhc, jjhc, lhc, khc,kkhc, phc, shc, ord,ord2, nhc, mhc  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
  real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
  integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  integer, dimension(1:nga+1) :: ens_rembeh, ens_rembeh_mes
  integer, dimension(1:5) :: kadhpos
  integer :: kadhcode, icode, jcode, kcode, ratio
  real*8, dimension(1:5) :: max_glim, min_glim
  real*8 :: intensa, intensb, intensc, intensd, intense, b_max
  character( len= 5 ) :: ithc, jthc, kthc
  
  
  
!******************************************************************************************************!
call getarg(1,input)    !! input file
call getarg(2,rangfile) !! file with the ranges "ranges.dat"

print*,trim(input)

!!HERE YOU READ THE RANGES FOR THE ACTIVATION OF PROPERTIES AND BEHAVIORS
call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh)  

!call fdate(when)
used="cp "//trim(rangfile)//" used_ranges.dat"!//trim(when)//".dat"
print*, used
call system(used)

b_max=max_glim(5)

ens_rembeh=0
ens_rembeh(1:nga)=rembeh(1:nga) !!use this vector that includes a final position for NO ACTIVATION
ens_rembeh(52)=1                !! WARNING!! this is only activated when nparam_per_node+8 is activated (see EXCEPTIONS)
ens_rembeh(1)=1                 !! WARNING!! this is the kadh, not using it yet 

ens_rembeh_mes=1
ens_rembeh_mes(1)=0
!ens_rembeh_mes(3)=0
!ens_rembeh_mes(4)=0
!ens_rembeh_mes(5)=0
!ens_rembeh_mes(6)=0
!ens_rembeh_mes(7)=0
ens_rembeh_mes(nparam_per_node+2)=0
ens_rembeh_mes(nparam_per_node+8)=0
ens_rembeh_mes(nga+1)=0
!!Activation intensities
intensa=0.0d0; intensb=0.0d0

!!MAIN LOOP

ord=0


      do iihc=1, nga+1                                  !! mesenchyme genes
         do jjhc=1, nga+1                               !! mesenchyme genes
            do kkhc=1, nga+1
              do mhc=1,1                                  !! intensity 1=medium, 2=high 
                 do nhc=1,1                               !! intensity 1=medium, 2=high
                    do shc=1,1
                       if(iihc>3.and.iihc.le.7)cycle
                       if(kkhc>3.and.kkhc.le.7)cycle
                       if(iihc==nga+1.and.jjhc==nga+1.and.kkhc==nga+1)cycle  !! do not run no activation + no activation
                       if(ens_rembeh_mes(iihc)==1)cycle      !! filter unused positions of e matrix 
                       if(ens_rembeh_mes(jjhc)==1)cycle      !! filter unused positions of e matrix 
                       if(ens_rembeh_mes(kkhc)==1)cycle      !! filter unused positions of e matrix 
                     
                       call iniread
                       call readsnap(input)                            !! Read inputfile     
            
                      if (mhc==1)then;  intensc=0.50d0; else; intensc=1.0d0; endif !! set up intensity of activation
                      if (nhc==1)then;  intensd=0.50d0; else; intensd=1.0d0; endif !! set up intensity of activation
                      if (shc==1)then;  intense=0.50d0; else; intense=1.0d0; endif !! set up intensity of activation
                  

                       if (iihc.le.nga.and.iihc>7)then
                          gen(4)%e(iihc)=max_elim(iihc)*intensc
                       endif

                       if (jjhc.le.nga.and.jjhc>7)then
                          gen(5)%e(jjhc)=max_elim(jjhc)*intensd
                       endif                  

                       if (kkhc.le.nga.and.kkhc>7)then
                          gen(6)%e(kkhc)=max_elim(kkhc)*intense
                       endif       
                           
                                  
                                 !!!!!!!!!!!!!EXCEPTIONS!!!!!!!!!!!!!!!!!!
                     
                                 !!PCP
                  
                       if (iihc==nparam_per_node+8)then            !! if nparam_per_node+8 (PCP) is activated, we have to activate nparam_per_node+2 (elongation)
                          gen(4)%e(iihc)=max_elim(iihc)
                          gen(3)%e(52)=max_elim(52)*intensc
                       endif
                  
                       if (jjhc==nparam_per_node+8)then
                           gen(5)%e(jjhc)=max_elim(jjhc)
                           gen(3)%e(52)=max_elim(52)*intensd
                       endif
                     
                       if (kkhc==nparam_per_node+8)then
                           gen(6)%e(kkhc)=max_elim(kkhc)
                           gen(3)%e(52)=max_elim(52)*intense
                       endif
                              

                       !! kadh
                       kadhpos=0; kadhcode=0
                       if (iihc==1)then; kadh(1,1)=b_max*intensc; kadhpos(1)=1; endif
                       if (iihc==2)then; kadh(1,2)=b_max*intensc; kadh(2,1)=b_max*intensc; kadhpos(4)=1; endif
                       if (iihc==3)then
                          kadh(1,2)=b_max*intensc; kadh(2,1)=b_max*intensc; kadh(1,1)=b_max*intensc
                          kadhpos(1)=1; kadhpos(4)=1
                       endif
                                 
                       if (kkhc==1)then; kadh(3,3)=b_max*intense; kadhpos(3)=1; endif
                       if (kkhc==2)then; kadh(3,2)=b_max*intense; kadh(2,3)=b_max*intense; kadhpos(5)=1; endif
                       if (kkhc==3)then
                          kadh(3,2)=b_max*intense; kadh(2,3)=b_max*intense; kadh(3,3)=b_max*intense
                          kadhpos(3)=1; kadhpos(5)=1
                       endif
                                 
                       if (jjhc==1)then; kadh(2,2)=b_max*intensd; kadhpos(2)=1; endif
                       if (jjhc==2)then
                          if (kadh(1,2)>0.0d0)cycle
                          kadh(1,2)=b_max*intensd; kadh(2,1)=b_max*intensd; kadhpos(4)=1
                       endif
                       if (jjhc==3)then
                          if (kadh(2,3)>0.0d0)cycle
                          kadh(2,3)=b_max*intensd; kadh(3,2)=b_max*intensd; kadhpos(5)=1
                       endif
                       if (jjhc==4)then
                          if (kadh(2,1)>0.0d0)cycle
                          kadh(2,1)=b_max*intensd; kadh(1,2)=b_max*intensd; kadh(2,2)=b_max*intensd
                          kadhpos(2)=1; kadhpos(4)=1
                       endif
                       if (jjhc==5)then
                          if (kadh(2,3)>0.0d0)cycle
                          kadh(2,3)=b_max*intensd; kadh(3,2)=b_max*intensd; kadh(2,2)=b_max*intensd
                          kadhpos(2)=1; kadhpos(5)=1
                       endif
                       if (jjhc==6)then
                          if (kadh(2,1)>0.0d0.or.kadh(2,3)>0.0d0)cycle
                          kadh(2,1)=b_max*intensd; kadh(1,2)=b_max*intensd; kadh(2,3)=b_max*intensd;
                          kadh(3,2)=b_max*intensd
                          kadhpos(4)=1; kadhpos(5)=1
                       endif  
                       if (jjhc==7)then
                          if (kadh(2,1)>0.0d0.or.kadh(2,3)>0.0d0)cycle
                          kadh(2,1)=b_max*intensd; kadh(1,2)=b_max*intensd
                          kadh(2,3)=b_max*intensd; kadh(3,2)=b_max*intensd; kadh(2,2)=b_max*intensd
                          kadhpos(2)=1; kadhpos(4)=1; kadhpos(5)=1
                       endif
                     
                       do phc=1,5
                          kadhcode=kadhcode+kadhpos(phc)*10**phc
                       enddo
                       kadhcode=kadhcode/10
                     
                       if(iihc.le.7)then; icode=1; else; icode=iihc; endif
                       if(jjhc.le.7)then; jcode=1; else; jcode=jjhc; endif
                       if(kkhc.le.7)then; kcode=1; else; kcode=kkhc; endif
                  
                       !!!!!!!!!!!!!OUTPUT!!!!!!!!!!!!!!!!!!   (uncomment this section  to make the output)
                       ord=ord+1
                     
                       write( ithc, '(I5.5)' )    iihc
                        write( jthc, '(I5.5)' )   jjhc
                         write( kthc, '(I5.5)' )  kkhc
                       print*, ord,  iihc, gen(4)%e(icode), jjhc, gen(5)%e(jcode), kkhc, gen(6)%e(kcode), ithc
                                 call iniio
                                 call writesnap
	                         !output is the name of the new file, which is ne name of the original one + jhc
                                 output=trim(input)//trim(ithc)//"_"//trim(jthc)//"_"//trim(kthc)//".modi.dat"

	                         !This actually makes the new file
                                 name_dat2=trim("cp fort.1 "//adjustl(output)) 
                                   call system(name_dat2)
                                 name_dat2=trim("gzip "//adjustl(output))  
                                  call system(name_dat2)
                                 !  print*, "THE NEW RECOMBINED INPUT FILE IS: ",output  
                  
                   enddo
                 enddo
              enddo
            enddo
         enddo
      enddo




end program


