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
  integer :: ihc, iihc, jhc, jjhc, lhc, khc, phc, shc, ord,ord2, nhc, mhc  !!>> HC 24-11-2020 counters (we use i for the gene that is going to mutate and j for the other gene in the interaction)
  real*8, dimension(1:nga) :: max_elim, min_elim  !!>> HC 27-11-2020 These vectors store the limits of e matrix 
  integer, dimension(1:nga) :: rembeh !>> HC 28-11-2020 Unused cellular behaviors/properties CHECK DIMENSION MANUALLY
  integer, dimension(1:nga+1) :: ens_rembeh, ens_rembeh_mes
  real*8, dimension(1:5) :: max_glim, min_glim
  real*8 :: intensa, intensb, intensc, intensd, b_max
  character( len= 5 ) :: ithc, jthc, kthc, lthc, mthc, nthc, pthc, qthc,rthc  
  
  
  
!******************************************************************************************************!
call getarg(1,input)    !! input file
call getarg(2,rangfile) !! file with the ranges "ranges.dat"

print*,trim(input)

!!HERE YOU READ THE RANGES FOR THE ACTIVATION OF PROPERTIES AND BEHAVIORS
call read_rang(rangfile, max_elim, min_elim, max_glim, min_glim, rembeh)  

call fdate(when)
used="cp "//trim(rangfile)//" used_ranges"//trim(when)//".dat"
call system(used)

b_max=max_glim(5)

ens_rembeh=0
ens_rembeh(1:nga)=rembeh(1:nga) !!use this vector that includes a final position for NO ACTIVATION
ens_rembeh(52)=1                !! WARNING!! this is only activated when nparam_per_node+8 is activated (see EXCEPTIONS)
ens_rembeh(1)=1                 !! WARNING!! this is the kadh, not using it yet 

ens_rembeh_mes=1
ens_rembeh_mes(1)=0
ens_rembeh_mes(nparam_per_node+2)=0
ens_rembeh_mes(nparam_per_node+4)=0
ens_rembeh_mes(nparam_per_node+8)=0
ens_rembeh_mes(nga+1)=0
!!Activation intensities
intensa=0.0d0; intensb=0.0d0

!!MAIN LOOP

ord=0
do ihc=1,nga+1                                          !! All spaces in matrix e + nothing (=no activation)                                                        
   do jhc =1, nga+1                                     !! All spaces in matrix e + nothing  
      do iihc=1, nga+1                                  !! mesenchyme genes
         do jjhc=1, nga+1                               !! mesenchyme genes
            do lhc=2,2                                        !! intensity 1=medium, 2=high
               do khc=2,2                                     !! intensity 1=medium, 2=high
!                  do mhc=1,2                                  !! intensity 1=medium, 2=high 
!                     do nhc=1,2                               !! intensity 1=medium, 2=high
                        do phc=1,2                                  !! sign 1=positive range 2=negative range
                           do shc=1,2                               !! sign 1=positive range 2=negative range
                              if(ihc==nga+1 .and. jhc==nga+1.and. iihc==nga+1.and.jjhc==nga+1)cycle  !! do not run no activation + no activation
                              if(ens_rembeh(ihc)==1)cycle           !! filter unused positions of e matrix     
                              if(ens_rembeh(jhc)==1)cycle           !! filter unused positions of e matrix 
                              if(ens_rembeh_mes(iihc)==1)cycle      !! filter unused positions of e matrix 
                              if(ens_rembeh_mes(jjhc)==1)cycle      !! filter unused positions of e matrix 
                              if (ihc.le.nga)then                             !! filter properties/behaviors with no negative range
                                 if(min_elim(ihc).ge.0.0d0.and.phc==2)cycle
                              endif
                              if (jhc.le.nga)then                             !! filter properties/behaviors with no negative range
                                 if(min_elim(jhc).ge.0.0d0.and.shc==2)cycle
                              endif
                  
                              call iniread
                              call readsnap(input)                            !! Read inputfile
               
            
                              if (lhc==1)then;  intensa=0.50d0; else; intensa=1.0d0; endif !! set up intensity of activation
                              if (khc==1)then;  intensb=0.50d0; else; intensb=1.0d0; endif !! set up intensity of activation
                           !   if (mhc==1)then;  intensc=0.50d0; else; intensc=1.0d0; endif !! set up intensity of activation
                          !    if (nhc==1)then;  intensd=0.50d0; else; intensd=1.0d0; endif !! set up intensity of activation
                  
                              if (ihc.le.nga)then                       
                                 if (phc==1)then                           !! give the activation to the e matrix with the right sign
                                    gen(1)%e(ihc)=max_elim(ihc)*intensa
                                 else
                                    gen(1)%e(ihc)=min_elim(ihc)*intensa
                                 endif
                              endif
                  
                              if (jhc.le.nga)then
                                 if (shc==1)then
                                    gen(2)%e(jhc)=max_elim(jhc)*intensb
                                 else
                                    gen(2)%e(jhc)=min_elim(jhc)*intensb
                                 endif
                              endif

                              if (iihc.le.nga.and.iihc.ne.1)then
                                 gen(4)%e(iihc)=max_elim(iihc)!*intensc
                              endif

                              if (jjhc.le.nga.and.jjhc.ne.1)then
                                 gen(5)%e(jjhc)=max_elim(jjhc)!*intensd
                              endif                  
                  
                  
                  
                              !!!!!!!!!!!!!EXCEPTIONS!!!!!!!!!!!!!!!!!!
                  
                              !!PCP
                  
                              if (ihc==nparam_per_node+8)then            !! if nparam_per_node+8 (PCP) is activated, we have to activate nparam_per_node+2 (elongation)
                                 gen(1)%e(ihc)=max_elim(ihc)
                                 gen(3)%e(52)=max_elim(52)!*intensa
                              endif
                  
                              if (jhc==nparam_per_node+8)then
                                 gen(2)%e(jhc)=max_elim(jhc)
                                 gen(3)%e(52)=max_elim(52)!*intensb
                              endif
                              
                              !! kadh
                              if(iihc==1.and.jjhc.ne.1)then
                                gen(4)%e(1)=1
                                if(allocated(kadh))deallocate(kadh)
                                allocate(kadh(1:1, 1:1))
                                ntipusadh=1
                                kadh(1,1)=b_max!*intensc
                              endif
                              
                              if(jjhc==1.and.iihc.ne.1)then
                                gen(5)%e(1)=1
                                if(allocated(kadh))deallocate(kadh)
                                allocate(kadh(1:1, 1:1))
                                ntipusadh=1
                                kadh(1,1)=b_max!*intensd
                              endif
                              
                              if(iihc==1.and.jjhc==1)then
                                gen(5)%e(1)=1
                                if(allocated(kadh))deallocate(kadh)
                                allocate(kadh(1:2,1:2))
                                ntipusadh=2
                                kadh(1,1)=b_max!*intensc
                                kadh(2,2)=b_max!*intensd
                              endif
                              
                  
                              !!!!!!!!!!!!!OUTPUT!!!!!!!!!!!!!!!!!!   (uncomment this section  to make the output)
                              ord=ord+1
                              !if (ihc.le.nga.and.jhc.le.nga)then
                              !   print*, ord, ihc, gen(1)%e(ihc), jhc, gen(2)%e(jhc)
                              !elseif(ihc.le.nga)then
                              !   print*, ord, ihc, gen(1)%e(ihc), jhc, 0.0d0
                              !elseif(jhc.le.nga)then
                               !  print*, ord, ihc, 0.0d0, jhc, gen(2)%e(jhc)
                              !endif
                              if(iihc==1.and.jjhc==1)then
                                if(iihc==1.and.jjhc.ne.1)then
                                  print*, ord, ihc, gen(1)%e(ihc), jhc, gen(2)%e(jhc), iihc, kadh(1,1), jjhc, gen(5)%e(jjhc)
                                endif
                                if(jjhc==1.and.iihc.ne.1)then
                                  print*, ord, ihc, gen(1)%e(ihc), jhc, gen(2)%e(jhc), iihc, gen(4)%e(iihc), jjhc, kadh(1,1)
                                endif
                                if(jjhc==1.and.iihc==1)then
                                  print*, ord, ihc, gen(1)%e(ihc), jhc, gen(2)%e(jhc), iihc, kadh(1,1), jjhc, kadh(2,2)
                                endif
                              else
                                 print*, ord, ihc, gen(1)%e(ihc), jhc, gen(2)%e(jhc), iihc, gen(4)%e(iihc), jjhc, gen(5)%e(jjhc)
                              endif
                              write( ithc, '(I5.5)' ) ihc
                              write( jthc, '(I5.5)' ) jhc
                              write( kthc, '(I5.5)' ) iihc
                              write( mthc, '(I5.5)' ) jjhc
                              write( nthc, '(I5.5)' ) phc
                              write( pthc, '(I5.5)' ) shc
                     
                              call iniio
                              call writesnap
	                      !output is the name of the new file, which is ne name of the original one + jhc
output=trim(input)//trim(ithc)//"_"//trim(nthc)//"_"//trim(jthc)//"_"//trim(pthc)//"_"//trim(kthc)//"_"//trim(mthc)//".modi.dat"

	                      !This actually makes the new file
                              name_dat2=trim("cp fort.1 "//adjustl(output)) 
                                call system(name_dat2)
                               name_dat2=trim("gzip "//adjustl(output))  
                                call system(name_dat2)
                              !  print*, "THE NEW RECOMBINED INPUT FILE IS: ",output   
                  
!                           enddo
!                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
enddo



end program


