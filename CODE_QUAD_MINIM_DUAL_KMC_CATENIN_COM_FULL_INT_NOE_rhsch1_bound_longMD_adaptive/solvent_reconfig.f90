subroutine solvent_reconfig(g_bool)

implicit none

 integer                                   :: i,j,k,l,m,n

 character(len=5),dimension(:),allocatable ::  atomtypesol,restypesol
 real(kind=8) , dimension(:),allocatable   ::  coord_xsol,coord_ysol,coord_zsol,coord_test_x,coord_test_y,coord_test_z
 integer      , dimension(:),allocatable   ::  resnumsol,atomnumsol

 real(kind=8)                              :: diff_x,diff_y,diff_z,diff_tot

 logical                                   :: changebool,changebool2,changeboolean3

 real(kind=8)                              :: diff_x1,diff_y1,diff_z1,diff_x2,diff_y2,diff_z2

 real(kind=8)                              :: box_x,box_y,box_z

 integer                                   :: ntotal,natoms

 integer                                   :: ierr

 real(kind=8)                              :: coord_x,coord_y,coord_z

 integer                                   :: resnum

 character(len=5)                          :: restype,atomtype

 real(kind=8)                              :: vel_x,vel_y,vel_z

 integer                                   :: changenumber,p
 
 real(kind=8)                              :: harvest3

 logical                                   :: g_bool

 integer                                   :: atomnum,out


open(unit=1,file='sol.gro')
read(1,*) ntotal
read(1,*) natoms

allocate(resnumsol(ntotal-natoms),stat=ierr)
allocate(atomtypesol(ntotal-natoms),stat=ierr)
allocate(restypesol(ntotal-natoms),stat=ierr)
allocate(atomnumsol(ntotal-natoms),stat=ierr)
allocate(coord_xsol(ntotal-natoms),stat=ierr)
allocate(coord_ysol(ntotal-natoms),stat=ierr)
allocate(coord_zsol(ntotal-natoms),stat=ierr)

write(*,*) 'h4'


do i=1,ntotal-natoms

           read(1,'(i5,2a5,i5,3f8.3,3f8.4)',end=11,err=11) resnumsol(i), restypesol(i), atomtypesol(i), atomnumsol(i), coord_xsol(i), &
                                             coord_ysol(i), coord_zsol(i), vel_x , vel_y , vel_z

enddo

11 continue

read(1,*) box_x,box_y,box_z


close(unit=1)

allocate(coord_test_x(natoms),stat=ierr)
allocate(coord_test_y(natoms),stat=ierr)
allocate(coord_test_z(natoms),stat=ierr)

if(g_bool == .false.) then

         open(unit=1111,file='minimized3.gro')

endif

if(g_bool == .true.) then

         open(unit=1111,file='minimized.gro')

endif

read(1111,*)
read(1111,*)

do m=1,NATOMS

           read(1111,'(i5,2a5,i5,3f8.3,3f8.4)',end=91111,err=91111) resnum, restype, atomtype, atomnum, coord_test_x(m), &
                                             coord_test_y(m), coord_test_z(m), vel_x , vel_y , vel_z

enddo

91111     continue

close(unit=1111)


if(g_bool == .false. ) then

         call system('/scratch/pee18323/GRO_bin/bin/editconf -f minimized3.gro -c -o minimized3.gro')

         open(unit=1111,file='minimized3.gro')

endif

if(g_bool == .true.) then

         call system('/scratch/pee18323/GRO_bin/bin/editconf -f minimized.gro -c -o minimized.gro')

         open(unit=1111,file='minimized.gro')

endif

         call system ('rm -f traj.trr mdout.mdp ener.edr md.log input.tpr \#minimized.gro.1#')

         call system('rm -f md_cheater2.gro')
!         call system('editconf -f minimized3.gro -o minimized_boxer.gro -d 0.75 -bt cubic')
!         call system('cp -f LOV2.top LOV2gen.top')
!         call system('genbox -cp minimized_boxer.gro -cs spc216.gro -o md_cheater2.gro -p LOV2gen.top')

         open(unit=1112,file='md_cheater2.gro')

         read(1111,*)
         read(1111,*)

         write(1112,*)
         write(1112,*) Ntotal

         do i=1,natoms
       
           read(1111,'(i5,2a5,i5,3f8.3,3f8.4)',end=92015,err=92015) resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z
          write(1112,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z

                                             changebool = .false.

                                             do n=1,ntotal-natoms

                                                if(atomtypesol(n) == '   OW') then

                                                diff_x = coord_x - coord_xsol(n)
                                                diff_y = coord_y - coord_ysol(n)
                                                diff_z = coord_z - coord_zsol(n)

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                

                                                endif

                                                

                                                diff_x = coord_x - coord_xsol(n) - box_x
                                                diff_y = coord_y - coord_ysol(n)
                                                diff_z = coord_z - coord_zsol(n)

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                
                                                endif

                                                
                                                diff_x = coord_x - coord_xsol(n) + box_x
                                                diff_y = coord_y - coord_ysol(n)
                                                diff_z = coord_z - coord_zsol(n)

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                endif

                                                
                                                diff_x = coord_x - coord_xsol(n) 
                                                diff_y = coord_y - coord_ysol(n) - box_y
                                                diff_z = coord_z - coord_zsol(n)

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                
                                                endif

                                                
                                                diff_x = coord_x - coord_xsol(n) 
                                                diff_y = coord_y - coord_ysol(n) + box_y
                                                diff_z = coord_z - coord_zsol(n)

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                
                                                endif

                                                
                                                diff_x = coord_x - coord_xsol(n) 
                                                diff_y = coord_y - coord_ysol(n) 
                                                diff_z = coord_z - coord_zsol(n) -box_z

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                
                                                endif

                                                diff_x = coord_x - coord_xsol(n) 
                                                diff_y = coord_y - coord_ysol(n) 
                                                diff_z = coord_z - coord_zsol(n) + box_z

                                                diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                if(diff_tot.le.0.15) then

                                                         write(*,*) diff_tot,'diff_tot'

                                                         changebool = .true.

                                                         changebool2 = .true.

                                                         changeboolean3 = .false.

                                                         diff_x1 = 0.017
                                                         diff_y1 = -0.09
                                                         diff_z1 = 0.039

                                                         diff_x2 = -0.078
                                                         diff_y2 = 0.043
                                                         diff_z2 = 0.047

                                                         changenumber = n

                                                         exit

                                                
                                                endif

                                             endif

                                             enddo


                                             if(changebool == .true.) then


                                                do

                                                         !       10601  continue

                                                                call random_number(harvest3)

                                                                coord_xsol(changenumber) = harvest3*box_x                   
                                                                                            
                                                         !       if(harvest3*box_x .gt. maxprotx .or. harvest3*box_x .lt. minprotx) then

                                                         !               continue

                                                         !               write(*,*) coord_xsol(changenumber)

                                                         !       else

                                                         !               goto 10601

                                                         !       endif

                                                do   l = 1,10 


                                                             !   10602  continue

                                                                call random_number(harvest3)

                                                                coord_ysol(changenumber) = harvest3*box_y                   
                                                                                            
                                                              !  if(harvest3*box_y .gt. maxproty .or. harvest3*box_y .lt. minproty) then

                                                              !         continue

                                                              !  else

                                                              !         goto 10602

                                                              ! endif  

                                                     do p = 1,10
 
                                                                changebool2 = .true.
                                                             !   10603  continue

                                                                call random_number(harvest3)

                                                                coord_zsol(changenumber) = harvest3*box_z                   
                                                                                            
                                                             !   if(harvest3*box_z .gt. maxprotz .or. harvest3*box_z .lt. minprotz) then

                                                             !           continue

                                                             !   else

                                                             !           goto 10603

                                                             !   endif                 


                                                                coord_xsol(changenumber+1) = coord_xsol(changenumber) - 0.017
                                                                coord_ysol(changenumber+1) = coord_ysol(changenumber) + 0.09
                                                                coord_zsol(changenumber+1) = coord_zsol(changenumber) -0.039

                                                                coord_xsol(changenumber+2) = coord_xsol(changenumber) + 0.078
                                                                coord_ysol(changenumber+2) = coord_ysol(changenumber) - 0.043
                                                                coord_zsol(changenumber+2) = coord_zsol(changenumber) - 0.047

  !                                                              if(coord_xsol(changenumber).gt.maxprotx .and. &
  !                                                                 coord_ysol(changenumber).gt.maxproty .and. &
  !                                                                 coord_zsol(changenumber).gt.maxprotz) then

  !                                                              else if(coord_xsol(changenumber).lt.minprotx .and. &
  !                                                                 coord_ysol(changenumber).lt.minproty      .and. &
  !                                                                 coord_zsol(changenumber).lt.minprotz ) then

                                      !                  write(*,*) 'solvent test'


                                                                          do m=1,NATOMS

                             
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif

                                                                             diff_x = coord_test_x(m)-box_x - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m)+box_x - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m)-box_y - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m)+box_y - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m)-box_z - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber)

                                                                             diff_z = coord_test_z(m)+box_z - coord_zsol(changenumber)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif

                                                                             diff_x = coord_test_x(m)-box_x - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m)+box_x - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m)-box_y - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m)+box_y - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m)-box_z - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+1)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+1)

                                                                             diff_z = coord_test_z(m)+box_z - coord_zsol(changenumber+1)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif

                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif

                                                                             diff_x = coord_test_x(m)-box_x - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m)+box_x - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m)-box_y - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m)+box_y - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m) - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m)-box_z - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                             diff_x = coord_test_x(m) - coord_xsol(changenumber+2)

                                                                             diff_y = coord_test_y(m) - coord_ysol(changenumber+2)

                                                                             diff_z = coord_test_z(m)+box_z - coord_zsol(changenumber+2)

                                                                             diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                                                                             if(diff_tot .lt. 0.3) then

                                                                                changebool2 = .false.


                                                                                exit

                                                                             endif
                                                                           enddo
                                                               do k=1,ntotal - natoms

                                                                  if(k.ne.changenumber.and.k.ne.changenumber+1.and.k.ne.changenumber+2) then

                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif
                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)-box_x
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif
                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)+box_x
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif
                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)-box_y
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)+box_y
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)-box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                       diff_x = coord_xsol(changenumber) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber) - coord_zsol(k)+box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.25) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                       diff_x = coord_xsol(changenumber+1) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif    
                                                                         
                                                                       diff_x = coord_xsol(changenumber+1) - coord_xsol(k)+box_x
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif 

                                                                        diff_x = coord_xsol(changenumber+1) - coord_xsol(k)-box_x
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif 

                                                                         diff_x = coord_xsol(changenumber+1) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)-box_y
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                         diff_x = coord_xsol(changenumber+1) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)+box_y
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                         diff_x = coord_xsol(changenumber+1) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)-box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                         diff_x = coord_xsol(changenumber+1) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+1) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+1) - coord_zsol(k)+box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif

                                                                       diff_x = coord_xsol(changenumber+2) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                 
                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)-box_x
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             
 
                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)+box_x
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                    
                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)-box_y
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                       
                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)+box_y
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                      
                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)-box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                      

                                                                        diff_x = coord_xsol(changenumber+2) - coord_xsol(k)
                                                                       diff_y = coord_ysol(changenumber+2) - coord_ysol(k)
                                                                       diff_z = coord_zsol(changenumber+2) - coord_zsol(k)+box_z

                                                                       diff_tot = sqrt(diff_x**2+diff_y**2+diff_z**2)

                               

                                                                          if(diff_tot .lt. 0.1) then

                                                                          changebool2 = .false.     
                                                                          exit
 
                                                                   
                                                                          endif                                                                             

                                                                      endif
                                                                 
                                                                      enddo

                                                                             if(changebool2 == .true.) then

                                                                          write(*,'(i5,2a5,i5,3f8.3,3f8.4)') resnumsol(changenumber), restypesol(changenumber), atomtypesol(changenumber), atomnumsol(changenumber), coord_xsol(changenumber), &
                                                                              coord_ysol(changenumber), coord_zsol(changenumber)
                                                                          write(*,'(i5,2a5,i5,3f8.3,3f8.4)') resnumsol(changenumber+1), restypesol(changenumber+1), atomtypesol(changenumber+1), atomnumsol(changenumber+1), coord_xsol(changenumber+1), &
                                                                              coord_ysol(changenumber+1), coord_zsol(changenumber+1)
                                                                          write(*,'(i5,2a5,i5,3f8.3,3f8.4)') resnumsol(changenumber+2), restypesol(changenumber+2), atomtypesol(changenumber+2), atomnumsol(changenumber), coord_xsol(changenumber+2), &
                                                                              coord_ysol(changenumber+2), coord_zsol(changenumber+2)
                                                                              write(*,*) 'exit success'

                                                                              changeboolean3 = .true.

                                                                                exit

                                                                              endif

!                                                                          endif

                                                                       enddo

                                                                           if(changeboolean3 == .true.) then

                                                                            write(*,*) 'exit2'

                                                                            exit

                                                                           endif

                                                             enddo

                                                         if(changeboolean3 == .true.) then

                                                             write(*,*) 'exit 3'

                                                             exit

                                                          endif

                                           
                                       enddo
                             
                     endif


        enddo

          
        do i=1,ntotal-natoms

           write(1112,'(i5,2a5,i5,3f8.3,3f8.4)') resnumsol(i), restypesol(i), atomtypesol(i), atomnumsol(i), coord_xsol(i), &
                                             coord_ysol(i), coord_zsol(i), vel_x , vel_y , vel_z
          
        enddo             

        write(1112,*) box_x , box_y , box_z


92015    continue

close(unit=1111)
close(unit=1112)



         call system('rm -f fullmd.tpr')


deallocate(coord_test_x,stat=ierr)
deallocate(coord_test_y,stat=ierr)
deallocate(coord_test_z,stat=ierr)

write(*,*) 'h7'

if(allocated(resnumsol)) then
deallocate(resnumsol,stat=ierr) 
if(ierr .ne. 0 ) then
write(*,*) 'dealloc resnumsol'
continue
endif
endif
if(allocated(restypesol)) then
deallocate(restypesol,stat=ierr)
if(ierr .ne. 0 ) then
write(*,*) 'dealloc restypesol'
continue
endif
endif
if(allocated(atomtypesol)) then
deallocate(atomtypesol,stat=ierr)
if(ierr.ne.0) then
write(*,*) 'dealloc atomtypesol'
continue
endif
endif
if(allocated(atomnumsol)) then
deallocate(atomnumsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc atomnumsol'
continue
endif
endif
if(allocated(coord_xsol)) then
deallocate(coord_xsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_xsol'
continue
endif
endif
if(allocated(coord_ysol)) then
deallocate(coord_ysol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_ysol'
continue
endif
endif
if(allocated(coord_zsol)) then
deallocate(coord_zsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_zsol'
continue
endif
endif

call system('/scratch/pee18323/GRO_bin/bin/grompp -f minim_final_A.mdp -c md_cheater2.gro -p LOV2gen.top -o t.tpr ')

call system('/scratch/pee18323/GRO_bin/bin/mdrun -v -s t.tpr -c md_cheater2.gro -nice 0')

call system('rm -f #*# out')

end subroutine solvent_reconfig