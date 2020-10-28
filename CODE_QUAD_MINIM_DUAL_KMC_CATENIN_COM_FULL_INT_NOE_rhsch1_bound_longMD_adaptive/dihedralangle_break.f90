subroutine BREAK_DIHEDRALANGLE_SCAN(aminoresnum,resnumber,scanbool,rot_change,scansteps,energy_scan,g_bool,intevil,intevil2,intevil3)

implicit none

double precision , parameter :: PI = 3.141592654

integer         :: i,k,n,scanstep,scannumber,select

integer         :: natoms

integer         :: count1


integer         :: resnumber

logical         :: scanbool

integer,dimension(:),allocatable      :: resnum,atomnum
real(kind=8),dimension(:),allocatable :: coord_x,coord_y,coord_z

character(len=5),dimension(:),allocatable :: restype,atomtype 

real(kind=8)                          :: vel_x,vel_y,vel_z 

real(kind=8),dimension(1000,4,100)     :: energy_scan

integer                               :: ierr,atomnum1,atomnum2

real(kind=8)                          :: rotangle

real(kind=8)                          :: atom_1_x,atom_1_y,atom_1_z

real(kind=8)                          :: atom_2_x,atom_2_y,atom_2_z

real(kind=8)                          :: atom_3_x,atom_3_y,atom_3_z

real(kind=8)                          :: atom_4_x,atom_4_y,atom_4_z

real(kind=8)                          :: diff_x1,diff_y1,diff_z1
real(kind=8)                          :: diff_x2,diff_y2,diff_z2
real(kind=8)                          :: diff_x3,diff_y3,diff_z3
real(kind=8)                          :: diff_x4,diff_y4,diff_z4

real(kind=8)                          :: kreuz_x1,kreuz_y1,kreuz_z1
real(kind=8)                          :: kreuz_x2,kreuz_y2,kreuz_z2

real(kind=8)                          :: angle_dih,diff_tot,angle_real,energy_val

integer                               :: scansteps

real(kind=8)                          :: rot_change,x,y,z,paramdist

real(kind=8)                          :: box_x,box_y,box_z,n1,n2,n3

integer       :: execsteps,kmc_step

real(kind=8)  :: time_end,rate

logical       :: execbool,g_bool,new_break

integer,dimension(1000) :: aminoresnum

integer                :: icounter,selectamino,intevil,intevil2,intevil3

logical                :: alert
! stepwise rotation

do i=1,scansteps


write(*,*) i,'step'

if(i==1) then
open(unit=111,file='test.ndx')
open(unit=87,file='start.gro')
open(unit=12,file='1.ndx')

write(12,'(a9)') '[ group ]'

read(87,*)
read(87,*) natoms

n = 1

allocate(restype(natoms),stat=ierr)
allocate(resnum(natoms),stat=ierr)
allocate(atomtype(natoms),stat=ierr)
allocate(atomnum(natoms),stat=ierr)
allocate(coord_x(natoms),stat=ierr)
allocate(coord_y(natoms),stat=ierr)
allocate(coord_z(natoms),stat=ierr)


do k=1,natoms

           read(87,'(i5,2a5,i5,3f8.3,3f8.4)',end=9000,err=9000) resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z             

           if(resnum(k) == aminoresnum(resnumber)) then

              write(12,*) atomnum(k)

           endif

enddo

9000 continue

read(87,*) box_x,box_y,box_z

close(unit=87)
close(unit=12)

call system ('cat 1.ndx system.ndx > 2.ndx')

endif

if(i.ge.2) then

open(unit=1,file='minimized.gro')

read(1,*)
read(1,*) natoms

do k=1,natoms

           read(1,'(i5,2a5,i5,3f8.3,3f8.4)',end=900,err=900) resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z             

enddo

900 continue

read(1,*) box_x,box_y,box_z

close(unit=1)

endif


open(unit=1,file='n1.gro')

 write(1,*)
 write(1,*) '12'

 count1 = 1

do k=1,natoms

  if(resnum(k) == aminoresnum(resnumber)) then

        write(1,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z

  endif

enddo

  write(1,*) box_x,box_y,box_z    

close(unit=1)
! rotation around 4 atoms


 count1 = 1

do k=1,natoms

if(resnum(k) == aminoresnum(resnumber) .and. count1 == 1) then

write(111,*) 'yes'

! normalization to z-axis

atom_1_x = coord_x(k+intevil)
atom_1_y = coord_y(k+intevil)
atom_1_z = coord_z(k+intevil)

atom_2_x = coord_x(k+intevil2)
atom_2_y = coord_y(k+intevil2)
atom_2_z = coord_z(k+intevil2)

atom_3_x = coord_x(k+intevil3)
atom_3_y = coord_y(k+intevil3) 
atom_3_z = coord_z(k+intevil3)


 count1 = count1 + 1

endif

enddo

diff_x1 = atom_2_x - atom_3_x
diff_y1 = atom_2_y - atom_3_y
diff_z1 = atom_2_z - atom_3_z

n1 = diff_x1 / sqrt(diff_x1**2+diff_y1**2+diff_z1**2)
n2 = diff_y1 / sqrt(diff_x1**2+diff_y1**2+diff_z1**2)
n3 = diff_z1 / sqrt(diff_x1**2+diff_y1**2+diff_z1**2)


write(111,*) n1,n2,n3


do k=1,natoms

 if(resnum(k) == aminoresnum(resnumber)) then

  coord_x(k) = coord_x(k) - atom_1_x
  coord_y(k) = coord_y(k) - atom_1_y 
  coord_z(k) = coord_z(k) - atom_1_z

 endif

enddo



do k=1,natoms

if(resnum(k) == aminoresnum(resnumber)) then

! normalization to z-axis

  coord_x(k) = coord_x(k) * (cos(rot_change*PI/180)+n1**2*(1-cos(rot_change*PI/180))) + coord_x(k)*((n2*n1)*(1-cos(rot_change*PI/180))+n3*sin(rot_change*PI/180)) &
               & + coord_x(k)*(n3*n1*(1-cos(rot_change*PI/180))-n2*sin(rot_change*PI/180)) 

  coord_y(k) = coord_y(k) * ( n1*n2*(1-cos(rot_change*PI/180)) - n3*(sin(rot_change*PI/180)) + cos(rot_change*PI/180) + n2**2*(1-cos(rot_change*PI/180)) &
               & + n3*n2*(1-cos(rot_change*PI/180))+n1*sin(rot_change*PI/180) )

  coord_z(k) = coord_z(k) * ( n1*n3*(1-cos(rot_change*PI/180)) + n2*(sin(rot_change*PI/180)) + n2*n3*(1-cos(rot_change*PI/180))-n1*sin(rot_change*PI/180) &
               & + cos(rot_change*PI/180) + n3**2*(1-cos(rot_change*PI/180))) 

endif



enddo

do k=1,natoms

 if(resnum(k) == aminoresnum(resnumber)) then

  coord_x(k) = coord_x(k) + atom_1_x
  coord_y(k) = coord_y(k) + atom_1_y 
  coord_z(k) = coord_z(k) + atom_1_z

 endif

enddo

open(unit=1,file='n.gro')

 write(1,*)
 write(1,*) '12'

 count1 = 1

do k=1,natoms

  if(resnum(k) == aminoresnum(resnumber)) then

        write(1,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z

  endif

enddo

  write(1,*) box_x,box_y,box_z    

close(unit=1)


open(unit=11,file='intermediate.gro')

write(11,*)
write(11,*) Natoms

do k=1,natoms

            write(11,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z    
enddo

write(11,*) box_x,box_y,box_z

close(unit=11)

call system('cp intermediate.gro minimized.gro')


! measure of dihedral angle

if(i==1) then

     call system('cp n1.gro first.gro')

endif

 count1 = 1

do k=1,natoms

  if(resnum(k) == aminoresnum(i) .and. count1 == 1) then

        atom_1_x = coord_x(k)
        atom_1_y = coord_y(k)
        atom_1_z = coord_z(k)

        atom_2_x = coord_x(k+3)
        atom_2_y = coord_y(k+3)
        atom_2_z = coord_z(k+3)        
                       
        atom_3_x = coord_x(k+4)
        atom_3_y = coord_y(k+4)
        atom_3_z = coord_z(k+4) 

        atom_4_x = coord_x(k+6)
        atom_4_y = coord_y(k+6)
        atom_4_z = coord_z(k+6)

        count1 = count1+1

  endif

enddo

diff_x1 = atom_1_x - atom_2_x
diff_y1 = atom_1_y - atom_2_y
diff_z1 = atom_1_z - atom_2_z

diff_x2 = atom_2_x - atom_3_x                        
diff_y2 = atom_2_y - atom_3_y
diff_z2 = atom_2_z - atom_3_z

kreuz_x1 = diff_y1*diff_z2 - diff_z1*diff_y2
kreuz_y1 = diff_z1*diff_x2 - diff_x1*diff_z2
kreuz_z1 = diff_x1*diff_y2 - diff_y1*diff_x2

diff_x3 = atom_2_x - atom_3_x
diff_y3 = atom_2_y - atom_3_y
diff_z3 = atom_2_z - atom_3_z

diff_x4 = atom_3_x - atom_4_x                        
diff_y4 = atom_3_y - atom_4_y
diff_z4 = atom_3_z - atom_4_z

kreuz_x2 = diff_y3*diff_z4 - diff_z3*diff_y4
kreuz_y2 = diff_z3*diff_x4 - diff_x3*diff_z4
kreuz_z2 = diff_x3*diff_y4 - diff_y3*diff_x4

                                                                 
diff_tot  = sqrt((kreuz_x1)**2+(kreuz_y1)**2+(kreuz_z1)**2)*sqrt((kreuz_x2)**2+(kreuz_y2)**2+(kreuz_z2)**2) 

angle_dih = 180/PI*acos((kreuz_x1*kreuz_x2 + kreuz_y1*kreuz_y2 + kreuz_z1*kreuz_z2)/diff_tot) 


call minimization_of_protein

if(scanbool == .true.) then

scanstep   = i
scannumber = 1


write(*,*) rot_change,'rot_change'
write(*,*) intevil   ,'intevil'
write(*,*) intevil2  ,'intevil2'
write(*,*) intevil3  ,'intevil3'

  call scan_energy_of_protein(energy_val,resnumber,scannumber,scanstep,alert)
 
  energy_scan(resnumber,1,i) = energy_val 

       if(alert == .true.) then

          exit

       endif
!  write(*,*) energy_scan(resnumber,1,i) ,'hier'

  if(i==scansteps) then

                 new_break = .true.

       else

                 new_break = .false.

   endif
          
   if(g_bool == .false.) then
      
           write(resnumber*1001,*) i,energy_val

           if(new_break == .true.) then

                write(resnumber*1001,*) ' '

           endif

   endif

   if(g_bool == .true.) then
      
           write(resnumber*1002,*) i,energy_val

           if(new_break == .true.) then

                write(resnumber*1002,*) ' '

           endif

   endif
!call write_file(scannumber,scanstep,energy_val,resnumber,selectamino,kmc_step,select,execsteps,g_bool,aminoresnum,time_end,rate,new_break)

endif


enddo



deallocate(restype,stat=ierr)
deallocate(resnum,stat=ierr)
deallocate(atomtype,stat=ierr)
deallocate(atomnum,stat=ierr)
deallocate(coord_x,stat=ierr)
deallocate(coord_y,stat=ierr)
deallocate(coord_z,stat=ierr)


end subroutine BREAK_DIHEDRALANGLE_SCAN 

subroutine minimization_of_protein

   call system('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_dih.mdp -c minimized.gro -p intermediate.top -n 2.ndx -o run.tpr')
   call system('/loctmp/pee18323/GRO_bin/bin/mdrun -v -s run.tpr -nice 0 -e md_ener.edr -c minimized.gro')
   call system('rm -f #*# out')

   write(*,*) 'minimized'
end subroutine minimization_of_protein

subroutine scan_energy_of_protein(energy_val,resnumber,scannumber,scanstep,alert)


character(len=1)                   :: finder

real(kind=8)                       :: a,b

integer                            :: i

real(kind=8)                       :: energy_val

logical                            :: alert
       
                call system('./energyfilereader.sh')

                call system('rm -f ./#*#')

                    ! energy read

                open(unit=4,file='t.xvg')

                do

                       read(4,'(a1)',end=2,err=2) finder

                       if(finder.ne.'@' .and. finder .ne. '#') then

                                        exit

                       endif

                enddo

2 continue

                do

                       read(4,*,end=1,err=1) a,b

                enddo

                1 continue

                close(unit=4)

                energy_val = b

                write(*,*) 'energy',energy_val

                if(energy_val .gt. 0) then

                   alert = .true.
 
                else

                   alert = .false.

                endif


end subroutine scan_energy_of_protein