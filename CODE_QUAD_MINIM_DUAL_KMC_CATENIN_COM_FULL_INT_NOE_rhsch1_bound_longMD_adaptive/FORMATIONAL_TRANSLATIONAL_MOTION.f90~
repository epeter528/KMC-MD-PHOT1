subroutine FORMATIONAL_TRANSLATIONAL_MOTION(aminoresnum,number_of_scans,resnumber,scanbool,energy_scan,g_bool,maxi,modulonumber)

implicit none

integer                 :: i,j,k,l,m,n

integer ,dimension(1000) :: aminoresnum

integer,dimension(:),allocatable      :: resnum,atomnum
real(kind=8),dimension(:),allocatable :: coord_x,coord_y,coord_z


character(len=5),dimension(:),allocatable :: restype,atomtype

integer                               :: natoms,ierr

integer                               :: scansteps,ende,resid,resatomnumber,modulonumber

integer,dimension(100)                :: hbmap_donor,hbmap_hydrogen,hbmap_acceptor

integer,dimension(100)                :: hbmap_donor2,hbmap_hydrogen2,hbmap_acceptor2,maximum

integer,dimension(100)                :: hbmap_donor3,hbmap_hydrogen3,hbmap_acceptor3

character(len=23)                     :: readchar

logical                               :: ex,bool,scanbool,g_bool,new_break

real(kind=8)                          :: acc_x,acc_y,acc_z
real(kind=8)                          :: don_x,don_y,don_z,maxim

integer                               :: number_of_scans

real(kind=8)                          :: expconst,energy_val

real(kind=8)                          :: d_x,d_y,d_z,d_a,box_x,box_y,box_z

real(kind=8) ,dimension(100)          :: scalekey

integer                               :: resnumber,maxi,p

real(kind=8),dimension(1000,4,100)     :: energy_scan

real(kind=8),dimension(1000,100)       :: energy_scan_break

integer                               :: scannumber,ende2

integer       :: execsteps,kmc_step

real(kind=8)  :: time_end,rate

logical       :: execbool,alert

real(kind=8)  :: vel_x,vel_y,vel_z,sc

real(kind=8),dimension(1000) :: new_x,new_y,new_z

real(kind=8),dimension(1000) :: don_x_res,don_y_res,don_z_res

integer                :: icounter,selectamino,scanstep,select

!aminoresnum(1) = 102

!number_of_scans = 5

!expconst = 0.2

!resnumber = 1

!open(unit=111,file='test.xvg')
open(unit=87,file='start.gro')
open(unit=12,file='1.ndx')
open(unit=13,file='3.ndx')

write(12,'(a10)')'[ group1 ]'
write(12,'(a9)') '[ group ]'
write(13,'(a10)') '[ system ]'

read(87,*)
read(87,*) natoms

n = 1

write(*,*) 'hey form 1'

allocate(restype(natoms),stat=ierr)
allocate(resnum(natoms),stat=ierr)
allocate(atomtype(natoms),stat=ierr)
allocate(atomnum(natoms),stat=ierr)
allocate(coord_x(natoms),stat=ierr)
allocate(coord_y(natoms),stat=ierr)
allocate(coord_z(natoms),stat=ierr)


write(*,*) 'hey form 2'

do k=1,natoms

           read(87,'(i5,2a5,i5,3f8.3,3f8.4)',end=9000,err=9000) resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z             


           if(resnum(k) == aminoresnum(resnumber)) then

              write(12,*) atomnum(k)

           endif

           if(resnum(k) .ne. aminoresnum(resnumber)) then

              write(13,*) atomnum(k)

           endif
enddo

9000 continue

read(87,*) box_x,box_y,box_z

close(unit=87)
close(unit=12)
close(unit=13)


call system ('cat 1.ndx 3.ndx > 2.ndx')


call system('./H_search_new2_int.sh')

inquire(file='table_hbond_B.ndx',exist=ex) 


open(unit=334,file='table_hbond_B.ndx')


k = 0

do

read(334,'(a23)',end=1,err=1) readchar

!write(*,'(a23)') readchar

if(readchar == '[ hbonds_group-system ]') then

      exit

endif

enddo

1 continue


do

     k=k+1
     read(334,*,end=111,err=111) hbmap_donor3(k) , hbmap_hydrogen3(k) , hbmap_acceptor3(k)

    
enddo



111 continue

ende=k-1
 
close(unit = 334)

call system('./H_search_new_int.sh')

inquire(file='table_hbond_B.ndx',exist=ex) 


open(unit=334,file='table_hbond_B.ndx')


k = 0

do

read(334,'(a23)',end=2,err=2) readchar

!write(*,'(a23)') readchar

if(readchar == '[ hbonds_group-system ]') then

      exit

endif

enddo

2 continue


do

     k=k+1
     read(334,*,end=112,err=112) hbmap_donor2(k) , hbmap_hydrogen2(k) , hbmap_acceptor2(k)

    
enddo



112 continue
 
close(unit = 334)
 
ende2=k-1

l = 1

do i=1,ende

                bool = .true.

       do k=1,ende2

              if(hbmap_donor2(k) == hbmap_donor3(i) .and. hbmap_hydrogen2(k) == hbmap_hydrogen3(i) .and. hbmap_acceptor2(k) ==  hbmap_acceptor3(i)) then

                 bool = .false.

              endif

       enddo

       if(bool == .true.) then

                hbmap_donor(l) = hbmap_donor3(i)
             hbmap_hydrogen(l) = hbmap_hydrogen3(i)
             hbmap_acceptor(l) = hbmap_acceptor3(i)

             l = l+1

       endif

enddo

ende = l-1


if(ex == .false.) then

        ende = 0

endif

   
do i=1,ende

         do k=1,natoms 

              if(atomnum(k) == hbmap_donor(i)) then

                     don_x = coord_x(k)
                     don_y = coord_y(k)
                     don_z = coord_z(k)

              endif

              if(atomnum(k) == hbmap_acceptor(i)) then

                     acc_x = coord_x(k)
                     acc_y = coord_y(k)
                     acc_z = coord_z(k)

              endif

         enddo

           D_x = acc_x - don_x
           D_y = acc_y - don_y
           D_z = acc_z - don_z

           D_A = sqrt(D_x**2 + D_y**2 + D_z**2)

           EXPCONST = D_A - 0.3 
          
           scalekey(i)  = -EXPCONST  / (number_of_scans*D_A)

           write(*,*) scalekey(i)
enddo

do i=1,ende

if(g_bool == .false.) then


open(unit=113,file='1.ndx')

write(113,'(a9)') '[ group ]'

write(113,*)                  hbmap_donor(i) , hbmap_hydrogen(i) , hbmap_acceptor(i)

close(unit=113)

call system('cat 1.ndx system.ndx > 2.ndx')

sc = 0
   
do k=1,natoms

    if(hbmap_donor(i) == atomnum(k)) then

       resid = resnum(k) 

    endif

enddo

l = 1

do k=1,natoms

     if(resnum(k) .ge. resid - modulonumber .and. resnum(k) .le. resid + modulonumber) then

       don_x_res(l) = coord_x(k)
       don_y_res(l) = coord_y(k) 
       don_z_res(l) = coord_z(k)


       l = l + 1

     endif

enddo

resatomnumber  = l - 1

sc = 0

       do k=1,natoms

              if(atomnum(k) == hbmap_donor(i)) then

                     don_x = coord_x(k)
                     don_y = coord_y(k)
                     don_z = coord_z(k)

              endif

              if(atomnum(k) == hbmap_acceptor(i)) then

                     acc_x = coord_x(k)
                     acc_y = coord_y(k)
                     acc_z = coord_z(k)

              endif

       enddo
 
   do n=1,number_of_scans

                  do p = 1,resatomnumber

                  new_x(p) = sc*(don_x-acc_x)+don_x_res(p) 
                  new_y(p) = sc*(don_y-acc_y)+don_y_res(p) 
                  new_z(p) = sc*(don_z-acc_z)+don_z_res(p)

                  enddo

                  sc = sc + scalekey(i)

l = 1

       do k=1,natoms

              if(resnum(k) .ge. resid - modulonumber .and. resnum(k) .le. resid + modulonumber .and. l == 1) then

                  do p = 1,resatomnumber                 

                     coord_x(k+p-1) = new_x(p)
                     coord_y(k+p-1) = new_y(p)
                     coord_z(k+p-1) = new_z(p)                    

                     l = l + 1

                  enddo

              endif

       enddo   
 


       open(unit=1,file='first.gro')

       write(1,*)
       write(1,*) natoms

       do k=1,natoms

           write(1,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z           

       enddo

       write(1,*) box_x,box_y,box_z

       close(unit=1)

       call system('cp first.gro minimized.gro')

       call minimization_of_protein

       scannumber = 3
       
       scanstep   = n


       call scan_energy_of_protein(energy_val,resnumber,scannumber,n,alert)

       energy_scan_break(i,n) = energy_val

       if(alert == .true.) then

          exit

       endif
!       write(*,*) energy_scan_break(i,n),'energy_break_scan(i,n)'
       if(n==number_of_scans) then

                 new_break = .true.

       else

                 new_break = .false.

       endif
        
!       call write_file(scannumber,scanstep,energy_val,resnumber,selectamino,kmc_step,select,execsteps,g_bool,aminoresnum,time_end,rate,new_break)

        if(g_bool == .false.) then
      
           write(resnumber*1012,*) n,energy_val

           if(new_break == .true.) then

                write(resnumber*1012,*) ' '

           endif
        endif

        if(g_bool == .true.) then
      
           write(resnumber*1004,*) n,energy_val

           if(new_break == .true.) then

                write(resnumber*1004,*) ' '

           endif

        endif

!       write(*,*) n,energy_scan(resnumber,scannumber,n)

       open(unit=1,file='minimized.gro')

       read(1,*)
       read(1,*)

       do k=1,natoms

           read(1,'(i5,2a5,i5,3f8.3,3f8.4)',end=9005,err=9005) resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z

       enddo

       9005 continue

    enddo                   

endif

enddo

do i=1,ende

if(g_bool == .true.) then

        if(i == maxi) then

open(unit=113,file='1.ndx')

write(113,'(a9)') '[ group ]'

write(113,*)                  hbmap_donor(i) , hbmap_hydrogen(i) , hbmap_acceptor(i)

close(unit=113)

call system('cat 1.ndx system.ndx > 2.ndx')

sc = 0
   
do k=1,natoms

    if(hbmap_donor(i) == atomnum(k)) then

       resid = resnum(k) 

    endif

enddo

l = 1

do k=1,natoms

     if(resnum(k) .ge. resid - modulonumber .and. resid .le. resid + modulonumber) then

       don_x_res(l) = coord_x(k)
       don_y_res(l) = coord_y(k) 
       don_z_res(l) = coord_z(k)


       l = l + 1

     endif

enddo

resatomnumber  = l - 1

sc = 0

       do k=1,natoms

              if(atomnum(k) == hbmap_donor(i)) then

                     don_x = coord_x(k)
                     don_y = coord_y(k)
                     don_z = coord_z(k)

              endif

              if(atomnum(k) == hbmap_acceptor(i)) then

                     acc_x = coord_x(k)
                     acc_y = coord_y(k)
                     acc_z = coord_z(k)

              endif

       enddo
 
   do n=1,number_of_scans

                  do p = 1,resatomnumber

                  new_x(p) = sc*(don_x-acc_x)+don_x_res(p) 
                  new_y(p) = sc*(don_y-acc_y)+don_y_res(p) 
                  new_z(p) = sc*(don_z-acc_z)+don_z_res(p)

                  enddo

                  sc = sc + scalekey(i)

l = 1

       do k=1,natoms

              if(resnum(k) .ge. resid - modulonumber .and. resnum(k) .le. resid + modulonumber .and. l == 1) then

                  do p = 1,resatomnumber                 

                     coord_x(k+p-1) = new_x(p)
                     coord_y(k+p-1) = new_y(p)
                     coord_z(k+p-1) = new_z(p)                    

                     l = l + 1

                  enddo

              endif

       enddo   
 


       open(unit=1,file='first.gro')

       write(1,*)
       write(1,*) natoms

       do k=1,natoms

           write(1,'(i5,2a5,i5,3f8.3,3f8.4)') resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z           

       enddo

       write(1,*) box_x,box_y,box_z

       close(unit=1)

       call system('cp first.gro minimized.gro')

       call minimization_of_protein

       scannumber = 3
       
       scanstep   = n

       if(scanbool == .true.) then

       call scan_energy_of_protein(energy_val,resnumber,scannumber,n,alert)

       if(alert == .true.) then

          exit

       endif

       energy_scan_break(i,n) = energy_val

!       write(*,*) energy_scan_break(i,n),'energy_break_scan(i,n)'
       if(n==number_of_scans) then

                 new_break = .true.

       else

                 new_break = .false.

       endif
        
!       call write_file(scannumber,scanstep,energy_val,resnumber,selectamino,kmc_step,select,execsteps,g_bool,aminoresnum,time_end,rate,new_break)

        if(g_bool == .false.) then
      
           write(resnumber*1012,*) n,energy_val

           if(new_break == .true.) then

                write(resnumber*1012,*) ' '

           endif

        endif

        if(g_bool == .true.) then
      
           write(resnumber*1004,*) n,energy_val

           if(new_break == .true.) then

                write(resnumber*1004,*) ' '

           endif
        endif

       endif

!       write(111,*) n,energy_scan(resnumber,scannumber,n)

       open(unit=1,file='minimized.gro')

       read(1,*)
       read(1,*)

       do k=1,natoms

           read(1,'(i5,2a5,i5,3f8.3,3f8.4)',end=9001,err=9001) resnum(k), restype(k), atomtype(k), atomnum(k), coord_x(k), &
                                             coord_y(k), coord_z(k), vel_x , vel_y , vel_z

       enddo

       9001 continue

    enddo                   

endif
endif

enddo

if(g_bool == .false.) then

do i=1,ende

          maximum(i) = energy_scan_break(i,1)

   do k=1,number_of_scans

          if(maximum(i) .lt. energy_scan_break(i,k)) then

             maximum(i)  =  energy_scan_break(i,k)


          endif

  enddo

enddo

maxi  = 1
maxim = maximum(1)

do i=1,ende

  if(maximum(i) .gt. maxim .and. maximum(i).lt. 0) then

     maxi = i

  endif

enddo


do k=1,number_of_scans

   energy_scan(resnumber,3,k) = energy_scan_break(maxi,k)

enddo


endif


deallocate(restype,stat=ierr)
deallocate(resnum,stat=ierr)
deallocate(atomtype,stat=ierr)
deallocate(atomnum,stat=ierr)
deallocate(coord_x,stat=ierr)
deallocate(coord_y,stat=ierr)
deallocate(coord_z,stat=ierr)



end subroutine FORMATIONAL_TRANSLATIONAL_MOTION 

