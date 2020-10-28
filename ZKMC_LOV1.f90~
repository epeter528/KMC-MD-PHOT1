program zkmc

implicit none

integer                :: i,j,k

integer                :: icounter

integer                :: KMC_simloops

integer                :: natoms

real(kind=8)           :: box_x,box_y,box_z,rot_change,expconst2

logical                :: timebool

integer,dimension(10000) :: aminoresnum

logical                :: scanbool , g_bool

integer                :: resnumber

real(kind=8),dimension(10000,4,100) :: energy_scan

integer                :: scansteps,scansteps2,sumo

real(kind=8)           :: min_delta_g,max_delta_g,noe_min

integer                :: execsteps,selectamino

real(kind=8)           :: sumrate,expconst,maximum_kmc_time

real(kind=8)           :: time,kmcbegintime,noe_max

real(kind=8)           :: maxtime,time_end,rate,total_time,mdquotient

integer                :: md_time,mod_noe

integer                :: select,scannumber,timer

integer                :: step,maxi,modulonumber,max_mdtime

integer                :: kmc_step,kmc_simnumber,dark_or_light
integer                :: after_KMC_MD_INT, AFTER_KMC_MD_TIMER

integer                :: intevil,intevil2,intevil3

integer                :: kmc_total_time,break

 character(len = 8)  :: date
 character(len = 5)  :: zone
 integer,dimension(8) :: values 
 integer,allocatable  :: iseed(:)
 integer              :: isize , idate(8),scandihnumber

 real(kind=8)         :: harvest3,rot_max,max_rmsd1

real(kind=8),dimension(10000) :: noe_val

character(len=4)      :: c11

integer,dimension(10000) :: sum_of_md_time



open(unit=1,file='KMCINPUT')

read(1,*)
read(1,*) KMC_simloops
read(1,*)
read(1,*) md_time
read(1,*) 
read(1,*) mdquotient
read(1,*) 
read(1,*) scansteps
read(1,*)
read(1,*) scandihnumber
read(1,*) 
read(1,*) rot_max
read(1,*)
read(1,*) max_rmsd1
read(1,*)
read(1,*) min_delta_g
read(1,*)
read(1,*) max_delta_g
read(1,*) 
read(1,*) maximum_kmc_time
read(1,*) 
read(1,*) dark_or_light ! 0== light , else == dark
read(1,*)
read(1,*) AFTER_KMC_MD_INT
read(1,*)
read(1,*) AFTER_KMC_MD_TIMER
read(1,*)
read(1,*) noe_max
read(1,*) noe_min
read(1,*)
read(1,*) modulonumber
read(1,*)
read(1,*) break
read(1,*)
read(1,*) sumo

mod_noe = modulonumber * 2

write(*,*) 'write rac1 or lov2 or phot (without NA+) or cate'

read(*,'(a4)') c11

write(*,*) 'continuation ? 1 = yes, 0 = no'

read(*,*) i

if(i==1) then

 write(*,*) 'total time for continuation ?'

 read(*,*) total_time

 else

  total_time = 0

endif

if(i==1) then

   inquire(file='TEST_RATE_INTERNAL',exist = g_bool)

   if(g_bool == .true.) then

      write(*,*) 'TEST_RATE_INTERNAL exists already !'

      stop

   endif

   inquire(file='md_trajx_internal.trr',exist = g_bool)

   if(g_bool == .true.) then

      write(*,*) 'md_trajx_internal.trr exists already !'

      stop

   endif

   inquire(file='md_trajx.trr',exist = g_bool)

   if(g_bool == .true.) then

      write(*,*) 'md_trajx.trr exists already !'

      stop

   endif
   inquire(file='TEST_RATE',exist = g_bool)

   if(g_bool == .true.) then

      write(*,*) 'TEST_RATE exists already !'

      stop

   endif

   inquire(file='minimized_water.gro.bak',exist = g_bool)

   if(g_bool == .true.) then

      write(*,*) 'minimized_water.gro.bak exists already !'

      stop

   endif

   call system('cp minimized_water.gro minimized_water.gro.bak')

   write(*,*) 'backup minimized_water.gro to minimized_water.gro.bak !'

   call system('cp md_cheater2.gro minimized_water.gro')

   write(*,*) 'cp md_cheater2.gro minimized_water.gro  !!!'


endif

kmc_total_time = KMC_simloops * md_time


call WRITE_MDP(md_time,dark_or_light,c11)

call system('rm *.xvg')

!call OPEN_FILES(aminoresnum,icounter)


open(unit=10000,file='TEST_RATE_INTERNAL')
open(unit=10001,file='TIME')

call date_and_time(values=idate)

call random_seed(size=isize)
allocate( iseed(isize))
call random_seed(get=iseed)
iseed = iseed* (idate(8)-500)
call random_seed(put=iseed)
deallocate( iseed )

do j=1,KMC_simloops

! mdrun

! MD/Hauptabteil

call random_number(harvest3) 

intevil = abs(harvest3*3)

intevil2 = abs(harvest3*3) + intevil

intevil3 = abs(harvest3*2) + intevil2 

if(intevil2 == intevil3) then

   intevil3 = intevil2 + 1

endif

call random_number(harvest3)

rot_change = harvest3*rot_max 


if(j==1) then 

  call system('cp minimized_water.gro start.gro')

endif

call MD_PROD

call SOLVENT_CONFIG_SEPARATOR

call noe_sep(noe_max,noe_min,mod_noe)

open(unit=144,file='index2.sh')

write(144,*) '#!/bin/bash'

write(144,'(a77)') '/loctmp/pee18323/GRO_bin/bin/make_ndx -f minimized3.gro -o system.ndx << EOF'

write(144,*) 'q'

write(144,*) 'EOF'

close(unit=144)

call system('chmod 744 index2.sh')

call system('./index2.sh')


call READ_SEARCHLIST(aminoresnum,noe_val,icounter)

!if(j==1) then

!call OPEN_FILES(aminoresnum,icounter)

!endif

scanbool = .true.

do i=1,icounter


g_bool = .false.

call BREAK_DIHEDRALANGLE_SCAN(aminoresnum,i,scanbool,rot_change,scandihnumber,energy_scan,g_bool,intevil,intevil2,intevil3)

write(*,*) 'start of break'

call BREAK_TRANSLATIONAL_MOTION(aminoresnum,scansteps,noe_val(i)*max_rmsd1,i,scanbool,energy_scan,g_bool,maxi,modulonumber,break)

open(unit=1,file='H_search_new2_int.sh')

write(1,*) '#!/bin/bash'
write(1,'(a97,f4.2,a7)') '/loctmp/pee18323/GRO_bin/bin/g_hbond -f start.gro -n 2.ndx -s run.tpr -hbn table_hbond_B.ndx -r ', 0.35+noe_val(i)*max_rmsd1,' << EOF'
write(1,*) '1'
write(1,*) '2'
write(1,*) 'EOF'

close(unit=1)

call FORMATIONAL_TRANSLATIONAL_MOTION(aminoresnum,scansteps,i,scanbool,energy_scan,g_bool,maxi,modulonumber)

enddo

scanbool = .true.

call DELTA_G_SEARCH(icounter,scansteps,min_delta_g,max_delta_g,energy_scan,execsteps,selectamino,select,sumrate)

call TIME_TEST(sumrate,time,g_bool,maximum_KMC_time)

if(g_bool == .true.) then

kmc_step = j

total_time = total_time + time

timer = timer + md_time

sum_of_md_time(j) = timer



write(10000,*) 'total_time',total_time
write(10000,*) 'time', time
write(10000,*) 'md_time', md_time
write(10000,*) select,'select'
write(10000,*) execsteps,'execsteps'
write(10000,*) selectamino,'selectamino'
write(10000,*) maxi,'maxi'
write(10000,*) j,'KMC_step'

if(select == 1) then

scansteps2 = execsteps

i = selectamino

write(*,*) 'selected break dihedral',selectamino

           call BREAK_DIHEDRALANGLE_SCAN(aminoresnum,i,scanbool,rot_change,scandihnumber,energy_scan,g_bool,intevil,intevil2,intevil3)
endif

if(select == 2) then

scansteps2 = execsteps
i = selectamino

expconst2 = noe_val(i)*max_rmsd1*execsteps/scansteps

write(*,*) 'selected break trans',selectamino


           call BREAK_TRANSLATIONAL_MOTION(aminoresnum,scansteps2,expconst2,i,scanbool,energy_scan,g_bool,maxi,modulonumber,break)

endif
 
if(select == 3) then

scansteps2 = execsteps
i = selectamino

write(*,*) 'selected form trans',selectamino

open(unit=1,file='H_search_new2_int.sh')

write(1,*) '#!/bin/bash'
write(1,'(a97,f4.2,a7)') '/loctmp/pee18323/GRO_bin/bin/g_hbond -f start.gro -n 2.ndx -s run.tpr -hbn table_hbond_B.ndx -r ', 0.35+noe_val(i)*max_rmsd1,' << EOF'
write(1,*) '1'
write(1,*) '2'
write(1,*) 'EOF'

close(unit=1)
           call FORMATIONAL_TRANSLATIONAL_MOTION(aminoresnum,scansteps2,i,scanbool,energy_scan,g_bool,maxi,modulonumber)

endif

endif


write(10001,*) total_time

if(scanbool == .false. .and. g_bool == .false.) then

   write(10000,*) '########## no event !!!'
   write(10000,*) kmc_step

endif

if(g_bool == .false.) then

   call system('cp minimized3.gro prot.gro')

endif

if(g_bool == .true.) then

   call system('cp minimized.gro prot.gro')

endif

call solvent

call system('cp first.gro md_cheater2.gro')

call EQUIL_RELAXATION

step = j

call TRAJ(sum_of_md_time,step)

if(total_time .ge. maximum_kmc_time) then

      timebool = .true.

else

      timebool = .false.

endif

if(timebool == .true.) then

                           kmc_simnumber = j

                           exit

endif


enddo ! END of internal KMC

if(timebool == .true.) then

call AFTER_KMC_MD(md_time,kmc_simnumber,dark_or_light,kmc_total_time)

endif

if(AFTER_KMC_MD_INT == 1) then

call AFTER_KMC_MD(AFTER_KMC_MD_TIMER,1,dark_or_light,kmc_total_time)

endif

! Start of external KMC

kmcbegintime = total_time

call Potential(kmcbegintime)

end program zkmc


