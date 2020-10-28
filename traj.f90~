subroutine TRAJ(sum_of_md_time,step)

implicit none

integer,dimension(10000) :: sum_of_md_time

integer      :: md_time

integer      :: step

if(step == 1) then

         call system('mv md_traj.xtc md_trajx_internal.xtc')
         call system('mv md_traj.trr md_trajx_internal.trr')

endif

if(step .ge. 2) then

!timer = (step-1)*md_time/1000
md_time  = sum_of_md_time(step-1) / 1000


         open(unit=453,file='trjcat.sh')

          write(453,*) '#!/bin/bash'
          write(453,*)
          write(453,'(a113)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx_internal.xtc md_traj.xtc -o md_trajx_internal.xtc -settime << EOF'
          write(453,*) '0'
          write(453,*) md_time
          write(453,*) 'EOF'

         close(unit=453)

         open(unit=455,file='trjcat2.sh')


          write(455,*) '#!/bin/bash'
          write(455,*)
          write(455,'(a113)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx_internal.trr md_traj.trr -o md_trajx_internal.trr -settime << EOF'
          write(455,*) '0'
          write(455,*) md_time
          write(455,*) 'EOF'          

          close(unit=455) 


         call system('chmod 744 trjcat.sh')         
         call system('chmod 744 trjcat2.sh')

         call system('./trjcat.sh')
         call system('./trjcat2.sh')      


endif

end subroutine TRAJ