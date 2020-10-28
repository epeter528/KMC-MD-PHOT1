subroutine TRAJ(md_time,step)

implicit none

integer :: timer

integer      :: md_time

integer      :: step

if(step == 1) then

         call system('mv md_traj.xtc md_trajx_internal.xtc')
         call system('mv md_traj.trr md_trajx_internal.trr')

endif

if(step .ge. 2) then

timer = (step-1)*md_time/1000

         open(unit=453,file='trjcat.sh')

          write(453,*) '#!/bin/bash'
          write(453,*)
          write(453,'(a113)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx_internal.xtc md_traj.xtc -o md_trajx_internal.xtc -settime << EOF'
          write(453,*) '0'
          write(453,*) timer
          write(453,*) 'EOF'

         close(unit=453)

         open(unit=455,file='trjcat2.sh')


          write(455,*) '#!/bin/bash'
          write(455,*)
          write(455,'(a113)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx_internal.trr md_traj.trr -o md_trajx_internal.trr -settime << EOF'
          write(455,*) '0'
          write(455,*) timer
          write(455,*) 'EOF'          

           write(455,*) timer
          
           close(unit=455) 
         call system('chmod 744 trjcat.sh')         
         call system('chmod 744 trjcat2.sh')

         call system('./trjcat.sh')
         call system('./trjcat2.sh')      


endif

end subroutine TRAJ