subroutine AFTER_KMC_MD(mdtime,kmccycles,dark_or_light,kmc_total_time)

implicit none

integer      :: relaxtime

integer      :: kmccycles

integer      :: mdtime

integer      :: timer

integer      :: dark_or_light

integer      :: kmc_total_time


timer = kmc_total_time / 1000

relaxtime = timer + mdtime

            call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f relax.mdp -c md_cheater2.gro -p LOV2gen.top -o md_relax.tpr')

            write(*,*) 'here is md relax time'

            call system ('/loctmp/pee18323/GRO_bin/bin/tpbconv -s md_relax.tpr -f md_trajx_internal.trr -o md_relax2.tpr -nice 0')

            call system ('mpirun -np 2 /loctmp/pee18323/GRO_bin/bin/mdrun -v -s md_relax2.tpr -nice 0 -c md_finaloutrelax.gro -o md_trajrelax.trr -x md_traj.xtc -e md_enerrelax.edr')


            open(unit=453,file='trjcat.sh')

             write(453,*) '#!/bin/bash'
             write(453,*)
             write(453,'(a129)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx_internal.trr md_trajrelax.trr -o md_trajx_internal_plus_relax.trr -settime << EOF'
             write(453,*) '0'
             write(453,*) timer
             write(453,*) 'EOF'

            close(unit=453)

            call system ('chmod 744 trjcat.sh')

            call system ('./trjcat.sh')

            write(*,*) 'relax finished'



end subroutine AFTER_KMC_MD