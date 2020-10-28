subroutine noe_sep(noe_max,noe_min,mod_noe)

implicit none


integer      :: i,noe_number,k,l,m,p,mod_noe

real(kind=8) :: noe_max,max_mode,max_rms,noe_min,min_mode

real(kind=8) , dimension(1000000) :: NOE_R3 

integer, dimension(1000000)       :: rnr,rnr2,noe_rnr,noe_rnr2,max_seg

integer :: ianr,anr,ianr2,anr2

character(len=5) :: anm,rnm,anm2,rnm2

character(len=1) :: char1


open(unit=1,file='noe.sh')

write(1,*) '#!/bin/bash'

write(1,'(a86)')'/loctmp/pee18323/GRO_bin/bin/g_rmsf -f md_traj.trr -s md.tpr -o rmsf.xvg -res << EOF'

write(1,*) '1'

write(1,*) 'EOF'

close(unit=1)

call system('chmod 744 noe.sh')

call system('./noe.sh')

open(unit=1,file='rmsf.xvg')

do

read(1,'(a1)',end=2,err=2) char1

if(char1 .ne. '#' .and. char1 .ne. '@') then

           exit

endif

enddo

2 continue

i = 1

do 

read(1,*,end=1,err=1) rnr(i),NOE_R3(i)

!write(*,*)  rnr(i),NOE_R3(i)

i = i + 1

enddo

1 continue

noe_number = i-1

close(unit=1)

open(unit=1,file='searchlist.ndx')

k = 1


!do l=1,rnr(noe_number)

max_mode = 0

do i=1,noe_number

       if(NOE_R3(i) .ge. max_mode) then

                     max_mode = NOE_R3(i)

       endif

enddo

k = 1

do i=1,noe_number

        if(NOE_R3(i) .ge. max_mode - noe_max*max_mode ) then

           max_seg(k) =  rnr(i)

!           write(*,*) max_seg(k),'max_seg'

           k = k + 1

        endif

enddo

min_mode = max_mode

do i=1,noe_number

       if(NOE_R3(i) .lt. min_mode) then

                     min_mode = NOE_R3(i)

       endif

enddo


max_rms = min_mode + noe_min*min_mode

write(*,*) max_rms,'maxrms'
!close(unit=1)



p = 1


do l=1,noe_number


do i=1,k-1

      if(rnr(l) .gt. max_seg(i) .and. rnr(l) .lt. max_seg(i+1) .and. (max_seg(i+1) - max_seg(i)).ge. 8 ) then

          m = 1

          

          exit

      else 

          m = 2

      endif

enddo

     if(m == 1) then

         p = p + 1
 
         if(noe_r3(l) .lt. max_rms) then

!            write(*,*) noe_r3(l),max_rms

            write(1,*) rnr(l),noe_r3(l)

         endif

     endif

enddo

close(unit=1)

end subroutine noe_sep
!end program test