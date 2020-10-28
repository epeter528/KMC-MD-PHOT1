subroutine read_searchlist(aminoresnum,noe_val,icounter)

implicit none

integer                 :: i

integer                 :: icounter

integer,dimension(10000)  :: aminoresnum

real(kind=8),dimension(10000)    :: noe_val

open(unit=1,file='searchlist.ndx')


i = 1

do 

read(1,*,end=1,err=1) aminoresnum(i),noe_val(i)

i = i + 1

enddo

1 continue

icounter = i - 1


end subroutine read_searchlist