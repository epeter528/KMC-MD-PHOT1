subroutine TIME_TEST(sumrate,timeout,g_bool,maxtime)

implicit none

 character(len = 8)  :: date
 character(len =10)  :: time
 character(len = 5)  :: zone
 integer,dimension(8) :: values 
 integer,allocatable  :: iseed(:)
 integer              :: isize , idate(8)

 real(kind=8)         :: harvest3

 real(kind=8)         :: sumrate

 real(kind=8)         :: timeout,maxtime,mdtime

 logical              :: g_bool
 
call date_and_time(values=idate)

call random_seed(size=isize)
allocate( iseed(isize))
call random_seed(get=iseed)
iseed = iseed* (idate(8)-500)
call random_seed(put=iseed)
deallocate( iseed )


call random_number(harvest3)

timeout = (-1)*log(harvest3)/ sumrate

if(timeout.ge.maxtime) then

                g_bool = .false.

else 

                g_bool = .true.

endif




end subroutine TIME_TEST