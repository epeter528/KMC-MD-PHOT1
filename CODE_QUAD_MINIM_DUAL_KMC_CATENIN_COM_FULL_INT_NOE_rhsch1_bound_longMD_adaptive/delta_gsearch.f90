subroutine DELTA_G_SEARCH(icounter,scansteps,min_delta_g,max_delta_g,energy_scan,execsteps,selectamino,select,sumrate)

implicit none

 double precision , parameter :: T = 300
 double precision , parameter :: R = 8.314
 double precision , parameter :: frequency = 6.250294745d9

integer     :: i,j,k,l,m,n

real(kind=8),dimension(1000,4,100) :: energy_scan

integer     :: select,secondminscan

integer     :: scansteps

integer     :: ende

real(kind=8)                    :: harvest3

real(kind=8)  , dimension(1000)  :: selectionrate
integer       , dimension(1000)  :: selectionaminoacid
integer       , dimension(1000)  :: selectionnumber
integer                         :: selectamino

integer                         :: execsteps
real(kind=8) , dimension(1000,3) :: delta_g,rate_g
real(kind=8)                    :: sumrate
integer                          :: p
real(kind=8) , dimension(100,3)  :: secondminimum
integer      , dimension(100,3)  :: sec_min_scan
integer  , dimension(100,3)     :: firstmax_scan

integer :: scanmin

real(kind=8) , dimension(100,3) :: firstmaximum

integer  :: maxbereich

real(kind=8)                      :: firstmax,secondmin,min_delta_g,max_delta_g

integer                           :: scanmax,icounter
real(kind=8),dimension(100,3)    :: min_step_en
integer     ,dimension(100,3)    :: exec_scanstep

real(kind=8),dimension(100)      :: firstmini
integer     ,dimension(100)      :: scanmini

 character(len = 8)  :: date
 character(len =10)  :: time
 character(len = 5)  :: zone
 integer,dimension(8) :: values 
 integer,allocatable  :: iseed(:)
 integer              :: isize , idate(8)

call date_and_time(values=idate)

call random_seed(size=isize)
allocate( iseed(isize))
call random_seed(get=iseed)
iseed = iseed* (idate(8)-500)
call random_seed(put=iseed)
deallocate( iseed )


do i=1,icounter

         scanmin         = 1
         select          = 1

do j=1,3

firstmax         =  energy_scan(i,j,1)
firstmaximum(i,j) = energy_scan(i,j,1)
scanmax         = 1

   do k=1,scansteps

m = 0
maxbereich = 0 

        if(energy_scan(i,j,k).gt.firstmax) then

                 firstmax     = energy_scan(i,j,k)
                 scanmax      = k
                 m            = k                     
  
        endif
     
  enddo

  firstmaximum(i,j)  = firstmax
  firstmax_scan(i,j) = scanmax

enddo

enddo

do i=1,icounter

do j=1,3

   
secondmin        = energy_scan(i,j,firstmax_scan(i,j))
secondminscan    = 1


    do k=firstmax_scan(i,j),scansteps

         if(energy_scan(i,j,k).lt.secondmin.and.k.gt.firstmax_scan(i,j)) then

                  secondmin    = energy_scan(i,j,k)
                secondminscan  = k     

         endif

    enddo


secondminimum(i,j) = secondmin
sec_min_scan(i,j)  = secondminscan

if(sec_min_scan(i,j) .lt. firstmax_scan(i,j)) then

   sec_min_scan(i,j) = firstmax_scan(i,j) 

endif

enddo

do m=1,icounter

     do l=1,3

          write(*,*) secondminimum(m,l),firstmaximum(m,l)

     enddo

enddo


do j=1,3       


min_step_en(i,j)   = energy_scan(i,j,1)

exec_scanstep(i,j) = sec_min_scan(i,j)

if(exec_scanstep(i,j) == 0 .and. firstmax_scan(i,j) == 0) then

                  firstmini(1) = energy_scan(i,j,1)
                   scanmini(1) = 1


                             p = 1

                             m = 0

                 do k=1,scansteps

                       if(k.gt.2 .and. energy_scan(i,j,k-1) .gt. energy_scan(i,j,k) .and. energy_scan(i,j,k) .lt. energy_scan(i,j,k+1) &
                          & .and. energy_scan(i,j,1) .gt. energy_scan(i,j,k)) then

                          firstmini(p) = energy_scan(i,j,k)

                          scanmini(p)  = k

                          p = p + 1 

                          m = k

                       endif

                 enddo
     
                 firstmax_scan(i,j) = scanmini(1)
                 firstmaximum(i,j)  = firstmini(1)

                 do k=m,scansteps

                         if(k .gt. firstmax_scan(i,j) .and. energy_scan(i,j,k) .gt. firstmaximum(i,j)) then

                               firstmax_scan(i,j) = k
                               firstmaximum(i,j)  = energy_scan(i,j,k)

                               n = k

                         endif

                 enddo

                 do k=n,scansteps

                                exec_scanstep(i,j) = n
                                secondminimum(i,j) = n

                                if(energy_scan(i,j,k) .lt. secondminimum(i,j)) then

                                                           secondminimum(i,j) = energy_scan(i,j,k)
                                                           exec_scanstep(i,j) = k     

                                endif

                 enddo               

          endif

enddo

do j=1,3

     if(exec_scanstep(i,j) .le. scansteps/2) then

                   exec_scanstep(i,j) = scansteps

                   firstmaximum(i,j)  = energy_scan(i,j,1)

                   secondminimum(i,j) = energy_scan(i,j,scansteps)

    endif         

enddo

do j=1,3 

     if(firstmaximum(i,j) == 0) then

            firstmaximum(i,j)  = energy_scan(i,j,scansteps)

            exec_scanstep(i,j) = scansteps
            
     endif

enddo

enddo

do i=1,icounter

   do j=1,3 

      rate_g(i,j) = 0

   enddo

enddo

sumrate = 0

do i=1,icounter

       do j=1,3

            delta_g(i,j) = firstmaximum(i,j) - energy_scan(i,j,1)

            write(*,*) delta_g(i,j),'delta_g'

            if(delta_g(i,j).gt.min_delta_g .and. delta_g(i,j).lt.max_delta_g ) then
     
             rate_g(i,j)  = frequency * exp(-delta_g(i,j)*1000/(R*T))
        
  
!            else

!             rate_g(i,j) =  0

            endif

            sumrate = sumrate + rate_g(i,j)   

       enddo

enddo

write(*,*) sumrate,'sumrate'

k = 2


selectionrate(1) = 0

do i=1,icounter

    do j=1,3

       selectionrate(k)   = selectionrate(k-1) + rate_g(i,j)

       selectionnumber(k) = j

       selectionaminoacid(k) = i

       k = k + 1

    enddo

enddo

k=1
do i=1,icounter

   do j=1,3

       write(*,*) selectionrate(k)
       k = k+1

   enddo

enddo

call random_number(harvest3)

do i = 1,icounter*3+1

    if(selectionrate(i-1) .lt. sumrate*harvest3 .and. sumrate*harvest3 .le. selectionrate(i)) then 
                           

                                select        = selectionnumber(i)

                                selectamino   = selectionaminoacid(i)

    endif

enddo

execsteps    = exec_scanstep(selectamino,select)

if(execsteps == 0) then

   execsteps = scansteps

endif

write(*,*) execsteps,'execsteps'
write(*,*) selectamino,'selectamino'

end subroutine DELTA_G_SEARCH