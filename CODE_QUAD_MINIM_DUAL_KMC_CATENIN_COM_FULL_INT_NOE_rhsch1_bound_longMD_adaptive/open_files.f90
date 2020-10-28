subroutine open_files(aminoacid_resnum,icounter)

implicit none

integer                :: icounter,i
integer,dimension(100) :: aminoacid_resnum

character(len=40),dimension(100) :: filename

do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'DIH_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'DIH_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

!   write(*,*) filename(i)
  open(unit=i*1001,file=trim(filename(i)))

enddo

do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'TRANS_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'TRANS_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

  open(unit=i*1011,file=trim(filename(i)))

enddo

do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'FORM_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'FORM_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

  open(unit=i*1012,file=trim(filename(i)))

enddo


do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'EXEC_DIH_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'EXEC_DIH_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

!   write(*,*) filename(i)
  open(unit=i*1002,file=trim(filename(i)))

enddo

do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'EXEC_TRANS_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'EXEC_TRANS_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

  open(unit=i*1003,file=trim(filename(i)))

enddo

do i=1,icounter

   if(aminoacid_resnum(i).lt.100) then

   write(filename(i),'(a,i2,a)') 'EXEC_FORM_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif
   if(aminoacid_resnum(i).ge.100) then

   write(filename(i),'(a,i3,a)') 'EXEC_FORM_CHANGE_',aminoacid_resnum(i),'.xvg'
 
   endif

enddo

do i=1,icounter

  open(unit=i*1004,file=trim(filename(i)))

enddo


end subroutine open_files