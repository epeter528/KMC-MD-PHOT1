!program solvent
subroutine solvent

implicit none

integer         :: m,k,i,p,l,n,ierr,j

integer         :: changenumber,nsol

integer         :: natoms,resnum,atomnum

real(kind=8) ,dimension(500000) :: coord_test_x,coord_test_y,coord_test_z

real(kind=8)                    :: vel_x,vel_y,vel_z

real(kind=8)                    :: coord_x,coord_y,coord_z,fact

character(len=5)                :: restype,atomtype

integer,     dimension(500000)  :: resnumsol,atomnumsol

real(kind=8) , dimension(500000) :: coord_xsol,coord_ysol,coord_zsol

character(len=5),dimension(500000) :: restypesol,atomtypesol

real(kind=8)                       :: diff_x,diff_y,diff_z,diff_tot,diff_x1,diff_x2

real(kind=8)                       :: diff_y1,diff_y2,diff_z1,diff_z2,ranger

real(kind=8)                       :: box_x,box_y,box_z,ff21,ff22,ff23,ff24

logical                            :: changebool,changebool2

logical                            :: changeboolean3

real(kind=8)                       :: harvest3,box_x1,box_y1,box_z1

character(len = 8)  :: date

character(len =10)  :: time

character(len = 5)  :: zone

integer             :: tot

integer,dimension(8) :: values 

integer,allocatable  :: iseed(:)

integer              :: isize , idate(8),nprotein

real(kind=8),dimension(:),allocatable      :: coord_x22,coord_y22,coord_z22

integer     ,dimension(:),allocatable      :: resnum22,atomnum22

character(len=5),dimension(:),allocatable  :: restype22,atomtype22

logical                                    :: exitboolean



write(*,*) ranger

open(unit=888,file='s11')


i = 1

do 

           read(888,'(i5,2a5,i5,3f8.3)',end=9716,err=9716) resnumsol(i), restypesol(i), atomtypesol(i), atomnumsol(i), coord_xsol(i), &
                                             coord_ysol(i), coord_zsol(i)

        !   write(*,'(i5,2a5,i5,3f8.3)') resnumsol(i), restypesol(i), atomtypesol(i), atomnumsol(i), coord_xsol(i), &
        !                                     coord_ysol(i), coord_zsol(i) 
           if(restypesol(i) == 'SOL  ') then

              i =  i + 1

           endif
enddo


9716 continue

close(unit=888)

nsol = i - 1

open(unit=1111,file='prot.gro')


read(1111,*)
read(1111,*) natoms

do m=1,NATOMS

           read(1111,'(i5,2a5,i5,3f8.3)',end=91111,err=91111) resnum, restype, atomtype, atomnum, coord_test_x(m), &
                                             coord_test_y(m), coord_test_z(m)
enddo

read(1111,*) box_x,box_y,box_z

91111     continue

close(unit=1111)

tot = natoms + nsol

write(*,*) nsol
write(*,*) tot

allocate(coord_x22(tot),stat=ierr)
allocate(coord_y22(tot),stat=ierr)
allocate(coord_z22(tot),stat=ierr)
allocate(resnum22(tot),stat=ierr)
allocate(atomnum22(tot),stat=ierr)
allocate(restype22(tot),stat=ierr)
allocate(atomtype22(tot),stat=ierr)

!call system('/scratch/pee18323/GRO_seq/bin/editconf -f prot.gro -c -o prot.gro')

nsol = nsol / 3

call system('/loctmp/pee18323/GRO_bin/bin/editconf -f prot.gro -c -o prot.gro -bt cubic ')

open(unit=1,file='solvent.sh')

write(1,*) '#!/bin/bash'
write(1,*)
write(1,'(a85,i7)') '/loctmp/pee18323/GRO_bin/bin/genbox -cp prot.gro -cs spc216.gro -o first.gro -maxsol ', nsol 

close(unit=1)

call system('chmod 744 solvent.sh')

call system('./solvent.sh')

call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_solequil.mdp -c first.gro -p LOV2gen.top -o equil.tpr')
call system('mpirun -np 2 /loctmp/pee18323/GRO_bin/bin/mdrun -v -s equil.tpr -c first.gro')

end subroutine solvent
!end program solvent
