subroutine Potential(kmcbegintime)

implicit none

 double precision , parameter :: R = 8.314
 double precision , parameter :: k_B = 1.380650524d-23
 double precision , parameter :: frequency = 6.250294745d9
 double precision , parameter :: T = 300
 double precision , parameter :: e_V = 1.602176487d-19
 double precision , parameter :: bohr = 0.052917721086d-9
 double precision , parameter :: epsilon_0 = 8.854187817d-12
 double precision , parameter :: PI = 3.141592654
 double precision , parameter :: COUL_CONST = 138.935485
 double precision , parameter :: RAD_CONST  = 1d-9
 double precision , parameter :: MOL        = 6.02214179d23
 double precision , parameter :: charge_e = 1.6021764870d-19
 double precision , parameter :: sigma_6_O =  0.00226195360
 double precision , parameter :: sigma_12_O = 7.41493210d-07
 double precision , parameter :: sigma_6_N = 0.00243640960
 double precision , parameter :: sigma_12_N = 1.6926010d-06
 double precision , parameter :: sigma_6_H = 0
 double precision , parameter :: sigma_12_H = 0
 double precision , parameter :: sigma_6_S  = 0.00998400640 
 double precision , parameter :: sigma_12_S  = 1.30754560d-05
 double precision , parameter :: sigma_6_P  = 0.014737960
 double precision , parameter :: sigma_12_P  = 2.21935210d-05
 double precision , parameter :: sigma_6_C  = 0.00234062440
 double precision , parameter :: sigma_12_C  = 3.3745690d-06
 
 real(kind=8)              :: confact,box_fact_x,box_fact_y,box_fact_z,angle_en
 real(kind=8)              :: scaling2c,scalekey2c,step_x,step_y,step_z,diffval,diffaver
 real(kind=8)              :: DIST_GRAPH_X, DIST_GRAPH_Y, DIST_GRAPH_Z
 real(kind=8)              :: factor_A,factor_B,factor_C,alpha,factor_D
 integer                   :: izhilf15,izhilf20,izhilf21,izhilf22,izhilf26,partnum,q
 integer                   :: o 
 integer                   :: y,i,j,k,n,p,mum,mum2,numbers, atom_nr,itplines,rtplines,countA,countB,a,l,m
 integer                   :: selectrate,selnumber
 integer                   :: counter , ende , endstep , endacc , m_acc2 ,b , mum3 ,endacc2 , m_don2,minimline
 integer                   :: FOOT,mini,maxi,mini2,maxi2,mini1,ATOMNUMBER,ende2,ende3,mam,RATER,RATER2,minu,kmccycles
 integer , dimension(100) :: hbmap_donor , hbmap_hydrogen , hbmap_acceptor , m_acc3 , m_don3
 integer , dimension(100) :: m_acc4,m_don4,endacc2b
 integer , dimension(100) :: Bhbmap_donor , Bhbmap_hydrogen , Bhbmap_acceptor 
 integer                   :: Bhbmap_donor2 , Bhbmap_hydrogen2 , Bhbmap_acceptor2
 integer                   :: NATOMS , resnum , atomnum , incount,endtable,endtable2,H_bond_number 
 integer                   :: m_donor,m_acc,minomax,minomax2,endmino,endmino2,matter,matter2,timer
 integer                   :: H_don_numspec,H_acc_numspec,hb_hyd_A,hb_acc_A,hb_don_A,mdtimestep
 integer                   :: inp,pcount_A,modnumber,mdtime,maxbereich,donorregionstart,donorregionend

 integer                   :: relaxtime

 integer                   :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15
 integer                   :: a16,a17,a18,a19,a20,a21,a22,a23,a24,ierr
 
 integer                   :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14
 
 integer                   :: deltamax,tu

 integer , dimension(50)   :: resnum_min_range, resnum_max_range 

 integer                   :: donormode,keepresmode1,keepresmode2,anglelines

 integer , dimension(100) :: hbmap_donor_resnum , hbmap_hyd_resnum , hbmap_acc_resnum
 integer , dimension(100) :: hbmap_donor_resnum_real, hbmap_acc_resnum_real
 integer , dimension(100) :: Bhbmap_donor_resnum ,Bhbmap_hyd_resnum , Bhbmap_acc_resnum
 integer , dimension(100) :: hbmap_donor_num , hbmap_hyd_num , hbmap_acc_num , count_printer
 integer , dimension(100) :: hbmap_donor_num_real,hbmap_acc_num_real 
 integer , dimension(100) :: Bhbmap_donor_num , Bhbmap_hyd_num , Bhbmap_acc_num
 integer , dimension(5000)   :: typenumberrange,typenumberacc
 integer , dimension(5000)   :: selectivenumber,selectivenumber2
 integer , dimension(5000,50) :: donor_resnum, donor_atomnum
 integer , dimension(5000,50) :: acc_resnum , acc_atomnum , acc_donornum , Bdonor_atomnum
 integer , dimension(5000,50) :: Bdonor_resnum , Bacc_resnum , Bacc_atomnum , Bacc_donornum 

 real(kind=8)                       :: AC_X , AC_Y , AC_Z , time2 , timeend , wrapper, EN_barrier, EN_barrierA
 real(kind=8)                       :: DO_X , DO_Y , DO_Z , zhilfB , sumrate , bsumrate,realsumrate,brealsumrate
 real(kind=8)                       :: DO_X2 , DO_Y2, DO_Z2
 
 real(kind=8)                       :: timekonst_delta,timefirst

 real(kind=8)                       :: AC_X2 , AC_Y2 , AC_Z2, damping 
 real(kind=8)                       :: DIST_X , DIST_Y , DIST_Z , TOTAL_Energy_Begin
 real(kind=8)                       :: DIST_X2 , DIST_Y2 , DIST_Z2, DIST_X3 , DIST_Y3 , DIST_Z3  
 real(kind=8)                       :: Dist_scale_x,Dist_scale_y,Dist_scale_z,Dist_scale_A

 real(kind=8)                       :: LJ_param_6,LJ_param_12,bmin_delta_g

 real(kind=8),dimension(100)        :: FIRSTMINI,DIST_FIRSTMINI 

 real(kind=8), dimension(100,100)      :: minix,maxix,distminix,distmaxix
 real(kind=8), dimension(100,100)      :: bminix,bmaxix,bdistminix,bdistmaxix
 
 real(kind=8), dimension(50)       :: EXPOSI,scalekeyB,scalekeyB2
 real(kind=8), dimension(50)       :: BEXPOSI,bEXPOSI2

 real(kind=8), dimension(100)        :: firstmax,dist_firstmax,secondmin,dist_secondmin

! real(kind=8), dimension(100,500)   ::  AC_XA

 real(kind=8), dimension(100)       :: MAXI_G,MINI_G,DIST_MAXI,MINI_A,rate_G,rate_H,DELTA_G,DELTA_H
 real(kind=8), dimension(100)       :: BMAXI_G,BMINI_G,BDIST_MAXI,BMINI_A,brate_G,brate_H,bDelta_G,bDelta_H
 real(kind=8), dimension(100)       :: BMAXI_G2,BMINI_G2,BDIST_MAXI2,BMINI_A2
 real(kind=8), dimension(100,100)    :: DISTANCE2,DISTANCE2B,DIST_GRAPH,DIST_GRAPH3
 real(kind=8)                       :: charge2,harvest3,Delta_G2,Delta_H2,Delta_time
 real(kind=8)                       :: cutoff , TOTALE , EXPCONST, EXPCONST2
 real(kind=8)                       :: scaling , atomcharge , range_betrag , scaling2
 real(kind=8)                       :: dist_acc_x,dist_acc_y,dist_acc_z
 real(kind=8)                       :: dist_acc_x2 , dist_acc_y2 , dist_acc_z2
 real(kind=8)                       :: coord_x , coord_y , coord_z , vel_x , vel_y , vel_z
 real(kind=8)                       :: LJ,TOTALCHARGE,ratemin,csumrate,kmcbegintime
 real(kind=8), dimension(5000)          :: LJ_6, LJ_12, LJ_pot
 real(kind=8), dimension(100)        :: selectionrate
 real(kind=8), dimension(5000)      :: LJ_6b, LJ_12b, LJ_potb 
 real(kind=8)                       :: COOR_ACCA_X,COOR_ACCA_Y,COOR_ACCA_Z,box_x,box_y,box_z
 real(kind=8)                       :: Chargeranger,Coulranger,LJranger,LJranger_pot
 real(kind=8)                       :: LJranger_12,LJranger_6,COOR_ACCA_BETRAG,CHARGE_ACCRES
 real(kind=8)                       :: accranger_charge,LJacc_6,LJacc_12,LJ_ranger_6,LJ_ranger_12,ow_solrange
 real(kind=8), dimension(100)       :: rate_1,rate_2,rate_1B,rate_2B,scalekey,scalekey2
 real(kind=8), dimension(5000)     :: COULOMB,hb_acc_charge,hb_acc_chargeB,COULOMB2
 real(kind=8), dimension(5000,50)   :: new_dir_acc_x , new_dir_acc_y , new_dir_acc_z
 real(kind=8), dimension(5000,50)   :: new_dir_acc_x2B , new_dir_acc_y2B, new_dir_acc_z2B
 real(kind=8), dimension(5000,50)   :: new_dir_acc_xA , new_dir_acc_yA , new_dir_acc_zA 
 real(kind=8), dimension(5000,50)   :: bnew_dir_acc_xA , bnew_dir_acc_yA , bnew_dir_acc_zA
 real(kind=8), dimension(100,200)    :: TOTAL_EN,TOTAL_EN2
 real(kind=8), dimension(5000)      :: charge_acc,charge,charge_atom,charge_atomB
 real(kind=8), dimension(5000)     :: Cut_range,Cut_range2
 real(kind=8), dimension(100)      :: hbmap_donor_x , hbmap_donor_y , hbmap_donor_z
 real(kind=8), dimension(100)      :: hbmap_donor_x_real , hbmap_donor_y_real , hbmap_donor_z_real
 real(kind=8), dimension(100)      :: Bhbmap_donor_x , Bhbmap_donor_y , Bhbmap_donor_z
 real(kind=8), dimension(100)      :: hbmap_donor_vx , hbmap_donor_vy , hbmap_donor_vz
 real(kind=8), dimension(100)      :: hbmap_donor_vx_real , hbmap_donor_vy_real , hbmap_donor_vz_real
 real(kind=8), dimension(100)      :: Bhbmap_donor_vx , Bhbmap_donor_vy , Bhbmap_donor_vz
 real(kind=8), dimension(100)      :: hbmap_hyd_x , hbmap_hyd_y , hbmap_hyd_z
 real(kind=8), dimension(100)      :: Bhbmap_hyd_x , Bhbmap_hyd_y , Bhbmap_hyd_z 
 real(kind=8), dimension(100)      :: hbmap_hyd_vx , hbmap_hyd_vy , hbmap_hyd_vz
 real(kind=8), dimension(100)      :: Bhbmap_hyd_vx , Bhbmap_hyd_vy , Bhbmap_hyd_vz
 real(kind=8), dimension(100)      :: hbmap_acc_x , hbmap_acc_y , hbmap_acc_z
 real(kind=8), dimension(100)      :: hbmap_acc_x_real , hbmap_acc_y_real , hbmap_acc_z_real 
 real(kind=8), dimension(100)      :: Bhbmap_acc_x , Bhbmap_acc_y , Bhbmap_acc_z
 real(kind=8), dimension(100)      :: hbmap_acc_vx , hbmap_acc_vy , hbmap_acc_vz
 real(kind=8), dimension(100)      :: hbmap_acc_vx_real , hbmap_acc_vy_real , hbmap_acc_vz_real
 real(kind=8), dimension(100)      :: Bhbmap_acc_vx , Bhbmap_acc_vy , Bhbmap_acc_vz
 real(kind=8), dimension(5000,50)  :: donor_coord_x , donor_coord_y , donor_coord_z
 real(kind=8), dimension(5000,50)  :: Bdonor_coord_x , Bdonor_coord_y , Bdonor_coord_z 
 
 real(kind=8), dimension(5000)   :: LJ_acc_6,LJ_acc_12,LJ_range_6,LJ_range_12,LJ_range_6b,LJ_range_12b
 real(kind=8), dimension(5000)   :: LJ_acc_6b,LJ_acc_12b
 

 character(len=4)                   :: RESSEARCH
 character(len=5) ,dimension(5000)  :: rangetype,rangetype2 
 character(len=5) ,dimension(5000)  :: RANGE_ZATOM,RANGE_ZATOM2
 character(len=5)                   :: restype,atomtype,nameatom,typeatom
 character(len=5)                   :: READFIRST,READFIRST2 , calc_type1 , calc_type2
 character(len=5) ,dimension(100)  :: hbmap_donor_type,hbmap_hyd_type,hbmap_acc_type
 character(len=5) ,dimension(100)  :: hbmap_donor_type_real
 character(len=5) ,dimension(100)  :: hbmap_acc_type_real,radius_donor
 character(len=5) ,dimension(100)  :: Bhbmap_donor_type,Bhbmap_hyd_type,Bhbmap_acc_type
 character(len=5) ,dimension(100)  :: hbmap_donor_restype,hbmap_hyd_restype,hbmap_acc_restype
 character(len=5) ,dimension(100)  :: hbmap_donor_restype2,hbmap_acc_restype2 
 character(len=5) ,dimension(100)  :: hbmap_donor_restype_real,hbmap_acc_restype_real 
 character(len=5) ,dimension(100)  :: Bhbmap_donor_restype,Bhbmap_hyd_restype,Bhbmap_acc_restype
 character(len=5) ,dimension(5000)  :: nameatomtypeb,rangetypeA2,radius_atom
 character(len=5)                   :: RANGER,ACCtype
 character(len=3)                   :: tester2
 character(len=4)                   :: tester 
 character(len=5)                   :: tester3
 character(len=2)                   :: tester4
 character(len=5) , dimension(100,100) :: acc_restype  , acc_atomtype 
 character(len=5) , dimension(5000,50) :: donor_restype, donor_atomtype, Bdonor_restype, Bdonor_atomtype
 character(len=5) , dimension(100,100) :: Bacc_restype  , Bacc_atomtype 
 character(len=5)                   :: Namedonor , Nameacc 
 character(len=11)                          :: enfilename

 character(len = 8)  :: date
 character(len =10)  :: time
 character(len=20)   :: hsearcher,hsearcher2,trjcat
 character(len = 5)  :: zone,H_bond_special,H_num_spec
 integer,dimension(8) :: values 
 integer,allocatable  :: iseed(:)
 integer              :: isize , idate(8),ntotal,nsol,changenumber

 logical              :: next,changebool,ex,ex2,changeboolean3,changebool2,exxbool,runbool

 logical              :: hazardbool,mdboolean,timeboolean

 character(len=5),dimension(:),allocatable ::  atomtypesol,restypesol
 real(kind=8) , dimension(:),allocatable   ::  coord_xsol,coord_ysol,coord_zsol,coord_test_x,coord_test_y,coord_test_z
 integer      , dimension(:),allocatable   ::  resnumsol,atomnumsol

 integer                           :: angleboolean1,angleboolean2,anglefirstbool,kmc_number

 real(kind=8)                      :: maxprotx,maxproty,maxprotz,minprotx,minproty,minprotz,diff_x1,diff_x2,diff_y1,diff_y2,diff_z1,diff_z2

 real(kind=8)                      :: diff_x,diff_y,diff_z,diff_tot,forceconstant,angle

 real(kind=8)                      :: atom_1_x,atom_2_x,atom_3_x,atom_4_x
  
 real(kind=8)                      :: atom_1_y,atom_2_y,atom_3_y,atom_4_y

 real(kind=8)                      :: atom_1_z,atom_2_z,atom_3_z,atom_4_z

 real(kind=8)                      :: maxtime

 real(kind=8)                      :: diff_x3,diff_x4

 real(kind=8)                      :: diff_y3,diff_y4

 real(kind=8)                      :: diff_z3,diff_z4 

 real(kind=8)                      :: kreuz_x1,kreuz_x2,kreuz_y1,kreuz_y2,kreuz_z1,kreuz_z2

 real(kind=8)                      :: angle_dih,angle_real

 real(kind=8)                      :: FORCE

 real(kind=8),dimension(20)        :: anglefirst1,anglefirst2

 character(len=6)                  :: ANGLESEARCH

 character(len=8)                  :: DIHANGLESEARCH

 character(len=5)                  :: atom_1,atom_2,atom_4,atom_3


 integer,dimension(:)         ,allocatable          :: atomnum_array,resnum_array
 character(len=5),dimension(:),allocatable          :: atomtype_array
 real(kind=8),dimension(:)    ,allocatable          :: coord_array_x,coord_array_y,coord_array_z
 
 real(kind=8)                                       :: N_x,N_y,N_z,O_x,O_y,O_z,H_x,H_y,H_z,hazardparamreal

 integer                                            :: count1,count2,count3,count4

 integer                                            :: resnum_O, resnum_N

 integer,dimension(50)                              :: segmentstart,segmentend,o_true,n_true

 integer                                            :: whole_res_bool,max_chain_B,mini_chain_B

 

open(unit=333,file='KMCinput9b')

read(333,*) 
read(333,*) 
read(333,*)
read(333,*) kmccycles
read(333,*)
read(333,*) timekonst_delta
read(333,*)
read(333,*) mdtime
read(333,*) 
read(333,*) 
read(333,*)
read(333,*) itplines
read(333,*)
read(333,*) rtplines
read(333,*)
read(333,*) minimline
read(333,*)
read(333,*) FOOT
read(333,*)
read(333,*) NATOMS
read(333,*)
read(333,*) cutoff
read(333,*)
read(333,*) EXPCONST
read(333,*)
read(333,*) EXPCONST2
read(333,*)
read(333,*) endacc
read(333,*)
read(333,*) endacc2
read(333,*) 
read(333,*) endtable
read(333,*)
read(333,*) endtable2
read(333,*)
read(333,*) box_x
read(333,*)
read(333,*) box_y
read(333,*)
read(333,*) box_z
read(333,*) 
read(333,*) 
read(333,*) 
read(333,*) ntotal
read(333,*)
read(333,*) donorregionstart
read(333,*)
read(333,*) alpha
read(333,*)
read(333,*) donormode
read(333,*)
read(333,*) keepresmode1
read(333,*) 
read(333,*) keepresmode2
read(333,*) 
read(333,*) Deltamax
read(333,*)
read(333,*) bmin_delta_g
read(333,*)
read(333,*) angleboolean1
read(333,*) 
read(333,*) angleboolean2
read(333,*)
read(333,*) anglefirstbool
read(333,*)
read(333,*) anglelines
read(333,*)
read(333,*) FORCE
read(333,*)
read(333,*) hazardparamreal
read(333,*)
read(333,*) maxtime
read(333,*) 
read(333,*) whole_res_bool
read(333,*) 
read(333,*) max_chain_B
read(333,*) 
read(333,*) mini_chain_B
  close(unit=333)


!     rtplines       = 10244
     hbmap_donor    = 0
     hbmap_hydrogen = 0
     hbmap_acceptor = 0
!     FOOT           = 95
!     NATOMS         = 1653
     counter        = 0
     n              = 0
     i              = 0 
     m_donor        = 0
     m_acc          = 0
     ende           = 0
     ende2          = 0
     scaling        = 0
     mum            = 0
     mum2           = 0
     hb_don_A       = 0
     hb_acc_A       = 0 
     hb_hyd_A       = 0
!     cutoff         = 3.0
     countA         = 0
     Namedonor      = 'DON'
     Nameacc        = 'ACC'
     matter         = 0 
     matter2        = 0


 call date_and_time(values=idate)

 call random_seed(size=isize)
 allocate( iseed(isize))
 call random_seed(get=iseed)
 iseed = iseed* (idate(8)-500)
 call random_seed(put=iseed)
 deallocate( iseed )


 open(unit=38,file='time')
    read(38,*) 
    write(38,*) '0'
    write(38,*) '0'
    write(38,*) '0'
       
 close(unit=38)

if(kmcbegintime==0) then

 call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_sol.mdp -c minimized_water.gro -p LOV2_bak.top -o fullmd.tpr')

endif

if(kmcbegintime.ne.0) then

 call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_sol.mdp -c md_finaloutrelax.gro -p LOV2gen.top -o fullmd.tpr')

 open(unit=38,file='time')
   
 write(38,*) '# time'
 write(38,*) kmcbegintime
 write(38,*) '0'
 write(38,*) '0'


 close(unit=38)

endif

 open(unit=111,file='TEST_RATE')

 do y=1,kmccycles 

 timer = (mdtime/1000)*(y-1) 

 RATER  = endacc
 RATER2 = endacc2

! close(unit=111)

 runbool = .true.

 call system('rm -f traj.trr mdout.mdp ener.edr md.log ')

 call system('rm -f md_final.gro') 

 call system('mpirun -np 2 /loctmp/pee18323/GRO_bin/bin/mdrun -nice 0 -v -s fullmd.tpr -o md_traj.trr -x md_traj.xtc -c md_final.gro -e md_ener.edr')
 
 call system('cp md_final.gro md_copy.gro')

 inquire(file='md_final.gro',err=55555,exist=runbool)

 55555 continue

 ex = .true.


 if(runbool==.true.) then

     open(unit=86,file='md_final.gro')

 else

     open(unit=86,file='md_copy.gro')

 endif

  hazardbool = .false.

        next = .true.

 timeboolean = .false.

! open(unit=87,file='md_start.gro')
 open(unit=87,file='minimized3.gro')
 open(unit=88888,file='s11')


 read(86,*,end=9000,err=9000)
 read(86,*,end=9000,err=9000)

 write(87,*)
 write(87,*) NATOMS
 
 do i=1,NATOMS
       
           read(86,'(i5,2a5,i5,3f8.3,3f8.4)',end=9000,err=9000) resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z
          write(87,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z
 enddo
          write(87,*) box_x,box_y,box_z

 do i=1,ntotal-natoms

           read(86,'(i5,2a5,i5,3f8.3,3f8.4)',end=9000,err=9000) resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z

         write(88888,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                             coord_y, coord_z, vel_x , vel_y , vel_z   
       
 enddo
 
 close(unit=87)
 close(unit=86)
 close(unit=88888)


9000 continue

 call system('/loctmp/pee18323/GRO_bin/bin/grompp -f minim.mdp -p LOV2.top -c minimized3.gro -o input.tpr')
 call system('/loctmp/pee18323/GRO_bin/bin/mdrun -v -s input.tpr -nice 0 -c minimized3.gro ')
 call system('rm -f ./#*# out')

! alpha-Helical parts

tu = 0

timefirst = 0

do 

 tu = tu + 1

 allocate(resnum_array(NATOMS),stat=ierr)
 allocate(atomtype_array(NATOMS),stat=ierr)
 allocate(coord_array_x(NATOMS),stat=ierr)
 allocate(coord_array_y(NATOMS),stat=ierr)
 allocate(coord_array_z(NATOMS),stat=ierr)

 allocate(atomnum_array(NATOMS),stat=ierr)


 count2  = 1


    open(unit=2,file='minimized3.gro')

    read(2,*,end=1212,err=1212)
    read(2,*,end=1212,err=1212)


    do k = 1,NATOMS
              
         read(2,'(i5,2a5,i5,3f8.3,3f8.4)',end=11212,err=11212) resnum_array(k), restype, atomtype_array(k), atomnum_array(k), coord_array_x(k), &
             &   coord_array_y(k), coord_array_z(k), vel_x , vel_y , vel_z

    enddo

    close(unit=2)

11212 continue

!write(*,*) 'Hier 1'

 count2  = 1

 count3  = 1

 o = 1

     do k = o,NATOMS                       

         if(atomnum_array(k) .gt. donorregionstart) then

                              if( atomtype_array(k) == '    O') then

                                            O_x   = coord_array_x(k)
                                            O_y   = coord_array_y(k)
                                            O_z   = coord_array_z(k)

!write(*,*) O_x,'O_x'
!write(*,*) O_y,'O_y'
!write(*,*) O_z,'O_z'

                                            i = k

!write(*,*) i,'i'

                                            resnum_O = resnum_array(k)

!write(*,*) resnum_O,'resnum_O'


                                            count2 = 0
                                      do 
                                     
                                            i = i + 1

                                            if(i.ge.NATOMS) then

                                                      
!write(*,*) 'exit'

                                                 exit

                                            endif


                                            if( atomtype_array(i) == '    N') then

                                              count2 = count2 + 1

!write(*,*)  count2,'count2'

                                                if(count2 == 4) then


                                                    N_x   = coord_array_x(i)
                                                    N_y   = coord_array_y(i)
                                                    N_z   = coord_array_z(i)


                                                    H_x   = coord_array_x(i+1)
                                                    H_y   = coord_array_y(i+1)
                                                    H_z   = coord_array_z(i+1)


                                                    resnum_N = resnum_array(i)
                                                    


                                                    exit

                                                 endif
                                                          

                                               endif

                                       enddo

                                    
           



                                     diff_x1 = H_x - N_x
                                     diff_y1 = H_y - N_y
                                     diff_z1 = H_z - N_z


                                     diff_x3 = O_x - N_x
                                     diff_y3 = O_y - N_y
                                     diff_z3 = O_z - N_z

                                     diff_tot  = sqrt((diff_x1)**2+(diff_y1)**2+(diff_z1)**2)*sqrt((diff_x3)**2+(diff_y3)**2+(diff_z3)**2)

                                     angle_real = 180/PI*acos((diff_x1*diff_x3 + diff_y1*diff_y3 + diff_z1*diff_z3)/diff_tot)    

!write(*,*)   'angle',angle_real

                                     diff_x1 = N_x - O_x
                                     diff_y1 = N_y - O_y
                                     diff_z1 = N_z - O_z

                                     diff_tot  = sqrt(diff_x1**2 + diff_y1**2 + diff_z1**2)


 
                              if(diff_tot .le. 0.35 .and. angle_real .le. 35) then

!write(*,*)  O_true(count3),'O_true(count3)',N_true(count3),'N_true(count3)'

!write(*,*)  count3,'count3'

                                     O_true(count3) = resnum_O
                                     N_true(count3) = resnum_N

                                     count3 = count3 + 1 


                              endif


                                     o = i
       

                            endif

                  endif

    enddo

    count4 = 1

    n = 1

 count3 = count3 - 1



exxbool = .false.

do 

  if(count3 .ge. 1) then

    do i  = n,count3



                                  count1 = O_true(i)

                                  count2 = N_true(i)


                                               


              if(O_true(i) + 1 .eq. O_true(i+1).and. 1.le.count3-i) then

                                  count1 = O_true(i) 


                                  count2   = N_true(i+1)



              endif
                                                          

              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2).and. 2 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+2)



              endif

              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3).and. 3 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+3)


              endif
              
              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3) &
                     & .and. O_true(i+3)+ 1 .eq. O_true(i+4) .and. 4 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+4)


              endif

              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3) &
                     & .and. O_true(i+3)+ 1 .eq. O_true(i+4) .and. O_true(i+4)+ 1 .eq. O_true(i+5) .and. 5 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+5)


              endif 

              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3) &
                     & .and. O_true(i+3)+ 1 .eq. O_true(i+4) .and. O_true(i+4)+ 1 .eq. O_true(i+5) .and. O_true(i+5) + 1 .eq. O_true(i+6) .and. 6 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+6)



              endif 

              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3) &
                     & .and. O_true(i+3)+ 1 .eq. O_true(i+4) .and. O_true(i+4)+ 1 .eq. O_true(i+5) .and. O_true(i+5) + 1 .eq. O_true(i+6) &
                     & .and. O_true(i+6)+ 1 .eq. O_true(i+7) .and. 7 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+7)




              endif 
              if(O_true(i) + 1 .eq. O_true(i+1) .and.O_true(i+1) + 1 .eq. O_true(i+2) .and. O_true(i+2) + 1 .eq. O_true(i+3) &
                     & .and. O_true(i+3)+ 1 .eq. O_true(i+4) .and. O_true(i+4)+ 1 .eq. O_true(i+5) .and. O_true(i+5) + 1 .eq. O_true(i+6) &
                     & .and. O_true(i+6)+ 1 .eq. O_true(i+7) .and. O_true(i+7)+ 1 .eq. O_true(i+8) .and. 8 .le. count3-i) then

                                  count1 = O_true(i)

                                  count2 = N_true(i+8)



              endif 

              if(count2==N_true(i)) then

                                    n = n + 1

                                    exit

              endif
              if(count2==N_true(i+1)) then

                                    n = n + 2

                                    exit

              endif
              if(count2==N_true(i+2)) then

                                    n = n + 3

                                    exit

              endif
              if(count2==N_true(i+3)) then

                                    n = n + 4
  
                                    exit

              endif
              if(count2==N_true(i+4)) then

                                    n = n + 5

                                    exit

              endif
              if(count2==N_true(i+5)) then

                                    n = n + 6

                                    exit

              endif
              if(count2==N_true(i+6)) then

                                    n = n + 7

                                    exit

              endif
              if(count2==N_true(i+7)) then

                                    n = n + 8

                                    exit

              endif
              if(count2==N_true(i+8)) then

                                    n = n + 9

                                    exit

              endif

              if(count2.ge.N_true(count3)) then

                                    exxbool = .true.

                                    exit

              endif

     enddo 
              segmentstart(count4) = count1

              segmentend(count4)   = count2


              if(segmentend(count4) == N_true(count3) .or. exxbool == .true. ) then

                                exit

              endif

              count4 = count4 + 1

     else

             exit

     endif

enddo

!do i=1,count4

!           write(*,*) segmentstart(i),'start'
!           write(*,*) segmentend(i),'end'

!enddo

               
 deallocate(resnum_array,stat=ierr)
 deallocate(atomtype_array,stat=ierr)
 deallocate(coord_array_x,stat=ierr)
 deallocate(coord_array_y,stat=ierr)
 deallocate(coord_array_z,stat=ierr)

 deallocate(atomnum_array,stat=ierr)

!write(*,*) 'Hier 2'
 

 call system('rm -f table_hbond_B.ndx')

 call system('./H_search_new.sh')

 inquire(file='table_hbond_B.ndx',exist=ex) 


 open(unit=334,file='table_hbond_B.ndx')

 do i=1,endtable

     read(334,*,end=111,err=111)
 
 enddo

k = 0

 do
     k=k+1
     read(334,*,end=111,err=111) hbmap_donor(k) , hbmap_hydrogen(k) , hbmap_acceptor(k)
    
 enddo

111 continue
 
  close(unit = 334)

 
ende=k-1


 if(ex ==.false.) then

    ende = 0

 endif


if(donormode==1) then

 do j=1,ende

           if(hbmap_donor(j) .le. donorregionstart ) then

              a = hbmap_donor(j)
              hbmap_donor(j) = hbmap_acceptor(j)
              hbmap_acceptor(j)    = a

           endif

 enddo 

endif
 
 do l=1,ende
 
     do k=1,endacc

            TOTAL_EN(l,k)=0
     enddo

 enddo

 a =  0

 do j=1,ende 

  open(unit=2,file='minimized3.gro')

  read(2,*,end=1212,err=1212)
  read(2,*,end=1212,err=1212)


    do k = 1,NATOMS
              
    read(2,'(i5,2a5,i5,3f8.3,3f8.4)',end=1212,err=1212) resnum, restype, atomtype, atomnum, coord_x, &
         &   coord_y, coord_z, vel_x , vel_y , vel_z
         
             if(atomnum.eq.hbmap_donor(j)) then
               hbmap_donor_resnum(j) = resnum
               hbmap_donor_restype(j) = restype
               hbmap_donor_type(j) = atomtype               
               hbmap_donor_num(j)  = atomnum
               hbmap_donor_x(j)    = coord_x
               hbmap_donor_y(j)    = coord_y
               hbmap_donor_z(j)    = coord_z
               hbmap_donor_vx(j)      = vel_x
               hbmap_donor_vy(j)      = vel_y
               hbmap_donor_vz(j)      = vel_z
             endif
             if(atomnum.eq.hbmap_hydrogen(j)) then
               hbmap_hyd_resnum(j) = resnum
               hbmap_hyd_restype(j) = restype
               hbmap_hyd_type(j) = atomtype
               hbmap_hyd_num(j)  = atomnum
               hbmap_hyd_x(j)    = coord_x
               hbmap_hyd_y(j)    = coord_y
               hbmap_hyd_z(j)    = coord_z
               hbmap_hyd_vx(j)      = vel_x
               hbmap_hyd_vy(j)      = vel_y
               hbmap_hyd_vz(j)      = vel_z
             endif
             if(atomnum.eq.hbmap_acceptor(j)) then
               hbmap_acc_resnum(j) = resnum
               hbmap_acc_restype(j) = restype
               hbmap_acc_type(j) = atomtype
               hbmap_acc_num(j)  = atomnum
               hbmap_acc_x(j)    = coord_x
               hbmap_acc_y(j)    = coord_y
               hbmap_acc_z(j)    = coord_z
               hbmap_acc_vx(j)      = vel_x
               hbmap_acc_vy(j)      = vel_y
               hbmap_acc_vz(j)      = vel_z
             endif               
                              
     enddo          

 close(unit=2)
  
 enddo


1212 continue

j = 0

 do i=1,ende
              if(hbmap_acc_resnum(i).ne.hbmap_donor_resnum(i).and.hbmap_acc_type(i).ne.'   XX' .and. &
                 & hbmap_acc_type(i).ne.'  XXX' .and. hbmap_donor_type(i) .ne. '   XX' .and. &
                 & hbmap_donor_type(i) .ne. '  XXX') then
               j = j + 1
               hbmap_acc_resnum_real(j) = hbmap_acc_resnum(i)
               hbmap_acc_restype_real(j)= hbmap_acc_restype(i)
               hbmap_acc_type_real(j) =   hbmap_acc_type(i)  
               hbmap_acc_num_real(j)  =   hbmap_acc_num(i)  
               hbmap_acc_x_real(j)    =   hbmap_acc_x(i)
               hbmap_acc_y_real(j)    =   hbmap_acc_y(i)
               hbmap_acc_z_real(j)    =    hbmap_acc_z(i) 
               hbmap_acc_vx_real(j)      =  hbmap_acc_vx(i)
               hbmap_acc_vy_real(j)      =  hbmap_acc_vy(i)
               hbmap_acc_vz_real(j)      =  hbmap_acc_vz(i)

               hbmap_donor_resnum_real(j) = hbmap_donor_resnum(i)
               hbmap_donor_restype_real(j) = hbmap_donor_restype(i)
               hbmap_donor_type_real(j) =  hbmap_donor_type(i)              
               hbmap_donor_num_real(j)  =  hbmap_donor_num(i)
               hbmap_donor_x_real(j)    =  hbmap_donor_x(i)
               hbmap_donor_y_real(j)    =  hbmap_donor_y(i)
               hbmap_donor_z_real(j)    =  hbmap_donor_z(i)
               hbmap_donor_vx_real(j)      = hbmap_donor_vx(i)
               hbmap_donor_vy_real(j)      = hbmap_donor_vy(i)
               hbmap_donor_vz_real(j)      = hbmap_donor_vz(i)
              endif
  enddo
              
ende = j                 

if(whole_res_bool == 0) then

 do i=1,ende

        resnum_min_range(i) = hbmap_donor_resnum_real(i) - 2

        resnum_max_range(i) = hbmap_donor_resnum_real(i) + 2

 enddo

 do i=1,ende

    do j=1,count4

       if(hbmap_donor_resnum_real(i) .ge. segmentstart(j) .and. hbmap_donor_resnum_real(i) .le. segmentend(j)) then

                                          resnum_min_range(i) = segmentstart(j)

 !                                         write(*,*) resnum_min_range(i),'resnum_min_range(i)'


                                          resnum_max_range(i) = segmentend(j)

 !                                         write(*,*) resnum_max_range(i),'resnum_max_range(i)'

       endif

    enddo

 enddo 

else

  do i=1,ende

       resnum_min_range(i) = mini_chain_B

       resnum_max_range(i) = max_chain_B

  enddo

endif


 do i=1,ende

        m_acc   = 1
        m_donor = 1

        open(unit=3,file='minimized3.gro')

        read(3,*,end=9001,err=9001)
        read(3,*,end=9001,err=9001)
    
        do j=1,NATOMS
        
        read(3,'(i5,2a5,i5,3f8.3,3f8.4)',end=9001,err=9001) resnum, restype, atomtype, atomnum, coord_x, &
                                           coord_y, coord_z, vel_x , vel_y , vel_z


           if (resnum.ge.resnum_min_range(i) .and. &
             & resnum.le.resnum_max_range(i)) then

              donor_resnum(m_donor,i)  = resnum
              donor_restype(m_donor,i) = restype
              donor_atomtype(m_donor,i)= atomtype
              donor_atomnum(m_donor,i) = atomnum
              donor_coord_x(m_donor,i) = coord_x
              donor_coord_y(m_donor,i) = coord_y
              donor_coord_z(m_donor,i) = coord_z

              m_donor= m_donor+1


          endif

        enddo

        close(unit=3)

9001 continue

 
        m_don3(i) = m_donor


        open(unit=4,file='minimized3.gro')
        read(4,*,end=9002,err=9002)
        read(4,*,end=9002,err=9002)
    
        do j=1,NATOMS
        
        read(4,'(i5,2a5,i5,3f8.3,3f8.4)',end=9002,err=9002) resnum, restype, atomtype, atomnum, coord_x, &
                                           coord_y, coord_z, vel_x , vel_y , vel_z    
           if (resnum.eq.hbmap_acc_resnum_real(i)) then 
               
              acc_resnum(m_acc,i)   = resnum
              acc_restype(m_acc,i)  = restype
              acc_atomnum(m_acc,i)  = atomnum
              acc_atomtype(m_acc,i) = atomtype
         
        m_acc=m_acc+1
        

          endif
           
        enddo
       close(unit=4)

9002 continue

       m_acc3(i) = m_acc
       
enddo
 


                                               open(unit=5,file='total_Energy1_x')
                                               open(unit=6,file='total_Energy2_x')
                                               open(unit=7,file='total_Energy3_x')
                                               open(unit=8,file='total_Energy4_x')
                                               open(unit=9,file='total_Energy5_x')
                                               open(unit=10,file='total_Energy6_x')
                                               open(unit=11,file='total_Energy7_x')
                                               open(unit=12,file='total_Energy8_x')
                                               open(unit=13,file='total_Energy9_x')
                                               open(unit=14,file='total_Energy10_x')
                                               open(unit=15,file='total_Energy11_x')
                                               open(unit=16,file='total_Energy12_x')
                                               open(unit=17,file='total_Energy13_x')
                                               open(unit=18,file='total_Energy14_x')
                                               open(unit=889,file='total_Energy15_x')
                                               open(unit=890,file='total_Energy16_x')
                                               open(unit=891,file='total_Energy17_x')
                                               open(unit=892,file='total_Energy18_x')
                                               open(unit=893,file='total_Energy19_x')
                                               open(unit=894,file='total_Energy20_x') 
                                               open(unit=895,file='total_Energy21_x')
                                               open(unit=896,file='total_Energy22_x')
                                               open(unit=897,file='total_Energy23_x')                                               
                                               open(unit=898,file='total_Energy24_x')
!                                               open(unit=888,file='testerCoul1.9')
                                               
DIST_scale_x = 0
DIST_scale_y = 0 
DIST_scale_z = 0
DIST_scale_A = 0

do l=1,ende
              
           DIST_scale_x = hbmap_acc_x_real(l) - hbmap_donor_x_real(l)
           DIST_scale_y = hbmap_acc_y_real(l) - hbmap_donor_y_real(l)
           DIST_scale_z = hbmap_acc_z_real(l) - hbmap_donor_z_real(l)

           DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
           scalekey(l)  = EXPCONST  / (endacc*DIST_scale_A)

!write(8881,*) scalekey(l),'scalekey'
!write(8881,*) DIST_scale_A,'DIST_scale_A'

enddo

m = 1



do l=1,ende

    m_acc2 = m_acc3(l)-1
    m_don2 = m_don3(l)-1
           
    scaling = 0 

     do k = 1,endacc

              if (k == 1) then
                  
               do p=1,m_don3(l)-1
                  
                  new_dir_acc_x(p,l) = scaling*(hbmap_donor_x_real(l)-hbmap_acc_x_real(l))+donor_coord_x(p,l) 
                  new_dir_acc_y(p,l) = scaling*(hbmap_donor_y_real(l)-hbmap_acc_y_real(l))+donor_coord_y(p,l) 
                  new_dir_acc_z(p,l) = scaling*(hbmap_donor_z_real(l)-hbmap_acc_z_real(l))+donor_coord_z(p,l)

               enddo
                           
             

              open(unit=19,file='minimized3.gro')
                
                       read(19,*,end=9003,err=9003)
                       read(19,*,end=9003,err=9003)

              open(unit=20,file='md_tester.gro')

                      write(20,*)
                      write(20,*) NATOMS
             
              p = 1

                       do m=1,NATOMS
                                

                       read(19,'(i5,2a5,i5,3f8.3,3f8.4)',end=9003,err=9003) resnum, restype, atomtype, atomnum, coord_x, &
                                         coord_y, coord_z, vel_x , vel_y , vel_z

                       
                       if(resnum .eq. donor_resnum(p,l) .and.p.le.m_don3(l)-1) then
                               
                                      

                            coord_x = new_dir_acc_x(p,l)
                            coord_y = new_dir_acc_y(p,l)
                            coord_z = new_dir_acc_z(p,l)

                            restype = Nameacc
                                       
                            p = p + 1
                                    
                        endif

                        if(resnum.eq.donor_resnum(p,l)) then

                            restype = Namedonor
                            
                     
      !                      write(*,*) m ,'m'
                        endif
            
              write(20,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                         coord_y, coord_z, vel_x , vel_y , vel_z
                                                                               
                        enddo

              close(unit=19)

9003 continue
                               

         write(20,*) box_x,box_y,box_z
              
              close(unit=20)
  
      endif
     

! wegen des Einsetzens in minimized.gro

      if(k.gt.1) then   


    
        do p=1,m_don3(l)-1
           
              
        new_dir_acc_x(p,l) = scaling*(hbmap_donor_x_real(l)-hbmap_acc_x_real(l))+donor_coord_x(p,l) 
        new_dir_acc_y(p,l) = scaling*(hbmap_donor_y_real(l)-hbmap_acc_y_real(l))+donor_coord_y(p,l) 
        new_dir_acc_z(p,l) = scaling*(hbmap_donor_z_real(l)-hbmap_acc_z_real(l))+donor_coord_z(p,l)
        
        enddo 
             
        open(unit=21,file='minimized.gro')

        read(21,*,end=9004,err=9004)
        read(21,*,end=9004,err=9004)
        
        open(unit=22,file='md_tester.gro')
     
        write(22,*)
        write(22,*) NATOMS

        p = 1
       
            do m=1,NATOMS
 

               read(21,'(i5,2a5,i5,3f8.3,3f8.4)',end=9004,err=9004) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

!               if(resnum.eq.hbmap_donor_resnum_real(l) .and. p.le.m_don3(l)-1) then
                       
                if(resnum.eq. donor_resnum(p,l) .and.p.le.m_don3(l)-1) then

                   coord_x = new_dir_acc_x(p,l)
                   coord_y = new_dir_acc_y(p,l)
                   coord_z = new_dir_acc_z(p,l)

                   p = p + 1


!               write(*,*) m ,'m'

               endif

               write(22,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(22,*)  box_x,box_y,box_z

               close(unit=21)
               close(unit=22)

9004  continue

       endif
          

       call system ('rm -f indexsim.ndx')
       
      
       open(unit=23,file='temp.ndx')
       
            write(23,'(a7)') '[ ACC ]'
        
          
           write(23,*)  hbmap_acc_num_real(l)         
  
  
           write(23,'(a7)') '[ DON ]'
      
           if(keepresmode1== 0 ) then

                write(23,*) hbmap_donor_num_real(l)

           else

              do p=1,m_don3(l)-1

                 write(23,*) donor_atomnum(p,l)

              enddo 
      
           endif           

       close(unit=23)
      
       m = m +1
      
       call system('cat temp.ndx index.ndx > indexsim.ndx')
       call system('rm -f temp.ndx')
        
       call system('cp -f indexsim.ndx testindex.ndx')
!99     continue
            
       call system ('rm -f minimized.gro ener.edr traj.trr mdout.mdp md.log')
       

       if(k==1) then
            
            call system ('cp -f minim.mdp minimtest.mdp') 
    
            open(unit=24,file='minimtest.mdp')

                 do j = 1,minimline
                      
                          read(24,*)
                 enddo
             write(24,*) 'freezegrps        =  ' , Nameacc , Namedonor
             write(24,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=24)

       endif
            
!         write(*,*) 'everything is ok'

         call system ('rm -f input.tpr minimized.gro traj.trr ener.edr md.log mdout.mdp')       
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest.mdp -c md_tester.gro -p LOV2.top -o input.tpr -n indexsim.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized_A.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp out')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro out')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized.gro -p LOV2.top -o input.tpr -n indexsim.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f ./#*# out')
      
         AC_X     = 0
         AC_Y     = 0
         AC_Z     = 0
         DO_X     = 0
         DO_Y     = 0
         DO_Z     = 0
         DIST_X   = 0
         DIST_Y   = 0
         DIST_Z   = 0
         DIST_GRAPH_X = 0
         DIST_GRAPH_Y = 0
         DIST_GRAPH_Z = 0

         open(unit=25,file='minimized.gro')

         read(25,*,end=9005,err=9005)
         read(25,*,end=9005,err=9005)

        do b=1,NATOMS

           read(25,'(i5,2a5,i5,3f8.3,3f8.4)',end=9005,err=9005) resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
           if(atomnum.eq.hbmap_acc_num_real(l)) then
                 AC_X = coord_x 
                 AC_Y = coord_y
                 AC_Z = coord_z
           endif
                   
           if(atomnum.eq.hbmap_donor_num_real(l)) then
                 DO_X = coord_x
                 DO_Y = coord_y
                 DO_Z = coord_z
           endif
!Distanz nicht auf Null skaliert!!
             
           DIST_X = (AC_X - DO_X) -  (hbmap_acc_x_real(l)-hbmap_donor_x_real(l))
           DIST_Y = (AC_Y - DO_Y) -  (hbmap_acc_y_real(l)-hbmap_donor_y_real(l))
           DIST_Z = (AC_Z - DO_Z) -  (hbmap_acc_z_real(l)-hbmap_donor_z_real(l))



           DISTANCE2(l,k) = sqrt(DIST_X**2+DIST_Y**2+DIST_Z**2)

           DIST_GRAPH_X = AC_X - DO_X
           DIST_GRAPH_Y = AC_Y - DO_Y
           DIST_GRAPH_Z = AC_Z - DO_Z

           DIST_GRAPH(l,k) = sqrt(DIST_GRAPH_X**2+DIST_GRAPH_Y**2+DIST_GRAPH_Z**2)                       
      
      enddo 
      close(unit=25)

9005  continue
        
      write(*,*)  hbmap_acc_num_real(l),'acc',hbmap_donor_num_real(l),'don',' A-cycle', k ,'ende',l
      write(*,*)   DISTANCE2(l,k),'distance2'
       


      scaling = scaling + scalekey(l)

                  
                       call system ('rm -f md_tester.gro')
                      
                               dist_acc_x = 0
                               dist_acc_y = 0
                               dist_acc_z = 0
                             range_betrag = 0 
                                     
                                     mum =  1


        open(unit =26,file='minimized.gro')

        read(26,*,end=9006,err=9006)
        read(26,*,end=9006,err=9006)
          
        do i=1,NATOMS

          read(26,'(i5,2a5,i5,3f8.3,3f8.4)',end=9006,err=9006) resnum, restype, atomtype, atomnum, coord_x, &
          coord_y, coord_z, vel_x , vel_y , vel_z

          if(resnum .ne. hbmap_donor_resnum_real(l) ) then

                dist_acc_x = coord_x - DO_X
                dist_acc_y = coord_y - DO_Y    
                dist_acc_z = coord_z - DO_Z

                range_betrag = sqrt(dist_acc_x**2 + dist_acc_y**2 + dist_acc_z**2)

               if(range_betrag.lt.cutoff .and. range_betrag .ne. 0.0) then
                            Cut_range(mum)       = range_betrag
                            rangetype(mum)       = restype
                            RANGE_ZATOM(mum)   = atomtype

           
                            mum = mum+1          
                             
               endif
         endif
         
       enddo

       close ( unit=26 )

9006   continue

matter = mum -1




n = 1

       do i=1,matter

           open(unit=27,file='ffG43a1B.rtp')

            do  j=1,rtplines
   
                read(27,'(2X,a4)') RESSEARCH
                
                          if(RESSEARCH.eq.rangetype(i)) then
                                read(27,*)
                                do  
                                  read(27,'(a5,1X,a5,4X,f8.5)',end=690,err=690) nameatom , typeatom,&
                                     & charge2

                                  if(nameatom.eq.RANGE_ZATOM(i)) then
                                       charge_atom(i) = charge2
!                                       write(888,*) charge_atom(n) , 'chargeatom'
                                       radius_atom(i) = typeatom
                    
                                       n = n + 1 
                                 
                                       exit
                                   endif
                                enddo
690                             exit
                           endif
             enddo

            close(unit=27)
       enddo 


if(n.ne.matter) then

!         write(*,*) 'matter',matter,'n',n   

endif

       
            open(unit=28,file='ffG43a1B.rtp')

             do j=1,rtplines

                read(28,'(2X,a4)') RESSEARCH                

                if(RESSEARCH.eq.hbmap_donor_restype_real(l)) then
                                read(28,*)
                                      do 
                                          read(28,'(a5,1X,a5,4X,f8.5)',end=1620,err=1620) nameatom , typeatom,&
                                          & charge2
                                          if(nameatom.eq.hbmap_donor_type_real(l)) then
                                             hb_acc_charge(l)  = charge2
                                             radius_donor(l)   = typeatom                                        

                                             exit
                                          endif 
                                       enddo

1620                  exit
                endif
                
            enddo
           
          close(unit=28)

n = 1

angle_en = 0

if(angleboolean1 == 1) then

            open(unit=28,file='angles')

             do j=1,anglelines

                read(28,'(2X,a4)') RESSEARCH                

                if(RESSEARCH.eq.hbmap_donor_restype_real(l)) then
                               
                                  do 
                                          read(28,'(3X,a6)',end=16251,err=16251) ANGLESEARCH
                                          
                                          if(ANGLESEARCH.eq.'angles') then

                                            
                                             do 

                                               read(28,'(a5,1X,a5,1X,a5,5X,f6.2,6X,f6.2)',end=16251,err=16251) atom_1 , atom_2 , atom_3 , angle , forceconstant                                       


                                               if(atom_1 == hbmap_donor_type_real(l) .or. &
                                                  atom_2 == hbmap_donor_type_real(l) .or. &
                                                  atom_3 == hbmap_donor_type_real(l) ) then


                                                            open(unit =26,file='minimized.gro')

                                                            read(26,*)
                                                            read(26,*)
          
                                                            do i=1,NATOMS

                                                               read(26,'(i5,2a5,i5,3f8.3,3f8.4)',end=10106,err=10106) resnum, restype, atomtype, atomnum, coord_x, &
                                                                                                                      coord_y, coord_z, vel_x , vel_y , vel_z
  
                                                               if(resnum == hbmap_donor_resnum_real(l)) then                                                                
                                     
                                                                             
                                                               if(atom_1 == atomtype) then

                                                                            atom_1_x = coord_x               
                                                                            atom_1_y = coord_y
                                                                            atom_1_z = coord_z

                                                               endif

                                                               if(atom_2 == atomtype) then

                                                                            atom_2_x = coord_x               
                                                                            atom_2_y = coord_y
                                                                            atom_2_z = coord_z

                                                               endif

                                                               if(atom_3 == atomtype) then

                                                                            atom_3_x = coord_x               
                                                                            atom_3_y = coord_y
                                                                            atom_3_z = coord_z

                                                               endif

                                                               endif

                                                            enddo

10106                                                          continue

                                                               close(unit=26)

                                                               diff_x1 = atom_1_x - atom_2_x
                                                               diff_y1 = atom_1_y - atom_2_y
                                                               diff_z1 = atom_1_z - atom_2_z


                                                               diff_x3 = atom_3_x - atom_2_x
                                                               diff_y3 = atom_3_y - atom_2_y
                                                               diff_z3 = atom_3_z - atom_2_z

                                                             diff_tot  = sqrt((diff_x1)**2+(diff_y1)**2+(diff_z1)**2)*sqrt((diff_x3)**2+(diff_y3)**2+(diff_z3)**2)

                                                            angle_real = 180/PI*acos((diff_x1*diff_x3 + diff_y1*diff_y3 + diff_z1*diff_z3)/diff_tot)

!                                                            write(*,*) angle_real,'angle'

                                                            if(k==1 .and. anglefirstbool ==1) then

                                                                       anglefirst1(n) = angle_real

                                                            endif

                                                            if(anglefirstbool == 1) then

                                                                       angle = anglefirst1(n)

                                                            endif 


                                                            angle_en = angle_en + 0.5*forceconstant*(cos(angle_real*PI/180)-cos(angle*PI/180))**2


                                                            n = n + 1
                                             
                                          endif 
                                       enddo

                                 endif
                            enddo
16251                  exit
                endif
                
            enddo
           
          close(unit=28)
          close(unit=26)

i = 1

            open(unit=28,file='angles')

             do j=1,anglelines

                read(28,'(2X,a4)') RESSEARCH                


                if(RESSEARCH.eq.hbmap_donor_restype_real(l)) then

                                read(28,*)

                                      do 
                                          read(28,'(3X,a8)',end=1624,err=1624) DIHANGLESEARCH
                                          
                                          if(DIHANGLESEARCH.eq.'dihangle') then
             
                             
                                                   exxbool = .false.

                                             do 

                                               read(28,'(a5,1X,a5,1X,a5,1X,a5,5X,f7.3,7X,f4.1,4x,i1)',end=1624,err=1624) atom_1 , atom_2 , atom_3 , atom_4, angle , forceconstant,o                                       


                                               if(angle == 0 .and. forceconstant == 0 .and. o == 0) then

                                                   exxbool = .true.

                                               endif

                                               if(atom_1 == hbmap_donor_type_real(l) .or. &
                                                  atom_2 == hbmap_donor_type_real(l) .or. &
                                                  atom_3 == hbmap_donor_type_real(l) .or. & 
                                                  atom_4 == hbmap_donor_type_real(l)) then

                                                            open(unit =26,file='minimized.gro')

                                                            read(26,*)
                                                            read(26,*)
          
                                                            do i=1,NATOMS

                                                               read(26,'(i5,2a5,i5,3f8.3,3f8.4)',end=10006,err=10006) resnum, restype, atomtype, atomnum, coord_x, &
                                                                                                                      coord_y, coord_z, vel_x , vel_y , vel_z
                                                                                                     
                                                               if(resnum == hbmap_donor_resnum_real(l)) then                                                                 
            
                                                               if(atom_1 == atomtype) then

                                                                            atom_1_x = coord_x               
                                                                            atom_1_y = coord_y
                                                                            atom_1_z = coord_z


                                                               endif

                                                               if(atom_2 == atomtype) then

                                                                            atom_2_x = coord_x               
                                                                            atom_2_y = coord_y
                                                                            atom_2_z = coord_z


                                                               endif

                                                               if(atom_3 == atomtype) then

                                                                            atom_3_x = coord_x               
                                                                            atom_3_y = coord_y
                                                                            atom_3_z = coord_z




                                                               endif

                                                               if(atom_4 == atomtype) then

                                                                            atom_4_x = coord_x               
                                                                            atom_4_y = coord_y
                                                                            atom_4_z = coord_z



                                                               endif

                                                               endif

                                                            enddo

                                                               10006 continue

                                                               close(unit=26)
                                                               
                                        !                      Bildung der Normalen

                                                               diff_x1 = atom_1_x - atom_2_x
                                                               diff_y1 = atom_1_y - atom_2_y
                                                               diff_z1 = atom_1_z - atom_2_z

                                                               diff_x2 = atom_2_x - atom_3_x                        
                                                               diff_y2 = atom_2_y - atom_3_y
                                                               diff_z2 = atom_2_z - atom_3_z

                                                               kreuz_x1 = diff_y1*diff_z2 - diff_z1*diff_y2
                                                               kreuz_y1 = diff_z1*diff_x2 - diff_x1*diff_z2
                                                               kreuz_z1 = diff_x1*diff_y2 - diff_y1*diff_x2

                                                               diff_x3 = atom_2_x - atom_3_x
                                                               diff_y3 = atom_2_y - atom_3_y
                                                               diff_z3 = atom_2_z - atom_3_z

                                                               diff_x4 = atom_3_x - atom_4_x                        
                                                               diff_y4 = atom_3_y - atom_4_y
                                                               diff_z4 = atom_3_z - atom_4_z

                                                               kreuz_x2 = diff_y3*diff_z4 - diff_z3*diff_y4
                                                               kreuz_y2 = diff_z3*diff_x4 - diff_x3*diff_z4
                                                               kreuz_z2 = diff_x3*diff_y4 - diff_y3*diff_x4

                                                                 
                                                              diff_tot  = sqrt((kreuz_x1)**2+(kreuz_y1)**2+(kreuz_z1)**2)*sqrt((kreuz_x2)**2+(kreuz_y2)**2+(kreuz_z2)**2) 

                                                              angle_dih = 180/PI*acos((kreuz_x1*kreuz_x2 + kreuz_y1*kreuz_y2 + kreuz_z1*kreuz_z2)/diff_tot) 

                                                               if(k==1 .and. anglefirstbool ==1) then

                                                                        anglefirst2(i) = angle_dih

                                                               endif

                                                               if(anglefirstbool == 1) then

                                                                       angle = anglefirst2(i)

                                                               endif 


                                                              angle_en = angle_en + forceconstant*(1+cos((o*angle_dih-angle)*PI/180))

                                                              i = i + 1

                                              endif 
                                       
                                     enddo

                               endif

                               if(exxbool == .true.) then

                                              exit

                               endif

                       enddo

1624             continue

                 exit
                endif
                
            enddo
           
          close(unit=28)
          close(unit=26)

endif
     
n = 1
                      open(unit=98989,file='ffG43a1nb.itp')

                      read(98989,*)
                      read(98989,*)

                      do p=1,itplines

                             read(98989,'(a5,32x,g13.10,g14.11)') tester3 , LJ_param_6, LJ_param_12

                             if(tester3 == radius_donor(l)) then

                                   LJ_acc_6(l) = LJ_param_6

                                   LJ_acc_12(l)= LJ_param_12

!                                   write(*,*) LJ_acc_6(l),'6acc'
!                                   write(*,*) LJ_acc_12(l),'12acc'
!                                   write(*,*) tester3,'t3',hbmap_donor_type_real(l),'rangeacc'
 
                                   n = n+1

                             endif

                      enddo

                      close(unit=98989)
 

    
        do i=1,matter       


                      open(unit=98989,file='ffG43a1nb.itp')

                      read(98989,*)
                      read(98989,*)

                             exxbool = .false.

                      do p=1,itplines

                             read(98989,'(a5,32x,g13.10,g14.11)') tester3 , LJ_param_6, LJ_param_12

                             if(tester3 == radius_atom(i)) then

                                   LJ_range_6(i) = LJ_param_6

                                   LJ_range_12(i)= LJ_param_12

                            !       write(*,*) LJ_range_6(i),'6'
                            !       write(*,*) LJ_range_12(i),'12'
                            !       write(*,*) tester3,'t3',range_zatom(i),'range'

                                   n = n+1

                                   exxbool = .true.


                             endif

                             if(p==itplines.and.exxbool==.false.) then

!                                   write(*,*) RANGE_ZATOM(i),'not taken'


                             endif

                      enddo

                      close(unit=98989)

           enddo
        
TOTALCHARGE = 0
   LJ    = 0 
  
                          factor_A = -((alpha+4)*cutoff - (alpha+1)*0)/(cutoff**(alpha+2)*(cutoff)**2)

                          factor_B =  ((alpha+3)*cutoff - (alpha+1)*0)/(cutoff**(alpha+2)*(cutoff)**3)

                          factor_C =  1/cutoff**alpha-factor_A/3*(cutoff)**3-factor_B/4*(cutoff)**4




      do i=1,matter 

                   if( Cut_range(i).gt.0.15 ) then

                         Coulomb(i) = COUL_CONST*hb_acc_charge(l)*charge_atom(i)*(1/Cut_range(i)**(alpha) - factor_A/3*(Cut_range(i))**3 &
                                       & - factor_B/4*(Cut_range(i))**4 -factor_C) 


                                          
                        TOTALCHARGE = TOTALCHARGE + Coulomb(i)

                   endif

                     
                   if(Cut_range(i).gt.0.15)then


                       LJ_6b(i)  = sqrt(abs(LJ_range_6(i)*LJ_acc_6(l)))

                       LJ_12b(i) = sqrt(abs(LJ_range_12(i)*LJ_acc_12(l)))

                       
                       LJ_potb(i) = LJ_12b(i)/Cut_range(i)**12 - LJ_6b(i)/Cut_range(i)**6



                       if(Cut_range(i).gt.cutoff) then

                          LJ_potb(i) = 0

                       endif
                                            
                       LJ = LJ + LJ_potb(i)

                   endif
                     
         enddo     

 
                       
                       TOTAL_EN(l,k) = TOTAL_EN(l,k) + LJ + TOTALCHARGE + angle_en - FORCE*k/endacc
                     
                   !   if(y == 1 .and. k==1 .and. l==1) then

                   !      open(unit=56565,file='Energytest')

                   !   endif

                   !   write(56565,*)  k,l,'k,l,y'
                   !   write(56565,*) TOTAL_En(l,k) , 'total_energy'
                   !   write(56565,*) '#################' 
                       
                      
  enddo
           call system ('rm -f minimtest.mdp')

enddo      


do l = 1,ende

m = 0
maxbereich = 0 

         firstmax(l)     = TOTAL_EN(l,1)
        dist_firstmax(l) = DISTANCE2(l,1)


     do k=1,endacc

        if(TOTAL_EN(l,k).gt.firstmax(l)) then

                 firstmax(l)  = TOTAL_EN(l,k)
             dist_firstmax(l) = Distance2(l,k)
                          
                           m  = k

        endif
     
     enddo

             secondmin(l)      = TOTAL_EN(l,m)
        dist_secondmin(l)      = DISTANCE2(l,m)


    do k=m,endacc


         if(TOTAL_EN(l,k).lt.secondmin(l).and.distance2(l,k).gt.dist_firstmax(l)) then

                  secondmin(l) = TOTAL_EN(l,k)
             dist_secondmin(l) = DISTANCE2(l,k)

         endif

    enddo
       
             maxi_g(l)         = firstmax(l)
             mini_a(l)         = total_EN(l,1)
             exposi(l)         = dist_secondmin(l)
             mini_g(l)         = secondmin(l)
          dist_maxi(l)         = dist_firstmax(l)

         if(dist_maxi(l)==0 .and. exposi(l) == 0) then

                    firstmini(1)    =  total_EN(l,1)
                 dist_firstmini(1)  =  distance2(l,1)
 
                        p        =   1

                        m        =   0

                    do k=1,endacc

                          if(k.gt.2 .and. total_EN(l,k-1) .gt. total_EN(l,k) .and. total_EN(l,k).lt.total_EN(l,k+1) &
                             & .and. total_EN(l,1).gt.total_EN(l,k)) then

                           firstmini(p)     = total_EN(l,k)

                          dist_firstmini(p) = distance2(l,k)

                          p = p + 1

                          m = k

                          endif
                    enddo 
    
                    dist_maxi(l) = dist_firstmini(1)
                    maxi_g(l)    = firstmini(1)

                    n            =    0

                    do k=m,endacc

                          if(distance2(l,k) .gt. dist_firstmini(1) .and. total_EN(l,k) .gt. maxi_g(l)) then

                                     dist_maxi(l) = distance2(l,k)

                                     maxi_g(l)    = total_EN(l,k)

                                     n = k

                           endif

                    enddo

                    do k=n,endacc

                                      exposi(l) = distance2(l,n)

                                      mini_g(l) = total_EN(l,n)  

                          if(total_EN(l,k).lt.mini_g(l)) then

                                      mini_g(l) = total_EN(l,k)

                                      exposi(l) = distance2(l,k)

                           endif

                     enddo
             endif    

     if(exposi(l) .lt. expconst/2 ) then

         exposi(l) = expconst

         maxi_g(l) = total_en(l,1)

         mini_a(L) = TOTAL_EN(L,ENDACC)

     endif
         
enddo



sumrate=0
Realsumrate=0

do l=1,ende

      rate_g(l) = 0

enddo

do l=1,ende

       DELTA_G(l) = MAXI_G(l) - MINI_A(l)
!       DELTA_H(l) = MINI_G(l) - MINI_A(l)

!     if(delta_g(l).gt.-150) then

       if(delta_g(l).gt.bmin_delta_g .and. delta_g(l).lt.deltamax .and. exposi(l) .gt. expconst/2 ) then
     
             rate_G(l)  = frequency * exp(-DELTA_G(l)*1000/(R*T))
        
  
       else

             rate_g(l) =  0

      endif

             sumrate = sumrate + rate_G(l)  

             Realsumrate = Realsumrate + rate_G(l) 

enddo

! Bubblesortverfahren

         

 a = 1

  
! open(unit = 111,file='TEST_RATE')

incount = 1

if(y==1) then

                          write(5,*) '@    legend on'
                          write(6,*) '@    legend on'
                          write(7,*) '@    legend on'
                          write(8,*) '@    legend on'
                          write(9,*) '@    legend on'
                          write(10,*) '@    legend on'                 
                          write(11,*) '@    legend on'
                          write(12,*) '@    legend on'
                          write(13,*) '@    legend on'
                          write(14,*) '@    legend on'
                          write(15,*) '@    legend on'
                          write(16,*) '@    legend on'
                          write(17,*) '@    legend on'
                          write(18,*) '@    legend on'
                          write(889,*) '@    legend on'
                          write(890,*) '@    legend on'
                          write(891,*) '@    legend on'
                          write(892,*) '@    legend on'
                          write(893,*) '@    legend on'       
                          write(894,*) '@    legend on'
                          write(895,*) '@    legend on'
                          write(896,*) '@    legend on'
                          write(897,*) '@    legend on'
                          write(898,*) '@    legend on'
                                  a1 = 0
                                  a2 = 0
                                  a3 = 0
                                  a4 = 0
                                  a5 = 0
                                  a6 = 0
                                  a7 = 0
                                  a8 = 0
                                  a9 = 0
                                 a10 = 0
                                 a11 = 0
                                 a12 = 0
                                 a13 = 0
                                 a14 = 0
                                 a15 = 0
                                 a16 = 0
                                 a17 = 0
                                 a18 = 0
                                 a19 = 0
                                 a20 = 0
                                 a21 = 0
                                 a22 = 0
                                 a23 = 0
                                 a24 = 0
endif



do l=1,ende
        do k=1,endacc
               
                         if(l==1) then
                           if(k == 1) then

                                     if(a1.lt.10) then
                           write(5,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a1,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                           write(5,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a1,'legend " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a1.ge.10) then                           
                           write(5,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a1,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                           write(5,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a1,'legend " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a1 = a1 + 1
                           endif  
                           write(5,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
                         if(l==2) then
                           if(k == 1) then
                                   
                                    if(a2.lt.10) then
                           write(6,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a2,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(6,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a2,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                                    if(a2.ge.10) then
                           write(6,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a2,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(6,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a2,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif               
                           a2 = a2 + 1
                           endif
                           write(6,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif              
                         if(l==3) then
                           if(k == 1) then
                                     if(a3.lt.10) then
                           write(7,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a3,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(7,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a3,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a3.ge.10) then
                           write(7,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a3,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(7,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a3,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a3 = a3 + 1
                           endif
                           write(7,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif  
                         if(l==4) then
                           if(k == 1) then
                                     if(a4 .lt. 10) then
                           write(8,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a4,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(8,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a4,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a4 .ge. 10) then
                           write(8,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a4,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(8,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a4,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a4 = a4 + 1
                           endif
                           write(8,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif  
                         if(l==5) then
                           if(k == 1) then
                                     if(a5 .lt. 10) then
                           write(9,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a5,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(9,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a5,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a5 .ge. 10) then
                           write(9,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a5,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(9,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a5,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a5 = a5 + 1
                           endif
                           write(9,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==6) then
                           if(k == 1) then
                                     if(a6 .lt. 10) then 
                           write(10,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a6,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(10,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a6,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a6 .ge. 10) then 
                           write(10,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a6,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(10,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a6,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a6 = a6 + 1
                           endif 
                           write(10,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==7) then
                           if(k == 1) then
                                     if(a7 .lt.10) then
                           write(11,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a7,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(11,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a7,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a7 .ge.10) then
                           write(11,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a7,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(11,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a7,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif 
                           a7 = a7 + 1
                           endif
                           write(11,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==8) then
                           if(k == 1) then
                                     if(a8 .lt.10) then
                           write(12,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a8,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(12,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a8,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a8 .ge.10) then
                           write(12,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a8,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(12,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a8,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a8 = a8 + 1
                           endif
                           write(12,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==9) then
                           if(k == 1) then
                                     if(a9 .lt. 10) then
                           write(13,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a9,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(13,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a9,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a9 .ge. 10) then
                           write(13,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a9,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(13,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a9,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a9 = a9 + 1
                           endif
                           write(13,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==10) then
                           if(k == 1) then
                                    if(a10 .lt. 10) then
                           write(14,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a10,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(14,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a10,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                                    if(a10 .ge. 10) then
                           write(14,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a10,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(14,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a10,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                           a10 = a10 + 1
                           endif  
                           write(14,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==11) then
                           if(k == 1) then
                                    if(a11.lt.10) then
                           write(15,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a11,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(15,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a11,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif  
                                    if(a11.ge.10) then
                           write(15,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a11,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(15,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a11,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif                       
                           a11 = a11 + 1
                           endif
                           write(15,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==12) then
                           if(k == 1) then
                                    if(a12.lt.10) then
                           write(16,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a12,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(16,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a12,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                                    if(a12.ge.10) then
                           write(16,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a12,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(16,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a12,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                           a12 = a12 + 1
                           endif
                           write(16,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif  
                         if(l==13) then
                           if(k == 1) then
                                    if(a13.lt.10) then
                           write(17,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a13,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(17,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a13,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                                    if(a13.ge.10) then
                           write(17,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a13,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(17,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a13,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                           a13 = a13 + 1
                           endif
                           write(17,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
                         if(l==14) then
                           if(k == 1) then
                                     if(a14.lt.10)then
                           write(18,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a14,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(18,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a14,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a14.ge.10)then
                           write(18,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a14,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(18,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a14,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a14 = a14 + 1
                           endif
                           write(18,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
                         if(l==15) then
                           if(k == 1) then
                                     if(a15.lt.10) then
                           write(889,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a15,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(889,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a15,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a15.ge.10) then
                           write(889,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a15,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(889,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a15,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a15 = a15 + 1
                           endif 
                           write(889,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==16) then
                           if(k == 1) then
                                     if(a16.lt.10) then
                           write(890,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a16,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(890,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a16,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a16.ge.10) then
                           write(890,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a16,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(890,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a16,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a16 = a16 + 1
                           endif
                           write(890,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==17) then
                           if(k == 1) then
                                     if(a17.lt.10) then
                           write(891,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a17,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(891,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a17,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif        
                                     if(a17.ge.10) then
                           write(891,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a17,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(891,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a17,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif 
                           a17 = a17 + 1
                           endif 
                           write(891,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==18) then
                           if(k == 1) then
                                     if(a18.lt.10) then
                           write(892,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a18,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(892,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a18,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a18.ge.10) then
                           write(892,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a18,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(892,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a18,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif                                     
                           a18 = a18 + 1
                           endif
                           write(892,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif  
                         if(l==19) then
                           if(k == 1) then
                                     if(a19.lt.10) then
                           write(893,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a19,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(893,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a19,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(a19.ge.10) then
                           write(893,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a19,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(893,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a19,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           a19 = a19 + 1
                           endif
                           write(893,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
                         if(l==20) then
                           if(k == 1) then
                                    if(a20.lt.10) then
                           write(894,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a20,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(894,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a20,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                                    if(a20.ge.10) then
                           write(894,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a20,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(894,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a20,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                    endif
                           a20 = a20 + 1
                           endif
                           write(894,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==21) then
                           if(k == 1) then
                                   if(a21.lt.10) then
                           write(895,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a21,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(895,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a21,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                   endif
                                   if(a21.ge.10) then
                           write(895,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a21,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(895,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a21,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                   endif
                           a21 = a21 + 1
                           endif
                           write(895,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif 
                         if(l==22) then
                           if(k == 1) then
                                  if(a22.lt.10) then
                           write(896,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a22,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(896,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a22,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif
                                   if(a22.ge.10) then
                           write(896,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a22,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(896,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a22,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif                                 
                           a22 = a22 + 1
                           endif
                           write(896,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif  
                         if(l==23) then 
                           if(k == 1) then
                                  if(a23.lt.10) then
                           write(897,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a23,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(897,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a23,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif
                                  if(a23.ge.10) then
                           write(897,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a23,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(897,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a23,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif
                          a23 = a23 + 1
                           endif
                           write(897,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
                         if(l==24) then
                           if(k == 1) then
                                  if(a24.lt.10) then
                           write(898,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a24,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(898,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',a24,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif 
                                  if(a24.ge.10) then
                           write(898,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a24,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(898,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',a24,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                  endif 
                           a24 = a24 + 1
                           endif
                           write(898,*) DIST_GRAPH(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
                         endif
         enddo
enddo


do l=1,ende

   if(l==1) then
      rate_G(l-1) = 0
   endif 

   write(111,*) '################'
   write(111,*) sumrate, ' Asumrate'
   write(111,*) rate_G(l)  , ' Arate_G(l)'
   write(111,*) '################' 
   write(111,*) '################'
   write(111,*) incount,'  incount'
   write(111,*) sumrate, ' sumrate'
   write(111,*) rate_G(l)  , ' rate_G(l)'
   write(111,*) MINI_A(l) , 'minimumanfang'
   write(111,*) MINI_G(l) , 'minimum_danach'
   write(111,*) MAXI_G(l) , 'maximum'
   write(111,*) EXPOSI(l) , 'exposi'
   write(111,*) DIST_MAXI(l) ,'dist_maxi'
   write(111,*) Delta_G(l) , 'deltaG'
   write(111,*) l,'l'  
   write(111,*) '###################111'
   write(111,*) hbmap_acc_num_real(l),'acc'
   write(111,*) hbmap_donor_num_real(l),'donor'
   write(111,*) 'ende'
   write(111,*) 'timekonst_delta',timekonst_delta

enddo            


 call system('rm -f table_hbond_B2.ndx table_hbondC2.ndx hbnum.xvg')

 call system('./H_search2b_new.sh')

 inquire(file='table_hbond_B2.ndx',exist=ex2)


 call system('rm -f input_H_search.tpr')
  
 open(unit=39,file='table_hbond_B2.ndx')

do i=1,endtable

    read(39,*,end=554,err=554)

enddo
 
k = 1

    
    do 

434  read(39,*,end = 554,err=554) Bhbmap_donor2 , Bhbmap_hydrogen2 , Bhbmap_acceptor2 
 
       do i=1,size(hbmap_donor)

          if(Bhbmap_donor2    .eq. hbmap_donor(i) .and. &
             Bhbmap_acceptor2 .eq. hbmap_acceptor(i)) then
!              write(*,*) 'case 1*'           
             goto 434
          endif
                          
       enddo

       do i=1,size(hbmap_donor)

          if(Bhbmap_donor2   .eq. hbmap_acceptor(i) .and.&
             Bhbmap_acceptor2 .eq. hbmap_donor(i)) then
!             write(*,*) 'case 2*'
             goto 434

          endif

       enddo

      Bhbmap_donor(k)    = Bhbmap_donor2
      Bhbmap_hydrogen(k) = Bhbmap_hydrogen2
      Bhbmap_acceptor(k) = Bhbmap_acceptor2
      
      
      k=k+1

    enddo

554 continue

 close(unit=39)
 

ende2 = k-1
 

!aa Hinzufügen von eben gebrochener H-Brueckenbindung

 if(k==0) then
      k=1
 endif



do l=1,ende2

   do k=1,endacc2

          TOTAL_EN2(l,k) = 0

   enddo

enddo


 if(ex2 == .false.) then

             ende2 = 0

 endif

!ende2 = 0
      if(h_bond_special.eq.'yes'.or.H_num_spec.eq.'yes') then
                      ende2 = 0
      endif

!do k=1,ende2

!      write(999,*) k , 'k2'
!      write(999,*) Bhbmap_donor(k) , Bhbmap_hydrogen(k) , Bhbmap_acceptor(k)
!      write(999,*) hb_don_A, hb_hyd_A, hb_acc_A
!      write(999,*) '##########' 
      
!enddo
if(donormode==1) then

  do j=1,ende2

               if(bhbmap_donor(j) .le. donorregionstart) then

                  a = bhbmap_donor(j)
                  bhbmap_donor(j) = bhbmap_acceptor(j)
                  bhbmap_acceptor(j) = a

               endif

   enddo

endif

! close(unit=999)

 do i=1,ende2

  open(unit=41,file='minimized3.gro')
  read(41,*,end=9008,err=9008)
  read(41,*,end=9008,err=9008)

    do k = 1,NATOMS
              
              read(41,'(i5,2a5,i5,3f8.3,3f8.4)',end=9008,err=9008) resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
         
             if(atomnum.eq.Bhbmap_donor(i)) then
               bhbmap_donor_resnum(i) = resnum
               bhbmap_donor_restype(i) = restype
               bhbmap_donor_type(i) = atomtype
               bhbmap_donor_num(i)  = atomnum
               bhbmap_donor_x(i)    = coord_x
               bhbmap_donor_y(i)    = coord_y
               bhbmap_donor_z(i)    = coord_z
               Bhbmap_donor_vx(i)      = vel_x
               Bhbmap_donor_vy(i)      = vel_y
               Bhbmap_donor_vz(i)      = vel_z
             endif
             if(atomnum.eq.Bhbmap_hydrogen(i)) then
               Bhbmap_hyd_resnum(i) = resnum
               Bhbmap_hyd_restype(i) = restype
               Bhbmap_hyd_type(i) = atomtype
               Bhbmap_hyd_num(i)  = atomnum
               Bhbmap_hyd_x(i)    = coord_x
               Bhbmap_hyd_y(i)    = coord_y
               Bhbmap_hyd_z(i)    = coord_z
               Bhbmap_hyd_vx(i)      = vel_x
               Bhbmap_hyd_vy(i)      = vel_y
               Bhbmap_hyd_vz(i)      = vel_z
             endif
             if(atomnum.eq.Bhbmap_acceptor(i)) then
               bhbmap_acc_resnum(i) = resnum
               bhbmap_acc_restype(i) = restype
               bhbmap_acc_type(i) = atomtype
               bhbmap_acc_num(i)  = atomnum
               bhbmap_acc_x(i)    = coord_x
               bhbmap_acc_y(i)    = coord_y
               bhbmap_acc_z(i)    = coord_z
               Bhbmap_acc_vx(i)      = vel_x
               Bhbmap_acc_vy(i)      = vel_y
               Bhbmap_acc_vz(i)      = vel_z
             endif
                              
     enddo          

 close(unit=41)

9008 continue
  
 enddo

j = 0

 do i=1,ende2

            if(bhbmap_donor_type(i).ne.'   XX'.and. bhbmap_donor_type(i).ne.'  XXX') then

               j=j+1
           
               bhbmap_donor_resnum(j) = bhbmap_donor_resnum(i)
               bhbmap_donor_restype(j) = bhbmap_donor_restype(i)
               bhbmap_donor_type(j) =  bhbmap_donor_type(i)
               bhbmap_donor_num(j)  =  bhbmap_donor_num(i)
               bhbmap_donor_x(j)    =  bhbmap_donor_x(i)
               bhbmap_donor_y(j)    =  bhbmap_donor_y(i)
               bhbmap_donor_z(j)    =  bhbmap_donor_z(i) 
               Bhbmap_donor_vx(j)   =  Bhbmap_donor_vx(i)
               Bhbmap_donor_vy(j)   =  Bhbmap_donor_vy(i)
               Bhbmap_donor_vz(j)   =  Bhbmap_donor_vz(i)
               Bhbmap_hyd_resnum(j) =  Bhbmap_hyd_resnum(i)             
               Bhbmap_hyd_restype(j) = Bhbmap_hyd_restype(i)
               Bhbmap_hyd_type(j)   =  Bhbmap_hyd_type(i)
               Bhbmap_hyd_num(j)    =  Bhbmap_hyd_num(i)
               Bhbmap_hyd_x(j)      =  Bhbmap_hyd_x(i)
               Bhbmap_hyd_y(j)      =  Bhbmap_hyd_y(i)
               Bhbmap_hyd_z(j)      =  Bhbmap_hyd_z(i)
               Bhbmap_hyd_vx(j)     =  Bhbmap_hyd_vx(i)
               Bhbmap_hyd_vy(j)     =  Bhbmap_hyd_vy(i)
               Bhbmap_hyd_vz(j)     =  Bhbmap_hyd_vz(i)
               bhbmap_acc_resnum(j) =  bhbmap_acc_resnum(i)
               bhbmap_acc_restype(j)=  bhbmap_acc_restype(i)
               bhbmap_acc_type(j)   =  bhbmap_acc_type(i)
               bhbmap_acc_num(j)    =  bhbmap_acc_num(i)
               bhbmap_acc_x(j)      =  bhbmap_acc_x(i)
               bhbmap_acc_y(j)      =  bhbmap_acc_y(i)
               bhbmap_acc_z(j)      =  bhbmap_acc_z(i)
               Bhbmap_acc_vx(j)     =  Bhbmap_acc_vx(i)
               Bhbmap_acc_vy(j)     =  Bhbmap_acc_vy(i)
               Bhbmap_acc_vz(j)     =  Bhbmap_acc_vz(i)

            endif

 enddo

 ende2 = j

 do i=1,ende2

        resnum_min_range(i) = bhbmap_donor_resnum(i) - 1

        resnum_max_range(i) = bhbmap_donor_resnum(i) + 1

 enddo

 do i=1,ende2

    do j=1,count4

       if(bhbmap_donor_resnum(i) .ge. segmentstart(j) .and. bhbmap_donor_resnum(i) .le. segmentend(j)) then

                                          resnum_min_range(i) = segmentstart(j)

!                                          write(*,*) resnum_min_range(i),'resnum_min_range(i)'


                                          resnum_max_range(i) = segmentend(j)

!                                          write(*,*) resnum_max_range(i),'resnum_max_range(i)'

       endif

    enddo

 enddo

 do i=1,ende2

        m_acc   = 1
        m_donor = 1

        open(unit=42,file='minimized3.gro')

        read(42,*,end=9009,err=9009)
        read(42,*,end=9009,err=9009)
    
        do j=1,NATOMS
        
        read(42,'(i5,2a5,i5,3f8.3,3f8.4)',end=9009,err=9009) resnum, restype, atomtype, atomnum, coord_x, &
                                           coord_y, coord_z, vel_x , vel_y , vel_z


!          CG-Element2
!          if (resnum.eq.bhbmap_donor_resnum(i)) then

           if (resnum.ge.resnum_min_range(i) .and. &
             & resnum.le.resnum_max_range(i)) then
          
              Bdonor_resnum(m_donor,i)  = resnum
              Bdonor_restype(m_donor,i) = restype
              Bdonor_atomtype(m_donor,i)= atomtype
              Bdonor_atomnum(m_donor,i) = atomnum
              Bdonor_coord_x(m_donor,i) = coord_x
              Bdonor_coord_y(m_donor,i) = coord_y
              Bdonor_coord_z(m_donor,i) = coord_z

              m_donor= m_donor+1

          endif

        enddo

        close(unit=42) 

9009   continue

        m_don4(i) = m_donor

       open(unit=43,file='minimized3.gro')
        read(43,*,end=9010,err=9010)
        read(43,*,end=9010,err=9010)
    
        do j=1,NATOMS
        
        read(43,'(i5,2a5,i5,3f8.3,3f8.4)',end=9010,err=9010) resnum, restype, atomtype, atomnum, coord_x, &
                                           coord_y, coord_z, vel_x , vel_y , vel_z    
          if (resnum.eq.bhbmap_acc_resnum(i)) then
              Bacc_resnum(m_acc,i)   = resnum
              Bacc_restype(m_acc,i)  = restype
              Bacc_atomnum(m_acc,i)  = atomnum
              Bacc_atomtype(m_acc,i) = atomtype
!              Bacc_coord_x(m_acc,i)  = coord_x
!              Bacc_coord_y(m_acc,i)  = coord_y
!              Bacc_coord_z(m_acc,i)  = coord_z
!              Bacc_vel_x(m_acc,i)    = vel_x
!              Bacc_vel_y(m_acc,i)    = vel_y
!              Bacc_vel_z(m_acc,i)    = vel_z

        m_acc=m_acc+1
        
          endif  
        enddo
       close(unit=43)

9010   continue

       m_acc4(i) = m_acc
       
enddo

 
DIST_scale_x = 0
DIST_scale_y = 0
DIST_scale_z = 0
DIST_scale_A = 0

do l=1,ende2
              
           DIST_scale_x = bhbmap_acc_x(l) - bhbmap_donor_x(l)
           DIST_scale_y = bhbmap_acc_y(l) - bhbmap_donor_y(l)
           DIST_scale_z = bhbmap_acc_z(l) - bhbmap_donor_z(l)

           DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
           scalekey2(l)  = (EXPCONST2-DIST_scale_A)  / (endacc2*DIST_scale_A)
enddo


do l=1,ende2
 

    m_acc2 = m_acc4(l)-1
    m_don2 = m_don4(l)-1
 
          
    scaling = 0 

     do k = 1,endacc2

              if (k == 1) then
                  
               do p=1,m_don4(l)-1
                  
                  new_dir_acc_x2B(p,l) = scaling*(bhbmap_donor_x(l)-bhbmap_acc_x(l))+Bdonor_coord_x(p,l) 
                  new_dir_acc_y2B(p,l) = scaling*(bhbmap_donor_y(l)-bhbmap_acc_y(l))+Bdonor_coord_y(p,l) 
                  new_dir_acc_z2B(p,l) = scaling*(bhbmap_donor_z(l)-bhbmap_acc_z(l))+Bdonor_coord_z(p,l)

               enddo

              open(unit=44,file='minimized3.gro')
                
                       read(44,*,end=9011,err=9011)
                       read(44,*,end=9011,err=9011)

              open(unit=45,file='md_tester2.gro')

                      write(45,*)
                      write(45,*) NATOMS
             
              p = 1

                       do m=1,NATOMS
                                

                       read(44,'(i5,2a5,i5,3f8.3,3f8.4)',end=9011,err=9011) resnum, restype, atomtype, atomnum, coord_x, &
                                         coord_y, coord_z, vel_x , vel_y , vel_z

                       
                       if(resnum.eq.Bdonor_resnum(p,l).and.p.le.m_don4(l)-1) then
                               
                                      

                            coord_x = new_dir_acc_x2B(p,l)
                            coord_y = new_dir_acc_y2B(p,l)
                            coord_z = new_dir_acc_z2B(p,l)
!                            restype = Nameacc
                                       
                            p = p + 1
                                    
                        endif

!                        if(resnum.eq.Bdonor_resnum(p,l)) then

!                            restype = Namedonor
                            
                     
!                            write(*,*) m ,'m'
!                        endif
            
              write(45,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                         coord_y, coord_z, vel_x , vel_y , vel_z
                                                                               
                        enddo

              close(unit=44)

9011          continue
                               

         write(45,*) box_x,box_y,box_z
              
              close(unit=45)
  
      endif


      if(k.gt.1) then   


   
        do p=1,m_don4(l)-1
           
              
                 new_dir_acc_x2B(p,l) = scaling*(bhbmap_donor_x(l)-bhbmap_acc_x(l))+Bdonor_coord_x(p,l) 
                 new_dir_acc_y2B(p,l) = scaling*(bhbmap_donor_y(l)-bhbmap_acc_y(l))+Bdonor_coord_y(p,l) 
                 new_dir_acc_z2B(p,l) = scaling*(bhbmap_donor_z(l)-bhbmap_acc_z(l))+Bdonor_coord_z(p,l)
        
        enddo 
        
!        DIST_TEST(l,k) = sqrt(new_dir_acc_z2B(2,l)**2+new_dir_acc_y2B(2,l)**2+new_dir_acc_z2B(2,l)**2)
             
        open(unit=46,file='minimized2.gro')

        read(46,*,end=9012,err=9012)
        read(46,*,end=9012,err=9012)
        
        open(unit=47,file='md_tester2.gro')
        
     
        write(47,*)
        write(47,*) NATOMS

        p = 1
       
            do m=1,NATOMS
 

               read(46,'(i5,2a5,i5,3f8.3,3f8.4)',end=9012,err=9012) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

               if(resnum.eq.Bdonor_resnum(p,l) .and. p.le.m_don4(l)-1) then

                   coord_x = new_dir_acc_x2B(p,l)
                   coord_y = new_dir_acc_y2B(p,l)
                   coord_z = new_dir_acc_z2B(p,l)   

                   p = p + 1

               write(*,*) m ,'m'

               endif

               write(47,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(47,*) box_x,box_y,box_z

               close(unit=47)
               close(unit=46)

9012   continue

       endif

      call system ('rm -f indexsim2.ndx')
       
      
       open(unit=48,file='temp.ndx')
       
            write(48,'(a7)') '[ ACC ]'
        
!           do p=1,m_acc4(l)-1

           write(48,*)  bhbmap_acc_num(l)
!           write(48,*) bacc_atomnum(p,l)

!           enddo
      
           write(48,'(a7)') '[ DON ]'

          if(keepresmode2 == 0) then

            write(48,*) bhbmap_donor_num(l)             
   
          else        
      
          do p=1,m_don4(l) -1
                    
           write(48,*)  bdonor_atomnum(p,l)
           
           enddo

          endif
    
       close(unit=48)
      
      
       call system('cat temp.ndx index.ndx > indexsim2.ndx')
       call system('rm -f temp.ndx')
      
       call system ('rm -f minimized2.gro traj.trr mdout.mdp md.log ener.edr')
       

       if(k==1) then
            
            call system ('cp minim3C.mdp minimtest2.mdp') 
    
            open(unit=49,file='minimtest2.mdp')

                 do j = 1,minimline
                      
                          read(49,*)
                 enddo

             write(49,*) 'freezegrps        =' , Nameacc , Namedonor
             write(49,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=49)

       endif

                write(*,*) 'everything is ok'
         
         call system ('rm -f input2.tpr minimized2.gro traj.trr mdout.mdp md.log ener.edr')
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest2.mdp -c md_tester2.gro -p LOV2.top -o input2.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input2.tpr -c minimized_A.gro ')
         call system ('rm -f md.log traj.trr mdout.mdp input2.tpr ener.edr out')
! steep --> lbfgs
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_C.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized2.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro out')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized2.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized2.gro ')
         call system ('rm -f ./#*# out')

!         if(modulo(k,modnumber)==0) then

!         call system('/loctmp/pee18323/GRO_bin/bin/grompp -f minim3C.mdp -c minimized2.gro -p LOV2.top -o input.tpr')
!         call system('/loctmp/pee18323/GRO_bin/bin/mdrun  -nice 0 -v -s input.tpr -c minimized2.gro')
!         call system('rm -f \#*# input.tpr md.log traj.trr ener.edr mdout.mdp')

!         endif

         AC_X2     = 0
         AC_Y2     = 0
         AC_Z2     = 0
         DO_X2     = 0
         DO_Y2     = 0
         DO_Z2     = 0
         DIST_X2   = 0
         DIST_Y2   = 0
         DIST_Z2   = 0
         DIST_X3   = 0
         DIST_Y3   = 0
         DIST_Z3   = 0

        open(unit=50,file='minimized2.gro')

         read(50,*,end=9013,err=9013)
         read(50,*,end=9013,err=9013)


        do b=1,NATOMS

           read(50,'(i5,2a5,i5,3f8.3,3f8.4)',end=9013,err=9013) resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
           if(atomnum.eq.bhbmap_acc_num(l)) then
                 AC_X2 = coord_x 
                 AC_Y2 = coord_y
                 AC_Z2 = coord_z
!                 AC_XA(l,k) = coord_x
            endif       
           if(atomnum.eq.bhbmap_donor_num(l)) then
                 DO_X2 = coord_x
                 DO_Y2 = coord_y
                 DO_Z2 = coord_z
            endif
           DIST_X2 = AC_X2 - DO_X2 -  (bhbmap_acc_x(l)-bhbmap_donor_x(l))
           DIST_Y2 = AC_Y2 - DO_Y2 -  (bhbmap_acc_y(l)-bhbmap_donor_y(l))
           DIST_Z2 = AC_Z2 - DO_Z2 -  (bhbmap_acc_z(l)-bhbmap_donor_z(l))

           DISTANCE2B(l,k) = sqrt(DIST_X2**2+DIST_Y2**2+DIST_Z2**2)

           DIST_X3 = AC_X2 - DO_X2
           DIST_Y3 = AC_Y2 - DO_Y2
           DIST_Z3 = AC_Z2 - DO_Z2

           DIST_GRAPH3(l,k) = sqrt(DIST_X3**2+DIST_Y3**2+DIST_Z3**2)
           
            
        enddo 
      close(unit=50)

9013  continue
     

      write(*,*)  'everything is ok B',k,'ende2',l
 
      scaling = scaling + scalekey2(l)

                             call system ('rm -f md_tester2.gro')
                      
                       call system ('rm -f md.log input2.tpr ener.edr traj.trr mdout.mdp')

                               dist_acc_x = 0
                               dist_acc_y = 0
                               dist_acc_z = 0
                             range_betrag = 0 
                                     
                                     mam =  1


        open(unit = 51,file='minimized2.gro')
        read(51,*,end=9014,err=9014)
        read(51,*,end=9014,err=9014)
          
        do i=1,NATOMS

          read(51,'(i5,2a5,i5,3f8.3,3f8.4)',end=9014,err=9014) resnum, restype, atomtype, atomnum, coord_x, &
          coord_y, coord_z, vel_x , vel_y , vel_z

          if(resnum .ne. bhbmap_donor_resnum(l) ) then

                dist_acc_x = coord_x - DO_X2
                dist_acc_y = coord_y - DO_Y2    
                dist_acc_z = coord_z - DO_Z2

                range_betrag = sqrt(dist_acc_x**2 + dist_acc_y**2 + dist_acc_z**2)

               if(range_betrag.lt.cutoff .and. range_betrag .ne. 0.0) then
                            Cut_range2(mam)       = range_betrag
                            rangetype2(mam)       = restype
                            RANGE_ZATOM2(mam)   = atomtype
!                            write(26,*) rangetype(mum)
           
                            mam = mam+1          
                             
               endif
         endif
         
       enddo

       close ( unit=51 )

9014   continue

matter2 = mam - 1
n = 1

       do i=1,matter2

           open(unit=52,file='ffG43a1B.rtp')

            do  j=1,rtplines
   
                read(52,'(2X,a4)') RESSEARCH
                
                          if(RESSEARCH.eq.rangetype2(i)) then
                                read(52,*)
                                do  
                                  read(52,'(a5,1X,a5,4X,f8.5)',end=5690,err=5690) nameatom , typeatom,&
                                     & charge2

                                  if(nameatom.eq.RANGE_ZATOM2(i)) then
                                        charge_atomB(i) = charge2

                                        radius_atom(i)  = typeatom
     
                                        n = n + 1 
                                        exit
                                   endif
                                enddo
5690                              exit
                           endif
             enddo

            close(unit=52)
       enddo    

            open(unit=53,file='ffG43a1B.rtp')

             do j=1,rtplines

                read(53,'(2X,a4)') RESSEARCH                

                if(RESSEARCH.eq.bhbmap_donor_restype(l)) then
                                read(53,*)
                                      do 
                                          read(53,'(a5,1X,a5,4X,f8.5)',end=5620,err=5620) nameatom , typeatom,&
                                          & charge2
                                          if(nameatom.eq.bhbmap_donor_type(l)) then
                                             hb_acc_chargeB(l)  = charge2

                                             radius_donor(l)    = typeatom

                                             exit
                                          endif 
                                       enddo

5620                  exit
                endif
                
            enddo
           
          close(unit=53)
                            
!AAA


angle_en = 0


if(angleboolean2 == 1) then

            open(unit=28,file='angles')

             do j=1,anglelines

                read(28,'(2X,a4)') RESSEARCH                

                if(RESSEARCH.eq.bhbmap_donor_restype(l)) then
                              

                                  do 
                                          read(28,'(3X,a6)',end=16451,err=16451) ANGLESEARCH
                                          
                                          if(ANGLESEARCH.eq.'angles') then

                                              n = 1
                                            
                                             do 

                                               read(28,'(a5,1X,a5,1X,a5,5X,f6.2,6X,f6.2)',end=16451,err=16451) atom_1 , atom_2 , atom_3 , angle , forceconstant                                       

                                               if(atom_1 == bhbmap_donor_type(l) .or. &
                                                  atom_2 == bhbmap_donor_type(l) .or. &
                                                  atom_3 == bhbmap_donor_type(l) ) then



                                                            open(unit =26,file='minimized2.gro')

                                                            read(26,*)
                                                            read(26,*)
          
                                                            do i=1,NATOMS

                                                               read(26,'(i5,2a5,i5,3f8.3,3f8.4)',end=14106,err=14106) resnum, restype, atomtype, atomnum, coord_x, &
                                                                                                                      coord_y, coord_z, vel_x , vel_y , vel_z
  
                                                               if(resnum == bhbmap_donor_resnum(l)) then                                                                
                                     
                                                                             
                                                               if(atom_1 == atomtype) then

                                                                            atom_1_x = coord_x               
                                                                            atom_1_y = coord_y
                                                                            atom_1_z = coord_z

                                                               endif

                                                               if(atom_2 == atomtype) then

                                                                            atom_2_x = coord_x               
                                                                            atom_2_y = coord_y
                                                                            atom_2_z = coord_z

                                                               endif

                                                               if(atom_3 == atomtype) then

                                                                            atom_3_x = coord_x               
                                                                            atom_3_y = coord_y
                                                                            atom_3_z = coord_z

                                                               endif

                                                               endif

                                                            enddo

14106                                                          continue

                                                               close(unit=26)

                                                               diff_x1 = atom_1_x - atom_2_x
                                                               diff_y1 = atom_1_y - atom_2_y
                                                               diff_z1 = atom_1_z - atom_2_z

                                                               diff_x2 = atom_3_x - atom_2_x                        
                                                               diff_y2 = atom_3_y - atom_2_y
                                                               diff_z2 = atom_3_z - atom_2_z

                                                             diff_tot  = sqrt((diff_x1)**2+(diff_y1)**2+(diff_z1)**2)*sqrt((diff_x2)**2+(diff_y2)**2+(diff_z2)**2)

                                                            angle_real = 180/PI*acos((diff_x1*diff_x2 + diff_y1*diff_y2 + diff_z1*diff_z2)/diff_tot)

                                                            write(*,*) angle_real,'angle'

                                                               if(k==1 .and. anglefirstbool ==1) then

                                                                        anglefirst1(n) = angle_real


                                                               endif

                                                               if(anglefirstbool == 1) then

                                                                       angle = anglefirst1(n)

                                                               endif 


                                                            angle_en = angle_en + 0.5*forceconstant*(cos(angle_real*PI/180)-cos(angle*PI/180))**2


                                                            n = n + 1
                                             
                                          endif 
                                       enddo

                                 endif
                            enddo
16451                  exit
                endif
                
            enddo
           
          close(unit=28)
          close(unit=26)

if(k==2) then


i = 1

endif

            open(unit=28,file='angles')

             do j=1,anglelines

                read(28,'(2X,a4)') RESSEARCH                


                if(RESSEARCH.eq.bhbmap_donor_restype(l)) then

                                read(28,*)

                                      do 
                                          read(28,'(3X,a8)',end=2424,err=2424) DIHANGLESEARCH
                                          
                                          if(DIHANGLESEARCH.eq.'dihangle') then
                            
                                                                              exxbool = .false.
                                         
                                             do 

                                               read(28,'(a5,1X,a5,1X,a5,1X,a5,5X,f7.3,7X,f4.1,4X,i1)',end=2424,err=2424) atom_1 , atom_2 , atom_3 , atom_4, angle , forceconstant , o                                       


                                               if(angle==0 .and. forceconstant == 0 .and. o==0 ) then

                                                                                exxbool = .true.
                                                                                exit
                                               endif

                                               if(atom_1 == bhbmap_donor_type(l) .or. &
                                                  atom_2 == bhbmap_donor_type(l) .or. &
                                                  atom_3 == bhbmap_donor_type(l) .or. & 
                                                  atom_4 == bhbmap_donor_type(l)) then


                                                            open(unit =26,file='minimized2.gro')

                                                            read(26,*)
                                                            read(26,*)
          
                                                            do i=1,NATOMS

                                                               read(26,'(i5,2a5,i5,3f8.3,3f8.4)',end=14006,err=14006) resnum, restype, atomtype, atomnum, coord_x, &
                                                                                                                      coord_y, coord_z, vel_x , vel_y , vel_z
                                                                                                     
                                                               if(resnum == bhbmap_donor_resnum(l)) then                                                                 
            
                                                               if(atom_1 == atomtype) then

                                                                            atom_1_x = coord_x               
                                                                            atom_1_y = coord_y
                                                                            atom_1_z = coord_z


                                                               endif

                                                               if(atom_2 == atomtype) then

                                                                            atom_2_x = coord_x               
                                                                            atom_2_y = coord_y
                                                                            atom_2_z = coord_z


                                                               endif

                                                               if(atom_3 == atomtype) then

                                                                            atom_3_x = coord_x               
                                                                            atom_3_y = coord_y
                                                                            atom_3_z = coord_z




                                                               endif

                                                               if(atom_4 == atomtype) then

                                                                            atom_4_x = coord_x               
                                                                            atom_4_y = coord_y
                                                                            atom_4_z = coord_z



                                                               endif

                                                               endif

                                                            enddo

                                                               14006 continue

                                                               close(unit=26)
                                                               
                                        !                      Bildung der Normalen

                                                               diff_x1 = atom_1_x - atom_2_x
                                                               diff_y1 = atom_1_y - atom_2_y
                                                               diff_z1 = atom_1_z - atom_2_z

                                                               diff_x2 = atom_2_x - atom_3_x                        
                                                               diff_y2 = atom_2_y - atom_3_y
                                                               diff_z2 = atom_2_z - atom_3_z

                                                               kreuz_x1 = diff_y1*diff_z2 - diff_z1*diff_y2
                                                               kreuz_y1 = diff_z1*diff_x2 - diff_x1*diff_z2
                                                               kreuz_z1 = diff_x1*diff_y2 - diff_y1*diff_x2

                                                               diff_x3 = atom_2_x - atom_3_x
                                                               diff_y3 = atom_2_y - atom_3_y
                                                               diff_z3 = atom_2_z - atom_3_z

                                                               diff_x4 = atom_3_x - atom_4_x                        
                                                               diff_y4 = atom_3_y - atom_4_y
                                                               diff_z4 = atom_3_z - atom_4_z

                                                               kreuz_x2 = diff_y3*diff_z4 - diff_z3*diff_y4
                                                               kreuz_y2 = diff_z3*diff_x4 - diff_x3*diff_z4
                                                               kreuz_z2 = diff_x3*diff_y4 - diff_y3*diff_x4

                                                                 
                                                              diff_tot  = sqrt((kreuz_x1)**2+(kreuz_y1)**2+(kreuz_z1)**2)*sqrt((kreuz_x2)**2+(kreuz_y2)**2+(kreuz_z2)**2) 

                                                              angle_dih = 180/PI*acos((kreuz_x1*kreuz_x2 + kreuz_y1*kreuz_y2 + kreuz_z1*kreuz_z2)/diff_tot) 

                                                               if(k==1 .and. anglefirstbool ==1) then

                                                                        anglefirst2(i) = angle_dih

                                                               endif

                                                               if(anglefirstbool == 1) then

                                                                       angle = anglefirst2(i)

                                                               endif 


                                                              angle_en = angle_en + forceconstant*(1 + cos((o*angle_dih - angle)*PI/180))

                                                              i = i+1
                                                endif 
                                       
                                     enddo

                               endif

                                    if(exxbool == .true.) then

                                                          exit

                                    endif

                       enddo

2424             continue

                 exit
                endif
                
            enddo
           
          close(unit=28)
          close(unit=26)
endif

                      open(unit=98989,file='ffG43a1nb.itp')

                      read(98989,*)
                      read(98989,*)

                      do p=1,itplines

                             read(98989,'(a5,32x,g13.10,g14.11)') tester3 , LJ_param_6, LJ_param_12

                             if(tester3 == radius_donor(l)) then

                                   LJ_acc_6b(l) = LJ_param_6

                                   LJ_acc_12b(l)= LJ_param_12

                                   write(*,*) LJ_acc_6b(l),'6acc'
                                   write(*,*) LJ_acc_12b(l),'12acc'
                                   write(*,*) tester3,'t3',hbmap_donor_type_real(l),'rangeacc'
 
                                   n = n+1

                             endif

                      enddo

                      close(unit=98989)
 

    
        do i=1,matter2       


                      open(unit=98989,file='ffG43a1nb.itp')

                      read(98989,*)
                      read(98989,*)

                             exxbool = .false.

                      do p=1,itplines

                             read(98989,'(a5,32x,g13.10,g14.11)') tester3 , LJ_param_6, LJ_param_12

                             if(tester3 == radius_atom(i)) then

                                   LJ_range_6b(i) = LJ_param_6

                                   LJ_range_12b(i)= LJ_param_12

                            !       write(*,*) LJ_range_6(i),'6'
                            !       write(*,*) LJ_range_12(i),'12'
                            !       write(*,*) tester3,'t3',range_zatom(i),'range'

                                   n = n+1

                                   exxbool = .true.


                             endif

                             if(p==itplines.and.exxbool==.false.) then

                                   write(*,*) RANGE_ZATOM2(i),'not taken'


                             endif

                      enddo

                      close(unit=98989)

           enddo
     



TOTALCHARGE = 0
   LJ     = 0 


           do i=1,matter2

! Cut_range von 0.25 auf 0.3 nm

  
                          factor_A = -((alpha+4)*cutoff - (alpha+1)*0)/(cutoff**(alpha+2)*(cutoff)**2)

                          factor_B =  ((alpha+3)*cutoff - (alpha+1)*0)/(cutoff**(alpha+2)*(cutoff)**3)

                          factor_C =  1/cutoff**alpha-factor_A/3*(cutoff)**3-factor_B/4*(cutoff)**4

                    
                     if(Cut_range2(i).gt.0.15) then

                        Coulomb2(i) = COUL_CONST*hb_acc_chargeB(l)*charge_atomB(i)*(1/Cut_range2(i)**(alpha) - factor_A/3*(Cut_range2(i))**3 &
                                       & - factor_B/4*(Cut_range2(i))**4 -factor_C) 
                     
                        TOTALCHARGE = TOTALCHARGE + Coulomb2(i)

                      endif
 
                    
                      if(Cut_range2(i).gt.0.15)then

                       LJ_6(i)  = sqrt(abs(LJ_range_6b(i)*LJ_acc_6b(l)))

                       LJ_12(i) = sqrt(abs(LJ_range_12b(i)*LJ_acc_12b(l)))

                       
                       LJ_pot(i) = LJ_12(i)/Cut_range2(i)**12 - LJ_6(i)/Cut_range2(i)**6
                      
                       if(Cut_range2(i).gt.cutoff) then

                          LJ_pot(i) = 0

                       endif
                       
                       LJ = LJ + LJ_pot(i)

                      endif
                       
         enddo     

                       
                  
                       TOTAL_EN2(l,k) = TOTAL_EN2(l,k) + TOTALCHARGE + LJ + angle_en + FORCE*k/endacc2 

                                  
  enddo

           call system ('rm -f minimtest2.mdp')
enddo


            
                                               open(unit=54,file='total_Energy1_n')
                                               open(unit=55,file='total_Energy2_n')
                                               open(unit=56,file='total_Energy3_n')
                                               open(unit=57,file='total_Energy4_n')
                                               open(unit=58,file='total_Energy5_n')
                                               open(unit=59,file='total_Energy6_n')
                                               open(unit=60,file='total_Energy7_n')
                                               open(unit=61,file='total_Energy8_n')
                                               open(unit=62,file='total_Energy9_n')
                                               open(unit=63,file='total_Energy10_n')
                                               open(unit=64,file='total_Energy11_n')
                                               open(unit=65,file='total_Energy12_n')
                                               open(unit=66,file='total_Energy13_n')
                                               open(unit=67,file='total_Energy14_n')
!                                               open(unit=69,file='total_check_27.6')
            
 

do l=1,ende2

             bexposi(l)    = distance2b(l,endacc2)

             bmini_g(l)    = total_EN2(l,1)

             bmini_a(l)    = total_EN2(l,1)
             
             bmaxi_g(l)    = total_EN2(l,endacc2)
             
             bdist_maxi(l) = distance2b(l,endacc2)

             do k=1,endacc2
 
                    if(total_EN2(l,endacc2-k+1).gt.bmaxi_g(l).and.total_EN2(l,endacc2-k+1).gt.bmini_a(l)) then

                       bmaxi_g(l)    = total_EN2(l,endacc2-k+1)
                       bdist_maxi(l) = distance2b(l,endacc2-k+1)

                    endif  

             enddo


enddo

!Verringerung des Schleifenparameters ende2 um 1 und damit Wegnahme des unerwünschten Potentialverlaufs
!Mit Bexposi = 0, d.h. kein Maximum
! falls glatte Kurve dann globales Minimum
! close(unit=70) 

bsumrate = 0
brealsumrate=0

do l=1,ende2

       brate_g(l) = 0

enddo

do l=1,ende2

    
       BDELTA_G(l) = BMAXI_G(l) - BMINI_A(l)

!       BDELTA_H(l) = total_En2(l,1) - total_En2(l,endacc2)

!       if(Bdelta_H(l) == bdelta_g(l)) then

!                         bdelta_h(l) =  0

!       endif

       if( bexposi(l).gt.0.05 .and. bdelta_g(l)-bdelta_H(l).gt.bmin_delta_g .and.bdelta_g(l)-bdelta_H(l).lt.deltamax ) then
     
           Brate_G(l)  = frequency  * exp(-BDELTA_G(l)*1000/(R*T))

       else
 
           Brate_G(l) = 0

       endif
!       Brate_H(l)  = Brate_G(l) * exp(BDELTA_H(l)*1000/(R*T))

        bsumrate = bsumrate + Brate_G(l)  
  
        brealsumrate=brealsumrate + Brate_G(l)

!      endif
      
enddo

 csumrate = 0

       csumrate = bsumrate + sumrate
! Bubblesort2
 

! Auswahl des events, Bildung der Summe der Raten

selectrate = 0
k = 2

selectionrate(1) = 0

do i=1,ende

       selectionrate(k) = selectionrate(k-1) + rate_g(i)

       selectivenumber(k) = i

       k = k + 1 

enddo


do i=1,ende2

       selectionrate(k) = selectionrate(k-1) + brate_g(i)

       selectivenumber2(k) = i

       k = k+1

enddo

call random_number(harvest3)

do i = 1,ende + ende2+1

    if(selectionrate(i-1) .lt. csumrate*harvest3 .and. csumrate*harvest3 .le. selectionrate(i)) then 

                           if(i.le.ende+1) then

                                selectrate = 1

                                selnumber     = selectivenumber(i)

                           endif

                           if(i.gt.ende+1) then

                                selectrate = 2

                                selnumber     = selectivenumber2(i)

                           endif

    endif

enddo

        write(111,*) '#######rates#######'

do i = 1,ende + ende2+1


        write(111,*) 'i'
        write(111,*) selectionrate(i),'selectionrate'

enddo

        write(111,*) selectrate,'selectrate'
        write(111,*) selnumber, 'selnumber'
        write(111,*) csumrate*harvest3,'csumrate*harvest3'
        write(111,*) csumrate,'csumrate'
        write(111,*) ende,'ende'
        write(111,*) ende2,'ende2'
        write(111,*) '######ende rates########'

write(*,*) 'sel1'


write(*,*) 'sel2'    

if(y==1) then

                          write(54,*) '@    legend on'
                          write(55,*) '@    legend on'
                          write(56,*) '@    legend on'
                          write(57,*) '@    legend on'
                          write(58,*) '@    legend on'
                          write(59,*) '@    legend on'                 
                          write(60,*) '@    legend on'
                          write(61,*) '@    legend on'
                          write(62,*) '@    legend on'
                          write(63,*) '@    legend on'
                          write(64,*) '@    legend on'
                          write(65,*) '@    legend on'
                          write(66,*) '@    legend on'
                          write(67,*) '@    legend on'
                          b1  = 0
                          b2  = 0
                          b3  = 0
                          b4  = 0
                          b5  = 0
                          b6  = 0
                          b7  = 0
                          b8  = 0
                          b9  = 0
                          b10 = 0
                          b11 = 0
                          b12 = 0
                          b13 = 0
                          b14 = 0 
endif



do l=1,ende2
        do k=1,endacc2
               
                         if(l==1) then
                           if(k == 1) then
                                        if(b1.lt.10) then
                           write(54,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b1,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                           write(54,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b1,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif     
                                        if(b1.ge.10) then
                           write(54,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b1,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                           write(54,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b1,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif 
                           b1 = b1 + 1  
                         endif  
                           write(54,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif
                         if(l==2) then
                           if(k == 1) then
                                        if(b2.lt.10) then
                           write(55,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b2,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(55,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b2,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                                        if(b2.ge.10) then
                           write(55,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b2,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(55,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b2,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                           b2 = b2 + 1  
                         endif
                           write(55,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif              
                         if(l==3) then
                           if(k == 1) then
                                        if(b3.lt.10) then
                           write(56,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b3,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(56,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b3,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                                        if(b3.ge.10) then
                           write(56,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b3,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(56,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b3,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                           b3 = b3 + 1
                           endif
                           write(56,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif  
                         if(l==4) then
                           if(k == 1) then
                                        if(b4.lt.10) then
                           write(57,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b4,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(57,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b4,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                                        if(b4.ge.10) then
                           write(57,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b4,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(57,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b4,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                        endif
                           b4 = b4 + 1
                           endif
                           write(57,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif  
                         if(l==5) then
                           if(k == 1) then
                                       if(b5.lt.10) then
                           write(58,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b5,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(58,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b5,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                       endif
                                       if(b5.ge.10) then
                           write(58,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b5,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(58,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b5,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                       endif
                           b5 = b5 + 1
                           endif
                           write(58,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==6) then
                           if(k == 1) then
                                      if(b6.lt.10) then
                           write(59,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b6,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(59,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b6,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                                      if(b6.ge.10) then
                           write(59,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b6,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(59,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b6,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                           b6 = b6 + 1
                           endif 
                           write(59,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==7) then
                           if(k == 1) then
                                      if(b7.lt.10) then
                           write(60,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b7,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(60,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b7,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                                      if(b7.ge.10) then
                           write(60,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b7,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(60,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b7,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                           b7 = b7 + 1
                           endif
                           write(60,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==8) then
                           if(k == 1) then
                                      if(b8 .lt.10) then
                           write(61,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b8,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(61,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b8,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                                      if(b8 .ge.10) then
                           write(61,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b8,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(61,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b8,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif

                           b8 = b8 + 1
                           endif
                           write(61,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==9) then
                           if(k == 1) then
                                      if(b9.lt.10) then
                           write(62,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b9,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(62,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b9,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                                      if(b9.ge.10) then
                           write(62,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b9,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(62,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b9,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                      endif
                           b9 = b9 + 1
                           endif
                           write(62,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==10) then
                           if(k == 1) then
                                     if(b10.lt.10) then
                           write(63,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b10,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(63,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b10,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(b10.ge.10) then
                           write(63,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b10,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(63,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b10,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           b10 = b10 + 1
                           endif  
                           write(63,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==11) then
                           if(k == 1) then
                                     if(b11.lt.10) then
                           write(64,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b11,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(64,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b11,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(b11.ge.10) then
                           write(64,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b11,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(64,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b11,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           b11 = b11 + 1
                           endif
                           write(64,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif 
                         if(l==12) then
                           if(k == 1) then
                                     if(b12.lt.10) then
                           write(65,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b12,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(65,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b12,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(b12.ge.10) then
                           write(65,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b12,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(65,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b12,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           b12 = b12 + 1
                           endif
                           write(65,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif  
                         if(l==13) then
                           if(k == 1) then
                                     if(b13.lt.10) then
                           write(66,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b13,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(66,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b13,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(b13.ge.10) then
                           write(66,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b13,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(66,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b13,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           b13 = b13 + 1
                           endif
                           write(66,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif
                         if(l==14) then
                           if(k == 1) then
                                     if(b14.lt.10) then
                           write(67,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b14,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(67,'(a6,i1,a15,i2,a5,i4,a1,i4,a1)') '@    s',b14,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                                     if(b14.ge.10) then
                           write(67,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b14,'comment " cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"' 
                           write(67,'(a6,i2,a15,i2,a5,i4,a1,i4,a1)') '@    s',b14,'legend" cycle',y,'type ',hbmap_acc_num_real(l),'-',hbmap_donor_num_real(l),'"'
                                     endif
                           b14 = b14 + 1
                           endif
                           write(67,*) DIST_GRAPH3(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
                         endif
                       
         enddo
enddo

do l=1,ende2

   if(l==1) then
      Brate_G(l-1) = 0
   endif 

incount = 1  

   write(111,*) Brate_G(l)  , ' Brate_G(l)'
   write(111,*) '######2##########nonevent2'      
   write(111,*) Brate_G(l-1), ' rate_G(l-1)'
   write(111,*) csumrate, ' csumrate'
   write(111,*) Brate_G(l)  , ' Brate_G(l)'
   write(111,*) Bdelta_h(l), 'bdelta_h(l)'
   write(111,*) BMINI_A(l) , 'minimumanfang'
   write(111,*) BMINI_G(l) , 'minimum_danach'
   write(111,*) BMAXI_G(l) , 'maximum'
   write(111,*) BEXPOSI(l) , 'exposi'
   write(111,*) BDIST_MAXI(l) ,'dist_maxi'
   write(111,*) BDelta_G(l)   ,'!bdeltaG!'
   write(111,*) bhbmap_donor_num(l),'!donornum!'
   write(111,*) bhbmap_acc_num(l),'accnum'
   write(111,*) '#################nonevent2'

enddo

      next = .true.

 if(y==1) then

    open(unit=125,file="xvg_rate_n")

 endif


  if(selectrate == 2) then
   next = .false.

   l = selnumber

   write(111,*) '#######2event########'
   write(111,*) incount,'  incount,form'
   write(111,*) csumrate, ' csumrate'
   write(111,*) Brate_G(l)  , ' Brate_G(l)'
   write(111,*) BMINI_A(l) , 'minimumanfang'
   write(111,*) Bdelta_g(l), 'bdelta_g(l)'
   write(111,*) bdelta_h(l), 'bdelta_h(l)'
   write(111,*) BMINI_G(l) , 'minimum_danach'
   write(111,*) BMAXI_G(l) , 'maximum'
   write(111,*) BEXPOSI(l) , 'exposi'
   write(111,*) BDIST_MAXI(l) ,'dist_maxi'
   write(111,*) l, 'l'
   write(111,*) bhbmap_donor_num(l),'bhbmap_donor_num(l)'
   write(111,*) bhbmap_hyd_num(l)  ,'bhbmap_hyd_num(l)'
   write(111,*) bhbmap_acc_num(l)  ,'bhbmap_acc_num(l)'
   write(111,*) '######2##########'

   incount  = incount + 1

DIST_scale_x = 0
DIST_scale_y = 0
DIST_scale_z = 0
DIST_scale_A = 0


                 write(125,*)         

         do k=1,endacc2
                 write(125,*) DISTANCE2B(l,k)  , TOTAL_EN2(l,k)-TOTAL_EN2(l,1)
         enddo
  
                 write(125,*)
  
!    bEXPOSI(l) = bEXPOSI(l) - damping
   
         
scalekey2c = 0
scaling2c  = 0


         do k=1,endacc2

         if(k==1) then

               DIST_scale_x = bhbmap_acc_x(l) - bhbmap_donor_x(l)
               DIST_scale_y = bhbmap_acc_y(l) - bhbmap_donor_y(l)
               DIST_scale_z = bhbmap_acc_z(l) - bhbmap_donor_z(l)

               DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
               scalekey2c   = - bEXPOSI(l) / (endacc2*DIST_scale_A)
               
      
              do p=1,m_don4(l)-1
                                        
                 bnew_dir_acc_xA(p,l) = scaling2c*(bhbmap_donor_x(l)-bhbmap_acc_x(l))+bdonor_coord_x(p,l) 
                 bnew_dir_acc_yA(p,l) = scaling2c*(bhbmap_donor_y(l)-bhbmap_acc_y(l))+bdonor_coord_y(p,l) 
                 bnew_dir_acc_zA(p,l) = scaling2c*(bhbmap_donor_z(l)-bhbmap_acc_z(l))+bdonor_coord_z(p,l)
        
              enddo 
             
        open(unit=71,file='minimized3.gro')

        read(71,*,end=9015,err=9015)
        read(71,*,end=9015,err=9015)
        
        open(unit=72,file='md_tester.gro')
     
        write(72,*)
        write(72,*) NATOMS

        p = 1
       
            do m=1,NATOMS
 

               read(71,'(i5,2a5,i5,3f8.3,3f8.4)',end=9015,err=9015) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

               if(resnum.eq.bdonor_resnum(p,l) .and. p.le.m_don4(l)-1) then

                   coord_x = bnew_dir_acc_xA(p,l)
                   coord_y = bnew_dir_acc_yA(p,l)
                   coord_z = bnew_dir_acc_zA(p,l)

                   p = p + 1

               write(*,*) m ,'m'

               endif

               write(72,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(72,*) box_x,box_y,box_z

               close(unit=72)
               close(unit=71)

9015       continue
               
           open(unit=73,file='temp.ndx')
       
           write(73,'(a7)') '[ ACC ]'
        
           
!           do  p=1,m_acc4(l)-1

           write(73,*)  bhbmap_acc_num(l)
!            write(73,*) bacc_atomnum(p,l)           

!           enddo

      
           write(73,'(a7)') '[ DON ]'

           if(keepresmode2==0) then

             write(73,*) bhbmap_donor_num(l)
  
           else

    
             do p=1,m_don4(l)-1     
        
                write(73,*)  bdonor_atomnum(p,l)
!           write(73,*) bhbmap_donor_num(l)
              enddo
    
           endif

           close(unit=73)

      
      
          call system('cat temp.ndx index.ndx > indexsim2.ndx')
          call system('rm -f temp.ndx')

          call system ('cp -f minim3C.mdp minimtest.mdp') 
    
          open(unit=74,file='minimtest.mdp')

          do j = 1,minimline
                      
             read(74,*)
          
          enddo
   
             write(74,*) 'freezegrps        =' , Nameacc , Namedonor
             write(74,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=74)             
                            
         call system ('rm -f input.tpr minimized.gro ener.edr md.log mdout.mdp traj.trr')
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest.mdp -c md_tester.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized_A.gro ')
         call system ('rm -f input.tpr md.log mdout.mdp traj.trr ener.edr')
! lbfgs --> steep        
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_C.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro out')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f ./#*# out')

        endif
       
         if(k.ge.2) then


!          DIST_scale_x = bhbmap_acc_x(l) - bhbmap_donor_x(l)
!          DIST_scale_y = bhbmap_acc_y(l) - bhbmap_donor_y(l)
!          DIST_scale_z = bhbmap_acc_z(l) - bhbmap_donor_z(l)

!          DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
!               scalekeyB2(l)  = -bEXPOSI(l)  / (RATER2*DIST_scale_A)
               
      
              do p=1,m_don4(l)-1
                                        
                 bnew_dir_acc_xA(p,l) = scaling2c*(bhbmap_donor_x(l)-bhbmap_acc_x(l))+bdonor_coord_x(p,l) 
                 bnew_dir_acc_yA(p,l) = scaling2c*(bhbmap_donor_y(l)-bhbmap_acc_y(l))+bdonor_coord_y(p,l) 
                 bnew_dir_acc_zA(p,l) = scaling2c*(bhbmap_donor_z(l)-bhbmap_acc_z(l))+bdonor_coord_z(p,l)
        
              enddo 
             
        open(unit=83,file='minimized.gro')

        read(83,*,end=9016,err=9016)
        read(83,*,end=9016,err=9016)
        
        open(unit=84,file='md_tester.gro')
     
        write(84,*)
        write(84,*) NATOMS

        p = 1
       
            do m=1,NATOMS
 

               read(83,'(i5,2a5,i5,3f8.3,3f8.4)',end=9016,err=9016) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

               if(resnum.eq.bdonor_resnum(p,l) .and. p.le.m_don4(l)-1) then

                   coord_x = bnew_dir_acc_xA(p,l)
                   coord_y = bnew_dir_acc_yA(p,l)
                   coord_z = bnew_dir_acc_zA(p,l)

                   p = p + 1

               write(*,*) m ,'m'

               endif

               write(84,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(84,*) box_x,box_y,box_z

               close(unit=83)
               close(unit=84)

9016       continue
               
           open(unit=73,file='temp.ndx')
       
           write(73,'(a7)') '[ ACC ]'
        
!           do p=1,m_acc4(l)-1

           write(73,*)  bhbmap_acc_num(l)
!           write(73,*) bacc_atomnum(p,l)

!           enddo
      
           write(73,'(a7)') '[ DON ]'

           if(keepresmode2==0) then

               write(73,*)  bhbmap_donor_num(l)            

           else


               do p=1,m_don4(l)-1
              
                  write(73,*)  bdonor_atomnum(p,l)
      
               enddo

           endif
    
           close(unit=73)

      
      
          call system('cat temp.ndx index.ndx > indexsim2.ndx')
          call system('rm -f temp.ndx')

          call system ('cp -f minim3C.mdp minimtest.mdp') 
    
          open(unit=74,file='minimtest.mdp')

          do j = 1,minimline
                      
             read(74,*)
          
          enddo
   
             write(74,*) 'freezegrps        =' , Nameacc , Namedonor
             write(74,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=74)             
        
         call system ('rm -f input.tpr minimized.gro ener.edr md.log traj.trr mdout.mdp')                    
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest.mdp -c md_tester.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized_A.gro ')
         call system ('rm -f traj.trr mdout.mdp ener.edr md.log input.tpr')
! bfgs --> steep         
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_C.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f ./#*# out')
        
!          if(modulo(k,modnumber)==0) then
          
!              call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim3C.mdp -c minimized_A.gro -p LOV2.top -o input.tpr')
!              call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0 -v -s input.tpr -c minimized.gro')
!              call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro')
         
!           endif
        endif 
               scaling2c = scaling2c + scalekey2c
      
        enddo 

       call system('rm -f input.tpr traj.trr ener.edr mdout.mdp md.log')
!       call system('/loctmp/pee18323/GRO_bin/bin/grompp -f minimsystem.mdp -c md_cheater2.gro -p LOV2_bak.top -o input.tpr')
!       call system('rm -f md_cheater3.gro')
!       call system('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0 -v -s input.tpr -c md_cheater3.gro')
!       call system('rm -f traj.trr ener.edr mdout.mdp md.log')

 call random_number(harvest3)
                     

 timeboolean = .false.

       Delta_time = - log(harvest3)/ csumrate

       write(111,*) Delta_time,'delta_time'

       if(Delta_time .ge. maxtime) then

                          write(111,*) Delta_time,'delta_time too large'

                             Delta_time = 0

                           timeboolean  = .true.

       endif

 
        
       open(unit=38,file='time')
       read(38,*,err=9999,end=9999) 
       read(38,*,err=9999,end=9999) time2


       
       timeend = time2 + Delta_time

       close(unit=38)

       open(unit=38,file='time')
       read(38,*)
       write(38,*) timeend 
       write(38,*) ende2 , 'ende'
       write(38,*) y

       close(unit=38)

9999   continue

       write(111,*) '###############'
       write(111,*) 'event_btime',timeend
       write(111,*) 'cycle',y
       write(111,*) '#####ende2#####'

!       exit
 
   endif
        
          
       call system('rm -f input.tpr mdout.mdp minimtest.mdp traj.trr md.log ener.edr')
       call system('rm -f indexsim2.ndx')


 if(y==1) then
 
     open(unit=128,file='xvg_rate_x')

 endif

   if(selectrate == 1) then

   l = selnumber

   next = .false.

   write(111,*) '################'
   write(111,*) incount,'  incount,break'
   write(111,*) rate_G(l-1), ' rate_G(l-1)'
   write(111,*) harvest3*sumrate, 'harvest3sum'
   write(111,*) harvest3 , 'harvest3'
   write(111,*) csumrate, ' csumrate'
   write(111,*) rate_G(l)  , ' rate_G(l)'
   write(111,*) MINI_A(l) , 'minimumanfang'
   write(111,*) MINI_G(l) , 'minimum_danach'
   write(111,*) MAXI_G(l) , 'maximum'
   write(111,*) EXPOSI(l) , 'exposi'
   write(111,*) DIST_MAXI(l) ,'dist_maxi'
   write(111,*) rate_g(l)   , 'rate'
   write(111,*) l           , 'l'
   write(111,*) hbmap_donor_num_real(l),'hbmap_donor_num_real(l)'
   write(111,*)  hbmap_hyd_num(l) ,'hbmap_hyd_num(l)'
   write(111,*) hbmap_acc_num_real(l)  ,'hbmap_acc_num_real(l)'
   write(111,*) '#################'

!   hb_don_A = hbmap_donor_num_real(l)
!   hb_hyd_A = hbmap_hyd_num(l)
!   hb_acc_A = hbmap_acc_num_real(l)

   incount  = incount + 1


         write(128,*)

   do k=1,endacc
         write(128,*) DISTANCE2(l,k)  , TOTAL_EN(l,k)-TOTAL_EN(l,1)
   enddo
 
         write(128,*)




scaling2c = 0
scalekey2c = 0


   do k=1,endacc

        if(k==1) then
              
               DIST_scale_x = hbmap_acc_x_real(l) - hbmap_donor_x_real(l)
               DIST_scale_y = hbmap_acc_y_real(l) - hbmap_donor_y_real(l)
               DIST_scale_z = hbmap_acc_z_real(l) - hbmap_donor_z_real(l)

               DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
               scalekey2c  = EXPOSI(l)  / (endacc*DIST_scale_A)
               
      
              do p=1,m_don3(l)-1
                                        
                 new_dir_acc_xA(p,l) = scaling2c*(hbmap_donor_x_real(l)-hbmap_acc_x_real(l))+donor_coord_x(p,l) 
                 new_dir_acc_yA(p,l) = scaling2c*(hbmap_donor_y_real(l)-hbmap_acc_y_real(l))+donor_coord_y(p,l) 
                 new_dir_acc_zA(p,l) = scaling2c*(hbmap_donor_z_real(l)-hbmap_acc_z_real(l))+donor_coord_z(p,l)
        
              enddo 
             
        open(unit=31,file='minimized3.gro')

        read(31,*,end=9017,err=9017)
        read(31,*,end=9017,err=9017)
        
        open(unit=32,file='md_tester.gro')
     
        write(32,*)
        write(32,*) NATOMS

        p = 1
       
            do m=1,NATOMS
 

               read(31,'(i5,2a5,i5,3f8.3,3f8.4)',end=9017,err=9017) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

               if(resnum.eq.donor_resnum(p,l) .and. p.le.m_don3(l)-1) then

                   coord_x = new_dir_acc_xA(p,l)
                   coord_y = new_dir_acc_yA(p,l)
                   coord_z = new_dir_acc_zA(p,l)

                   p = p + 1

               write(*,*) m ,'m'

               endif

               write(32,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(32,*) box_x,box_y,box_z

               close(unit=32)
               close(unit=31)

9017       continue
               
           open(unit=33,file='temp.ndx')
       
           write(33,'(a7)') '[ ACC ]'
        
           
!           do p=1,m_acc3(l)-1  

           write(33,*)  hbmap_acc_num_real(l)
!            write(33,*)  acc_atomnum(p,l)

!           enddo
      
           write(33,'(a7)') '[ DON ]'

           if(keepresmode1==0) then

              write(33,*)  hbmap_donor_num_real(l) 
      
           else

             do p=1,m_don3(l)-1
        
!           write(33,*)  hbmap_donor_num_real(l)
               write(33,*)  donor_atomnum(p,l)
           
             enddo
    
           endif

           close(unit=33)

           
          call system('cat temp.ndx index.ndx > indexsim2.ndx')
          call system('rm -f temp.ndx')

          call system('cp -f indexsim2.ndx testindex2.ndx')

          call system ('cp -f minim.mdp minimtest.mdp') 
    
          open(unit=34,file='minimtest.mdp')

          do j = 1,minimline
                      
             read(34,*)
          
          enddo
   
             write(34,*) 'freezegrps        =' , Nameacc , Namedonor
             write(34,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=34)               

         call system ('rm -f input.tpr minimized.gro traj.trr ener.edr mdout.mdp md.log')                            
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest.mdp -c md_tester.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized_A.gro ')
         call system ('rm -f input.tpr ener.edr mdout.mdp traj.trr md.log out')
! bfgs --> steep
        
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro out')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f ./#*# out')         


         endif

    if(k.ge.2)  then


!               DIST_scale_x = hbmap_acc_x_real(l) - hbmap_donor_x_real(l)
!               DIST_scale_y = hbmap_acc_y_real(l) - hbmap_donor_y_real(l)
!               DIST_scale_z = hbmap_acc_z_real(l) - hbmap_donor_z_real(l)

!               DIST_scale_A = sqrt(DIST_scale_x**2 + DIST_scale_y**2 + DIST_scale_z**2) 
          
!               scalekeyB(l)  = EXPOSI(l)  / (RATER*DIST_scale_A)

              do p=1,m_don3(l)-1
                                        
                 new_dir_acc_xA(p,l) = scaling2c*(hbmap_donor_x_real(l)-hbmap_acc_x_real(l))+donor_coord_x(p,l) 
                 new_dir_acc_yA(p,l) = scaling2c*(hbmap_donor_y_real(l)-hbmap_acc_y_real(l))+donor_coord_y(p,l) 
                 new_dir_acc_zA(p,l) = scaling2c*(hbmap_donor_z_real(l)-hbmap_acc_z_real(l))+donor_coord_z(p,l)
        
              enddo 

         open(unit=80,file='minimized.gro')
         open(unit=81,file='md_tester.gro')
          
         read(80,*,end=9018,err=9018)
         read(80,*,end=9018,err=9018)
         write(81,*)
         write(81,*)   NATOMS
               
         p = 1
       
            do m=1,NATOMS
 

               read(80,'(i5,2a5,i5,3f8.3,3f8.4)',end=9018,err=9018) resnum, restype, atomtype, atomnum, coord_x, &
                                     coord_y, coord_z, vel_x , vel_y , vel_z

               if(resnum.eq.donor_resnum(p,l) .and. p.le.m_don3(l)-1) then

                   coord_x = new_dir_acc_xA(p,l)
                   coord_y = new_dir_acc_yA(p,l)
                   coord_z = new_dir_acc_zA(p,l)

                   p = p + 1

               write(*,*) m ,'m'

               endif

               write(81,'(i5,2a5,i5,3f8.3,3f8.4)') resnum, restype, atomtype, atomnum, coord_x, &
                                                coord_y, coord_z, vel_x , vel_y , vel_z
            enddo
              
              write(81,*) box_x,box_y,box_z

               close(unit=80)
               close(unit=81)

9018     continue
      
         open(unit=33,file='temp.ndx')
       
           write(33,'(a7)') '[ ACC ]'
        
!           do p=1,m_acc3(l)-1           
           
           write(33,*)  hbmap_acc_num_real(l)
!           write(33,*)  acc_atomnum(p,l)

!           enddo
      
           write(33,'(a7)') '[ DON ]'

           if(keepresmode1==0) then

             write(33,*) hbmap_donor_num_real(l)
 
           else
     
             do p=1,m_don3(l)-1
        

              write(33,*) donor_atomnum(p,l)
           
             enddo

           endif
    
           close(unit=33)

      
      
          call system('cat temp.ndx index.ndx > indexsim2.ndx')
          call system('rm -f temp.ndx')

          call system ('cp -f minim.mdp minimtest.mdp') 
    
          open(unit=34,file='minimtest.mdp')

          do j = 1,minimline
                      
             read(34,*)
          
          enddo
   
             write(34,*) 'freezegrps        =' , Nameacc , Namedonor
             write(34,*) 'freezedim         =    Y Y Y  Y Y Y    '
            

            close(unit=34)               
                            
         call system ('rm -f input.tpr minimized.gro md.log mdout.mdp  traj.trr')
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minimtest.mdp -c md_tester.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized_A.gro ')
         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp out')
! A--> bfgs
! B--> steep
         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B.mdp -c minimized_A.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_B2.mdp -c minimized.gro -p LOV2.top -o input.tpr -n indexsim2.ndx ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0  -s input.tpr -c minimized.gro ')
         call system ('rm -f ./#*#') 

         call system ('rm -f input.tpr md.log traj.trr ener.edr mdout.mdp minimized_A.gro out')

         endif
   
!         if(modulo(k,modnumber)==0) then

!         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim.mdp -p LOV2.top -o input.tpr -c minimized.gro')
!         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -nice 0 -v -s input.tpr -c minimized.gro')
!         call system ('rm -f \#*# input.tpr md.log traj.trr ener.edr mdout.mdp')

!         endif
!         RATER  = RATER-1
         scaling2c = scaling2c + scalekey2c


enddo


       timeboolean = .false.

 call random_number(harvest3)

       
       Delta_time = - log(harvest3)/ csumrate
 
                      write(111,*) 'delta_time',delta_time


       if(Delta_time .ge. maxtime) then

                      write(111,*) 'delta_time too large !!!',delta_time

                      delta_time = 0


                      timeboolean = .true.

       endif
        
       open(unit=38,file='time')
       read(38,*,err=9991,end=9991) 
       read(38,*,err=9991,end=9991) time2
       
       timeend = time2 + Delta_time

       close(unit=38)

       open(unit=38,file='time')
       read(38,*)
       write(38,*) timeend
       write(38,*) ende , 'ende'
       write(38,*) y

       close(unit=38)

9991   continue

       write(111,*) '###############'
       write(111,*) 'event_time',timeend
       write(111,*) 'cycle',y
       write(111,*) '#####ende1#####' 
  
!       exit 
 
   endif
           
       call system('rm -f input.tpr mdout.mdp minimtest.mdp traj.trr md.log ener.edr')
       call system('rm -f indexsim2.ndx')

!enddo
hazardbool = .false.

   

            allocate(resnum_array(NATOMS),stat=ierr)
            allocate(atomtype_array(NATOMS),stat=ierr)
            allocate(coord_array_x(NATOMS),stat=ierr)
            allocate(coord_array_y(NATOMS),stat=ierr)
            allocate(coord_array_z(NATOMS),stat=ierr)

            allocate(atomnum_array(NATOMS),stat=ierr)


 
           open(unit=2,file='minimized.gro')

            read(2,*,end=12128,err=12128)
            read(2,*,end=12128,err=12128)


            do k = 1,NATOMS
              
             read(2,'(i5,2a5,i5,3f8.3,3f8.4)',end=12128,err=12128) resnum_array(k), restype, atomtype_array(k), atomnum_array(k), coord_array_x(k), &
             &   coord_array_y(k), coord_array_z(k), vel_x , vel_y , vel_z

            enddo

            close(unit=2)

            12128 continue

            do p=1,NATOMS

                   do m=1,NATOMS

                          if(atomnum_array(p) .ne. atomnum_array(m)) then

                          diff_x = coord_array_x(p) - coord_array_x(m)
                          diff_y = coord_array_y(p) - coord_array_y(m)
                          diff_z = coord_array_z(p) - coord_array_z(m)

                          diff_tot = sqrt(diff_x**2 + diff_y**2 + diff_z**2)



                          if(diff_tot .lt. hazardparamreal) then


                             write(*,*) diff_tot


                             write(111,*) 'diff_tot',diff_tot
                             write(111,*) atomnum_array(p),'atomnum_array(p)',atomnum_array(m),'atomnum_array(m)'

                             hazardbool = .true.

                             exit
                          endif 

                          endif

                   enddo

                          if(hazardbool == .true.) then

                             exit

                          endif

             enddo

       deallocate(resnum_array,stat=ierr)
       deallocate(atomtype_array,stat=ierr)
       deallocate(coord_array_x,stat=ierr)
       deallocate(coord_array_y,stat=ierr)
       deallocate(coord_array_z,stat=ierr)

       deallocate(atomnum_array,stat=ierr)


maxprotx = -100
maxproty = -100
maxprotz = -100

minprotx = 100
minproty = 100
minprotz = 100


timefirst = timefirst + Delta_time


if(next == .true. .or. hazardbool == .true. .or. timeboolean == .true.) then

write(111,'(a38,f6.4)') 'HAZARD OCCURED, DIFF_TOT in PROTEIN <<',hazardparamreal

write(111,*) 'kmc_number, hazard',tu

write(111,*) timeend,'timeend'

    if(timeboolean == .true.) then

              write(111,*) 'time exceeded !!!',maxtime

    endif

exit

endif

if(next == .false. .and. hazardbool == .false. .and. timeboolean == .false.) then

              write(111,*) 'step occured !!!'

              write(111,*) tu,'kmc_number'

              write(111,*) timeend ,'total_time'

              write(111,*) '#############'

              if(tu.ge.2) then

                          call system('cp minimized.gro minimized3.gro')

              endif

endif

if(timefirst .ge. timekonst_delta ) then

                   write(111,*) 'timekonst_delta reached !' , timefirst

                   exit

endif

if(ende == 0) then

                   write(111,*) 'no hydrogen bonds found !'

                   exit

endif

enddo

write(111,*) next,'next',hazardbool,'hazardbool',timeboolean,'timeboolean'


write(*,*) maxprotx,maxproty,maxprotz,'max'
write(*,*) minprotx,minproty,minprotz,'min'






allocate(coord_test_x(NATOMS),stat=ierr)
allocate(coord_test_y(NATOMS),stat=ierr)
allocate(coord_test_z(NATOMS),stat=ierr)

if(next == .false..and. hazardbool == .false. .and. timeboolean == .false.) then

         call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f minim.mdp -c minimized.gro -p LOV2.top -o finalminim.tpr -nice 0 ')
         call system ('/loctmp/pee18323/GRO_bin/bin/mdrun -v -s finalminim.tpr -nice 0 -c minimized.gro')
         call system ('rm -f ./#*#')

endif

if(next ==.true.  .or. hazardbool == .true. .or. timeboolean == .true.) then

         open(unit=1111,file='minimized3.gro')

endif

if(next == .false..and. hazardbool == .false. .and. timeboolean == .false.) then

         open(unit=1111,file='minimized.gro')

endif

read(1111,*)
read(1111,*)

do m=1,NATOMS

           read(1111,'(i5,2a5,i5,3f8.3,3f8.4)',end=91111,err=91111) resnum, restype, atomtype, atomnum, coord_test_x(m), &
                                             coord_test_y(m), coord_test_z(m), vel_x , vel_y , vel_z

enddo

91111     continue

close(unit=1111)


if(next ==.true. .or. hazardbool == .true. .or. timeboolean == .true. ) then

         call system('/loctmp/pee18323/GRO_bin/bin/editconf -f minimized3.gro -c -o minimized3.gro')

         call system('cp minimized3.gro prot.gro')

endif

if(next == .false..and. hazardbool == .false. .and. timeboolean == .false.) then

         call system('/loctmp/pee18323/GRO_bin/bin/editconf -f minimized.gro -c -o minimized.gro')

         call system('cp minimized.gro prot.gro')
endif


call solvent

call system('cp first.gro md_cheater2.gro')


         call system('rm -f fullmd.tpr')

if(next == .true.) then

         write(111,*)'#########nothing happened############'
         write(111,*)'##########delta_G exceeded###########'
         write(111,*)'nextstep             '
         write(111,*)'#####################################'
         write(111,*) y

endif

deallocate(coord_test_x,stat=ierr)
deallocate(coord_test_y,stat=ierr)
deallocate(coord_test_z,stat=ierr)

!write(*,*) 'h7'

if(allocated(resnumsol)) then
deallocate(resnumsol,stat=ierr) 
if(ierr .ne. 0 ) then
write(*,*) 'dealloc resnumsol'
continue
endif
endif
if(allocated(restypesol)) then
deallocate(restypesol,stat=ierr)
if(ierr .ne. 0 ) then
write(*,*) 'dealloc restypesol'
continue
endif
endif
if(allocated(atomtypesol)) then
deallocate(atomtypesol,stat=ierr)
if(ierr.ne.0) then
write(*,*) 'dealloc atomtypesol'
continue
endif
endif
if(allocated(atomnumsol)) then
deallocate(atomnumsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc atomnumsol'
continue
endif
endif
if(allocated(coord_xsol)) then
deallocate(coord_xsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_xsol'
continue
endif
endif
if(allocated(coord_ysol)) then
deallocate(coord_ysol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_ysol'
continue
endif
endif
if(allocated(coord_zsol)) then
deallocate(coord_zsol,stat=ierr)
if(ierr .ne. 0) then
write(*,*) 'dealloc coord_zsol'
continue
endif
endif

write(*,*) 'h8'

          
       a = a + 1
       
       if(y.eq.1) then

         call system('mv md_traj.xtc md_trajx.xtc')

         call system('mv md_traj.trr md_trajx.trr')

       endif
        
       if(y.ge.2.and.runbool==.true.) then

!         call system('mv \#md_traj.trr.1# \md_trajy.trr')
      
         open(unit=453,file='trjcat.sh')

          write(453,*) '#!/bin/bash'
          write(453,*)
          write(453,'(a96)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx.xtc md_traj.xtc -o md_trajx.xtc -settime << EOF'
          write(453,*) '0'
          write(453,*) timer
          write(453,*) 'EOF'

         close(unit=453)

         open(unit=455,file='trjcat2.sh')


          write(455,*) '#!/bin/bash'
          write(455,*)
          write(455,'(a96)') '/loctmp/pee18323/GRO_bin/bin/trjcat -f md_trajx.trr md_traj.trr -o md_trajx.trr -settime << EOF'
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

         open(unit=456,file='time')
         
        read(456,*)
        read(456,*)
        read(456,*) 
       write(456,*) y
        

       close(unit=456)
            
!       call system('./Ionizer.sh')

       call system('/loctmp/pee18323/GRO_bin/bin/grompp -f minim_finalA.mdp -c md_cheater2.gro -p LOV2gen.top -o input.tpr')

       call system('/loctmp/pee18323/GRO_bin/bin/mdrun -v  -s input.tpr -c md_cheater2.gro')

       call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_solequil.mdp -c md_cheater2.gro -p LOV2gen.top -o relax.tpr')

       call system('mpirun -np 2 /loctmp/pee18323/GRO_bin/bin/mdrun -v -s relax.tpr -nice 0 -c md_cheater2.gro')

       mdboolean = .false.


       if(timeend .ge. maxtime) then

                              mdboolean = .true.

                              exit

       endif


       call system('rm -f ./#*#')

       call system('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_sol.mdp -c md_cheater2.gro -p LOV2gen.top -o fullmd.tpr')

       call system('rm -f table_hbond_B.ndx table_hbond_B2.ndx table_hbondC.ndx table_hbondC2.ndx')
enddo

       if(mdboolean == .true.) then

            open(unit=456,file='time')
       
             read(456,*)
             read(456,*)
             read(456,*)
             read(456,*) 
             write(456,*) y,'maxtime reached!!!'
             write(456,*) maxtime,'maxtime'
        

       close(unit=456)
            

            relaxtime = kmccycles * mdtime

 open(unit=432,file='fullmd_sol.mdp')

 write(432,*) 'title                    = MD simulation LOV2-Jalpha'

 write(432,*) ';Preprocessor'
 write(432,*) 'cpp			 = /lib/cpp'
 write(432,*) ';Directories to include in the topology format'
 write(432,*) 'include 		 = -I../top'
 write(432,*) ';Run control: A leap-frog algorithm for integrating Newtons equations.' 
 write(432,*) 'integrator		 = md'
 write(432,*) ';Total simulation time: 100 ps'
 write(432,*) ':time step in femtoseconds'  
 write(432,*) 'dt			 = 0.001'
 write(432,*) ';number of steps'
 write(432,'(1x,a12,i7)') 'nsteps  		 =',relaxtime
 write(432,*) ';frequency to write coordinates to output trajectory file'
 write(432,*) 'nstxout 		 = 1000'
 write(432,*) ';frequency to write velocities to output trajectory file'
 write(432,*) 'nstvout 		 = 1000'
 write(432,*) ';frequency to write energies to log file'
 write(432,*) 'nstlog  		 = 1000'
 write(432,*) ';frequency to write energies to energy file'
 write(432,*) 'nstenergy		 = 1000'
 write(432,*) ';frequency to write coordinates to xtc trajectory' 
 write(432,*) 'nstxtcout		 = 1000'
 write(432,*) ';group(s) to write to xtc trajectory'
 write(432,*) 'xtc_grps		 = protein CFP'  
 write(432,*) 'group(s) to write to energy file' 
 write(432,*) 'energygrps		 = protein sol CFP'
 write(432,*) ';Frequency to update the neighbor list (and the long-range forces,' 
 write(432,*) ';when using twin-range cut-offs).' 
 write(432,*) 'nstlist 		 = 10'
 write(432,*) ';Make a grid in the box and only check atoms in neighboring grid cells'
 write(432,*) ';when constructing a new neighbor list every nstlist steps.' 
 write(432,*) 'ns_type 		 = grid'
 write(432,*) ';cut-off distance for the short-range neighbor list'
 write(432,*) 'rlist		 = 1.4 '
 write(432,*) ';treatment of electrostatic interactions'
 write(432,*) 'coulombtype		 = pme'
 write(432,*) 'pme_order            = 4'
 write(432,*) 'fourierspacing       = 0.12' 
 write(432,*) 'rcoulomb		 = 1.4' 
 write(432,*) ';treatment of van der waals interactions'
 write(432,*) 'rvdw			 = 1.4'
 write(432,*) '; Periodic boudary conditions in all the directions' 
 write(432,*) 'pbc                      = xyz'
 write(432,*) ';Temperature coupling'
 write(432,*) 'tcoupl  		 = nose-hoover'
 write(432,*) 'tc-grps 		 = protein  CFP SOL NA+ GTP'   
 write(432,*) 'tau_t	        = 0   0    0.1 0.1 0' 
 write(432,*) 'ref_t		 = 300 300 300 300 300' 
 write(432,*) ';Pressure coupling'
 write(432,*) 'Pcoupl  		 = Parrinello-Rahman'
 write(432,*) 'Pcoupltype           = isotropic'
 write(432,*) 'tau_p		 = 1.0'
 write(432,*) 'compressibility 	 = 4.5e-5'
 write(432,*) 'ref_p		 = 1.0'
 write(432,*) ';Velocity generation'
 write(432,*) 'gen_vel 		 = yes' 
 write(432,*) 'gen_temp		 = 300'
 write(432,*) 'gen_seed		 = 173529'
 write(432,*) ';Constrain all bonds'
 write(432,*) 'constraints		 = none'

 close(unit=432)

            call system ('/loctmp/pee18323/GRO_bin/bin/grompp -f fullmd_sol.mdp -c md_cheater2.gro -p LOV2gen.top -o md_relax.tpr')

            write(*,*) 'here is md relax time'

            call system ('/loctmp/pee18323/GRO_bin/bin/tpbconv -s md_relax.tpr -f md_trajx.trr -o md_relax2.tpr -nice 0')

            call system ('mpirun -np 2 /loctmp/pee18323/GRO_bin/bin/mdrun -v -s md_relax2.tpr -nice 0 -c md_finaloutrelax.gro -o md_trajrelax.trr -x md_traj.xtc -e md_enerrelax.edr')

            call system ('trjcat -f md_trajx.trr md_trajrelax.trr -o md_trajtotalx.trr -nice 0')

            write(*,*) 'relax finished'

endif

end subroutine Potential                 
