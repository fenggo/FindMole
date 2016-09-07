!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@hotmail.com                      =
!=========================================================================
!
!-------------------------------------------------------------------------
!  ------             MODULE                    ------                   -
!-------------------------------------------------------------------------
!
    module mode
    implicit none
!
    integer,  parameter:: mbond=8000
    integer,  parameter:: matom=10
    integer,       save:: nbond
    integer,       save:: ico   !coodinate foramte 1: lammpstrj, 2: xyz
!
    integer,       save:: npart,step    
    integer,       save:: nstep = -1  
!   real(8),       save:: cube(3),boxlo(3),boxhi(3)
    real(8),       save:: rcut(5,5),rcutsq(5,5)
    real(8),       save:: mass(5)
    real(8),       save:: dt = 0.1 ! fs

    integer,allocatable,   save:: types(:)
    integer,       save:: ktype(5)
    character(2),  save:: ctype(5) 
!
    real(8),allocatable,   save:: x(:,:)
    real(8),allocatable,   save:: v(:,:)
!
    real(8),allocatable,   save:: q(:)      ! charge
    real(8),allocatable,   save:: totq(:)   ! total charge of a molecular fragment

    integer,  allocatable, save:: atable(:,:)   ! atom-table within the cut-off
    integer,               save:: bondlist(2,mbond)
!
!   mlist(), molelist() : molecule lookup table!
!   mlist() is the 'table of contents'
!   the 'i' molecule is from mlist(i-1)+1 to mlist(i) in the molelist
!
!
!   cflag signs which atom has already considered
!   molecue, mlist, and molelist establelish a list
!   that we can find atom that belong to one molecule
!
!   n_mol    :  Total molecule number
!   molecule :  Atom i belong to molecue 'molecule(i)'
!
!   Atoms in molecue i are stored in molelist() from-
!         mlist(i-1)+1 to mlist(i)
!
    integer,allocatable,   save:: cflag(:)
    integer,allocatable,   save:: mlist(:)
    integer,allocatable,   save:: molecule(:)
    integer,               save:: n_mol
    integer,allocatable,   save:: molelist(:)

    integer,               save:: istep   ! the step number to be analysised in trajectories
!
!
!    data rcut /1.7892,1.3266,1.7200,2.1200,&   !!!! Set in setpara.f90 Now !!!!
!               1.3266,0.8640,1.2574,1.2267,&
!               1.7200,1.2574,1.6509,1.6201,&
!               2.1200,1.2267,1.6201,1.5894/

    data ctype /'C','H','O','N','Al'/
    data ktype /6,1,8,7,13/
    data mass  /12.0000,1.0080,15.9990,14.0000,26.9820/
!
    end module mode
!
!-----------------------------华丽的分割线---------------------------------
!    
    module mod_info
    implicit none

    integer            ,save:: typ(5)
    integer,     allocatable:: n_1(:)
    integer,     allocatable:: n_2(:)
    integer,     allocatable:: n_3(:)
    integer,     allocatable:: n_4(:)
    integer,     allocatable:: n_5(:)

    end module mod_info
!
!-----------------------------华丽的分割线---------------------------------
!
    module mod_pbc
    implicit none
!
    real(8)            ,save:: a(3)   ! lattice parameter
    real(8),            save:: b(3)
    real(8),            save:: c(3)
    real(8),            save:: r1,r2,r3
    real(8),            save:: alpha = 90.0
    real(8),            save:: beta  = 90.0
    real(8),            save:: gamma = 90.0
!
    end module mod_pbc
!
!-----------------------------华丽的分割线---------------------------------
!
    module mod_mass
    implicit none

    real(8), save::      mass_cdust
    real(8), save::      mass_oxid
    real(8), save::      mass_final
    real(8), save::      mass_all

    end module mod_mass
!
!-----------------------------华丽的分割线---------------------------------
!
    module mod_pke
    implicit none
!
    real(8),     save:: pkflag=-1
!
    end module mod_pke
!
!-----------------------------华丽的分割线---------------------------------
!
