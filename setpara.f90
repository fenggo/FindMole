!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                       email: gfeng.alan@qq.com                         =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!
    Subroutine setpara
    use mode
    use mod_pke
    implicit none
!
!
    integer i,j
    real(8) bondcutoffp(5,5)
    real(8) dep(5,5),r0p(5,5)
    real(8) rc(5)
!
!
!   data rc /1.3742, 0.6867, 1.3142, 1.2456/
    data rc /1.3742, 0.80, 1.3142, 1.25,2.1892/
! 
!    data ((bondcutoffp(i,j),i=1,4),j=1,4) /1.87,1.46,1.78,1.76, &
!                                           1.46,1.06,1.37,1.36, &
!                                           1.78,1.37,1.68,1.67, &
!                                           1.76,1.36,1.67,1.25/
!
!  original sigma bond of N: 1.29945
!  modified it to 1.4 because of the C-N bond
!  vibration

   !data r0p /1.37630,1.02045,1.32310,1.29945,&   ! Not used !!!!!
            ! 1.02045,0.66460,0.96725,0.94360,&
             !1.32310,0.96725,1.26990,1.24625,&
             !1.29945,0.94360,1.24625,1.22260/
!
!
    do i=1,5
       do j=1,5
          rcut(i,j)=1.35*(0.5*rc(i)+0.5*rc(j))
          rcutsq(i,j)=rcut(i,j)*rcut(i,j)
          !de(i,j)=dep(i,j)
          !r0(i,j)=r0p(i,j)
        enddo
    enddo

! For C-N
    rcut(1,4)=1.35*(0.5*1.3742+0.5*1.60)
    rcut(4,1)=1.35*(0.5*1.3742+0.5*1.60)
    rcutsq(1,4)=rcut(1,4)*rcut(1,4)
    rcutsq(4,1)=rcut(4,1)*rcut(4,1)

    ktype(1)=6
    ktype(2)=1
    ktype(3)=8
    ktype(4)=7
    ktype(5)=13
!
    mass(1)=12.000
    mass(2)=1.0080
    mass(3)=15.9990
    mass(4)=14.0000
    mass(5)=26.9820
!
    return
    end
