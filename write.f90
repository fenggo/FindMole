!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@hotmail.com                      =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!
!
!
    subroutine write_xyz(fin)
    use mode,only:npart,nstep,ctype,x,types
    implicit none
!
!   Local Variables
!
    integer i,fin
!
    write(fin,*)npart
    write(fin,*)'nstep: ',nstep
    do i=1,npart
       write(fin,500)ctype(types(i)),x(1,i),x(2,i),x(3,i)
    enddo
!
500 FORMAT(A3,3F20.6)
!
    return
    end
