!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfengs@hotmail.com                         =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!
!
    subroutine setmemory
    use mode
    use mod_info
    implicit none
!
    integer       i

    allocate(cflag(npart))
    allocate(molecule(npart))

!
    allocate(x(3,npart))
    allocate(v(3,npart))
    allocate(q(npart))
    allocate(totq(npart))
    allocate(types(npart))
    allocate(atable(1+matom,npart))
!
!
    allocate(molelist(npart))
    allocate(mlist(npart))
!
    do i=1,npart
       q(i)   = 0.0d0
       totq(i)= 0.0d0
    enddo
!
    return
    end 
!
!
    subroutine release
    use mode
    use mod_info
    implicit none

!
    deallocate(cflag)
    deallocate(molecule)
    deallocate(n_1)
    deallocate(n_2)
    deallocate(n_3)
    deallocate(n_4)
    deallocate(n_5)
    deallocate(x)
    deallocate(v)
    deallocate(types)
    deallocate(atable)
    deallocate(mlist)
    deallocate(molelist)
!
!
    return
    end
!
!
