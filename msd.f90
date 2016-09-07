!=========================================================================
!                       WRITED BY GUO F.                                 =
!                                                                        =
!          Institute of Atomic & Molecular Physics of Si-                =
!                  chuan University of China.                            =
!                                                                        =
!                 email: gfeng.alan@gmail.com                            =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!

subroutine msd
use mode     ,only: npart,x,nstep,mlist,molelist,types,ctype,n_mol,dt,mass,istep,ico
use mod_pbc
implicit none
      
    integer                 i,j,jj
    integer                 rflag,iatom
    integer                 ina,inb,inc

    real(8)                 tx,ty,tz,xct,yct,zct,masses,rt
    real(8),allocatable::   rt0(:,:)
    real(8)                 r_a,r_b,r_c

    rflag=0
    nstep=0

    allocate(rt0(3,npart))
    rt0(:,:) = 0.0 
    open(27,file='msd.txt',status='unknown')
    write(27,'(A9)')'#time msd'

    do   ! main cycle loop

       nstep=nstep+1

       if    (ico==1)then 
             call read_tra(rflag)
       elseif(ico==2)then
             call read_xyz(rflag)
       endif
       if(rflag==0)exit
       if(nstep<istep)cycle

       call neighbor
       call findmole
!
       do i=1, n_mol    ! cycle from molecules
          if(i>1)then
             jj=mlist(i-1) + 1
          else
             jj=1
          endif

!         determine center of mass of the molecule
          xct = 0.0
          yct = 0.0
          zct = 0.0
          masses = 0.0
          rt = 0.0

          do j=jj, mlist(i)
             iatom = molelist(j)
             xct = xct + x(1,iatom)*mass(types(iatom))
             yct = yct + x(2,iatom)*mass(types(iatom))
             zct = zct + x(3,iatom)*mass(types(iatom))
             masses = masses + mass(types(iatom))
          enddo
          xct = xct/masses
          yct = yct/masses
          zct = zct/masses    ! mass of center
          if(nstep==istep)then
              rt0(1,i)=xct
              rt0(2,i)=yct
              rt0(3,i)=zct
          elseif(nstep>istep)then
              tx = xct-rt0(1,i)
              ty = yct-rt0(2,i)
              tz = zct-rt0(3,i)
!   apply PBC
        
              r_c = (a(1)*b(2)-a(2)*b(1))*(b(2)*tz-b(3)*ty) - (a(3)*b(2)-a(2)*b(3))*(b(2)*tx-b(1)*ty)
              r_c = r_c/((a(1)*b(2)-a(2)*b(1))*(c(3)*b(2)-c(2)*b(3)) - (a(3)*b(2)-a(2)*b(3))*(c(1)*b(2)-c(2)*b(1)))

              r_a = b(2)*tx-b(1)*ty - (c(1)*b(2)-c(2)*b(1))*r_c
              r_a = r_a/(a(1)*b(2)-a(2)*b(1))
              r_b = ty- r_a*a(2) - r_c*c(2)
              r_b = r_b/b(2)

              ina = anint(r_a)
              inb = anint(r_b)
              inc = anint(r_c)

              tx = tx - ina*a(1) - inb*b(1) - inc*c(1)
              ty = ty - ina*a(2) - inb*b(2) - inc*c(2)
              tz = tz - ina*a(3) - inb*b(3) - inc*c(3)
              rt=rt+tx**2+ty**2+tz**2
          endif
       enddo

       rt = rt/n_mol
       write(27,*)(nstep-istep)*dt,rt
       print *,nstep-istep,rt
    enddo   

    deallocate(rt0)
    close(27)

return
end

