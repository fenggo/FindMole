!=========================================================================
!                                                            WRITED BY GUO F.                                                     =
!                                                                                                                                                =
!                                       Institute of Atomic & Molecular Physics of Si-                              =
!                                                     chuan University of China.                                               =
!                                                                                                                                                =
!                                                 email: gfeng.alan@foxmail.com                                           =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------             FIND MOLECULE             ------                   -
!-------------------------------------------------------------------------
!
!
     subroutine findmole
     use mode,only: molelist,mlist,molecule,cflag,n_mol,npart,matom,atable,totq,q
     implicit none    
!
!
     integer        i
     integer        ii
     integer        j
     integer        jj
     integer        cc
     integer        l
     real(8)        tq
!
!
     n_mol=1
!
!   cflag signs which atom has already considered
!   molecue, mlist, and molelist establelish a list
!   that we can find atom that belong to one molecule
!   n_mol    :  Total molecule number
!   molecule :  Atom i belong to molecue 'molecule(i)'
!   Atoms in molecue i are stored in molelist from-
!         mlist(i-1)+1 to mlist(i)
!
     do i=1,npart
        cflag(i)=1
     enddo
! 
!
     do i=1, npart
!
       ii=i
       if(cflag(ii)/=0)then
          cflag(ii)=0
          call trackbond(ii)
          n_mol=n_mol+1
       endif
!
     enddo

!
     n_mol=n_mol-1
     cc=0
!     print 100,'Total Molecule Number:',n_mol
!
      do i=1, n_mol
          tq = 0.0
          do j=1, npart
             if(molecule(j)==i)then
                tq= tq + q(j)
                cc=cc+1
                molelist(cc)=j
             endif
          enddo
          mlist(i)=cc
          totq(i)=tq
     enddo
!
100  format(2X,A,I6)
!
     return
     end
!
!-------------------------------------------------------------------------
!  ------             RECURSIVE TRACBOND             ------              -
!-------------------------------------------------------------------------
!
   
    recursive subroutine trackbond(m)
    use mode, only: atable,cflag,n_mol,molecule,matom
    implicit none
!
    integer      mmm(matom)
    integer      m, mm,i
!
!
!    if(cflag(m)/=0)return
!
    mm=m
    do i=1,matom
       mmm(i)=atable(i+1,mm)
    enddo
!
    molecule(mm)=n_mol
    cflag(mm)=0
!
    do i=1,matom
!
       if(mmm(i)/=0.and.cflag(mmm(i))/=0)then
           cflag(mmm(i))=0
           molecule(mmm(i))=n_mol
           call trackbond(mmm(i))
      endif
!
    enddo
!
    return
    end
!
!-------------------------------------------------------------------------
!

