!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@hotmail.com                     =
!=========================================================================
!
!-------------------------------------------------------------------------
!
      subroutine read_tra(iflag)
      use mode,only:step,npart,types,q,x,ctype
      use mod_pbc
      implicit none
!
!
!
      character             head
      integer               i,l,m, iflag
      integer               id
      character(80)         line
      real(8)               lo,hi
      real(8)               xlo_bound, xhi_bound, xy, xlo, xhi
      real(8)               ylo_bound, yhi_bound, xz, ylo, yhi
      real(8)               zlo_bound, zhi_bound, yz, zlo, zhi 
      real(8)               cos_alpha,cos_beta,cos_gamma
!
      iflag=0
!
      read(10,*,end=100)head
!
      read(10,*)step
      read(10,*)
      read(10,*)npart
      read(10,'(A80)')line
!
!
      if(index(line,'ITEM: BOX BOUNDS')/=0)then
           if(index(line,'xy xz yz')/=0)then      
                read(10,*)xlo_bound, xhi_bound, xy
                read(10,*)ylo_bound, yhi_bound, xz
                read(10,*)zlo_bound, zhi_bound, yz
                xlo= xlo_bound - min(0.0,xy,xz,xy+xz)
                xhi= xhi_bound - MAX(0.0,xy,xz,xy+xz)
                ylo= ylo_bound - MIN(0.0,yz)
                yhi= yhi_bound - MAX(0.0,yz)
                zlo= zlo_bound 
                zhi= zhi_bound
                a = (/xhi-xlo, 0.0d0, 0.0d0/)
                b = (/xy, yhi-ylo, 0.0d0/)
                c = (/xz, yz, zhi-zlo/)

                r1 = a(1)*a(1)+ a(2)*a(2)+ a(3)*a(3)
                r1 = dsqrt(r1)
                r2 = b(1)*b(1)+ b(2)*b(2)+ b(3)*b(3)
                r2 = dsqrt(r2)
                r3 = c(1)*c(1)+ c(2)*c(2)+ c(3)*c(3)
                r3 = dsqrt(r3)

                cos_alpha = (b(1)*c(1) + b(2)*c(2) + b(3)*c(3))/(r2*r3)
                cos_beta  = (c(1)*a(1) + c(2)*a(2) + a(3)*c(3))/(r3*r1)
                cos_gamma = (a(1)*b(1) + a(2)*b(2) + a(3)*b(3))/(r1*r2)

                alpha = dacos(cos_alpha)
                beta  = dacos(cos_beta)
                gamma = dacos(cos_gamma)
           else
               read(10,*)lo,hi
               a(1)= hi-lo
               r1= a(1)
               read(10,*)lo,hi
               b(2)= hi-lo
               r2= b(2)
               read(10,*)lo,hi
               c(3)= hi-lo
               r3= c(3)
           endif
      endif
!
      read(10,*)
!
       do i=1,npart
          read(10,*)id,types(id),(x(l,id),l=1,3),q(id)!,(v(m,i),m=1,3)
       enddo
!
      iflag=1
!
100   continue
!
      return
      end
!
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
    subroutine read_mol(fin)
    use mode,only:npart,x,types
    implicit none
!
!   Local Variables
!
    integer fin,i
    character cc
!
    read(fin,*)
    read(fin,*)
    read(fin,*)
    read(fin,*)npart
!
    do i=1,npart
       read(fin,*)x(1,i),x(2,i),x(3,i),cc
       if(cc=='C')types(i)=1
       if(cc=='H')types(i)=2
       if(cc=='N')types(i)=4
       if(cc=='O')types(i)=3
    enddo
!
   return
   end
!
!
!-------------------------------------------------------------------------
!
      subroutine read_xyz(iflag)
      use mode,only:step,npart,types,q,x,ctype,v
      implicit none
!
      integer            nn
      integer            i,l,m, iflag
      character(2)       cty
      
!      real(8)     q
!
      iflag=0
!
      read(10,*,end=100)nn
      read(10,*)
!
       do i=1,npart
          read(10,*)cty,(x(l,i),l=1,3)!,(v(m,i),m=1,3)
          if(cty==ctype(1))types(i)=1
          if(cty==ctype(2))types(i)=2
          if(cty==ctype(3))types(i)=3
          if(cty==ctype(4))types(i)=4
       enddo
!
      iflag=1
!
100   continue
!
      return
      end
!
!-------------------------------------------------------------------------
!
