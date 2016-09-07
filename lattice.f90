
subroutine lattice
    use mode, only: npart,types,mass,ctype,ico
    use mod_mass
    use mod_info
    use mod_pbc
    implicit none

    integer                  i,ii
    integer                  id
    integer                  frame, l, m, n
    real(8)                  den,den_cm
    real(8)                  vol
    real(8)                  cos_alpha,cos_beta,cos_gamma
    real(8) ::               pi = 3.1415927
    real(8)                  lo,hi
    real(8)                  xlo_bound, xhi_bound, xy, xlo, xhi
    real(8)                  ylo_bound, yhi_bound, xz, ylo, yhi
    real(8)                  zlo_bound, zhi_bound, yz, zlo, zhi 
    real(8)                 alp, bet, gam, aa,ab,ac  
    character(5)            cty
    character(80)           fin   
    character(80)           line   
  
   

!      allocate memory 

    !allocate(n_1(typ(1)))
    !allocate(n_2(typ(2)))
    !allocate(n_3(typ(3)))
    !allocate(n_4(typ(4)))
    !allocate(n_5(typ(5)))

    print '(A23)','Input l,m,n: (la*mb*nc)'
    read(*,*)l,m,n

    print '(A23)','Input equlibrium steps:'
    read(*,*) frame

    open(59,file='lattice.txt',status='unknown')

    ii = 0

    do
       if(ico==1)then
          read(10,*,end=100)
          read(10,*,end=100)
          read(10,*,end=100)
          read(10,*,end=100)npart
          read(10,'(A80)',end=100)line
!
!
          if(index(line,'ITEM: BOX BOUNDS')/=0)then
           if(index(line,'xy xz yz')/=0)then      
                read(10,*,end=100)xlo_bound, xhi_bound, xy
                read(10,*)ylo_bound, yhi_bound, xz
                read(10,*)zlo_bound, zhi_bound, yz
                xlo= xlo_bound - min(0.0,xy,xz,xy+xz)
                xhi= xhi_bound - MAX(0.0,xy,xz,xy+xz)
                ylo= ylo_bound - MIN(0.0,yz)
                yhi= yhi_bound - MAX(0.0,yz)
                zlo= zlo_bound 
                zhi= zhi_bound
                a = (/xhi-xlo,0.0d0 ,0.0d0/)
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
               read(10,*,end=100)lo,hi
               a(1)= hi-lo
               r1= a(1)
               read(10,*)lo,hi
               b(2)= hi-lo
               r2= b(2)
               read(10,*)lo,hi
               c(3)= hi-lo
               r3= c(3)

               alpha = dacos(0.0d0)
               beta  = dacos(0.0d0)
               gamma = dacos(0.0d0)

           endif
          endif
       elseif(ico==2)then
          read(10,*,end=100)npart
          open(11,file='pbc.txt',status='unknown')
          read(11,*)(a(i),i=1,3)
          read(11,*)(b(i),i=1,3)
          read(11,*)(c(i),i=1,3)
          close(11)

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
       endif
 
       read(10,*,end=100) 
!
!       print '(A70)','***********************************************************************'!华丽的分割线
!       print '(2X,A40)','Cell parameters (Angstroms/Degrees):'
!       print '(A70)','***********************************************************************'
!       print '(2(A8,F12.4))','a=',r1,'alpha=',alpha*180.0d0/pi
!       print '(2(A8,F12.4))','b=',r2,'beta =',beta*180.0d0/pi
!       print '(2(A8,F12.4))','c=',r3,'gamma=',gamma*180.0d0/pi


       do i=1,npart
          if    (ico==1)then
                read(10,*,end=100)id,types(id)!,(x(l,i),l=1,3)!,q,(v(m,i),m=1,3)
                mass_all = mass_all + mass(types(id))
                typ(types(id))=typ(types(id))+1
          elseif(ico==2)then
                read(10,*,end=100)cty
                !cty= trim(adjustl(cty))    
!                print *,cty            
                !if(cty(1:2)==trim(adjustl(ctype(1))))types(i)=1
                !if(cty(1:2)==trim(adjustl(ctype(2))))types(i)=2
                !if(cty(1:2)==trim(adjustl(ctype(3))))types(i)=3
                !if(cty(1:2)==trim(adjustl(ctype(4))))types(i)=4
                !if(cty(1:2)==trim(adjustl(ctype(5))))types(i)=5
                !mass_all = mass_all + mass(types(i))
                !typ(types(i))=typ(types(i))+1
          endif
       enddo



       vol   = c(1)*b(2)*a(3) - c(1)*b(3)*a(2) + c(2)*b(3)*a(1) - c(2)*b(1)*a(3) + c(3)*a(2)*b(1)&
               - c(3)*b(2)*a(1)

       den   = mass_all/abs(vol)
       den_cm= den*10.0000/6.02253
       write(59,*) ii,r1/l,r2/m,r3/n,alpha*180.0d0/pi,beta*180.0d0/pi,gamma*180.0d0/pi,den_cm

       !print '(A70)','***********************************************************************'
       !print '(A32,I6,31X)','Totally read in partical number:',npart
       !do i=1,4
          !print '(A15,1X,I5,A3,1X,A9,35X)','Totally read in:',typ(i),ctype(i),'atoms...'
       !enddo

       !print '(A1,1X,A37,30X)','*','The total mass in the simulation box:'
       !print '(A1,1X,F20.8,2X,A6,39X)','*',mass_all,'g/mole'
       !print '(A1,1X,A31,36X)','*','The total density of the system:'
       !print '(A1,1X,F20.8,2X,A10,35X)','*',den,'g/mole/A^3'
       !print '(A1,1X,F20.8,2X,A10,35X)','*',den_cm,'g/cm^3'
       !print '(A70)','***********************************************************************'

       ii = ii + 1
    enddo

100 continue
    close(59)
    !close(58)

    return
end
