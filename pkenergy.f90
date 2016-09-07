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
!
     subroutine pkenergy(br2,i,j,ty)
     use mode
     use mod_pke
     use mod_uspl
     implicit none
!
!    Local variables
!
     integer i,j
     integer ty
     real(8) br2,rr,de_en,ki_en
     real(8) uconver
     real(8), external:: uspl

!    br2 represent bond radius^2
     rr=sqrt(br2)

     !if(rr>r0(types(i),types(j)))then
                de_en= uspl(ty,rr)-d0(ty)   !    bond energy 
     !else
     !           de_en= de(ty)-d0(ty)
     !endif
!
!        
!        Potential energy well
!
!        Compute Kinetic energy 
!
!        uconver=2.3901D+3
        uconver=0.5*1.0D+4/4.184 
!
!    
        ki_en=(v(1,i)-v(1,j))**2 + (v(2,i)-v(2,j))**2 + (v(3,i)-v(3,j))**2
        ki_en=0.5*(mass(types(i))+mass(types(j)))*(ki_en)*uconver
!       
!   Unit convertion 
!   1m/s=1*10E-5A/fs
!   Atomic number to Kg
!   (unified)atomic mass unit
!   1amu = 1.66053886 E -27 Kg
!   unconver= ==> m/s ==> Kg/mol ==> J/mol ==> Kcal/mol 


!    Comparation between DE and Ki Energy 
!
     pkflag=de_en*23.0605 + ki_en*1.0e-6
     !print *,de_en*23.0605,ki_en*1.0e-6,rr

     return
     end

function uspl(ty,rd)
    use mod_uspl
    implicit none

    integer    ty,ind
    real(8)       rd
    real(8)       rr
    real(8)       uspl

    if (ty==2)then 
              if  (rd<=cnst(1))then
                   ind = 1
              elseif(rd>=cned(cns))then
                   ind = cns
              else
                   ind = int((rd-cnst(1))/(cned(1)-cnst(1)))+1
              endif
              !if (rd<=cnst(1).or.rd>=cned(cns)) print '(A36)','Warning: Radius out of spline range!'
              rr = rd-cnst(ind)
              uspl = cn0(ind) + cn1(ind)*rr + cn2(ind)*rr*rr + cn3(ind)*rr*rr*rr
              !print *,rd,ecn
     elseif(ty==1)then
              if  (rd<=chst(1))then
                   ind = 1
              elseif(rd>=ched(chs))then
                   ind = chs
              else
                   ind = int((rd-chst(1))/(ched(1)-chst(1)))+1
              endif
              !if (rd<=chst(1).or.rd>=ched(chs)) print '(A36)','Warning: Radius out of spline range!'
              rr = rd-chst(ind) 
              uspl = ch0(ind) + ch1(ind)*rr + ch2(ind)*rr*rr + ch3(ind)*rr*rr*rr
              !print *,rd,ecn
    endif
    
    return
end

