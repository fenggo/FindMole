
module mod_uspl
    implicit none

    ! for spline potential
    integer,save :: cns 
    integer,save :: chs 
    real(8) , allocatable::  cnst(:),cned(:)                ! parameters for CN pair
    real(8) , allocatable::  cn0(:),cn1(:),cn2(:),cn3(:)
    real(8) , allocatable::  chst(:),ched(:)                ! parameters for CH pair
    real(8) , allocatable::  ch0(:),ch1(:),ch2(:),ch3(:)
    real(8)                              d0(2), de(2)

    data d0  /-313.766517386, -309.947859867/ 
    data de  /-321.908320106, -321.919115819/

end module

! ------------------------------------------------------------------------------------------------------------------

    subroutine initialsp
    use mod_uspl
    implicit none

    integer i,j
    real (8)   c(6)
    !chs
    open(110,file='CH.sp',status='old')
    i = 1
    do 
          read(110,*,end=998)(c(j),j=1,6)
          i= i+1
    enddo
998  continue
    chs = i-1

    allocate(chst(chs),ched(chs))
    allocate(ch0(chs),ch1(chs),ch2(chs),ch3(chs))

    rewind(110)
    do i=1,chs
          read(110,*)chst(i),ched(i),ch0(i),ch1(i),ch2(i),ch3(i)
    enddo
    close(110)

    !cns
    open(111,file='CN.sp',status='old')
    i= 1
    do 
          read(111,*,end=998)(c(j),j=1,6)
          i= i+1
    enddo
999  continue
    cns = i-1

    allocate(cnst(cns),cned(cns))
    allocate(cn0(cns),cn1(cns),cn2(cns),cn3(cns))

    rewind(111)
    do i=1,cns
          read(111,*)cnst(i),cned(i),cn0(i),cn1(i),cn2(i),cn3(i)
    enddo
    close(111)

    return
    end
