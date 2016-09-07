!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@foxmail.com                      =
!=========================================================================


    subroutine products
    implicit none

    real                   ctime
    integer                i
    integer                nmol
    integer                imol
    integer                iline
    integer            ::  new= 1
    integer,parameter  ::  max_mol = 5000
    character(80)          line
    character(20)          t1,t2
    character(20)          prod(max_mol)
!    character(40)          file_name
    logical                equal

!    print '(A20)','Input in file name:'
!    read(*,*)file_name

    open(10,file='molecular_structure.txt',status='old')
    open(20,file='product.txt',status='unknown')

  
    iline = 0
    imol  = 0
    do
       read(10,'(A80)',end=999)line
       iline = iline+1
       if(index(line,'Time(fs)')/=0)then
               read(line,'(8X,F20.4)')ctime
!               print *,t1,ctime
                
       elseif(index(line,'Molecular structure -vs- number')/=0)then
               read(line,*)t1
       else
               read(line,'(A20,1X,A1,1X,I8)')t1,t2,nmol
!               print *,t1
               if(imol==0)then
                       prod(1)= adjustl(t1)
                       imol=1
                       write(20,'(A20,f10.2)')prod(1),ctime
               else
                       new= 1
                       do i=1,imol
                          t1= trim(adjustl(t1))
                          t2= trim(prod(i))
                          if(t1==t2)then
!                                 print *,t1,t2
                                  new=0
                                  exit
                          endif
                       enddo
                       if(new==1)then
                               imol= imol + 1
                               if(imol>max_mol)print '(A16)','Max_mol overfow!'
                               prod(imol)=trim(adjustl(t1))
!                               print *,prod(imol)
                               write(20,'(A20,f10.2)')prod(imol),ctime
                       endif
               endif
       endif        

    enddo
999 continue

    stop
    end


