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

    subroutine component
    use mode     ,only: npart,nstep,mlist,molelist,types,ctype,n_mol,dt,ico,istep
    implicit none

    integer                     rflag,i,j,jj,l,m

    integer                     n_str
    integer,allocatable::       str_num(:)  
!   e.g. n.u. of CHON <==> mol_table(2,2,2,2) because 
!                NO can be mol_table(1,1,2,2)
    real(8)                     xold(3)   ! for Larger moleculars split by PBC
    character(20)            :: mol_name, ml
    character(4)                str1
    integer,allocatable      :: str_name(:,:)
!   character(20),allocatable:: str_name(:,:)

    integer                     leng,lenn,lens
    integer                     nty(5)

    integer                     nof1,nof2,nof3,nof4,nof5
    integer                     c_cluster, other


    open(36,file='mol_str.txt',status='unknown')
!    allocate(mol_tab(nty1,nty2,nty3,nty4))   ! memory overflow!!!

    allocate(str_num(npart))
    allocate(str_name(5,npart))

!   ctype(1-4) = C, H, O, N
!   str_name 表示： 
!   分子 i 的分子式为 C_str_name(1:i)H_str_name(2:i)N_str_name(3:i)O_str_name(4:i)Al_str_name(5:i)


    rflag=0
    nstep=0

    do   ! main cycle loop

       nstep=nstep+1

       if    (ico==1)then 
             call read_tra(rflag)
       elseif(ico==2)then
             call read_xyz(rflag)
       endif
       if(rflag==0) print *,'* The last frame!'
       if(nstep>istep)  exit
       if(nstep/=istep .and. rflag /=0) cycle
!
       call neighbor
       call findmole
!
!
       str_num(:) = 0
       str_name(:,:)= 0
       n_str = 1

       do i=1, n_mol    ! cycle from molecules
          if(i>1)then
             jj=mlist(i-1) + 1
          else
             jj=1
          endif

          nof1= 0
          nof2= 0
          nof3= 0
          nof4= 0
          nof5= 0

          do j=jj, mlist(i)
             if        (types(molelist(j))==1)then
                   nof1=nof1 + 1
             elseif(types(molelist(j))==2)then
                   nof2=nof2 + 1
             elseif(types(molelist(j))==3)then
                   nof3=nof3 + 1
             elseif(types(molelist(j))==4)then
                   nof4=nof4 + 1
             elseif(types(molelist(j))==5)then
                   nof5=nof5 + 1
             endif
          enddo

!  molecular structure
!          mol_tab(n_1(i)+1,n_2(i)+1,n_3(i)+1,n_4(i)+1) = mol_tab(n_1(i)+1,n_2(i)+1,n_3(i)+1,n_4(i)+1) + 1
          
          if(i>1)then
                  do l=1,n_str
!                     read(str_name(l),'(4I5)')(nty(j),j=1,4) 
                     nty(:)= str_name(:,l)
                     if(nof1==nty(1).and.nof2==nty(2).and.nof3==nty(3).and.nof4==nty(4).and.nof5==nty(5))then
                             str_num(l)= str_num(l) + 1
                             goto 55
                     endif
                  enddo
!         new structure
                  n_str= n_str + 1
                  str_num(n_str)= 1
                  str_name(1,n_str)= nof1
                  str_name(2,n_str)= nof2
                  str_name(3,n_str)= nof3
                  str_name(4,n_str)= nof4
                  str_name(5,n_str)= nof5   
!                 write(str_name(n_str),'(4I5)')n_1(i),n_2(i),n_3(i),n_4(i)
          else
!                 write(str_name(n_str),'(4I5)')n_1(i),n_2(i),n_3(i),n_4(i)
                  str_name(1,n_str)= nof1
                  str_name(2,n_str)= nof2
                  str_name(3,n_str)= nof3
                  str_name(4,n_str)= nof4  
                  str_name(5,n_str)= nof5   
                  str_num(1)=1
          endif        
55        continue

       enddo
 
       write(35,'(A10,F20.4)')'Time(ps): ',nstep*dt
       write(35,'(A31)')'Molecular structure -vs- number'
!   ~~~~~~~~
       do l=1,n_str  
!         print *, str_name(:,l),l,nstep
!         read(str_name(l),'(4I5)')(nty(j),j=1,4)
          nty(:)= str_name(:,l)
          mol_name='                    '
          if (str_name(1,l)>=5) cycle   ! filter out the larger cluster ~~
          do j=1,5
             if(nty(j)>0)then
                  leng = len_trim(mol_name)
                  lenn = len_trim(ctype(j)) + leng
                  mol_name(leng+1:lenn)=ctype(j)
                  if(nty(j)>1)then
                         write(str1,'(I4)')nty(j)
                         str1=adjustl(str1)
                         lens=len_trim(str1)+lenn
                         mol_name(lenn+1:lens)=str1
                  endif
             endif
          enddo  
          write(35,'(A20,1X,A1,1X,I8)')mol_name,' ',str_num(l)      
       enddo
!   ~~~~~~~~
       do l=1,n_str
!         print *, str_name(:,l),l,nstep
!         read(str_name(l),'(4I5)')(nty(j),j=1,4)
          nty(:)= str_name(:,l)
          mol_name='                    '
          if (str_name(1,l)<=4) cycle  ! larger cluster ~~
          do j=1,5
             if(nty(j)>0)then
                  leng = len_trim(mol_name)
                  lenn = len_trim(ctype(j)) + leng
                  mol_name(leng+1:lenn)=ctype(j)
                  if(nty(j)>1)then
                         write(str1,'(I4)')nty(j)
                         str1=adjustl(str1)
                         lens=len_trim(str1)+lenn
                         mol_name(lenn+1:lens)=str1
                  endif
             endif
          enddo  
          write(35,'(A20,1X,A1,1X,I8)')mol_name,' ',str_num(l)      
       enddo

!   ~~~~~~~~


!   ~~~~~~~~
       do l=1,n_str  
!         print *, str_name(:,l),l,nstep
!         read(str_name(l),'(4I5)')(nty(j),j=1,4)
          nty(:)= str_name(:,l)
          mol_name='                    '
          if (str_name(1,l)<5.and.str_num(l)>12) then   ! filter out the larger cluster ~~
            do j=1,5
               if(nty(j)>0)then
                  leng = len_trim(mol_name)
                  lenn = len_trim(ctype(j)) + leng
                  mol_name(leng+1:lenn)=ctype(j)
                  if(nty(j)>1)then
                         write(str1,'(I4)')nty(j)
                         str1=adjustl(str1)
                         lens=len_trim(str1)+lenn
                         mol_name(lenn+1:lens)=str1
                  endif
               endif
            enddo  
            write(36,'(A20,1X,A1,1X,I8)')mol_name,' ',str_num(l) 
          endif     
       enddo

!   ~~~~~~~~
       c_cluster = 0
       do l=1,n_str
          nty(:)= str_name(:,l)
          mol_name='                    '
          if (str_name(1,l)<=4) cycle  ! larger cluster ~~
          c_cluster =  c_cluster + str_num(l)    
       enddo
       write(36,'(A20,1X,A1,1X,I8)')'C cluster (C>=5)',' ',c_cluster  
!   ~~~~~~~~
       other = 0
       do l=1,n_str
          nty(:)= str_name(:,l)
          mol_name='                    '
          if (str_name(1,l)<=4.and.str_num(l)<=12) then  ! other (num < 12) ~~
              other =  other + str_num(l)      
          endif
       enddo
       write(36,'(A20,1X,A1,1X,I8)')'Others)',' ',other
!   ~~~~~~~~

       if(rflag==0)then
         print *,'* exit now!'
         exit
       endif
       
    enddo     !  main cyle loop


    deallocate(str_name)
    deallocate(str_num)

    close(35)
    close(36)

    return
    end
