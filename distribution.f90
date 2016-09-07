!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================


    subroutine  distribution
    implicit none

    integer             i
    integer             nstep
    integer             iline
    integer             lent(30)
    real(8)                ctime,ctime_old
!    character(30)       file_name
    character(80)       line

    integer             n_species(30),n_species_old(30),n_species_old_old(30)


    real(8)             life_time

    character(12)        spec_name(30)
    character(22)       spec_index(30)
    character(22)       spec_id

!    print '(A20)','Input in file name:'
!    read(*,*)file_name

    open(110,file='molecular_structure.txt',status='old')
    open(120,file='distribution.txt',status='unknown')


    nstep = 0
    iline = 0
    n_species(:) = 0
    n_species_old(:) = 0    
    ctime = 0.0


    data spec_name(1:22) /'CH3O2N','H','HO','H2O','O',&
                          'ON','O2N','H2N','H3N','CO',&
                          'CO2','CH2O','C3O2','CH2O2N','CH4O2N',&
                          'CH3','CH4','CON','HON','C9O7',&
                          'C4O5','C6O3N7'/

    write(120,'(A10,22A6)')'Time      ',(spec_name(i),i=1,22)


    do i=1,22
       spec_name(i) = adjustl(spec_name(i))
       lent(i)= len_trim(spec_name(i))
       spec_id = '                      '
       spec_id(1:lent(i))= spec_name(i)
       spec_id(22:22)= 'x'
       spec_index(i)=spec_id
    enddo

!    open(12,file='debug.txt',status='unknown')
!    do i=1,22
!       write(12,'(A22)')spec_index(i)
!    enddo
!    close(12)
    
    do
       read(110,'(A80)',end=999)line
       iline = iline+1
       if    (index(line,'Time(fs)')/=0)then
               ctime_old = ctime
               read(line,'(8X,F20.4)')ctime
               nstep = nstep + 1
               if(nstep>1)then
                    write(120,'(F10.2,22I6)')ctime_old,(n_species(i),i=1,22)
                    n_species(:) = 0
               endif
       elseif(index(line,spec_index(1))/=0)then
               read(line,'(23X,I8)')n_species(1)
       elseif(index(line,spec_index(2))/=0)then
               read(line,'(23X,I8)')n_species(2)
       elseif(index(line,spec_index(3))/=0)then
               read(line,'(23X,I8)')n_species(3)
       elseif(index(line,spec_index(4))/=0)then
               read(line,'(23X,I8)')n_species(4)
       elseif(index(line,spec_index(5))/=0)then
               read(line,'(23X,I8)')n_species(5) 
       elseif(index(line,spec_index(6))/=0)then
               read(line,'(23X,I8)')n_species(6)
       elseif(index(line,spec_index(7))/=0)then
               read(line,'(23X,I8)')n_species(7) 
       elseif(index(line,spec_index(8))/=0)then
               read(line,'(23X,I8)')n_species(8)
       elseif(index(line,spec_index(9))/=0)then
               read(line,'(23X,I8)')n_species(9)
       elseif(index(line,spec_index(10))/=0)then
               read(line,'(23X,I8)')n_species(10)    
       elseif(index(line,spec_index(11))/=0)then
               read(line,'(23X,I8)')n_species(11)                
       elseif(index(line,spec_index(12))/=0)then
               read(line,'(23X,I8)')n_species(12) 
       elseif(index(line,spec_index(13))/=0)then
               read(line,'(23X,I8)')n_species(13)  
       elseif(index(line,spec_index(14))/=0)then
               read(line,'(23X,I8)')n_species(14) 
       elseif(index(line,spec_index(15))/=0)then
               read(line,'(23X,I8)')n_species(15)                 
       elseif(index(line,spec_index(16))/=0)then
               read(line,'(23X,I8)')n_species(16)   
       elseif(index(line,spec_index(17))/=0)then
               read(line,'(23X,I8)')n_species(17)                 
       elseif(index(line,spec_index(18))/=0)then
               read(line,'(23X,I8)')n_species(18)
       elseif(index(line,spec_index(19))/=0)then
               read(line,'(23X,I8)')n_species(19)  
       elseif(index(line,spec_index(20))/=0)then
               read(line,'(23X,I8)')n_species(20)  
       elseif(index(line,spec_index(21))/=0)then
               read(line,'(23X,I8)')n_species(21)  
       elseif(index(line,spec_index(22))/=0)then
               read(line,'(23X,I8)')n_species(22)     
       endif
    enddo
999 continue
    write(120,'(F10.2,22I6)')ctime,(n_species(i),i=1,22)


    close(110)
    close(120)
    
    return
    end  


