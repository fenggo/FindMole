!=========================================================================
!                                                       WRITED BY GUO F.                                                          =
!                                                                                                                                                =
!                                      Institute of Atomic & Molecular Physics of Si-                               =
!                                                  chuan University of China.                                                  =
!                                                                                                                                                =
!                                                  email: gfeng.alan@gmail.com                                             =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!

    program main
    use mode
    use mod_mass
    implicit none    

    integer       ianal
    character(80) fin                      
!    character(80) fout

    print *,"* Input in file name:"       
    read(*,'(A80)')fin                   
!    print *,"Input out file name:"        
!    read(*,'(A80)')fout                    

    print *,'* Coordinate Format:'
    print *,'* 1. Lammmpstrj'
    print *,'* 2. xyz'
    read(*,*)ico
    print *,'* In put time interval between trajectory frames (ps):'
    read(*,*)dt
!    fin='tt.lammpstrj'
!    fout='tt.out'

    fin =trim(fin)                          
!    fout=trim(fout)

    open(10,file=fin,status='old')                                                   
    open(35,file='molecular_structure.txt',status='unknown')

    call setpara
    call info
    !call initialsp
!
    print *,'1. analysis MS among all the trajectorys'
    print *,'2. analysis MS from one of the trajectorys'
    print *,'3. analysis the Mean squared displacement of the trajectorys'
    print *,'4. compute the Mean lattice parameters of the system'
    read(*,*) ianal


    if(ianal==1) then
          call find_reac
          close(35)
!         call distribution
          call products
    elseif(ianal==2)then
          print *,'* please input the step number of the trajectorys to analysis:'
          read(*,*) istep
          call component
          close(35)
    elseif(ianal==3)then
          print *,'* please input the step number of the trajectorys to analysis:'
          read(*,*) istep
          print *,'* Mean squared displacement ...'
          call msd
    elseif(ianal==4)then
         call lattice
    endif
!
!
    call release
    close(10)



    stop
    end
!
