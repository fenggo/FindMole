!============================================================
!                    WRITED BY GUO F.                       =
!                                                           =
!       Institute of Atomic & Molecular Physics of Si-      =
!                chuan University of China.                 =
!                                                           =
!                email: gfengs@foxmail.com                  =
!============================================================ 
!
!------------------------------------------------------------ 
!
! -----------------------------------------------------------------------
! string and number routines

! -----------------------------------------------------------------------
! returns TRUE if str1 matches 1st chars of str2, FALSE otherwise
! also returns m = loc of next char in str2
! could test if next char is a space or tab

      logical function match(str1,str2,m)
      implicit none

      character*(*) str1,str2
      integer m

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      if (str1.eq.str2(1:len(str1))) match = .TRUE.

      return
      end
!
!    When Error Accurred!
! 
     subroutine err(str)
     implicit none

     character*(*) str

     Write(*,'(4x,A)')str

     stop
     end
