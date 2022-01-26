!********************************************************************!
!********************************************************************!
!                                                                    !
!   precision -- A module to set the precision of variables          !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 30Oct17                - raw54 !
!                    Renamed parameters: 21Nov20             - raw54 !
!                    Variables added: 21Nov20                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!
!
!>  Doxygen Section:
!
!>  @author
!>  Rob Watson
!>
!>  @brief 
!>  This module sets up the precision of the calculation and some
!>  useful parameters to avoid having to define them regularly.
!>
!*******************************************************************
!*******************************************************************

module precision

   implicit none

   ! Declare variables

   integer, parameter :: WI = kind(0)
   integer, parameter :: WP = kind(0d0)
   integer, parameter :: CP = 80

   real(kind=WP), parameter :: zero = 0.0_wp
   real(kind=WP), parameter :: half = 0.5_wp
   real(kind=WP), parameter :: one  = 1.0_wp

   real(kind=WP), parameter :: pi  = 4.0_wp * atan(one)

end module precision

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module precision                                          !
!                                                                    !
!********************************************************************!
!********************************************************************!

