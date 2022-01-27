!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_setIBM -- This subroutine sets the IBM locations             !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Module created: 27Jan22                 - raw54 !
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
!>  This subroutine turns on and off the IBM in different places,
!>  depending on either a wall distance real, a matlab TIFF, or a
!>  simple geometrical set
!>
!*******************************************************************
!*******************************************************************

module mod_setIBM

   ! Declare modules

   use precision
   use mod_globalParameters
   use mod_ReadAlloc

   ! Turn off implicit typing

   implicit none

   ! Set everything private

   private

   ! Turn on what should be publicly visible

   public :: ibm_setSolid

contains

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine ibm_setSolid()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Set up the IBM solid locations based on a simple geometrical if
   !>  test
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine ibm_setSolid()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      real   (kind=WP) :: r

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         if ( sqrt((xc(1,i,j,k) - 0.030d0)**2 + xc(2,i,j,k)**2) .lt. 0.0025 ) then
            iS(1,i,j,k) = 1
           
         else

            iS(1,i,j,k) = 0

         end if

      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine ibm_setSolid

   !*******************************************************************
   !*******************************************************************


end module mod_setIBM

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_setIBM                                         !
!                                                                    !
!********************************************************************!
!********************************************************************!

