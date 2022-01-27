!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_update -- Update the vector u based on residual r            !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Module created: 22Jan22                 - raw54 !
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
!>  This module performs some kind of update to the conserved 
!>  variable vector, based on the residual, r.
!>  
!>
!*******************************************************************
!*******************************************************************

module mod_update

   ! Declare modules

   use precision
   use mod_globalParameters
   use mod_ReadAlloc
   use mod_Equations
   
   ! Turn off implicit typing

   implicit none

   ! Set everything private

   private

   ! Turn on what should be publicly visible

   public :: update_explicit

contains

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine update_explicit()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Perform an explicit update step to determine the value of u at
   !>  some future distance. Determine the time step to achieve this
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine update_explicit()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k, iPde

      real   (kind=WP) :: dx

      ! Set the maximum residual to a very large negative number

      maxRes = -1000000.0_wp

      ! First, get the time step needed to perform the update

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         dx = sqrt(vol(1,i,j,k))

         dt(:,i,j,k) = CFL * lambda(dx, q(:,i,j,k))

      end do
      end do
      end do

      ! Then perform the explicit update step

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         do iPde = 1, nPdes
         
            u(iPde,i,j,k) = u(iPde,i,j,k) + dt(iPde,i,j,k) * res(iPde,i,j,k)

            maxRes(iPde) = max(maxRes(iPde), abs(res(iPde,i,j,k)))

         end do

      end do
      end do
      end do
      
      ! Return to calling program

      return

   end subroutine update_explicit

   !*******************************************************************
   !*******************************************************************

end module mod_update

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_update                                         !
!                                                                    !
!********************************************************************!
!********************************************************************!

