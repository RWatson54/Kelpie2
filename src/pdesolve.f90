!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_PDESolve -- A module containing PDE solution routines        !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Module created: 25Jan22                 - raw54 !
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
!>  This module contains general routines for solving a PDE.
!>  Ultimately, we want a single "getResidual" which just gets how
!>  far from A(q) q = 0 we are. This should later let us tie in to
!>  better solvers (JFNK, or something like that?)
!>
!*******************************************************************
!*******************************************************************

module mod_PDESolve

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

   public :: pde_getResiduals

contains

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_getResiduals()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  This subroutine calculates A(u) u = R(u). This residual vector,
   !>  R(u) can then be used to solve the system for u, using, e.g. 
   !>  JFNK, or, far more easily but probably less efficiently,
   !>  time marching. Have a guess at the one we use here...
   !>
   !*******************************************************************
   !*******************************************************************

   subroutine pde_getResiduals()

      implicit none

      ! Compute the primitive variables from the conserved variables

      call pde_getPrimitive()

      ! Zero the residuals

      call pde_zeroResiduals()

      ! Get the gradients at each of the cells

      call pde_getGradients()

      ! Accumulate the fluxes

      call pde_sumFluxes()

      ! Accumulate any volume sources, including DTS, IBM

      call pde_applySources()

      ! Gather all of the terms together into a single scaled residual

      call pde_sumResidual()

      ! Return to calling program

      return

   end subroutine pde_getResiduals

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_sumResidual()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Gathers the residuals into a convenient single vector
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_sumResidual()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      ! Gather the residuals together

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         res(:,i,j,k) = ires(:,i,j,k) + vres(:,i,j,k) + sts(:,i,j,k)
      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_sumResidual

   !*******************************************************************
   !*******************************************************************


   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_zeroResiduals()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Zeros the residuals at the start of the loop to gather the flow
   !>  into.
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_zeroResiduals()

      implicit none

      ! Zero the residual arrays

      ires = zero
      vres = zero

      ! Return to calling program

      return

   end subroutine pde_zeroResiduals

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_getConserved()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Gets the conserved variables from the primitive variables
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_getConserved()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      ! Loop through, getting the conserved variables from the primitive

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         u(:,i,j,k) = q2u(q(:,i,j,k))
      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_getConserved

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_getPrimitive()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Gets the primitive variables from the conserved variables
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_getPrimitive()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      ! Loop through, getting the conserved variables from the primitive

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         q(:,i,j,k) = u2q(u(:,i,j,k))
      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_getPrimitive

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_sumFluxes()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Accumulate the fluxes into the residuals
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_sumFluxes()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      real   (kind=WP) :: fni(nPdes), fnv(nPdes)

      ! Loop over all the cells, accumulating i->i+1, j->j+1, k->k+1 faces

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc - 1

         ! First, in the i direction
         fni = iflux(q(:,i,j,k), q(:,i+1,j,k), nwt(:,2,i,j,k))
         fnv = vflux(q(:,i,j,k), q(:,i+1,j,k), qp(:,:,i,j,k), qp(:,:,i+1,j,k), nwt(:,2,i,j,k))

         ! And accumulate
         ires(:,i  ,j,k) = ires(:,i  ,j,k) - fni
         ires(:,i+1,j,k) = ires(:,i+1,j,k) + fni
         vres(:,i  ,j,k) = vres(:,i  ,j,k) - fnv
         vres(:,i+1,j,k) = vres(:,i+1,j,k) + fnv

      end do
      end do
      end do

      do k = 1, nKc
      do j = 1, nJc - 1
      do i = 1, nIc         

         ! Next, in the j direction
         fni = iflux(q(:,i,j,k), q(:,i,j+1,k), nwt(:,4,i,j,k))
         fnv = vflux(q(:,i,j,k), q(:,i,j+1,k), qp(:,:,i,j,k), qp(:,:,i,j+1,k), nwt(:,4,i,j,k))

         ! And accumulate
         ires(:,i,j  ,k) = ires(:,i,j  ,k) - fni
         ires(:,i,j+1,k) = ires(:,i,j+1,k) + fni
         vres(:,i,j  ,k) = vres(:,i,j  ,k) - fnv
         vres(:,i,j+1,k) = vres(:,i,j+1,k) + fnv

      end do
      end do
      end do

      do k = 1, nKc - 1
      do j = 1, nJc
      do i = 1, nIc         

         ! Finally, in the k direction
         fni = iflux(q(:,i,j,k), q(:,i,j,k+1), nwt(:,6,i,j,k))
         fnv = vflux(q(:,i,j,k), q(:,i,j,k+1), qp(:,:,i,j,k), qp(:,:,i,j,k+1), nwt(:,6,i,j,k))

         ! And accumulate
         ires(:,i,j,k  ) = ires(:,i,j,k  ) - fni
         ires(:,i,j,k+1) = ires(:,i,j,k+1) + fni
         vres(:,i,j,k  ) = vres(:,i,j,k  ) - fnv
         vres(:,i,j,k+1) = vres(:,i,j,k+1) + fnv

      end do
      end do
      end do

      ! Now loop over the boundaries

      ! i = 1, inflow
      ! i = nI, outflow
      ! j = 1, invwall
      ! j = nJ, invwall
      ! k = 1, invwall
      ! k = nK, invwall

      ! i = 1, inflow surface

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, 1

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), inflow(q(:,i,j,k), nwt(:,1,i,j,k)), nwt(:,1,i,j,k))
         fnv = vflux(q(:,i,j,k), inflow(q(:,i,j,k), nwt(:,1,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,1,i,j,k))
         
         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
                 
      end do
      end do
      end do

      ! i = nIc, outflow surface

      do k = 1, nKc
      do j = 1, nJc
      do i = nIc, nIc

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), outflow(q(:,i,j,k), nwt(:,2,i,j,k)), nwt(:,2,i,j,k))
         fnv = vflux(q(:,i,j,k), outflow(q(:,i,j,k), nwt(:,2,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,2,i,j,k))
         
         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
         
      end do
      end do
      end do

      ! j = 1, invwall surface

      do k = 1, nKc
      do j = 1, 1
      do i = 1, nIc

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,3,i,j,k)), nwt(:,3,i,j,k))
         fnv = vflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,3,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,3,i,j,k))
         
         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
         
      end do
      end do
      end do

      ! j = nJc, invwall surface

      do k = 1, nKc
      do j = nJc, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,4,i,j,k)), nwt(:,4,i,j,k))
         fnv = vflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,4,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,4,i,j,k))
         
         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
         
      end do
      end do
      end do

      ! k = 1, invwall surface

      do k = 1, 1
      do j = 1, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,5,i,j,k)), nwt(:,5,i,j,k))
         fnv = vflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,5,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,5,i,j,k))

         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
         
      end do
      end do
      end do

      ! k = nKc, invwall surface

      do k = nKc, nKc
      do j = 1, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fni = iflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,6,i,j,k)), nwt(:,6,i,j,k))
         fnv = vflux(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,6,i,j,k)), qp(:,:,i,j,k), qp(:,:,i,j,k), nwt(:,6,i,j,k))
         
         ! And accumulate
         ires(:,i,j,k) = ires(:,i,j,k) - fni
         vres(:,i,j,k) = vres(:,i,j,k) - fnv
         
      end do
      end do
      end do

      ! Rescale residuals by the volume

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         ires(:,i,j,k) = ires(:,i,j,k) / vol(1,i,j,k)
         vres(:,i,j,k) = vres(:,i,j,k) / vol(1,i,j,k)

      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_sumFluxes

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_getGradients()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Calculate the gradients needed for the viscous calculations
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_getGradients()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      real   (kind=WP) :: fN(nPdes,3)
      
      ! Zero qp

      qp = zero

      ! Loop over the internal faces, gathering the gradient contributions

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc - 1

         ! First, in the i direction
         fN = grad(q(:,i,j,k), q(:,i+1,j,k), nwt(:,2,i,j,k))

         ! And accumulate
         qp(:,:,i  ,j,k) = qp(:,:,i  ,j,k) + fN
         qp(:,:,i+1,j,k) = qp(:,:,i+1,j,k) - fN

      end do
      end do
      end do

      do k = 1, nKc
      do j = 1, nJc - 1
      do i = 1, nIc

         ! Next, in the j direction
         fN = grad(q(:,i,j,k), q(:,i,j+1,k), nwt(:,4,i,j,k))

         ! And accumulate
         qp(:,:,i,j  ,k) = qp(:,:,i,j  ,k) + fN
         qp(:,:,i,j+1,k) = qp(:,:,i,j+1,k) - fN

      end do
      end do
      end do

      do k = 1, nKc - 1
      do j = 1, nJc
      do i = 1, nIc

         ! Finally, in the k direction
         fN = grad(q(:,i,j,k), q(:,i,j,k+1), nwt(:,6,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k  ) = qp(:,:,i,j,k  ) + fN
         qp(:,:,i,j,k+1) = qp(:,:,i,j,k+1) - fN

      end do
      end do
      end do


      ! Now loop over the boundaries

      ! i = 1

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, 1

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), inflow(q(:,i,j,k), nwt(:,1,i,j,k)), nwt(:,1,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN
                 
      end do
      end do
      end do

      ! i = nIc, outflow surface

      do k = 1, nKc
      do j = 1, nJc
      do i = nIc, nIc

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), outflow(q(:,i,j,k), nwt(:,2,i,j,k)), nwt(:,2,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN

      end do
      end do
      end do

      ! j = 1, invwall surface

      do k = 1, nKc
      do j = 1, 1
      do i = 1, nIc

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,3,i,j,k)), nwt(:,3,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN
         
      end do
      end do
      end do

      ! j = nJc, invwall surface

      do k = 1, nKc
      do j = nJc, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,4,i,j,k)), nwt(:,4,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN
         
      end do
      end do
      end do

      ! k = 1, invwall surface

      do k = 1, 1
      do j = 1, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,5,i,j,k)), nwt(:,5,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN
         
      end do
      end do
      end do

      ! k = nKc, invwall surface

      do k = nKc, nKc
      do j = 1, nJc
      do i = 1, nIc

         ! Compute the boundary fluxes
         fN = grad(q(:,i,j,k), invwall(q(:,i,j,k), nwt(:,6,i,j,k)), nwt(:,6,i,j,k))

         ! And accumulate
         qp(:,:,i,j,k) = qp(:,:,i,j,k) + fN
         
      end do
      end do
      end do

      ! Rescale gradients by the volume

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         qp(:,:,i,j,k) = qp(:,:,i,j,k) / vol(1,i,j,k)

      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_getGradients

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine pde_applySources()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  This subroutine applies any sources to the equations, including
   !>  the IBM momentum sources, and any DTS unsteady source terms
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine pde_applySources()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k

      ! Loop over the elements, applying the source

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         sts(:,i,j,k) = source(q(:,i,j,k), iS(1,i,j,k), e(:,i,j,k))
      
      end do
      end do
      end do

      ! Return to calling program

      return

   end subroutine pde_applySources

   !*******************************************************************
   !*******************************************************************


end module mod_PDESolve

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_PDESolve                                       !
!                                                                    !
!********************************************************************!
!********************************************************************!

