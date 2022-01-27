!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_Equations -- Contains functions pertaining to the actual eqs !
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
!>  This subroutine contains the information pertaining to the
!>  equations which are being solved, including defining the PDEs
!>  
!>
!*******************************************************************
!*******************************************************************

module mod_Equations

   ! Declare modules

   use precision
   use mod_globalParameters

   ! Turn off implicit typing

   implicit none

   ! Declare module variables

   integer(kind=WI), parameter :: nPdes = 4

   real   (kind=WP), parameter :: gma = 1.4_wp, rgas = 287.4_wp

   real   (kind=WP), parameter :: ACM_beta = 2.00_wp, rho0 = 1.226_wp

   real   (kind=WP), parameter :: nuL = 0.00001_wp

   real   (kind=WP), parameter :: ic(nPdes) = (/ zero, one, one, one /)

   ! Set everything private

   private

   ! Turn on what should be publicly visible

   public :: nPdes, lambda, q2u, u2q, iflux, vflux, source
   public :: invWall, visWall, outFlow, inFlow, grad

contains

   !*******************************************************************!
   !*******************************************************************!
   !  Function lambda
   !> @brief Get the timestep stability condition
   !> @param[in]  dx - measure of cell size
   !> @param[in]  qL - vector of primitive variables
   !> @param[out] dt - approximate timestep
   !*******************************************************************!
   !*******************************************************************!

   function lambda(dx, qL) result(dt)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP) :: dt, dx, qL(npdes)

      ! Declare internal variables

      real   (kind=WP) :: lmb

      ! Set the timestep scaling

      lmb  = sqrt(qL(2)**2 + qL(3)**2 + qL(4)**2) + &
             sqrt(qL(2)**2 + qL(3)**2 + qL(4)**2 + ACM_beta)

      dt = min(dx*dx/nuL, dx/lmb)

      ! Return to calling subprogram 

      return 

   end function lambda

   !*******************************************************************!
   !*******************************************************************!
   !  Function u2q
   !> @brief Determine the primitive variables from the conserved
   !> @param[in]  uL - vector of conserved variables
   !> @param[out] qL - vector of primitive variables
   !*******************************************************************!
   !*******************************************************************!

   function u2q(uL) result(qL)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP) :: uL(npdes), qL(npdes)

      ! Set the primitive variables from the conserved

      qL(1) = uL(1)
      qL(2) = uL(2)
      qL(3) = uL(3)
      qL(4) = uL(4)

      ! Return to calling subprogram 

      return 

   end function u2q

   !*******************************************************************!
   !*******************************************************************!
   !  Function q2u
   !> @brief Determine the conserved variables from the primitive
   !> @param[in]  qL - vector of primitive variables
   !> @param[out] uL - vector of conserved variables
   !*******************************************************************!
   !*******************************************************************!

   function q2u(qL) result(uL)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP) :: qL(npdes), uL(npdes)

      ! Set the conserved variables from the conserved

      uL(1) = qL(1)
      uL(2) = qL(2)
      uL(3) = qL(3)
      uL(4) = qL(4)

      ! Return to calling subprogram 

      return 

   end function q2u

   !*******************************************************************!
   !*******************************************************************!
   !  Function iflux
   !> @brief Compute the inviscid intercell flux from the equations
   !> @param[in]  qL - the left side primitives
   !> @param[in]  qR - the right side primitives
   !> @param[in]  n  - the left-to-right pointing normals
   !> @param[out] Fn - the normal flux
   !*******************************************************************!
   !*******************************************************************!

   function iflux(qL, qR, n) result(Fn)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: qL(npdes), qR(npdes), n(3)

      real   (kind=WP) :: Fn(npdes)

      ! Declare local variables

      real   (kind=WP) :: nx, ny, nz, aa
      real   (kind=WP) :: vxL, vxR, vyL, vyR, vzL, vzR
      real   (kind=WP) :: pL, pR             
      real   (kind=WP) :: vnL, vnR
      real   (kind=WP) :: fL(nPdes), fR(nPdes)
      real   (kind=WP) :: cR, cL, smax

      ! Get flow area 

      aa = norm2(n)

      ! Normal vector

      nx = n(1) / aa
      ny = n(2) / aa
      nz = n(3) / aa

      ! Primitive and other variables

      !  Left state

      pL   = qL(1)
      vxL  = qL(2)
      vyL  = qL(3)
      vzL  = qL(4)
      vnL  = vxL*nx+vyL*ny+vzL*nz

      !  Right state

      pR   = qR(1)
      vxR  = qR(2)
      vyR  = qR(3)
      vzR  = qR(4)
      vnR  = vxR*nx+vyR*ny+vzR*nz

      ! Compute the max eigenvalues for the ACM Rusanov flux

      cL  = sqrt(vnL**2) + sqrt(ACM_beta + vnL**2)
      cR  = sqrt(vnR**2) + sqrt(ACM_beta + vnR**2)

      smax = max(cL,cR)

      ! Compute the NS and turbulent model fluxes - Multiplier for dissipation added 

      fL(1) = ACM_beta*vnL
      fL(2) = vnL * vxL + nx*pL/rho0
      fL(3) = vnL * vyL + ny*pL/rho0
      fL(4) = vnL * vzL + nz*pL/rho0

      fR(1) = ACM_beta*vnR
      fR(2) = vnR * vxR + nx*pR/rho0
      fR(3) = vnR * vyR + ny*pR/rho0
      fR(4) = vnR * vzR + nz*pR/rho0

      Fn(1:4) = aa * half * (fL(1:4) + fR(1:4) - 0.50_wp * smax * (qR(1:4) - qL(1:4)))

      ! Return to calling subprogram 

      return 

   end function iflux

   !*******************************************************************!
   !*******************************************************************!
   !  Function vflux
   !> @brief Compute the viscous intercell flux from the equations
   !> @param[in]  qL - the left side primitives
   !> @param[in]  qR - the right side primitives
   !> @param[in]  n  - the left-to-right pointing normals
   !> @param[out] Fn - the normal flux
   !*******************************************************************!
   !*******************************************************************!

   function vflux(qL, qR, qpL, qpR, n) result(Fn)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: qL(npdes), qR(npdes), qpL(npdes,3), qpR(npdes,3), n(3)

      real   (kind=WP) :: Fn(npdes)

      ! Declare local variables

      real   (kind=WP) :: q(npdes), qP(npdes,3)
      real   (kind=WP) :: aa, nx, ny, nz, tau(3,3)
      real   (kind=WP) :: sigxn, sigyn, sigzn

      ! Get flow area 

      aa = norm2(n)

      ! Normal vector

      nx = n(1) / aa
      ny = n(2) / aa
      nz = n(3) / aa

      ! Get the centred variables and gradients (could do some fancy Bassi-Rebay stuff here)

      q  = half * ( qL  + qR )
      qp = half * ( qpL + qpR)

      ! Compute the viscous stresses - ignore compressibility correction

      tau(1,1) = 2.0_wp * nuL * ( half * (qp(2,1) + qp(2,1)))
      tau(1,2) = 2.0_wp * nuL * ( half * (qp(3,1) + qp(2,2)))
      tau(1,3) = 2.0_wp * nuL * ( half * (qp(4,1) + qp(2,3)))

      tau(2,1) = 2.0_wp * nuL * ( half * (qp(2,2) + qp(3,1)))
      tau(2,2) = 2.0_wp * nuL * ( half * (qp(3,2) + qp(3,2)))
      tau(2,3) = 2.0_wp * nuL * ( half * (qp(4,2) + qp(3,3)))

      tau(3,1) = 2.0_wp * nuL * ( half * (qp(2,3) + qp(4,1)))
      tau(3,2) = 2.0_wp * nuL * ( half * (qp(3,3) + qp(4,2)))
      tau(3,3) = 2.0_wp * nuL * ( half * (qp(4,3) + qp(4,3)))

      ! Compute the normal components

      sigxn = tau(1,1)*nx + tau(1,2)*ny + tau(1,3)*nz
      sigyn = tau(2,1)*nx + tau(2,2)*ny + tau(2,3)*nz
      sigzn = tau(3,1)*nx + tau(3,2)*ny + tau(3,3)*nz

      ! Return the fluxes

      Fn(1) =  zero
      Fn(2) = -sigxn
      Fn(3) = -sigyn
      Fn(4) = -sigzn

      Fn(1:4) = aa * Fn(1:4)

      ! Return to calling subprogram 

      return 

   end function vflux

   !*******************************************************************!
   !*******************************************************************!
   !  Function source
   !> @brief Get the volume sources to be applied
   !> @param[in]  q - the primitive variables
   !> @param[in]  ibm - the solid (1) / fluid (0) switch
   !> @param[in]  error - the previous value of the error, for integrating
   !> @param[out] s - the volume sources
   !*******************************************************************!
   !*******************************************************************!

   function source(q, ibm, error) result(s)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP) :: q(npdes), s(npdes), error(npdes)

      integer(kind=WI) :: ibm(1)

      ! Declare internal variables

      real   (kind=WP) :: newerror(npdes)

      ! Avoid a logical test - set IBM switch to 1 for wall and 0 for fluid

      newerror(1) = zero
      newerror(2) = zero - ibm(1)*q(2)
      newerror(3) = zero - ibm(1)*q(3)
      newerror(4) = zero - ibm(1)*q(4)

      ! Integrate the current error for integral controller

      error = error + newerror

      ! Set the volume source vector

      s(1) = zero
      s(2) = zero*newerror(2) + 100.0_wp*error(2)
      s(3) = zero*newerror(3) + 100.0_wp*error(3)
      s(4) = zero*newerror(4) + 100.0_wp*error(4)

      ! Return to calling variables

      return 

   end function source

   !*******************************************************************!
   !*******************************************************************!
   !  Function invWall
   !> @brief Compute the boundary primitives for an inviscid wall
   !> @param[in]  q1 - the primitive variables inside the domain
   !> @param[in]  n  - the normal pointing from 1 to 2
   !> @param[out] q2 - the primitive variables in the "ghost" cell
   !*******************************************************************!
   !*******************************************************************!

   function invWall(q1, n) result(q2)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: n(3), q1(npdes)

      real   (kind=WP) :: q2(npdes), nb(3)

      ! Declare internal variables

      real   (kind=WP) :: vxL, vyL, vzL

      ! Normalise the boundary weights for computing the new velocities

      nb = n / norm2(n)

      ! Calculate the primitive variables in the "ghost cell"

      vxL = q1(2) - 2*(q1(2)*nb(1) + q1(3)*nb(2) + q1(4)*nb(3))*nb(1)
      vyL = q1(3) - 2*(q1(2)*nb(1) + q1(3)*nb(2) + q1(4)*nb(3))*nb(2)
      vzL = q1(4) - 2*(q1(2)*nb(1) + q1(3)*nb(2) + q1(4)*nb(3))*nb(3)

      ! Set the LHS variables from the flow field and BCs

      q2(1) = q1(1)
      q2(2) = vxL
      q2(3) = vyL
      q2(4) = vzL

      ! Return to calling subprogram 

      return 

   end function invWall

   !*******************************************************************!
   !*******************************************************************!
   !  Function visWall
   !> @brief Compute the boundary ghost cells for a viscous wall
   !> @param[in]  q1 - the primitive variables inside the domain
   !> @param[in]  n  - the normal pointing from 1 to 2
   !> @param[out] q2 - the primitive variables in the "ghost" cell
   !*******************************************************************!
   !*******************************************************************!

   function visWall(q1, n) result(q2)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: n(3), q1(npdes)

      real   (kind=WP) :: q2(npdes)

      ! Declare internal variables

      real   (kind=WP) :: vxL, vyL, vzL

      ! Calculate the primitive variables in the "ghost cell"

      vxL = -q1(2) 
      vyL = -q1(3) 
      vzL = -q1(4) 

      ! Set the LHS variables from the flow field and BCs

      q2(1) = q1(1)
      q2(2) = vxL
      q2(3) = vyL
      q2(4) = vzL

      ! Return to calling subprogram 

      return 

   end function visWall

   !*******************************************************************!
   !*******************************************************************!
   !  Function outFlow
   !> @brief Compute the ghost cell primitives for a static pressure outlet
   !> @param[in]  q1 - the primitive variables inside the domain
   !> @param[in]  n  - the normal pointing from 1 to 2
   !> @param[out] q2 - the primitive variables in the "ghost" cell
   !*******************************************************************!
   !*******************************************************************!

   function outFlow(q1, n) result(q2)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: n(3), q1(npdes)

      real   (kind=WP) :: q2(npdes)

      ! Declare internal variables

      real   (kind=WP) :: pB

      ! Set the back pressure

      pB = 100000.0_wp

      ! Set the LHS variables from the flow field and BCs

      q2(1) = pB
      q2(2) = q1(2)
      q2(3) = q1(3)
      q2(4) = q1(4)

      ! Return to calling subprogram 

      return 

   end function outFlow

   !*******************************************************************!
   !*******************************************************************!
   !  Function inFlow
   !> @brief Compute the ghost cell values for a velocity inlet
   !> @param[in]  q1 - the primitive variables inside the domain
   !> @param[in]  n  - the normal pointing from 1 to 2
   !> @param[out] q2 - the primitive variables in the "ghost" cell
   !*******************************************************************!
   !*******************************************************************!

   function inFlow(q1, n) result(q2)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: n(3), q1(npdes)

      real   (kind=WP) :: q2(npdes)

      ! Declare internal variables

      real   (kind=WP) :: ub, vb, wb

      ! Set the boundary conditions

      ! Set the LHS variables from the flow field and BCs

      ub = inflowVel
      vb = 0.0_wp
      wb = 0.0_wp
      
      q2(1) = q1(1)
      q2(2) = ub
      q2(3) = vb
      q2(4) = wb

      ! Return to calling subprogram 

      return 

   end function inFlow

   !*******************************************************************!
   !*******************************************************************!
   !  Function grad
   !> @brief Compute the gradient contribution between two cells
   !> @param[in]  q1 - the primitive variables in cell 1
   !> @param[in]  q2 - the primitive variables in cell 2
   !> @param[in]  n  - the normal pointing from 1 to 2
   !> @param[out] qp - the gradient contribution
   !*******************************************************************!
   !*******************************************************************!

   function grad(q1, q2, n) result(qp)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: n(3), q1(npdes), q2(npdes)

      real   (kind=WP) :: qp(npdes,3)

      ! Set the gradient contribution

      qp(:,1) = half * (q1(:) + q2(:)) * n(1)
      qp(:,2) = half * (q1(:) + q2(:)) * n(2)
      qp(:,3) = half * (q1(:) + q2(:)) * n(3)

      ! Return to calling subprogram 

      return 

   end function grad

end module mod_Equations

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_Equations                                      !
!                                                                    !
!********************************************************************!
!********************************************************************!

