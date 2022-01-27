!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_globalParameters -- A module which stores key variables      !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 21May21                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!
! 
!   Doxygen section: 
! 
!>  @author 
!>  Rob Watson
! 
!>  @brief This module stores variables for easy access by other
!>         modules, such as filenames, etc.
! 
!********************************************************************!
!********************************************************************!

module mod_globalParameters

   ! Declare modules

   use precision

   implicit none

   ! First, the variables which are important for the overall calculation

   ! Filenames various (read from input.dat)
   character(len=CP) :: meshFile, flowRestartFile, flowOutputFile, flowInputFile

   ! Main iteration loop numbers (read from input.dat)
   integer(kind=WI) :: nInnerIt, nOuterIt

   ! A switch for if we're steady or unsteady (read from input.dat)
   character(len=CP) :: steadySwitch

   ! Frequency of write out and number of files (read from input.dat)
   integer(kind=WI) :: nWriteFreq, nOutputFiles

   ! CFL number and real timestep size (read from input.dat)
   real(kind=WP) :: CFL, dtR

   ! Set velocity at the inflow from the input file
   real(kind=WP) :: inflowVel

   ! Set the number of PDEs
   integer(kind=WI) :: nPdes

contains

end module mod_globalParameters

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_globalParameters                               !
!                                                                    !
!********************************************************************!
!********************************************************************!
