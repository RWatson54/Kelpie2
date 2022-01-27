!********************************************************************!
!********************************************************************!
!                                                                    !
!   Kelpie -- A refactored Kelpie IBM CFD code                       !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 20Jan22                - raw54 !
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
!>  This version of Kelpie is a structured FV CFD solver, designed to
!>  be used to carry out immersed boundary investigation of flow 
!>  through porous media.
!>
!*******************************************************************
!*******************************************************************

program Kelpie

   ! Declare modules

   use precision
   use mod_HelloGoodbye
   use mod_GlobalParameters
   use mod_Input
   use mod_ReadAlloc
   use mod_PDESolve
   use mod_Update
   use mod_Output
   use mod_setIBM

   ! Turn off implicit typing

   implicit none

   ! Say hello, welcoming everyone to the execution of the run

   call hg_Hello()

   ! Read the input file

   write(6,*) ' *** '
   write(6,*) ' *** Reading in the input file...'

   call input_getInput()

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '

   ! Read in the mesh and flow files, and allocate solution memory

   write(6,*) ' *** '
   write(6,*) ' *** Reading in and preparing the mesh and flow...'

   call ra_readalloc()

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '   

   ! Set up the IBM locations

   write(6,*) ' *** '
   write(6,*) ' *** Setting up the IBM locations '

   call ibm_setSolid()

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '   

   ! Print out the solution

   write(6,*) ' *** '
   write(6,*) ' *** Printing out a VTU file for examination '

   call output_VTU()

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '   

   ! Open a file to store the residual history (overwrites and existing file)

   open(file='history.dat', unit=15)

   ! ********************************************
   ! Beginning of real timestepping loop
   ! ********************************************
   do iRTime = 1, nOuterIt

      ! Accumulate the "total residual counter"

      iIteration = iIteration + 1

      ! Get the residuals

      call pde_getResiduals()

      ! Update the solution

      call update_explicit()

      ! Print progress to screen

      call output_hist()

      if ( mod(iRTime, 500) .eq. 0 ) then
         call output_VTU()
      end if
      
   end do
   ! ********************************************
   ! End of real timestepping loop
   ! ********************************************

   ! Close the history file

   close(15)

   ! Say goodbye, letting everyone know the dream is over

   call hg_Goodbye()

contains

end program Kelpie

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of program Kelpie                                            !
!                                                                    !
!********************************************************************!
!********************************************************************!

