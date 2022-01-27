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

   ! Say goodbye, letting everyone know the dream is over

   call hg_Goodbye()

   ! Stop and exit cleanly

   stop

contains

end program Kelpie

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of program Kelpie                                            !
!                                                                    !
!********************************************************************!
!********************************************************************!

