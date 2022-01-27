!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_HelloGoodbye -- A module for saying hello and goodbye        !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Module created: 23Jan22                 - raw54 !
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
!>  This is a simple module which just prints a welcome message to 
!>  the screen, and a termination message at completion.
!>  
!>  
!>
!*******************************************************************
!*******************************************************************

module mod_HelloGoodbye

   ! Declare modules

   use precision

   ! Turn off implicit typing

   implicit none

   ! Set everything private

   private   

   ! Turn on what should be publicly visible

   public :: hg_Hello, hg_Goodbye

contains

   !*******************************************************************
   !*******************************************************************

   subroutine hg_Hello()

      implicit none

      ! Declare variables

      ! Welcome screen

      write(6,*)
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' *****        Kelpie 2.0 CFD IBM Code        ***** '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' *** '

      ! Return to calling program

      return

   end subroutine hg_Hello

   !*******************************************************************
   !*******************************************************************

   subroutine hg_Goodbye()

      implicit none

      ! Declare variables

      ! Leaving screen

      write(6,*) ' *** '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' *****           End of execution.           ***** '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) 

      ! Return to calling program

      return

   end subroutine hg_Goodbye
   
   !*******************************************************************
   !*******************************************************************

end module mod_HelloGoodbye

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_HelloGoodbye                                   !
!                                                                    !
!********************************************************************!
!********************************************************************!

