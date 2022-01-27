!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_Input -- A module for reading the input file                 !
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
!>  @brief This module reads the input file and stores global input
!>         variables, such as filenames, etc.
! 
!********************************************************************!
!********************************************************************!

module mod_Input

   ! Declare modules

   use precision
   use mod_globalParameters

   implicit none

contains

   !*******************************************************************!
   !*******************************************************************!
   !  Subroutine input_getInput
   !
   !> @author 
   !> Rob Watson
   !
   !> @brief Reads the input file
   !> @param[in]  
   !> @param[out] 
   !*******************************************************************!
   !*******************************************************************!

   subroutine input_getInput()

      ! Turn off implicit typing 

      implicit none

      ! Declare internal variables

      integer(kind=WI)  :: ioErr, nLine, iLine, cLine

      character(len=CP) :: inputs(100), testStr

      logical           :: fileExists

      ! Set up the number of partial differential equations (from file later)
      nPdes = 3

      ! Test if the input file exists

      inquire(file="input.dat", exist=fileExists)

      if ( .not. fileExists) then
         write(6,*) ' Can not find input.dat file...exiting.'
         stop
      end if

      ! Open the input file

      open(file='input.dat', unit=101)

      ! Get the contents of the file

      nLine = 1

      do

         ! Read the contents, with end of file checking

         read(101,'(A70)',iostat=ioErr) inputs(nLine)

         if ( ioErr .eq. 0 ) then
            nLine = nLine + 1
         else
            exit
         end if

      end do
      
      ! Close the input file

      close(101)

      ! Set the current line pointer to the top of the file

      cLine = 0

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - input file name
    
      read(inputs(cLine), *) meshFile

        ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - flow input name
      
      read(inputs(cLine), *) flowInputFile
 

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - flow output names
      
      read(inputs(cLine), *) flowRestartFile, flowOutputFile

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - output and restart write frequency
      
      read(inputs(cLine), *) nWriteFreq, nOutputFiles

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - steady or unsteady calculation
    
      read(inputs(cLine), *) steadySwitch

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - number of inner iterations
    
      read(inputs(cLine), *) nInnerIt

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - CFL number
    
      read(inputs(cLine), *) CFL

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - number of outer iterations for unsteady
    
      read(inputs(cLine), *) nOuterIt

      ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - real time step, for unsteady
    
      read(inputs(cLine), *) dtR

           ! Find the next line with ###s, move read cursor to after section

      iLine = cLine
      do 
         iLine = iLine + 1
         read(inputs(iLine),'(A5)') testStr
         if ( trim(testStr) .eq. '#####') then
            cLine = iLine + 3
            exit
         end if
      end do

      ! Read the next piece of data - flow input name
      
      read(inputs(cLine), *) inflowVel
      
      ! Return to calling subprogram 

      return 

   end subroutine input_getInput

end module mod_Input

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_Input                                          !
!                                                                    !
!********************************************************************!
!********************************************************************!
