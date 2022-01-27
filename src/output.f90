!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_output -- A module for outputting data                       !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 03Jan21                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!

module mod_Output

   ! Declare modules

   use precision
   use mod_globalParameters
   use mod_ReadAlloc

   implicit none

contains

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine output_restart()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  This subroutine prints a "plot3d-like" file for restarting the 
   !>  simulation.
   !>
   !*******************************************************************!
   !*******************************************************************!

   subroutine output_restart()

      ! Turn off implicit typing 

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: i, j, k, ipde

      ! Write out the restart file
      
      open(file=flowRestartFile, unit=102, form='unformatted')

      write(102) 1
      write(102) nIc, nJc, nKc
      write(102) ((((q(ipde,i,j,k), i=1,nIc), j=1, nJc), k=1, nKc), iPde=1,nPdes), &
                 ((((e(ipde,i,j,k), i=1,nIc), j=1, nJc), k=1, nKc), iPde=1,nPdes)

      close(102)

      ! Return to calling subprogram 

      return 

   end subroutine output_restart
  
   !*******************************************************************!
   !*******************************************************************!
   !>
   !>  Subroutine output_VTU
   !>
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  This subroutine prints a legacy VTU file for viewing the results
   !>  of the simulation.
   !>
   !*******************************************************************!
   !*******************************************************************!

   subroutine output_VTU()

      ! Turn off implicit typing 

      implicit none

      ! Declare internal variables

      character(len=CP) :: filename

      character(len=5)  :: countStr

      integer(kind=WI) :: i, j, k, ma, iE, iN

      integer(kind=WI) :: Cell(8,nIc*nJc*nKc)

      integer(kind=WI) :: nFile

      ! For convenience, get a list of the cells

      ma = 0

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         ma = ma + 1

         Cell(1,ma) = (k-1)*nJp*nIp + (j-1)*nIp + (i-1) + 1
         Cell(2,ma) = (k-1)*nJp*nIp + (j-1)*nIp + (i  ) + 1
         Cell(3,ma) = (k-1)*nJp*nIp + (j  )*nIp + (i  ) + 1
         Cell(4,ma) = (k-1)*nJp*nIp + (j  )*nIp + (i-1) + 1
         Cell(5,ma) = (k  )*nJp*nIp + (j-1)*nIp + (i-1) + 1
         Cell(6,ma) = (k  )*nJp*nIp + (j-1)*nIp + (i  ) + 1
         Cell(7,ma) = (k  )*nJp*nIp + (j  )*nIp + (i  ) + 1
         Cell(8,ma) = (k  )*nJp*nIp + (j  )*nIp + (i-1) + 1

      end do
      end do
      end do

      ! Work out the file name in question

      nFile = mod(iIteration/nWriteFreq - 1, nOutputFiles) + 1
      
      write(countStr(1:1),'(A1)') '_'
      write(countStr(2:5),'(I0.4)') nFile
      
      filename = trim(flowOutputFile)//countStr//'.vtu'

      ! And print out to a VTU format

      open(file=trim(filename),unit=12)

      write(12,'(A80)') '<VTKFile type="UnstructuredGrid" byte_order="LittleEndian">                     '
      write(12,'(A80)') '<UnstructuredGrid>                                                              '
      write(12,'(A23,I10,A17,I10,A2)') '<Piece NumberOfPoints="',nIp*nJp*nKp,'" NumberOfCells="',nIc*nJc*nKc,'">'
      write(12,'(A80)') '<CellData>                                                                      '
      write(12,'(A80)') '<DataArray type="Float32" Name="Pressure" format="ascii">                        '
      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         write(12,*) q(1,i,j,k)
      end do
      end do
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         write(12,*) q(2,i,j,k), q(3,i,j,k), q(4,i,j,k)
      end do
      end do
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '<DataArray type="Float32" Name="IBM Switch" format="ascii">                     '
      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc
         write(12,*) is(1,i,j,k)
      end do
      end do
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '</CellData>                                                                     '
      write(12,'(A80)') '<Points>                                                                        '
      write(12,'(A80)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      do k = 1, nKp
      do j = 1, nJp
      do i = 1, nIp
         write(12,*) x(1,i,j,k), x(2,i,j,k), x(3,i,j,k)
      end do
      end do
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '</Points>                                                                       '
      write(12,'(A80)') '<Cells>                                                                         '
      write(12,'(A80)') '<DataArray type="Int32" Name="connectivity" format="ascii">                     '
      do iE = 1, nIc*nJc*nKc
         write(12,'(8(I10,2X))') (Cell(in,iE)-1, in=1,8)
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '<DataArray type="Int32" Name="offsets" format="ascii">                          '
      ma = 0
      do iE = 1, nIc*nJc*nKc
         ma = ma + 8
         write(12,'(I10)') ma 
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '<DataArray type="UInt8" Name="types" format="ascii">'
      do iE = 1, nIc*nJc*nKc
         write(12,'(A2)') '12'
      end do
      write(12,'(A80)') '</DataArray>                                                                    '
      write(12,'(A80)') '</Cells>                                                                        '  
      write(12,'(A80)') '</Piece>                                                                        '  
      write(12,'(A80)') '</UnstructuredGrid>                                                             '  
      write(12,'(A80)') '</VTKFile>                                                                      '  
      write(12,'(A80)') '                                                                                '

      ! Close the file

      close(12)

      ! Return to calling subprogram 

      return 

   end subroutine output_VTU

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine output_hist()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  This subroutine prints a convergence history to file and to 
   !>  screen.
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine output_hist()

      implicit none

      ! Declare internal variables

      integer(kind=WI) :: iPde

      ! Print residuals to screen and file

      if ( mod(iIteration,10) .eq. 0) then

         write(6,*) iRTime, iPTime, (log10(MaxRes(iPde)), iPde = 1, nPdes)

      end if

      write(15,*) iRTime, iPTime, (log10(MaxRes(iPde)), iPde = 1, nPdes)

      ! Return to calling subprogram 

      return 

   end subroutine output_hist

end module mod_Output

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module mod_Output                                         !
!                                                                    !
!********************************************************************!
!********************************************************************!
