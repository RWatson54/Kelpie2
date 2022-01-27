!********************************************************************!
!********************************************************************!
!                                                                    !
!   Write_P3D -- Writes a simple plot3D grid                         !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 17Oct20                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!

program Write_P3D

   ! Declare modules

   use precision

   implicit none

   ! Declare variables

   ! Mesh parameters

   integer(kind=WI) :: nIp, nJp, nKp
   integer(kind=WI) :: nIc, nJc, nKc

   real   (kind=WP), allocatable :: x(:,:,:,:)

   ! Flow parameters

   integer(kind=WI), parameter :: nPdes = 4

   real   (kind=WP), allocatable :: q(:,:,:,:), e(:,:,:,:)

   ! Internal variables
   
   integer(kind=WI) :: i, j, k, ipde

   real   (kind=WP) :: xa, ya, za

   ! Say hello

   call main_Hello()

   ! Here, look at the input.dat file to set up the mesh as required

   nIp = 600
   nJp = 300
   nKp = 5

   nIc = nIp - 1
   nJc = nJp - 1
   nKc = nKp - 1

   ! Allocate memory for the coordinates

   allocate(x(3,nIp,nJp,nKp))
   
   ! Set up the mesh

   write(6,*) ' *** '
   write(6,*) ' *** Generating mesh...'

   do k = 1, nKp
   do j = 1, nJp
   do i = 1, nIp

      xa = (0.120_wp/(one*(nIp-1))) * (i-1)
      ya = (0.020_wp/(one*(nJp-1))) * (j-1) - 0.01_wp
      za = (0.005_wp/(one*(nKp-1))) * (k-1) - 0.0025_wp

      x(1,i,j,k) = xa
      x(2,i,j,k) = ya
      x(3,i,j,k) = za

   end do
   end do
   end do

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '

   ! Allocate memory for the initial flow

   allocate(q(nPdes,nIc,nJc,nKc))
   allocate(e(nPdes,nIc,nJc,nKc))

   ! Set up the initial flow, read in from input.dat file if required

   write(6,*) ' *** '
   write(6,*) ' *** Generating initial flow guess...'

   q(1,:,:,:) = 100000d0
   q(2,:,:,:) = 5d0
   q(3,:,:,:) = 0d0
   q(4,:,:,:) = 0d0

   e = 0d0

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '

   ! Write out in P3D format (as binary, for size)

   write(6,*) ' *** '
   write(6,*) ' *** Writing the mesh file...'

   open(file='mesh.xyz',unit=101, form='unformatted')

   write(101) 1
   write(101) nIp, nJp, nKp
   write(101) (((x(1,i,j,k), i=1,nIp), j=1,nJp), k=1,nKp), &
              (((x(2,i,j,k), i=1,nIp), j=1,nJp), k=1,nKp), &
              (((x(3,i,j,k), i=1,nIp), j=1,nJp), k=1,nKp) 

   close(101)

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '

   ! Write out an initial flow field

   write(6,*) ' *** '
   write(6,*) ' *** Writing the initial flow file...'

   open(file='flow.q', unit=102, form='unformatted')

   write(102) 1
   write(102) nIc, nJc, nKc
   write(102) ((((q(ipde,i,j,k), i=1,nIc), j=1, nJc), k=1, nKc), iPde=1,nPdes), &
              ((((e(ipde,i,j,k), i=1,nIc), j=1, nJc), k=1, nKc), iPde=1,nPdes)
   close(102)

   write(6,*) ' *** ...done.'
   write(6,*) ' *** '

   ! Say bye bye
   
   call main_Goodbye()

contains

   !*******************************************************************
   !*******************************************************************

   subroutine main_Hello()

      implicit none

      ! Declare variables

      ! Welcome screen

      write(6,*)
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' *****   Kelpie Meshing and Initialisation   ***** '
      write(6,*) ' *****                                       ***** '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' ************************************************* '
      write(6,*) ' *** '


      ! Return to calling program

      return

   end subroutine main_Hello

   !*******************************************************************
   !*******************************************************************

   subroutine main_Goodbye()

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

   end subroutine main_Goodbye

   !*******************************************************************
   !*******************************************************************

end program Write_P3D

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of program Write_P3D                                         !
!                                                                    !
!********************************************************************!
!********************************************************************!


