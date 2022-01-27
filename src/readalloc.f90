!********************************************************************!
!********************************************************************!
!                                                                    !
!   mod_ReadAlloc -- Read mesh and flow, set up mesh, sort memory    !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Module created: 20Jan22                 - raw54 !
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
!>  This module does a lot of work reading in the mesh and flow 
!>  files, and setting up the mesh for future computation. It also
!>  allocates the memory necessary for the eventual computation
!>
!*******************************************************************
!*******************************************************************

module mod_ReadAlloc

   ! Declare modules

   use precision
   use mod_globalParameters

   ! Turn off implicit typing

   implicit none

   ! Declare module variables

   ! Mesh variables

   integer(kind=WI) :: nIp, nJp, nKp, nIc, nJc, nKc
   
   real   (kind=WP), allocatable :: x(:,:,:,:), vol(:,:,:,:), xc(:,:,:,:)

   real   (kind=WP), allocatable :: nwt(:,:,:,:,:)

   ! Solution variables

   real   (kind=WP), allocatable :: q(:,:,:,:)

   real   (kind=WP), allocatable :: error(:,:,:,:)

   ! Set everything public

   public 

   ! Hide the helper routines

   private :: Cross, Triple, HexVolume, HexCentre, HexNormal

contains

   !*******************************************************************
   !*******************************************************************
   !
   !>  Subroutine ra_readalloc()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Read a plot3D mesh and flow file, allocate memory, and set up
   !>  the solution process.
   !>  
   !*******************************************************************
   !*******************************************************************

   subroutine ra_readalloc()

      implicit none

      ! Declare external variables

      character(len=CP) :: filename

      ! Declare internal variables

      integer(kind=WI) :: fileState, endState

      integer(kind=WI) :: i, j, k, iPde

      ! Open the plot3D file

      open(unit=191, file=trim(meshFile), status='old', iostat=fileState, form='unformatted')

      if ( fileState .ne. 0 ) then

         write(6,*) ' P3D mesh input - fatal error! ' 
         write(6,*) '    Count not open mesh input file ' // trim(filename)
         stop

      end if
      
      ! Get the data sizes

      read(191,iostat=endState) 

      read(191,iostat=endState) nIp, nJp, nKp

      if ( endState .lt. 0 ) then

         write(6,*) ' P3D data input - fatal error! ' 
         write(6,*) '    Unexpected end of file detected - sizing data ' // trim(filename)
         stop
         
      end if

      nIc = nIp - 1
      nJc = nJp - 1
      nKc = nKp - 1

      ! Allocate space for the mesh point coordinates

      allocate(x(3,nIp,nJp,nKp))
      
      ! Get the mesh data

      read(191) (((x(1,i,j,k), i=1, nIp), j=1, nJp), k=1, nKp), &
                (((x(2,i,j,k), i=1, nIp), j=1, nJp), k=1, nKp), &
                (((x(3,i,j,k), i=1, nIp), j=1, nJp), k=1, nKp)

      if ( endState .lt. 0 ) then

         write(6,*) ' P3D data input - fatal error! ' 
         write(6,*) '    Unexpected end of file detected - coordinate data ' // trim(filename)
         stop
         
      end if

      ! Close the plot3D file

      close(191)

      ! Allocate memory for the cell volumes and for the cell outward facing normals

      allocate(vol(1,nIc,nJc,nKc))
      allocate(xc(3,nIc,nJc,nKc))
      allocate(nwt(3,6,nIc,nJc,nKc))

      ! Calculate the cell volumes and centres

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         vol(1,i,j,k) = HexVolume(x(:,i  ,j  ,k  ),  x(:,i+1,j  ,k  ), &
                                  x(:,i  ,j+1,k  ),  x(:,i+1,j+1,k  ), &
                                  x(:,i  ,j  ,k+1),  x(:,i+1,j  ,k+1), &
                                  x(:,i  ,j+1,k+1),  x(:,i+1,j+1,k+1) )
         xc(:,i,j,k)  = HexCentre(x(:,i  ,j  ,k  ),  x(:,i+1,j  ,k  ), &
                                  x(:,i  ,j+1,k  ),  x(:,i+1,j+1,k  ), &
                                  x(:,i  ,j  ,k+1),  x(:,i+1,j  ,k+1), &
                                  x(:,i  ,j+1,k+1),  x(:,i+1,j+1,k+1) )

         write(60,*) xc(1,i,j,k), xc(2,i,j,k), xc(3,i,j,k)
         write(61,*) vol(1,i,j,k)


      end do
      end do
      end do

      ! Calculate the outward facing normals for each cell (1,2 are i and i+1, 3,4 j and j+1, 5,6 k and k+1 faces)

      do k = 1, nKc
      do j = 1, nJc
      do i = 1, nIc

         nwt(:,:,i,j,k) = HexNormal(x(:,i  ,j  ,k  ),  x(:,i+1,j  ,k  ), &
                                    x(:,i  ,j+1,k  ),  x(:,i+1,j+1,k  ), &
                                    x(:,i  ,j  ,k+1),  x(:,i+1,j  ,k+1), &
                                    x(:,i  ,j+1,k+1),  x(:,i+1,j+1,k+1) )

      end do
      end do
      end do

      ! Allocate memory for the flow primitive variables

      allocate(    q(nPdes,nIc,nJc,nKc))
      allocate(error(nPdes,nIc,nJc,nKc))

      ! Open the plot3D-like flow file

      open(unit=191, file=trim(flowInputFile), status='old', iostat=fileState, form='unformatted')

      if ( fileState .ne. 0 ) then

         write(6,*) ' P3D flow input - fatal error! ' 
         write(6,*) '    Count not open flow input file ' // trim(filename)
         stop

      end if
      
      ! Read the data

      read(191,iostat=endState) 

      read(191,iostat=endState) nIc, nJc, nKc

      if ( endState .lt. 0 ) then

         write(6,*) ' P3D flow input - fatal error! ' 
         write(6,*) '    Unexpected end of file detected - sizing data ' // trim(filename)
         stop
         
      end if
      
      ! Get the mesh data

      read(191) ((((    q(iPde,i,j,k), i=1, nIc), j=1, nJc), k=1, nKc), iPde = 1, nPdes), &
                ((((error(iPde,i,j,k), i=1, nIc), j=1, nJc), k=1, nKc), iPde = 1, nPdes) 

      if ( endState .lt. 0 ) then

         write(6,*) ' P3D flow input - fatal error! ' 
         write(6,*) '    Unexpected end of file detected - flow data ' // trim(filename)
         stop
         
      end if

      ! Close the plot3D file

      close(191)

      ! Return to calling subprogram 

      return 

   end subroutine Ra_Readalloc

   !*******************************************************************
   !*******************************************************************
   !
   !>  Function Cross()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Compute the cross product of two 3D vectors
   !>
   !>  @param[in]  x - the first input vector
   !>  @param[in]  y - the second input vector
   !>  @param[out] z - the output vector, x cross y
   !>
   !*******************************************************************
   !*******************************************************************

   function Cross(x, y) result(z)

      ! Turn off implicit typing 

      implicit none

      ! Declare external variables

      real   (kind=WP), intent(in) :: x(3), y(3)
     
      real   (kind=WP) :: z(3)

      ! Calculate the cross product

      z(1) = x(2) * y(3) - x(3) * y(2)
      z(2) = x(3) * y(1) - x(1) * y(3)
      z(3) = x(1) * y(2) - x(2) * y(1)

      ! Return to calling subprogram 

      return 

   end function Cross

   !*******************************************************************
   !*******************************************************************
   !
   !>  Function Triple()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Computes the vector triple product from two 3D vectors
   !>  
   !>  @param[in]  x - the first input vector
   !>  @param[in]  y - the second input vector
   !>  @param[in]  z - the third input vector
   !>  @param[out] a - the VTP, x dot (y cross z)
   !>
   !*******************************************************************
   !*******************************************************************

   function Triple(x, y, z) result(a)

      implicit none

      ! Declare external variables

      real   (kind=WP) :: x(3), y(3), z(3), a

      ! Compute the triple product
      
      a = dot_product(x, cross(y, z))
      
      ! Return to calling program

      return

   end function Triple

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Function HexVolume()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Computes the volume of a hexahedra using the outer point 
   !>  coordinates, as discussed in "Efficient Computation of Volume of
   !>  Hexahedral Cells", J. Grandy, 1997.
   !>  
   !>  @param[in]  x0 - the first hexahedral vertex coordinate
   !>  @param[in]  x1 - the second hexahedral vertex coordinate
   !>  @param[in]  x2 - the third hexahedral vertex coordinate
   !>  @param[in]  x3 - the fourth hexahedral vertex coordinate
   !>  @param[in]  x4 - the fifth hexahedral vertex coordinate
   !>  @param[in]  x5 - the sixth hexahedral vertex coordinate
   !>  @param[in]  x6 - the seventh hexahedral vertex coordinate
   !>  @param[in]  x7 - the eigth first hexahedral vertex coordinate
   !>  @param[out] V  - the volume of the hexahedron
   !>
   !*******************************************************************
   !*******************************************************************

   function HexVolume(x0, x1, x2, x3, x4, x5, x6, x7) result(V)

      implicit none

      ! Declare external variables

      real   (kind=WP) :: x0(3), x1(3), x2(3), x3(3), x4(3), x5(3), x6(3), x7(3)

      real   (kind=WP) :: V

      ! Compute the volume
      
      V = (1.0_wp / 12.0_wp) * & 
             Triple((x7-x1) + (x6-x0), (x7-x2), (x3-x0)) + &
             Triple((x6-x0), (x7-x2) + (x5-x0), (x7-x4)) + &
             Triple((x7-x1), (x5-x0), (x7-x4) + (x3-x0)) 

      ! Return to calling program

      return

   end function HexVolume

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Function HexCentre()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Computes the centre of a hexahedra using the outer point 
   !>  coordinates.
   !>  
   !>  @param[in]  x0 - the first hexahedral vertex coordinate
   !>  @param[in]  x1 - the second hexahedral vertex coordinate
   !>  @param[in]  x2 - the third hexahedral vertex coordinate
   !>  @param[in]  x3 - the fourth hexahedral vertex coordinate
   !>  @param[in]  x4 - the fifth hexahedral vertex coordinate
   !>  @param[in]  x5 - the sixth hexahedral vertex coordinate
   !>  @param[in]  x6 - the seventh hexahedral vertex coordinate
   !>  @param[in]  x7 - the eigth first hexahedral vertex coordinate
   !>  @param[out] C  - the centre of the hexahedron
   !>
   !*******************************************************************
   !*******************************************************************

   function HexCentre(x0, x1, x2, x3, x4, x5, x6, x7) result(C)

      implicit none

      ! Declare external variables

      real   (kind=WP) :: x0(3), x1(3), x2(3), x3(3), x4(3), x5(3), x6(3), x7(3)

      real   (kind=WP) :: C(3)

      ! Compute the centre

      C = (x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7) * (1.0_wp/8.0_wp)
      
      ! Return to calling program

      return

   end function HexCentre

   !*******************************************************************
   !*******************************************************************

   !*******************************************************************
   !*******************************************************************
   !
   !>  Function HexNormal()
   !
   !>  @author
   !>  Rob Watson
   !>
   !>  @brief 
   !>  Computes the face normals of a hexahedra using the outer point 
   !>  coordinates.
   !>  
   !>  @param[in]  x0 - the first hexahedral vertex coordinate
   !>  @param[in]  x1 - the second hexahedral vertex coordinate
   !>  @param[in]  x2 - the third hexahedral vertex coordinate
   !>  @param[in]  x3 - the fourth hexahedral vertex coordinate
   !>  @param[in]  x4 - the fifth hexahedral vertex coordinate
   !>  @param[in]  x5 - the sixth hexahedral vertex coordinate
   !>  @param[in]  x6 - the seventh hexahedral vertex coordinate
   !>  @param[in]  x7 - the eigth first hexahedral vertex coordinate
   !>  @param[out] N  - the outward facing normals of the hexahedron
   !>
   !*******************************************************************
   !*******************************************************************

   function HexNormal(x0, x1, x2, x3, x4, x5, x6, x7) result(N)

      implicit none

      ! Declare external variables

      real   (kind=WP) :: x0(3), x1(3), x2(3), x3(3), x4(3), x5(3), x6(3), x7(3)

      real   (kind=WP) :: N(3,6)

      ! Compute the normals

      N(:,1) = Cross(x4-x2, x6-x0)
      N(:,2) = Cross(x7-x1, x5-x3)
      N(:,3) = Cross(x5-x0, x4-x1)
      N(:,4) = Cross(x6-x3, x7-x2)
      N(:,5) = Cross(x2-x1, x3-x0)
      N(:,6) = Cross(x7-x4, x6-x5)
      
      ! Return to calling program

      return

   end function HexNormal

   !*******************************************************************
   !*******************************************************************

end module mod_ReadAlloc

!********************************************************************!
!********************************************************************!
!                                                                    !
!   End of module ra_ReadAlloc                                       !
!                                                                    !
!********************************************************************!
!********************************************************************!

