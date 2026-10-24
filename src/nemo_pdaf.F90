!> Include NEMO variables and initialize NEMO grid information
!!
!! This module includes variables from NEMO. For PDAF in its
!! online coupling it is the single point which directly links
!! to NEMO. All other routines include from this module.
!! For the offline case the variables are declared here and
!! separately initialized.
!!
!! Next to including information from NEMO, the routine
!! set_nemo_grid initializes arrays holding the grid information
!! for use with state vectors in the PDAF user code for NEMO.
!!
module nemo_pdaf

   use mod_kind_pdaf
   implicit none

   ! *** NEMO domain variables ***
   integer :: jpiglo, jpjglo, jpk        ! Global NEMO grid dimensions
   integer :: numcat                     ! Number of categories in the sea ice model
   integer :: halo0(2), halo1(2)         ! The halo information for restart files
   integer :: time_counter(1)            ! Time counter from restart file
   integer  :: nn_time0                 ! initial time of day in hhmm
   real(pwp) :: ndastp                   ! NEMO time string
                                         ! spcified in NEMO namelist namrun

   real(pwp), allocatable :: glamt(:,:), glamu(:,:), glamv(:,:)       ! Longitudes
   real(pwp), allocatable :: gphit(:,:), gphiu(:,:), gphiv(:,:)       ! Latitudes
   real(4), allocatable   :: nav_lev(:)      ! Depths
   real(4), allocatable   :: nav_lon(:,:), nav_lat(:,:)   ! Restart file grid
   real(pwp), allocatable :: tmask(:,:,:)     ! Temperature mask array
   ! wet points for state vectors
   integer :: use_wet_state=2               ! 1: State vector contains full columns where surface grid point is wet
                                            ! 2: State vector only contains wet grid points
                                            ! other: State vector contains 2d/3d grid box
   integer :: nwet                          ! Number of surface wet grid points
   integer :: nwet3d                        ! Number of 3d wet grid points
   integer, allocatable :: wet_pts(:,:)     ! Index array for wet grid points
                                            ! (1) Global latitude index,
                                            ! (2) Global longitude index,
                                            ! (3) number wet layers at given latlon
                                            ! (4) index in 2d grid box
                                            ! (5) starting index in all wet points for vertical column
                                            ! (6) local longitude index in subdomain
                                            ! (7) local latitude index in subdomain
   integer, allocatable :: idx_wet_2d(:,:)  ! Index array for wet_pts row index in 2d box
   integer, allocatable :: idx_nwet(:,:)    ! Index array for wet_pts row index in wet surface grid points
   integer, allocatable :: nlev_wet_2d(:,:) ! Number of wet layers for ij position in 2d box

   integer :: ni_p, nj_p, nk_p               ! Size of decomposed grid
   integer :: i0, j0                         ! Start indices for internal local domain
   integer :: dim_2d_p, dim_3d_p             ! Dimension of 2d/3d grid box of sub-domain
   integer :: sdim2d, sdim3d                 ! 2D/3D dimension of field in state vector
   ! Constants for coordinate calculations
   real(pwp), parameter :: pi = 3.14159265358979323846_pwp
   real(pwp), parameter :: deg2rad = pi / 180.0_pwp       ! Conversion from degrees to radian


contains
   !> Initialize grid information for the DA
   !!
   !! This routine initializes grid information for the data assimilation
   !! In particular index information is initialize to map in between
   !! the model grid and the state vector
   !!
   subroutine set_nemo_grid()
      use PDAF, only: PDAFomi_set_domain_limits
      use config_pdaf, only: screen
      use mod_kind_pdaf
      use parallel_pdaf, only: mype_model, npes_model, comm_model, &
                               MPI_INT, MPI_SUM, MPIerr
      implicit none
      ! *** Local variables ***
      integer :: i, j, k                    ! Counters
      integer :: cnt, cnt_all, cnt_layers   ! Counters
      integer :: nwet_g, nwet3d_g           ! Global sums of wet grid point
      real(pwp) :: lim_coords(2,2)          ! Limiting coordinates of sub-domain
      ! *** set dimension of 2d and 3d fields in state vector ***
      ! Size of 2d/3d boxes without halo
      dim_2d_p = ni_p * nj_p
      dim_3d_p = ni_p * nj_p * nk_p
      ! Count number of surface points
      cnt = 0
      do k = 1, nk_p
         do j = 1, nj_p
            do i = 1, ni_p
               cnt = cnt + 1
               if (tmask(i, j, k) == 1.0_pwp) then
                  if (k==1) nwet = nwet + 1
                  nwet3d = nwet3d + 1
               endif
            enddo
         enddo
      enddo
      ! Initialize index arrays
      ! - for mapping from nx*ny grid to vector of wet points
      ! - mask for wet points
      if (nwet > 0) then
         allocate(wet_pts(7, nwet))
         cnt = 0
         cnt_all = 0
         do j = 1, nj_p
            do i = 1, ni_p
               cnt_all = cnt_all + 1
               if (tmask(i, j, 1) == 1.0_pwp) then
                  cnt = cnt + 1
                  wet_pts(1,cnt) = i + i0 - 1   ! Global longitude index
                  wet_pts(2,cnt) = j + j0 - 1   ! Global latitude index
                  wet_pts(6,cnt) = i            ! Longitude index in subdomain
                  wet_pts(7,cnt) = j            ! Latitidue index in subdomain
                  ! Determine number of wet layers
                  cnt_layers = 0
                  do k = 1, nk_p
                     if (tmask(i, j, k) == 1.0_pwp) cnt_layers = cnt_layers + 1
                  end do
                  wet_pts(3,cnt) = cnt_layers
                  wet_pts(4,cnt) = cnt_all
               end if
            end do
         end do

         ! row 5 stores wet_pts index for vertical column storage in state vector
         wet_pts(5,1) = 1
         do i = 2 , nwet
            wet_pts(5,i) = wet_pts(5,i-1) + wet_pts(3,i-1)
         end do

      else
         write (*, '(8x,a,i3)') 'WARNING: No valid local domains, PE=', mype_model
         nwet = 0
         allocate(wet_pts(3, 1))
         wet_pts(:, :) = 0
      end if

      ! Initialize index arrays
      ! - for mapping from vector of wet points to 2d box
      ! - for mapping from vector model grid point coordinats to wet point index
      ! - for retrieving number of wet layers at a grid point coordinate
      ! these arrays also serve as mask arrays (0 indicates land point)

      allocate(idx_wet_2d(ni_p, nj_p))
      allocate(idx_nwet(ni_p, nj_p))
      allocate(nlev_wet_2d(ni_p, nj_p))

      idx_wet_2d = 0
      idx_nwet = 0
      nlev_wet_2d = 0
      do i = 1 , nwet
         idx_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(4,i)
         nlev_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(3,i)
      end do
      if (use_wet_state/=2) then
         do i = 1 , nwet
            idx_nwet(wet_pts(6,i), wet_pts(7,i)) = i
         end do
      else
         do i = 1 , nwet
            idx_nwet(wet_pts(6,i), wet_pts(7,i)) = wet_pts(5,i)
         end do
      end if

      ! ********************************************************
      ! *** Set dimension of 2D and 3D field in state vector ***
      ! ********************************************************
      if (use_wet_state==1) then
         ! State vector contains full columns when surface grid point is wet
         sdim3d = abs(nwet)*jpk
         sdim2d = abs(nwet)
      elseif (use_wet_state==2) then
         ! State vector only contains wet grid points
         sdim3d = abs(nwet3d)
         sdim2d = abs(nwet)
      else
         ! State vector contains 2d/3d grid box
         sdim3d = dim_3d_p
         sdim2d = dim_2d_p
      end if

      ! *** Screen output ***
      if (npes_model==1) then
         write(*,'(a,5x,a,3x,i11)') 'NEMO-PDAF', 'Number of wet surface points', nwet
         write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', 'Number of 3D wet points', nwet3d
         write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', '2D wet points * nlayers', nwet*nk_p
      else
         if (screen==1 .or. screen==2) then
            ! Get global sums
            call MPI_Reduce (nwet, nwet_g, 1, MPI_INT, MPI_SUM, &
                  0, COMM_model, MPIerr)
            call MPI_Reduce (nwet3d, nwet3d_g, 1, MPI_INT, MPI_SUM, &
                  0, COMM_model, MPIerr)

            if (mype_model==0) then
               write(*,'(a,5x,a,3x,i11)') &
                     'NEMO-PDAF', 'Number of global wet surface points', nwet_g
               write(*,'(a,5x,a,8x,i11)') &
                     'NEMO-PDAF', 'Number of global 3D wet points', nwet3d_g
               write(*,'(a,5x,a,8x,i11)') &
                     'NEMO-PDAF', 'global 2D wet points * nlayers', nwet_g*nk_p
            end if
         elseif (screen>2) then
            write(*,'(a,2x,a,1x,i4,2x,a,3x,i11)') &
                  'NEMO-PDAF', 'PE', mype_model, 'Number of wet surface points', nwet
            write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
                  'NEMO-PDAF', 'PE', mype_model, 'Number of 3D wet points', nwet3d
            write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
                  'NEMO-PDAF', 'PE', mype_model, '2D wet points * nlayers', nwet*nk_p
         end if
      end if
      ! ******************************************************************
      ! *** Specify domain limits to limit observations to sub-domains ***
      ! ******************************************************************
      lim_coords(1,1) = minval(glamt(:, :)) * deg2rad
      lim_coords(1,2) = maxval(glamt(:, :)) * deg2rad
      lim_coords(2,1) = maxval(gphit(:, :)) * deg2rad
      lim_coords(2,2) = minval(gphit(:, :)) * deg2rad

      call PDAFomi_set_domain_limits(lim_coords)

   end subroutine set_nemo_grid

end module nemo_pdaf
