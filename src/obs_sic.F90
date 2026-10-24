!> PDAF-OMI observation module for ssh observations (on model grid)
!!
!! Observation type: SSH on model grid
!!
!! The subroutines in this module are for the particular handling of
!! ssh observations available on the model grid. The observation module
!! also allows to perform a twin experiment in which model output
!! is read and used as observations after adding noise to the values.
!!
!! The routines are called by the different call-back routines of PDAF.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initialized by the generic routines from `PDAFomi`.
!!
!! Author: Nicholas Byrne, NCEO & University of Reading, UK
!!
module obs_sic
   use mod_kind_pdaf
   use parallel_pdaf, only: mype_filter, abort_parallel
   use PDAF, only: obs_f, obs_l
   use nemo_pdaf, only: i0, j0
   use netcdf
   implicit none
   save

   !> Whether to assimilate this data type
   logical :: assim_sic = .false.
   !> (1) use global obs.; (0) use domain-reduced full obs.
   integer :: use_global_obs = 1
   !> index in sfields array
   integer :: id_sfields = 1
   !> Observation error standard deviation (for constant errors)
   real(pwp) :: rms = 0.1 !1
   !> Localization cut-off radius
   real(pwp) :: lradius = 1.0
   !> Support radius for weight function
   real(pwp) :: sradius = 1.0
   !> Standard deviation for Gaussian noise in twin experiment
   real(pwp) :: noise_amp = 1
   !> path to obs file
   character(lc) :: file_sic = 'sic_obs.nc'
   !> Name of SIC variable in the observation file
   character(lc) :: varname_sic = 'sic'
   !> Time step index to read from the observation file
   integer :: nc_step = 10 !WILL NEED TO BE ABLE TO LINE THESE UP
   ! the path to restart files is constructed by the following variables
   ! path_sic_op_root//ens_prefix[1-N]//path_sic_op_suffix
   logical :: do_obs_op_from_model_output = .false. ! whether to use model output as obs for twin experiment
   character(len=lc) :: path_op_root
   character(len=lc) :: ens_prefix = 'ens_'
   character(len=lc) :: path_op_suffix
   integer :: op_time_index = 1
   !> Variable name of SIC in the model output file for observed 15 day-averaged SIC
   character(lc) :: varname_obs_op = 'siconc'

   namelist /obs_sic_nml/ assim_sic, use_global_obs, id_sfields, rms, &
                            lradius, sradius, noise_amp, file_sic, &
                            varname_sic, nc_step, &
                            do_obs_op_from_model_output, &
                            path_op_root, ens_prefix, path_op_suffix, &
                            varname_obs_op, op_time_index

   !> Instance of full observation data type - see `PDAFomi` for details.
   type(obs_f), target, public :: thisobs
   !> Instance of local observation data type - see `PDAFomi` for details.
   type(obs_l), target, public :: thisobs_l

!$OMP THREADPRIVATE(thisobs_l)

contains
   !> print configuration of this observation type to screen
   subroutine print_sic_configuration()
      write (*, '(a,3x,a)') 'NEMO-PDAF', '[obs_sic_nml]:'
      write (*, '(a,5x,a,5x,l)') 'NEMO-PDAF',      'assim_sic', assim_sic
      if (assim_sic) then
          write (*, '(a,5x,a,f12.4)')  'NEMO-PDAF','rms            ', rms
          write (*, '(a,5x,a,i0)')     'NEMO-PDAF','use_global_obs ', use_global_obs
          write (*, '(a,5x,a,f12.4)') 'NEMO-PDAF', 'noise_amp      ', noise_amp
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','lradius        ', lradius
          write (*, '(a,5x,a,es12.4)') 'NEMO-PDAF','sradius        ', sradius
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'obs_file   ', trim(file_sic)
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'varname', trim(varname_sic)
          write (*, '(a,5x,a,i0)')   'NEMO-PDAF', 'nc_step          ', nc_step
          write (*, '(a,5x,a,l)')    'NEMO-PDAF', 'do_obs_op_from_model_output', do_obs_op_from_model_output
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'obs_op_file', trim(path_op_root)//trim(ens_prefix &
                                                                 )//'_XXX'//trim(path_op_suffix)
          write (*, '(a,5x,a,a)')    'NEMO-PDAF', 'obs_op_varname', trim(varname_obs_op)
          write (*, '(a,5x,a,i0)')   'NEMO-PDAF', 'op_time_index', op_time_index
      end if
   end subroutine print_sic_configuration

   !>Initialize information on the observation
   !!
   !! The routine is called by each filter process.
   !! at the beginning of the analysis step before
   !! the loop through all local analysis domains.
   !!
   !! It has to count the number of process-local and full
   !! observations, initialize the vector of observations
   !! and their inverse variances, initialize the coordinate
   !! array and index array for indices of observed elements
   !! of the state vector.
   !!
   !! The following four variables have to be initialized in this routine:
   !!
   !! - **thisobs%doassim** - Whether to assimilate ssh
   !! - **thisobs%disttype** - type of distance computation for localization
   !! with ssh
   !! - **thisobs%ncoord** - number of coordinates used for distance
   !! computation
   !! - **thisobs%id_obs_p** - index of module-type observation in PE-local state
   !! vector
   !!
   !!
   !! Optional is the use of:
   !! - **thisobs%icoeff_p**       - Interpolation coefficients for obs. operator (only if interpolation is used)
   !! - **thisobs%domainsize**     - Size of domain for periodicity for *disttype=1* (<0 for no periodicity)
   !! - **thisobs%obs_err_type**   - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
   !! - **thisobs%use_global_obs** - Whether to use global observations or restrict the observations
   !!                               to the relevant ones (default: *.true.* i.e use global full observations)
   !!
   !! Further variables are set when the routine PDAFomi_gather_obs is called.
   !!
   subroutine init_dim_obs_sic(step, dim_obs)
      use PDAF, only: PDAFomi_gather_obs
      use statevector_pdaf, only: sfields
      use nemo_pdaf, only: nwet, wet_pts, glamt, gphit, i0, j0, ni_p, nj_p
      use io_pdaf, only: check
      use, intrinsic :: ieee_arithmetic ! to detect Nans

      integer, intent(in)    :: step    !< Current time step
      integer, intent(inout) :: dim_obs !< Dimension of full observation vector

      integer :: i, j, i_obs, i_w             !> Counters
      integer :: dim_obs_p                     !> Number of process-local observations
      integer :: ncid, id_var
      integer :: pos(3), cnt(3)
      real(pwp), allocatable :: obs_p(:)       !> PE-local observation vector
      real(pwp), allocatable :: ivar_obs_p(:)  !> PE-local inverse observation error variance
      real(pwp), allocatable :: ocoord_p(:, :) !> PE-local observation coordinates
      real(pwp), allocatable :: obs_global(:,:) !< Global obs field read from file
      integer :: ndims, natts, xtype
      integer :: dimids(NF90_MAX_VAR_DIMS)
      integer :: ik, dimlen
      character(len=NF90_MAX_NAME) :: varname, dimname
      logical :: maytosept = .true.

      real(pwp) :: rad_conv = 3.141592653589793/180.0 !> Degree to radian conversion

      ! *****************************
      ! *** Global setting config ***
      ! *****************************
      if (mype_filter == 0) &
         write (*, '(a,4x,a)') 'NEMO-PDAF', 'Assimilate observations - obs_sic'
      ! Store whether to assimilate this observation type (used in routines
      ! below)
      if (assim_sic) thisobs%doassim = 1
      ! Specify type of distance computation
      thisobs%disttype = 3   ! 3=Haversine
      ! Number of coordinates used for distance computation.
      ! The distance compution starts from the first row
      thisobs%ncoord = 2
      ! Set to use limited full observations
      thisobs%use_global_obs = use_global_obs
      ! Omit large innovation to omit nwet = 0 points
      thisobs%inno_omit = 3!1e6    ! nwet= 0 points will have obs = 0? but no corresponding model state(?)

      ! **********************************
      ! *** Read observations from file ***
      ! **********************************
      if (mype_filter == 0) then
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- obs file: ', trim(file_sic)
         write (*, '(a,4x,a,a)')  'NEMO-PDAF', '--- variable: ', trim(varname_sic)
         write (*, '(a,4x,a,i0)') 'NEMO-PDAF', '--- reading nc_step: ', nc_step
      end if

      allocate(obs_global(ni_p, nj_p)) !180 X 148
      obs_global = -1
      call check( nf90_open(trim(file_sic), NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, trim(varname_sic), id_var) )

      call check( nf90_inquire_variable(ncid, id_var, name=varname, xtype=xtype, &
                                 ndims=ndims, dimids=dimids, natts=natts) )

      print *, 'Variable:', trim(varname)
      print *, 'Number of dimensions:', ndims

      ! Loop over dimensions in the order Fortran sees them to check same as ens
      do ik = 1, ndims
         call check( nf90_inquire_dimension(ncid, dimids(ik), name=dimname, len=dimlen) )
         print *, 'Dim', ik, ':', trim(dimname), ' size =', dimlen
      end do

      pos = (/i0, j0, nc_step/)
      cnt = (/ni_p, nj_p, 1/)
      call check( nf90_get_var(ncid, id_var, obs_global, start=pos, count=cnt) )
      call check( nf90_close(ncid))

      where (ieee_is_nan(obs_global))
          obs_global = -100.0            !set nan values to negative no., then filter out
      end where

      obs_global = obs_global / 100.0          !divide by 10 to convert from % to fraction
      ! ***********************************************************
      ! *** Count available observations for the process domain ***
      ! *** and initialize index and coordinate arrays.         ***
      ! ***********************************************************
      ! Set number of local observations
      dim_obs_p = 0
      do i_w = 1, nwet
         i = wet_pts(6,i_w)
         j = wet_pts(7,i_w)
         ! skip negative values (fill values)
         if (obs_global(i, j) < 0.0) cycle
         dim_obs_p = dim_obs_p + 1
      end do

      ! set observation vector
      if (dim_obs_p > 0) then
         ! Vector of observations on the process sub-domain
         allocate (obs_p(dim_obs_p))
         ! Coordinate array of observations on the process sub-domain
         allocate (ocoord_p(2, dim_obs_p))
         ! Coordinate array for observation operator
         allocate (thisobs%id_obs_p(1, dim_obs_p))
         allocate (ivar_obs_p(dim_obs_p))
         i_obs = 0
         do i_w = 1, nwet
            i = wet_pts(6,i_w)
            j = wet_pts(7,i_w)
            if (obs_global(i, j) < 0.0_pwp) cycle

            i_obs = i_obs + 1
            obs_p(i_obs) = obs_global(i, j)
            ocoord_p(1,i_obs) = glamt(i,j)*rad_conv
            ocoord_p(2,i_obs) = gphit(i,j)*rad_conv
            thisobs%id_obs_p(1,i_obs) = i_w + sfields(id_sfields)%off
         end do
         ! specify observation error variances
         if (maytosept) then
            do i=  1, dim_obs_p
               !ivar_obs_p(i) = 1.0/(0.2*obs_p(i)*0.2*obs_p(i))
               ivar_obs_p(i) = 1.0/(0.2*0.2)!(0.2*0.2)
               if (obs_p(i) == 0.0) ivar_obs_p(i) = 1.0/(0.2*0.02)
            enddo
         else
         ! if oct - apr inc
            do i=1,dim_obs_p
               if (obs_p(i) > 0.8) then
                  ivar_obs_p(i) = 1.0/(0.1*0.1)
                  !write(*,*) ivar_obs_p(i), 'std:', 0.1*obs_p(i)
               else
                  ivar_obs_p(i) = 1.0/(0.15*0.15)
                  !write(*,*) ivar_obs_p(i), 'std::', 0.15*obs_p(i)
               endif
            enddo
         endif
      else
         ! Set empty arrays if no local observations
         allocate (obs_p(1))
         obs_p(1) = -1.0! so that gets LARGE obs error applied thisobs%inno_omit
         allocate (ocoord_p(2, 1))
         allocate (thisobs%id_obs_p(1, 1))
         allocate (ivar_obs_p(1))
      end if

      deallocate(obs_global)

      ! ****************************************
      ! *** Gather global observation arrays ***
      ! ****************************************
      call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
                              thisobs%ncoord, lradius, dim_obs)

      ! ********************
      ! *** Finishing up ***
      ! ********************
      ! Deallocate all local arrays
      deallocate (obs_p, ocoord_p, ivar_obs_p)
   end subroutine init_dim_obs_sic

   !>###Implementation of observation operator from state vector
   !>
   !>This routine applies the full observation operator
   !>for the sic observations.
   !>
   !>The routine is called by all filter processes.
   !>
   subroutine obs_op_sic(dim_p, dim_obs_p, state_p, ostate)
      !> PE-local state dimension
      integer, intent(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      integer, intent(in) :: dim_obs_p
      !> PE-local model state
      real(pwp), intent(in) :: state_p(dim_p)
      !> Full observed state
      real(pwp), intent(inout) :: ostate(dim_obs_p)

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      !Just need to pick out the grid cells and then sum over cats
      if (thisobs%doassim == 1) then
         if (do_obs_op_from_model_output) then
            call obs_op_sic_output(dim_p, dim_obs_p, state_p, ostate)
         else
            call obs_op_sic_state(dim_p, dim_obs_p, state_p, ostate)
         end if
      end if
   end subroutine obs_op_sic

   !>###Implementation of observation operator from state vector
   !>
   !>This routine applies the full observation operator
   !>for the sic observations.
   !>
   !>The routine is called by all filter processes.
   !>
   subroutine obs_op_sic_state(dim_p, dim_obs_p, state_p, ostate)
      use statevector_pdaf, only: n_fields, sfields
      use PDAF, only: PDAFomi_obs_op_extern
      use, intrinsic :: ieee_arithmetic
      !> PE-local state dimension
      integer, intent(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      integer, intent(in) :: dim_obs_p
      !> PE-local model state
      real(pwp), intent(in) :: state_p(dim_p)
      !> Full observed state
      real(pwp), intent(inout) :: ostate(dim_obs_p)
      ! local variables
      integer :: i, j, i_off
      real(pwp), allocatable :: ostate_thisobs(:) !local part of state vector

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      !Just need to pick out the grid cells and then sum over cats
      allocate(ostate_thisobs(thisobs%dim_obs_p))
      ostate_thisobs = 0.
      do j = 1, thisobs%dim_obs_p
         i_off = thisobs%id_obs_p(1,j) - sfields(id_sfields)%off
         do i = 1, n_fields
            if (sfields(i)%name_rest_n == 'a_i') then
               ostate_thisobs(j) = ostate_thisobs(j) + state_p(i_off + sfields(i)%off)
            endif
         end do
      end do
      CALL PDAFomi_obs_op_extern(thisobs, ostate_thisobs, ostate)
      deallocate(ostate_thisobs)
   end subroutine obs_op_sic_state

   !>###Implementation of observation operator from model output
   !>
   !>This routine reads model output from 14 / 15 day average
   !>from model output consistent with observation data.
   !>
   !>The routine is called by all filter processes.
   !>
   subroutine obs_op_sic_output(dim_p, dim_obs_p, state_p, ostate)
      use statevector_pdaf, only: n_fields, sfields
      use PDAF, only: PDAFomi_obs_op_extern, PDAF_get_obsmemberid
      use nemo_pdaf, only: dim_2d_p
      use io_pdaf, only: read_model_output
      use config_pdaf, only: dim_ens
      use, intrinsic :: ieee_arithmetic
      !> PE-local state dimension
      integer, intent(in) :: dim_p
      !> Dimension of full observed state (all observed fields)
      integer, intent(in) :: dim_obs_p
      !> PE-local model state
      real(pwp), intent(in) :: state_p(dim_p)
      !> Full observed state
      real(pwp), intent(inout) :: ostate(dim_obs_p)
      ! *** local variables ***
      integer :: j, i_off
      integer :: ens_member ! index of ensemble member
      real(pwp), allocatable :: state_p_local(:) !local part of state vector
      real(pwp), allocatable :: ostate_thisobs(:) ! local part of observed state vector
                                            ! for this observation type

      ! ******************************************************
      ! *** Apply observation operator H on a state vector ***
      ! ******************************************************
      allocate(state_p_local(dim_2d_p))
      allocate(ostate_thisobs(thisobs%dim_obs_p))
      ! get ensemble member
      call PDAF_get_obsmemberid(ens_member)
      ! construct
      if (ens_member /= 0) then
         write(*,*) 'NEMO-PDAF --- obs_op_sic_output: ens_member=', ens_member
         call read_model_output(path_op_root, ens_prefix, path_op_suffix, &
            ens_member, varname_obs_op, op_time_index, 0, '', 0, state_p_local)
            do j = 1, thisobs%dim_obs_p
               i_off = thisobs%id_obs_p(1,j) - sfields(id_sfields)%off
               ostate_thisobs(j) = state_p_local(i_off)
            end do
      else
         write(*,*) 'NEMO-PDAF --- obs_op_sic_output: getting ensemble mean'
         ostate_thisobs = 0.0
         do ens_member = 1, dim_ens
            call read_model_output(path_op_root, ens_prefix, path_op_suffix, &
               ens_member, varname_obs_op, op_time_index, 0, '', 0, state_p_local)
               do j = 1, thisobs%dim_obs_p
                  i_off = thisobs%id_obs_p(1,j) - sfields(id_sfields)%off
               ostate_thisobs(j) = ostate_thisobs(j) + state_p_local(i_off)/dim_ens
            end do
         end do
      end if
      CALL PDAFomi_obs_op_extern(thisobs, ostate_thisobs, ostate)
      deallocate(state_p_local, ostate_thisobs)
   end subroutine obs_op_sic_output


   subroutine init_dim_obs_l_sic(domain_p, step, dim_obs, dim_obs_l)

      use PDAF, only: PDAFomi_init_dim_obs_l
      use assimilation_pdaf, only: domain_coords
      use config_pdaf, only: locweight

      !> Index of current local analysis domain
      integer, intent(in)  :: domain_p
      !> Current time step
      integer, intent(in)  :: step
      !> Full dimension of observation vector
      integer, intent(in)  :: dim_obs
      !> Local dimension of observation vector
      integer, intent(out) :: dim_obs_l

      ! **********************************************
      ! *** Initialize local observation dimension ***
      ! **********************************************

      call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
           locweight, lradius, sradius, dim_obs_l)

   end subroutine init_dim_obs_l_sic

 end module obs_sic
