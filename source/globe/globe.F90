#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe.F90
!> \brief Main routines of GLOBE

!> \file globe.F90
!> This file includes the main program which calls the included routines
!! of initialization and the main loop. The main loop calls the included
!! routines which execute the modules of the model and the model output.

!     ******************************************************************
!
!     G L O B E
!
!     ******************************************************************

!> \brief The main program

!> \detail The main program which initializes GLOBE by calling
!! globe_init(). Then,
!! it starts the main loop by calling globe_loop().
!! After the main loop is finished, globe_stop() is
!! called and the main program ends.\n

!> GLOBE accepts the option \b -V, to show how it was compiled.

!     ******************************************************************

      program globe
      use globe_mod
      implicit none

      character(len=32) :: arg
      integer :: i

!     * check for the -V option (output of the version string)

      do i = 1, command_argument_count()
        call get_command_argument(i, arg)

        select case (arg)
          case ('-v', '--version')
            call globe_show_version
            stop
          case default
            print '(a,a,/)', 'Unrecognized command-line option: ', arg
            stop
        end select
      end do

!     * initialization
      call globe_init

!     * main loop
      call globe_loop

!     * stop model
      call globe_stop

      end

!     ******************************************************************
!     GLOBE_LOOP
!     ******************************************************************

!> \brief The standard main loop

!> \detail This routine sets the time variables and calls the
!! corresponding time step subroutines. At the end of every month it
!! calls the diagnostic module and the output routine.\n

!> \b Variables \b set: kglobe_year, kglobe_month, kglobe_day, kglobe_ts

!     ******************************************************************

      subroutine globe_loop ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_loop: main loop')

!     * yearly loop
      do while ((kglobe_year + kglobe_firstyear - 1) .le. kglobe_lastyear)

        if (kglobe_status == 1) then
          __diag_num(kglobe_fdiag,'globe_loop: yearly loop ',kglobe_year+kglobe_firstyear-1)
        endif

!       * run pre-yearly steps
        call globe_step(GLOBE_STEP_PREYEAR)

!       * monthly loop
        do kglobe_month = 1, 12

!         * run pre-monthly steps
          call globe_step(GLOBE_STEP_PREMONTH)

          call globe_status

!         * daily loop
          do kglobe_day = 1, kglobe_dpm(kglobe_month) + kglobe_leapday

!           * run pre-daily steps
            call globe_step(GLOBE_STEP_PREDAY)

!           * timestep loop
            do kglobe_ts = 1, kglobe_tspd
              call globe_step(GLOBE_STEP_TIMESTEP)
            enddo

!           * run post-daily steps
            call globe_step(GLOBE_STEP_DAY)

          enddo

!         * run post-monthly steps
          call globe_step(GLOBE_STEP_MONTH)

!         * run diagnostic modules
          call globe_diag

!         * write out monthly averages
          call globe_output

        enddo

!       * run post-yearly steps
        call globe_step(GLOBE_STEP_YEAR)

!       * yearly loop variable
        kglobe_year = kglobe_year + 1
      enddo

      return
      end subroutine globe_loop

!     ******************************************************************
!     GLOBE_INIT
!     ******************************************************************

!> \brief Main initialization routine

!> \detail This main initialization routine starts MPI, calls the
!! reading of the GLOBE namelist, calls the initialization of the
!! surface parameters, calls the setup of fields for all CPUs, calls the
!! allocation of variables and calls the initialization steps of all
!! models used.\n

!     ******************************************************************

      subroutine globe_init ()
      use globe_mod
      implicit none

!     * function declaration
      integer      :: unlink

!     * variable definitions
      logical      :: fexist
      integer      :: ioerr
      character*40 :: sfile

!     * start mpi
      call globe_mpstart

      if (kglobe_mypid == kglobe_kroot) then
        write(*,*) '--------------------------------------------------'
        write(*,*) '                    G L O B E                     '
        write(*,*) '--------------------------------------------------'
        write(*,*)
      endif

!     * namelist parameters
      call globe_read_namelist

!     * get list of grid points

      if (kglobe_restart .gt. 0) then
        call globe_read_restart
      else
        inquire(FILE=sglobe_fgridpoint, EXIST=fexist)
        if (fexist) then
          call globe_read_landpoint
        else

!         * read surface parameter file (if not existing, create it)

          inquire(FILE=sglobe_fsurface, EXIST=fexist)
          if (.not. fexist) then
            call globe_create_surface_file
          endif
!         -> needs exclusive file access
          call globe_read_surface_file

        endif

!       * set up fields to the number of CPUs used
        call globe_npro_fields

!       * allocate state and flux fields
        call globe_alloc

!       * initialize fields
        call globe_fields_init
      endif

      __diag_num(kglobe_fdiag,'globe_init: nglobe_landpts = ',nglobe_landpts)

!     * initialize modules

      __diag(kglobe_fdiag,'globe_init: atmos_init')
      call atmos_init

      __diag(kglobe_fdiag,'globe_init: soil_init')
      call soil_init

      __diag(kglobe_fdiag,'globe_init: vegetation_init')
      call vegetation_init

!     * end
      __diag(kglobe_fdiag,'globe_init: done')

      return
      end subroutine globe_init

!     ******************************************************************
!     GLOBE_STEP
!     ******************************************************************

!> \brief Calls the time step routines of all models used

!> \detail This routine calls the time step routines of all models used.
!! The parameter \c globe_step_id defines the called time step.
!! At certain time steps additional actions are triggered, like the
!!  modification of some counters.\n

!> \b Variables \b modified: kglobe_diy, kglobe_dim, kglobe_timestep

!> \param globe_step_id Time step marker

!     ******************************************************************

      subroutine globe_step (globe_step_id)
      use globe_mod
      implicit none

      integer, INTENT(in) :: globe_step_id

!     ==================================================================
!     GLOBE_STEP_PREYEAR
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREYEAR) then

!       * initialize day in year

        kglobe_diy = 1

!       * run pre-yearly model steps

        call atmos_step(GLOBE_STEP_PREYEAR)
        call soil_step(GLOBE_STEP_PREYEAR)
        call vegetation_step(GLOBE_STEP_PREYEAR)

      endif

!     ==================================================================
!     GLOBE_STEP_PREMONTH
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREMONTH) then

!       * run pre-monthly model steps

!       * atmos_step needs to run first (date handling)
        call atmos_step(GLOBE_STEP_PREMONTH)
        call soil_step(GLOBE_STEP_PREMONTH)
        call vegetation_step(GLOBE_STEP_PREMONTH)

      endif

!     ==================================================================
!     GLOBE_STEP_PREDAY
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREDAY) then

!       * run pre-daily model steps

        call atmos_step(GLOBE_STEP_PREDAY)
        call soil_step(GLOBE_STEP_PREDAY)
        call vegetation_step(GLOBE_STEP_PREDAY)

      endif

!     ==================================================================
!     GLOBE_STEP_TIMESTEP
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_TIMESTEP) then

!       * run single model steps

        call atmos_step(GLOBE_STEP_TIMESTEP)
        call soil_step(GLOBE_STEP_TIMESTEP)
        call vegetation_step(GLOBE_STEP_TIMESTEP)
!       call globe_vegetation_flux(GLOBE_STEP_TIMESTEP)

!       * set counter

        kglobe_timestep = kglobe_timestep + 1

      endif

!     ==================================================================
!     GLOBE_STEP_DAY
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_DAY) then

!       * run post-daily model steps

        call atmos_step(GLOBE_STEP_DAY)
        call soil_step(GLOBE_STEP_DAY)
        call vegetation_step(GLOBE_STEP_DAY)

!       * set counter

        kglobe_diy = kglobe_diy + 1
        kglobe_dim = kglobe_day

      endif

!     ==================================================================
!     GLOBE_STEP_MONTH
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_MONTH) then

!       * run post-monthly model steps

        call atmos_step(GLOBE_STEP_MONTH)
        call soil_step(GLOBE_STEP_MONTH)
        call vegetation_step(GLOBE_STEP_MONTH)

      endif

!     ==================================================================
!     GLOBE_STEP_YEAR
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_YEAR) then

!       * run post-yearly model steps

        call atmos_step(GLOBE_STEP_YEAR)
        call soil_step(GLOBE_STEP_YEAR)
        call vegetation_step(GLOBE_STEP_YEAR)

      endif

      return
      end subroutine globe_step

!     ******************************************************************
!     GLOBE_STOP
!     ******************************************************************

!> \brief Calls the post-execution steps

!> \detail This routine calls the post-execution steps of each model,
!! creates the restart file of GLOBE, stops MPI,
!! deallocates all GLOBE variables, and closes the diagnostic output.

!     ******************************************************************

      subroutine globe_stop ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_stop')

!     * stop models

      call atmos_stop
      call soil_stop
      call vegetation_stop

!     * write GLOBE restart

      if (kglobe_restart .ge. 0) call globe_write_restart

      __diag(kglobe_fdiag,'globe_stop: globe_mpstop')

!     * stop mpi

      call globe_mpstop

!     * deallocation of all fields to clean free the memory

      call globe_dealloc

!     * end

      __diag(kglobe_fdiag,'globe_stop: done')

!     * close the diagnostic output file

      call globe_close_diag(kglobe_fdiag)

      return
      end subroutine globe_stop

!     ******************************************************************
!     GLOBE_OUTPUT
!     ******************************************************************

!> \brief Calls the model output

!> \detail This routine calls the output routines of each model.\n

!> \b Variables \b modified: kglobe_day

!     ******************************************************************

      subroutine globe_output ()
      use globe_mod
      implicit none

!     * set variable kglobe_day to the last loop day of the last month
      kglobe_day = kglobe_dim

!     * run model output
      call atmos_output
      call soil_output
      call vegetation_output

      return
      end subroutine globe_output

!     ******************************************************************
!     GLOBE_DIAG
!     ******************************************************************

!> \brief Calls the model diagnostics

!> \detail This routine calls the diagnostic routines of each model.

!     ******************************************************************

      subroutine globe_diag ()
      use globe_mod
      implicit none

!     * create additional diagnostic output
      call atmos_diag
      call soil_diag
      call vegetation_diag

      return
      end subroutine globe_diag
