#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam.F90
!> \brief Main routines of JAM

!> \file jam.F90
!> This file includes the main routines of JAM, the atmospheric module
!! of GLOBE, which handles the climate forcing and initializes the
!! pre-defined calendar and time step of GLOBE.

!     ******************************************************************
!     JAM_INIT
!     ******************************************************************

!> \brief Main initialization routine

!> \detail This routine calls all the steps needed for the
!! initialization of the JAM module.

!     ******************************************************************

      subroutine jam_init ()
      use jam_mod
      implicit none

      call globe_open_diag(sfile_diag, kfile_diag)

      __diag(kfile_diag,'jam_init')

!     ------------------------------------------------------------------
!     * read namelist
!     ------------------------------------------------------------------

      call jam_read_namelist

!     ------------------------------------------------------------------
!     * initialize input files
!     ------------------------------------------------------------------

      call jam_input_init

!     ------------------------------------------------------------------
!     * allocate all fields
!     ------------------------------------------------------------------

      call jam_alloc

!     ------------------------------------------------------------------
!     * initialize variables
!     ------------------------------------------------------------------

      if (ksingle .ne. 0 .and. kRandYear .eq. 0) then
        call jam_spt_preread_climate
      endif

      if (kRandYear .eq. 1) call RANDOM_SEED()

      __diag(kfile_diag,'jam_init: end')

      return
      end subroutine jam_init

!     ******************************************************************
!     JAM_STEP_S
!     ******************************************************************

!> \brief Calculations at every time step

!     ******************************************************************

      subroutine jam_step_s ()
      use jam_mod
      implicit none

!     nothing here for now

      return
      end subroutine jam_step_s

!     ******************************************************************
!     JAM_STEP_M
!     ******************************************************************

!> \brief Main routine to read climate forcing at every monthly step

!> \detail This routine manages the monthly calendar and calls the
!! reading of the appropriate climate forcing data from files. The
!! monthly data are temporary saved in memory to speed up the model. The
!! routine atmos_step() (jam_globe.F90) is used to extract the data at
!! every time step.\n

!> \b Variables \b modified: nleapday, kdpy

!     ******************************************************************

      subroutine jam_step_m ()
      use jam_mod
      implicit none

      integer :: zdatayear         ! actual year in the data set
      integer :: zryears           ! rest of data years by mod with actual year
      integer :: zlyr              ! indicator for leap year

!     * check for leap year

      zryears = MOD(kyear - 1, n_data_years) + 1
      zdatayear = zryears + k_first_dyear - 1
      zlyr = 0
      if ((kdpm(13) .ne. 0) .and. (MOD(zdatayear - kdpm(13),4) .eq. 0)   &
     &  .and. (zdatayear .ne. kdpm(14))) zlyr = 1

!     * set leap day for actual month

      nleapday = 0
      if ((zlyr .eq. 1) .and. (kmonth .eq. 2)) nleapday = 1

!     * setup the number of days of the actual year

      if (kmonth .eq. 1) kdpy = SUM(kdpm(1:12)) + zlyr

!     * read and scatter the data sets for this month

      call jam_read_climate_month

      return
      end subroutine jam_step_m

!     ******************************************************************
!     JAM_STOP
!     ******************************************************************

!> \brief Finish of JAM

!> \detail This routine finishes JAM at the end of the model run,
!! including closing open file handlers, and
!! freeing up used memory.

!     ******************************************************************

      subroutine jam_stop ()
      use jam_mod
      implicit none

      __diag(kfile_diag,'jam_stop')

      call jam_input_stop

      call jam_dealloc

      __diag(kfile_diag,'jam_stop: done')

!     * close diagnostic output file

      call globe_close_diag(kfile_diag)

      return
      end subroutine jam_stop
