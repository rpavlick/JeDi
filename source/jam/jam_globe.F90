#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam_globe.F90
!> \brief GLOBE interface of JAM

!> \file jam_globe.F90
!> This file includes the interface between JAM and GLOBE. All data
!! exchange is done by these routines, and all subroutine calls from
!! GLOBE are routed there.

!     ==================================================================
!     ------------------------------------------------------------------
!     GLOBE INTERFACE FOR JAM
!     ------------------------------------------------------------------
!     ==================================================================

!     ******************************************************************
!     ATMOS_INIT
!     ******************************************************************

!> \brief Initialization of JAM with GLOBE variables

!> \detail This routine initializes some variables with values from
!! GLOBE, calls jam_init, and sends some calendar variables to GLOBE.\n

!> \b JAM \b variables \b set: krestart, kdiag, ksingle,
!! nCPts, nLPts, nLPts2, kyear, nlon, nlat, kroot, kmypid,
!! kLPID, pdt\n

!> \b GLOBE \b variables \b set: kglobe_spts, kglobe_tspd, kglobe_dpm

!     ******************************************************************

      subroutine atmos_init ()
      use globe_mod
      use jam_mod
      implicit none

!     * copy GLOBE variables into JAM

      krestart           = kglobe_restart
      kdiag              = kglobe_diag
      ksingle            = kglobe_spt

      nCPts              = nglobe_cgpt
      nLPts              = nglobe_landpts
      nLPts2             = nglobe_gpts2

      kyear              = kglobe_year

      nlon               = nglobe_nlon
      nlat               = nglobe_nlat

!     * copy MPI variables to JAM

      kroot              = kglobe_kroot
      kmypid             = kglobe_mypid

!     * allocate and copy grid list

      if (kmypid == kroot) then
        allocate(kLPID(nLPts))
        kLPID(:)         = kglobe_gLPIDs(:)
      endif

!     * do initialization

      call jam_init

!     * copy variables back to globe

      pdt = pdt / ktspd ! seconds per time step
      kglobe_spts        = pdt

      kglobe_tspd        = ktspd
      kglobe_dpm(:)      = kdpm(:)

      return
      end subroutine atmos_init

!     ******************************************************************
!     ATMOS_STEP
!     ******************************************************************

!> \brief Transfer of variables at the different time steps

!> \detail This routine handles the different time steps called by
!! GLOBE.\n

!> At the \a pre-year \a time \a step, the optional randomization of the
!! year of climate forcing, and the CO<sub>2</sub> data input from a
!! text file is called.\n

!> At the \a pre-month \a time \a step, some calendar variables are
!! transferred between GLOBE and JAM, and the monthly step routine of JAM
!! is called.\n

!> \b Variables \b set: kdpm, kmonth, kyear, kglobe_leapday, kglobe_dpy\n

!> At the \a pre-day \a time \a step, the actual day is set from GLOBE.
!> \b Variables \b set: kday\n

!> At the call of the \a actual \a time \a step, some calendar variables
!! and water fluxes are transfered to JAM, the time step routine of JAM
!! is called, the daily climate forcing is transferred to GLOBE, and some
!! related atmospheric GLOBE variables are calculated.\n

!> \b JAM \b variables \b set:
!> kts, ktsim\n

!> \b Climate \b variables \b set: xglobe_T_a, fglobe_RADs_as,
!! fglobe_RADl_as, fglobe_H2Ol_as, xglobe_H2Og_a,
!> xglobe_CO2g_a\n

!> \b GLOBE \b variables \b modified: fglobe_H2Ol_as, fglobe_H2Os_as,
!! fglobe_CO2d_as

!> \param globe_step_id Time step marker

!     ******************************************************************

      subroutine atmos_step (globe_step_id)
      use globe_mod
      use jam_mod
      implicit none

      integer, INTENT(in) :: globe_step_id

      integer :: i
      real    :: zcCO2d_sat
      real    :: tmp_clim

      integer :: flag1, flag2, flag3

!     ==================================================================
!     GLOBE_STEP_PREYEAR
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREYEAR) then

        if (kmypid == kroot) then

          if (kRandYear .eq. 1) then

!           * rewind all file pointers

            if (input_pTEMP)     call globe_rewind(kfile_temp2m)
            if (input_fRADs_as)  call globe_rewind(kfile_solrad)
            if (input_fRADl_as)  call globe_rewind(kfile_terrad)
            if (input_fH2Ol_as)  call globe_rewind(kfile_precip)
            if (input_cH2Og_a)   call globe_rewind(kfile_relhum)
            if (input_xCO2g_a)   call globe_rewind(kfile_co2a)

!           * seek the climate files to a random year

            call jam_seek_clim_file

          endif

        endif

      endif

!     ==================================================================
!     GLOBE_STEP_PREMONTH
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREMONTH) then

!       * copy variables from GLOBE

        kdpm(:)            = kglobe_dpm(:)
        kmonth             = kglobe_month
        kyear              = kglobe_year

!       * do monthly step

        call jam_step_m

!       * read yearly time-varying CO2 concentrations from text file

        if (kCO2_dyn == 1) then
          call jam_step_carbon_y(kglobe_year + kglobe_firstyear - 1)
        endif

!       * copy variables back to GLOBE

        if (kRandYear .eq. 0) kglobe_leapday = nleapday
        kglobe_dpy         = kdpy

      endif

!     ==================================================================
!     GLOBE_STEP_PREDAY
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_PREDAY) then

        kday               = kglobe_day

      endif

!     ==================================================================
!     GLOBE_STEP_TIMESTEP
!     ==================================================================

      if (globe_step_id .eq. GLOBE_STEP_TIMESTEP) then

!       ----------------------------------------------------------------
!       * copy variables from GLOBE
!       ----------------------------------------------------------------

        kts                = kglobe_ts
        ktsim              = (kday - 1) * ktspd + kts

!       ------------------------------------------------------------------
!       * do daily step
!       ------------------------------------------------------------------

!        call jam_step_s

!       ------------------------------------------------------------------
!       * copy variables back to GLOBE arrays
!       ------------------------------------------------------------------

!       * energy variables

        if (input_pTEMP) then
          xglobe_T_a(:) = pTEMP(:,ktsim)
        else
            xglobe_T_a(:) = pTEMP_init
        endif
        if (input_fRADs_as) then
          fglobe_RADs_as(:) = fRADs_as(:,ktsim)
        else
          fglobe_RADs_as(:) = fRADs_as_init
        endif
        if (input_fRADl_as) then
          fglobe_RADl_as(:) = fRADl_as(:,ktsim)
        else
          fglobe_RADl_as(:) = fRADl_as_init
        endif

!       * water variables

        if (input_fH2Ol_as) then
          fglobe_H2Ol_as(:) = fH2Ol_as(:,ktsim)
        else
          fglobe_H2Ol_as(:) = fH2Ol_as_init
        endif
        if (input_cH2Og_a) then
          xglobe_H2Og_a(:) = cH2Og_a(:,ktsim)
        else
          xglobe_H2Og_a(:) = cH2Og_a_init
        endif

!       * carbon variable

        if (input_xCO2g_a) then
          xglobe_CO2g_a(:) = xCO2g_a(:,ktsim)
        else
          xglobe_CO2g_a(:) = xCO2g_a_init
        endif

!       * calculate snowfall from precipitation rate and temperature

        fglobe_H2Ol_as(:) = max(0.0, fglobe_H2Ol_as(:))
        fglobe_H2Os_as(:) = max(0.0, (3.3 - (xglobe_T_a(:) - cTMelt)) / 4.4)
        fglobe_H2Os_as(:) = min(1.0, fglobe_H2Os_as(:)) * fglobe_H2Ol_as(:)
        fglobe_H2Ol_as(:) = fglobe_H2Ol_as(:) - fglobe_H2Os_as(:)

      endif

      return
      end subroutine atmos_step

!     ******************************************************************
!     ATMOS_STOP
!     ******************************************************************

!> \brief Call to finish JAM

!> \detail This routine calls a subroutine to finish JAM.

!     ******************************************************************

      subroutine atmos_stop ()
      implicit none

      call jam_stop

      return
      end subroutine atmos_stop

!     ******************************************************************
!     ATMOS_OUTPUT
!     ******************************************************************

!> \brief Call of JAM output handling

!> \detail This routine calls the JAM output handling subroutine.

!     ******************************************************************

      subroutine atmos_output ()
      use globe_mod
      use jam_mod
      implicit none

      return
      end subroutine atmos_output

!     ******************************************************************
!     ATMOS_DIAG
!     ******************************************************************

!> \brief Call of JAM diagnostics

!> \detail This routine calls the JAM diagnostics.

!     ******************************************************************

      subroutine atmos_diag ()
      implicit none

!      call jam_diag

      return
      end subroutine atmos_diag
