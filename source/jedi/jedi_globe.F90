#define __ACTIVATE
#include "../globe/globe_macros.f90"

!     ------------------------------------------------------------------
!     GLOBE INTERFACE FOR JEDI
!     ------------------------------------------------------------------

!     ******************************************************************
!     VEGETATION_INIT
!     ******************************************************************

      subroutine vegetation_init ()
      use globe_mod
      use jedi_mod
      implicit none

      if (kglobe_mypid == kglobe_kroot) write(kglobe_fdiag,*) 'vegetation_init'

!     copy GLOBE flags to JEDI

      nonepoint = kglobe_spt
      krestart  = kglobe_restart
      kdiag     = kglobe_diag

!     * copy GLOBE dimensions to JEDI

      kNumGPts  = nglobe_landpts
      nlon      = nglobe_nlon
      nlat      = nglobe_nlat

!     * copy MPI variables to jedi

      kroot     = kglobe_kroot
      kmypid    = kglobe_mypid

      NHOR      = nglobe_cgpt
      gNHOR     = nglobe_gpts2

!     * copy grid

      if (kmypid == kroot) then
        allocate(gkGPID(nglobe_landpts))
        gkGPID(:) = kglobe_gLPIDs(:)
      endif

!     * scatter grid list from GLOBE to jedi

      allocate (pjedi_lat(nglobe_cgpt))
      __globe_mpsc_from_to(pglobe_gLP_lat2,pjedi_lat)

!     * scatter area from GLOBE to jedi

      allocate (pjedi_area(nglobe_cgpt))
      __globe_mpsc_from_to(pglobe_gLP_area2,pjedi_area)

!     * soil texture

      allocate(pjedi_paw(nglobe_cgpt))
      __globe_mpsc_from_to(pglobe_gLP_paw2,pjedi_paw)

!     * initialize jedi

      call jedi_init

!     * land-index of nearest neighbors

      if (kmypid == kroot) klpt_nn(:,:) = kglobe_lpt_nn(:,:)

      return
      end subroutine vegetation_init

!     ******************************************************************
!     VEGETATION_STEP
!     ******************************************************************

      subroutine vegetation_step (globe_step_id)
      use globe_mod
      use jedi_mod
      implicit none

      integer :: globe_step_id
      integer :: flag1, flag2, flag3, flag4
      integer :: flag5, flag6, flag7, flag8
      real    :: save1, save2, save3, save4, save5
      integer :: save6

      select case(globe_step_id)

!     ==================================================================
!     GLOBE_STEP_PREDAY
!     ==================================================================

!      case (GLOBE_STEP_PREDAY)

!     ==================================================================
!     GLOBE_STEP_PREYEAR
!     ==================================================================

      case (GLOBE_STEP_PREYEAR)

        YEAR     = kglobe_year
        OUTYEAR  = kglobe_year + kglobe_firstyear - 1

!       * optimizing version

        if (kOpti == 1) then
          if (mod(kglobe_year, kOptiYears) == 0) then
            call jedi_opti_filter
          endif
        endif

!       * set flag for speciation running

        if (kSpec_dyn .eq. 1) then
          if ( ((YEAR .gt. 50) .or. (kspinup .eq. 0)) .and.            &
     &         (OUTYEAR .lt. (kglobe_lastyear - 10)) ) then
            jedi_specruns = 1
          endif
        endif

!       * set flag for migration running

        if (kMig_dyn .eq. 1) then
          if ( ((YEAR .gt. 20) .or. (kspinup .eq. 0)) .and.            &
     &         (OUTYEAR .lt. (kglobe_lastyear - 10)) )then
            jedi_migruns = 1
          endif
        endif

!     ==================================================================
!     GLOBE_STEP_DAY
!     ==================================================================

      case (GLOBE_STEP_DAY)

!       * copy GLOBE date variables into jedi

        MONTH = kglobe_month
        kdiy  = kglobe_diy
        kdpy  = kglobe_dpy

!       * copy GLOBE climate variables into jedi

        pjedi_temp(:)   = xglobe_T_a(:) - cT0  ! K -> degree C
        pjedi_rain(:)   = fglobe_H2Ol_as(:) * 86400.0  ! mm/s -> mm/day
        pjedi_snow(:)   = fglobe_H2Os_as(:) * 86400.0  ! mm/s -> mm/day
        pjedi_swdown(:) = fglobe_RADs_as(:) ! W/m2
        pjedi_lwnet(:)  = fglobe_RADl_as(:) ! W/m2

!       * copy GLOBE CO2 concentration into jedi

        pPCO2 = xglobe_CO2g_a(1)

!       * step jedi

        call jedi_step

!     ==================================================================
!     GLOBE_STEP_MONTH
!     ==================================================================

      case (GLOBE_STEP_MONTH)

        if (kSpec_dyn .eq. 1) call jedi_extinction

        if (jedi_migruns .eq. 1) call jedi_migration

        if (jedi_specruns .eq. 1) call jedi_speciation

        if (((YEAR .gt. 50) .or. (kspinup .eq. 0)) .and. kDist_dyn .eq. 1 ) call jedi_disturbance

      end select

      return
      end subroutine vegetation_step

!     ******************************************************************
!     VEGETATION_STOP
!     ******************************************************************

      subroutine vegetation_stop ()
      use globe_mod
      use jedi_mod
      implicit none

      call jedi_stop

      return
      end subroutine vegetation_stop

!     ******************************************************************
!     VEGETATION_OUTPUT
!     ******************************************************************

      subroutine vegetation_output ()
      use globe_mod
      use jedi_mod
      implicit none

      integer :: zyr_modulo, zperiod

      zperiod    = kglobe_yrsskip + kglobe_yrsout
      zyr_modulo = modulo(kglobe_year, zperiod)

      if ((zyr_modulo .gt. kglobe_yrsskip) .or. (zyr_modulo == 0)) then
        if (koutput == 1) then
          call jedi_output
        else
          call jedi_output_reset
        endif
        if (kspec_yrout == 1 .and. kglobe_month == 12) call jedi_output_sppgrid
      else
        call jedi_output_reset
        if (kspec_yrout == 1 .and. kglobe_month == 12) call jedi_output_species_reset
      endif

      return
      end subroutine vegetation_output

!     ******************************************************************
!     VEGETATION_DIAG
!     ******************************************************************

      subroutine vegetation_diag ()
      use globe_mod
      use jedi_mod
      implicit none

      YEAR  = kglobe_year
      OUTYEAR = kglobe_year + kglobe_firstyear - 1
      MONTH = kglobe_month
      kdiy  = kglobe_diy

      if (kdiag .eq. 1) call jedi_diag

      return
      end subroutine vegetation_diag
