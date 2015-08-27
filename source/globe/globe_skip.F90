#define __ACTIVATE
#include "../globe/globe_macros.f90"

!     ******************************************************************
!     GLOBE_SKIP_INIT
!     ******************************************************************

      subroutine globe_skip_init ()
      use globe_mod
      implicit none

      integer ioerr

      __diag(kglobe_fdiag,'globe_skip_init')

!     * call globe_skip_readnamelist

      allocate(pglobe_skip_d2dt2(nglobe_cgpt))
      allocate(pglobe_skip_tau(nglobe_cgpt))
      allocate(pglobe_skip_tau2(nglobe_gpts2))

!     * open output file for skip

      open(kglobe_fskip, FILE=sglobe_fskip, STATUS='replace',          &
     &                   FORM='unformatted', IOSTAT=ioerr)
      call error_open_new(ioerr, sglobe_fskip)

!     * initialize counters

      nglobe_skip_lastyrs  = 0
      nglobe_skip_firstrun = 0

!     * done

      __diag(kglobe_fdiag,'globe_skip_init: done')

      return
      end subroutine globe_skip_init

!     ******************************************************************
!     GLOBE_SKIP_LOOP
!     ******************************************************************

      subroutine globe_skip_loop ()
      use globe_mod
      implicit none

      integer   kyr_run_start, kyr_run_end
      integer   nyrs_skip
      integer   istep

      real      pDTskip

      __diag(kglobe_fdiag,'globe_skip_loop: main loop')

!     * initialization

      istep = 0

!     ------------------------------------------------------------------
!     * main loop
!     ------------------------------------------------------------------

      do while ((kglobe_year + kglobe_firstyear - 1) .le. kglobe_lastyear)

!       * reset accumulation

        call globe_skip_reset

!       ----------------------------------------------------------------
!       * full run loop
!       ----------------------------------------------------------------

        kyr_run_start = kglobe_year
        kyr_run_end   = kyr_run_start + nglobe_skip_runyrs - 1

        do kglobe_year = kyr_run_start, kyr_run_end
          istep       = istep + 1

          if (kglobe_mypid == kglobe_kroot) then
            write(*,*) 'globe_skip_loop: current year = ', kglobe_year, ' full step = ', istep
          endif

!         * run yearly steps
          call globe_step(GLOBE_STEP_PREYEAR)

!         * monthly loop
          do kglobe_month = 1, 12

!           * run monthly steps
            call globe_step(GLOBE_STEP_PREMONTH)

            call globe_status

!           * daily loop
            do kglobe_day = 1, kglobe_dpm(kglobe_month) + kglobe_leapday

!             * run daily steps
              call globe_step(GLOBE_STEP_PREDAY)

            enddo

!           * run diagnostic modules
            call globe_diag

!           * write out monthly averages
            call globe_output

          enddo

        enddo

!       ----------------------------------------------------------------
!       * test for skipping
!       ----------------------------------------------------------------

        if (nglobe_skip_firstrun .eq. 0.) then

!         * calculate maximum skip length

          pDTskip = real(nglobe_skip_maxyrs)

!         call atmos_skip_max(pDTskip)
!         call soil_skip_max(pDTskip)
          call vegetation_skip_max(pDTskip)

          nyrs_skip = MAX(0, FLOOR(pglobe_skip_fractau * pDTskip))
          nyrs_skip = MIN(nyrs_skip, kglobe_lastyear - (kglobe_year + kglobe_firstyear - 1) - nglobe_skip_runyrs + 1)

          if (kglobe_mypid == kglobe_kroot) then
            write(*,*) 'globe_skip_loop: nyr_skip = ', nyrs_skip
          endif

!         * skip years

          if (nyrs_skip .ge. 1) then
            kyr_run_start = kglobe_year
            kyr_run_end   = kglobe_year + nyrs_skip - 1

            do kglobe_year = kyr_run_start, kyr_run_end
!             call atmos_skip(1)
!             call soil_skip(1)
              call vegetation_skip(1)
              kglobe_timestep = kglobe_timestep + 365
              call globe_status
            enddo

            nglobe_skip_lastyrs = nyrs_skip
          else
            nglobe_skip_lastyrs = 0
          endif
        else

          nglobe_skip_firstrun  = 0

          if (kglobe_mypid == kglobe_kroot) then
            write(*,*) 'globe_skip_loop: first calculation of derivatives, no second derivative'
          endif

        endif

      enddo

      return
      end subroutine globe_skip_loop

!     ******************************************************************
!     GLOBE_SKIP_RESET
!     ******************************************************************

      subroutine globe_skip_reset ()
      use globe_mod
      implicit none

!     call atmos_skip_reset
!     call soil_skip_reset
      call vegetation_skip_reset

      return
      end subroutine globe_skip_reset

!     ******************************************************************
!     GLOBE_SKIP_STOP
!     ******************************************************************

      subroutine globe_skip_stop ()
      use globe_mod
      implicit none

      __diag(kglobe_fdiag,'globe_skip_stop')

      __deallocate(pglobe_skip_d2dt2)
      __deallocate(pglobe_skip_tau)
      __deallocate(pglobe_skip_tau2)

      close(kglobe_fskip)

      return
      end subroutine globe_skip_stop

!     ******************************************************************
!     GLOBE_SKIP_MAX
!     ******************************************************************

!     * subroutine to calculate maximum length of interval that is safe
!       to skip for given mean aX, first derivative adXdt, and second
!       derivative ad2Xdt2.

      subroutine globe_skip_max (aX, adXdt, ad2Xdt2, pDTskip)
      use globe_mod
      implicit none

      real            :: aX(nglobe_cgpt)
      real            :: adXdt(nglobe_cgpt)
      real            :: ad2Xdt2(nglobe_cgpt)
      real            :: pDTskip

      real            :: cSecsYr
      parameter       (cSecsYr = 86400. * 365.)   ! seconds per year

      __diag(kglobe_fdiag,'globe_skip_max')

      pglobe_skip_tau(:)    = real(nglobe_skip_maxyrs)

!     * check if skip interval results in sign change of variable

      where (aX(:) * (aX(:) + adXdt(:) * pglobe_skip_tau(:) * cSecsYr) .lt. 0.0)
        pglobe_skip_tau(:)  = MIN(pglobe_skip_tau(:), ABS(aX(:) / adXdt(:) / cSecsYr))
      end where

      __globe_mpga_to_from(pglobe_skip_tau2,pglobe_skip_tau)
      pDTskip               = MINVAL(pglobe_skip_tau2(:))
      if (kglobe_mypid == kglobe_kroot) write(*,*) 'globe_skip_max: pDTskip = ', pDTskip

!     * check if skip interval results in sign change of derivative -- does not work well

!      where (abs(ad2Xdt2(:) * pglobe_skip_tau(:) * cSecsYr) .gt. pglobe_skip_fchange * aX(:))
!      where (adXdt(:) * (adXdt(:) + ad2Xdt2(:) * pglobe_skip_tau(:) * cSecsYr) .lt. 0.0)
!        pglobe_skip_tau(:)  = MIN(pglobe_skip_tau(:), abs(adXdt(:)/ad2Xdt2(:)/cSecsYr))
!      end where
!      end where

!     * reduce skipping interval if large change in presence of small mean

      where (abs(adXdt(:) * pglobe_skip_tau(:) * cSecsYr) .gt. ABS(0.5 * aX(:)) .and. aX(:) .ne. 0. )
        pglobe_skip_tau(:)  = MIN(pglobe_skip_tau(:), 0.5 * pglobe_skip_tau(:) * ABS(aX(:) / adXdt(:) / cSecsYr))
      end where

      __globe_mpga_to_from(pglobe_skip_tau2,pglobe_skip_tau)
      pDTskip               = MINVAL(pglobe_skip_tau2(:))
      if (kglobe_mypid == kglobe_kroot) write(*,*) 'globe_skip_max: pDTskip = ', pDTskip

!     * done

      __diag(kglobe_fdiag,'globe_skip_max: done')

      return
      end subroutine globe_skip_max
