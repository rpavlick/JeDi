#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam_sub.F90
!> \brief Other subroutines of the JAM module

!> \file jam_sub.F90
!> This file includes subroutines of the JAM module for initialization,
!! to get and scatter climate forcing, and derive the calendar from it.

!     ******************************************************************
!     JAM_ALLOC
!     ******************************************************************

!> \brief Allocation of all JAM fields

!> \detail This routine allocates all fields of JAM.

!     ******************************************************************

      subroutine jam_alloc ()
      use jam_mod
      implicit none

      __diag(kfile_diag,'jam_alloc')

!     ------------------------------------------------------------------
!     * allocate fields
!     ------------------------------------------------------------------

!     * arrays for the full grid (monthly)

      if (input_pTEMP)    then; __allocate(pTEMP2,(nLPts2,mdim)); endif
      if (input_fRADs_as) then; __allocate(fRADs_as2,(nLPts2,mdim)); endif
      if (input_fRADl_as) then; __allocate(fRADl_as2,(nLPts2,mdim)); endif
      if (input_fH2Ol_as) then; __allocate(fH2Ol_as2,(nLPts2,mdim)); endif
      if (input_cH2Og_a)  then; __allocate(cH2Og_a2,(nLPts2,mdim)); endif
      if (input_xCO2g_a)  then; __allocate(xCO2g_a2,(nLPts2,mdim)); endif

!     * arrays of CPU grid

      if (input_pTEMP)    then; __allocate(pTEMP,(nCPts,mdim)); endif
      if (input_fRADs_as) then; __allocate(fRADs_as,(nCPts,mdim)); endif
      if (input_fRADl_as) then; __allocate(fRADl_as,(nCPts,mdim)); endif
      if (input_fH2Ol_as) then; __allocate(fH2Ol_as,(nCPts,mdim)); endif
      if (input_cH2Og_a)  then; __allocate(cH2Og_a,(nCPts,mdim)); endif
      if (input_xCO2g_a)  then; __allocate(xCO2g_a,(nCPts,mdim)); endif

!     * array of single point pre-read climate

      if (ksingle .ne. 0 .and. kRandYear .eq. 0 .and. ndatafields .gt. 0) then
        __allocate(preclim,(ndatafields,mdim*12*n_data_years))
      endif

      __diag(kfile_diag,'jam_alloc: done')

      return
      end subroutine jam_alloc

!     ******************************************************************
!     JAM_DEALLOC
!     ******************************************************************

!> \brief Deallocation of all JAM fields

!> \detail This routine deallocates all fields of JAM.

!     ******************************************************************

      subroutine jam_dealloc ()
      use jam_mod
      implicit none

      __diag(kfile_diag,'jam_dealloc')

!     ------------------------------------------------------------------
!     * deallocate fields
!     ------------------------------------------------------------------

      __deallocate(kLPID)

      if (input_pTEMP)    then; __deallocate(pTEMP2); endif
      if (input_fRADs_as) then; __deallocate(fRADs_as2); endif
      if (input_fRADl_as) then; __deallocate(fRADl_as2); endif
      if (input_fH2Ol_as) then; __deallocate(fH2Ol_as2); endif
      if (input_cH2Og_a)  then; __deallocate(cH2Og_a2); endif
      if (input_xCO2g_a)  then; __deallocate(xCO2g_a2); endif



      if (input_pTEMP)    then; __deallocate(pTEMP); endif
      if (input_fRADs_as) then; __deallocate(fRADs_as); endif
      if (input_fRADl_as) then; __deallocate(fRADl_as); endif
      if (input_fH2Ol_as) then; __deallocate(fH2Ol_as); endif
      if (input_cH2Og_a)  then; __deallocate(cH2Og_a); endif
      if (input_xCO2g_a)  then; __deallocate(xCO2g_a); endif

!     * arrays of single point pre-read climate

      if (ksingle .ne. 0 .and. kRandYear .eq. 0 .and. ndatafields .gt. 0) then
        __deallocate(preclim)
      endif

      __diag(kfile_diag,'jam_dealloc: done')

      return
      end subroutine jam_dealloc

!     ******************************************************************
!     JAM_INPUT_INIT
!     ******************************************************************

!> \brief Initialization of JAM input and date

!> \detail This routine opens all climate input files, positions the file
!! pointers, sets up the calendar, and scatters the related variables.\n

!> \b Variables \b set: input_pTEMP, input_fRADs_as, input_fRADl_as,
!! input_fH2Ol_as, input_cH2Og_a, input_xCO2g_a,
!> input_pco2list, mdim, kmonmeans\n

!> \b Variables \b modified: kCO2_dyn, kdpm

!     ******************************************************************

      subroutine jam_input_init ()
      use jam_mod
      implicit none

      integer :: i, ioerr

      __diag(kfile_diag,'jam_input_init')

      if (kmypid == kroot) then
!       * count the number of available data files
        ndatafields = 0

!       * try to open input files

        call jam_input_open(kfile_temp2m, sfile_temp2m, input_pTEMP)
        call jam_input_open(kfile_solrad, sfile_solrad, input_fRADs_as)
        call jam_input_open(kfile_terrad, sfile_terrad, input_fRADl_as)
        call jam_input_open(kfile_precip, sfile_precip, input_fH2Ol_as)
        call jam_input_open(kfile_relhum, sfile_relhum, input_cH2Og_a)
        call jam_input_open(kfile_co2a,   sfile_co2a,   input_xCO2g_a)

        input_pco2list = .false.
        if (kCO2_dyn == 1) then
!         * open CO2 concentration list
          inquire(FILE=sfile_pco2list, EXIST=input_pco2list)
          if (input_pco2list) then
            open(kfile_pco2list, FILE=sfile_pco2list, STATUS='old',    &
     &                           FORM='FORMATTED', IOSTAT=ioerr)
            call error_open_old(ioerr, TRIM(sfile_pco2list))
          else
            __diag(kfile_diag, 'CO2 concentration list file not found.')
            __diag_num(kfile_diag, 'Using namelist value: ',xCO2g_a_init)
            kCO2_dyn = 0
          endif
        endif

!       * get date from climate data

        call jam_get_clim_date

!       * rewind all file pointers

        if (input_pTEMP)     call globe_rewind(kfile_temp2m)
        if (input_fRADs_as)  call globe_rewind(kfile_solrad)
        if (input_fRADl_as)  call globe_rewind(kfile_terrad)
        if (input_fH2Ol_as)  call globe_rewind(kfile_precip)
        if (input_cH2Og_a)   call globe_rewind(kfile_relhum)
        if (input_xCO2g_a)   call globe_rewind(kfile_co2a)

!       * seek the climate files to the right date (in case of restart)

        if (krestart .gt. 0) call jam_seek_clim_file

!       * count the maximum number of days per month for allocation of fields
        mdim = MAXVAL(kdpm(1:12)) * ktspd
        if (mdim .eq. 1) kmonmeans = 1       ! case of monmeans in data set

!       * set 30 days per month with monthly means (Must be done after jam_seek_clim_file!)

        if (kmonmeans .eq. 1) mdim = 30
        if (kmonmeans .eq. 1) kdpm(1:12) = 30

      endif

!     * scatter the global variables

      __globe_mpbc(kCO2_dyn)
      __globe_mpbc(mdim)
      __globe_mpbc(n_data_years)
      __globe_mpbc(k_first_dyear)
      __globe_mpbc(kdpm)
      __globe_mpbc(ktspd)

      __globe_mpbc(input_pTEMP)
      __globe_mpbc(input_fRADs_as)
      __globe_mpbc(input_fRADl_as)
      __globe_mpbc(input_fH2Ol_as)
      __globe_mpbc(input_cH2Og_a)
      __globe_mpbc(input_xCO2g_a)

      __diag(kfile_diag,'jam_input_init: done')

      return
      end subroutine jam_input_init

!     ******************************************************************
!     JAM_READ_CLIMATE_MONTH
!     ******************************************************************

!> \brief Reads the monthly climate forcing

!> \detail This routine calls the monthly reading of all climate data
!! sets, and fills up all the day fields, if the climate forcing
!! includes only monthly means.\n

!> \b Variables \b set: pTEMP2, fRADs_as2, fRADl_as2, fH2Ol_as2,
!! cH2Og_a2, xCO2g_a2

!     ******************************************************************

      subroutine jam_read_climate_month ()
      use jam_mod
      implicit none

      integer :: i
      integer :: zmean_days = 0

      integer :: zryears           ! rest of data years by mod with actual year
      integer :: zdatayear         ! computed actual year in the data set

      if (kmypid == kroot) then

!       * compute data year to check it

        zryears = MOD(kyear-1, n_data_years) + 1
        zdatayear = zryears + k_first_dyear - 1

        if (kmonmeans .eq. 1) zmean_days = 29
        do i=1, (kdpm(kmonth) - zmean_days + nleapday) * ktspd
          if (input_pTEMP) then
            call jam_read_climate_file(pTEMP2(1:nLPts,i)   ,nLPts,kfile_temp2m,sfile_temp2m,zdatayear)
          endif
          if (input_fRADs_as) then
            call jam_read_climate_file(fRADs_as2(1:nLPts,i),nLPts,kfile_solrad,sfile_solrad,zdatayear)
          endif
          if (input_fRADl_as) then
            call jam_read_climate_file(fRADl_as2(1:nLPts,i),nLPts,kfile_terrad,sfile_terrad,zdatayear)
          endif
          if (input_fH2Ol_as) then
            call jam_read_climate_file(fH2Ol_as2(1:nLPts,i),nLPts,kfile_precip,sfile_precip,zdatayear)
          endif
          if (input_cH2Og_a) then
            call jam_read_climate_file(cH2Og_a2(1:nLPts,i),nLPts,kfile_relhum,sfile_relhum,zdatayear)
          endif
          if (input_xCO2g_a) then
            call jam_read_climate_file(xCO2g_a2(1:nLPts,i),nLPts,kfile_co2a,sfile_co2a,zdatayear)
          endif
        enddo

!       * fill up 30 days of the month with the monthly mean

        if (kmonmeans .eq. 1) then
          do i = 2, 30
            if (input_pTEMP)     pTEMP2(1:nLPts,i)     = pTEMP2(1:nLPts,1)
            if (input_fRADs_as)  fRADs_as2(1:nLPts,i)  = fRADs_as2(1:nLPts,1)
            if (input_fRADl_as)  fRADl_as2(1:nLPts,i)  = fRADl_as2(1:nLPts,1)
            if (input_fH2Ol_as)  fH2Ol_as2(1:nLPts,i)  = fH2Ol_as2(1:nLPts,1)
            if (input_cH2Og_a)   cH2Og_a2(1:nLPts,i)   = cH2Og_a2(1:nLPts,1)
            if (input_xCO2g_a)   xCO2g_a2(1:nLPts,i)   = xCO2g_a2(1:nLPts,1)
          enddo
        endif

      endif

      call jam_scatter_climate_month

      return
      end subroutine jam_read_climate_month

!     ******************************************************************
!     JAM_SCATTER_CLIMATE_MONTH
!     ******************************************************************

!> \brief Scatter of monthly climate forcing

!> \detail This routine fills the additional fields in the variable
!! arrays with sizes of multiple of the number of CPUs, adds an anomaly
!! to temperature and precipitation, and scatters the climate forcing
!! fields.\n

!> \b Variables \b modified: pTEMP2, fRADs_as2, fRADl_as2, fH2Ol_as2,
!! cH2Og_a2, xCO2g_a2\n

!> \b Variables \b set: pTEMP, fRADs_as, fRADl_as, fH2Ol_as, cH2Og_a, xCO2g_a

!     ******************************************************************

      subroutine jam_scatter_climate_month ()
      use jam_mod
      implicit none

      integer :: i, j

!     * extend climate data

      if (kmypid == kroot) then

        do i = 1, mdim
          do j = nLPts + 1, nLPts2
            if (input_pTEMP)     pTEMP2(j,i)     = pTEMP2(nLPts,i)
            if (input_fRADs_as)  fRADs_as2(j,i)  = fRADs_as2(nLPts,i)
            if (input_fRADl_as)  fRADl_as2(j,i)  = fRADl_as2(nLPts,i)
            if (input_fH2Ol_as)  fH2Ol_as2(j,i)  = fH2Ol_as2(nLPts,i)
            if (input_cH2Og_a)   cH2Og_a2(j,i)   = cH2Og_a2(nLPts,i)
            if (input_xCO2g_a)   xCO2g_a2(j,i)   = xCO2g_a2(nLPts,i)
          enddo
        enddo
      endif

!     * add temperature anomaly

      if (input_pTEMP) pTEMP2(:,:) = pTEMP2(:,:) + pDeltaT

!     * add precipitation anomaly

      if (input_fH2Ol_as) then
        fH2Ol_as2(:,:) = fH2Ol_as2(:,:) * (1.0 + pDeltaP / 100.0)
      endif

!     * scatter data

      if (input_pTEMP)    __globe_mpsc_from_to(pTEMP2,pTEMP)
      if (input_fRADs_as) __globe_mpsc_from_to(fRADs_as2,fRADs_as)
      if (input_fRADl_as) __globe_mpsc_from_to(fRADl_as2,fRADl_as)
      if (input_fH2Ol_as) __globe_mpsc_from_to(fH2Ol_as2,fH2Ol_as)
      if (input_cH2Og_a)  __globe_mpsc_from_to(cH2Og_a2,cH2Og_a)
      if (input_xCO2g_a)  __globe_mpsc_from_to(xCO2g_a2,xCO2g_a)

      return
      end subroutine jam_scatter_climate_month

!     ******************************************************************
!     JAM_GET_CLIM_DATE
!     ******************************************************************

!> \brief Get calendar from climate forcing

!> \detail This routine reads a climate file and sets the calendar, the
!! first and last used year, and checks the climate file for
!! continuity.\n

!> \b Variables \b set: k_first_dyear, kdpm, ktspd, n_data_years

!     ******************************************************************

      subroutine jam_get_clim_date ()
      use jam_mod
      implicit none

      integer :: zlyear       ! last leap yera
      integer :: zts          ! actual time step
      integer :: zyear,zmonth ! actual year and month
      integer :: zndays       ! number of days of february of the year before
      integer :: zhead(8), ioerr
      integer :: tmp_fileID
      character*80 :: sfile
      character*160 :: stemp
      character*10 :: ttemp
      real,allocatable :: ztemp(:)

      __diag(kfile_diag,'jam_get_clim_date')

      allocate(ztemp(nlon*nlat))

!     * find existing file ID

      if (input_pTEMP) then
        tmp_fileID = kfile_temp2m
      elseif (input_fRADs_as) then
        tmp_fileID = kfile_solrad
      elseif (input_fRADl_as) then
        tmp_fileID = kfile_terrad
      elseif (input_fH2Ol_as) then
        tmp_fileID = kfile_precip
      elseif (input_cH2Og_a) then
        tmp_fileID = kfile_relhum
      elseif (input_xCO2g_a) then
        tmp_fileID = kfile_co2a

      else
        write(*,*) "ERROR: jam_get_clim_date(): No climate data file found"
        stop
      endif

!     * read and check field size

      call globe_rewind(tmp_fileID)
      __globe_read_srvnc(tmp_fileID,zhead,iostat=ioerr)
      inquire(tmp_fileID, NAME=sfile)
      __error_read(ioerr,sfile)
      if ((nlat .ne. zhead(6)) .or. (nlon .ne. zhead(5))) then
        write(*,*) 'ERROR: The grid-size of the climate data differ.'
        stop
      endif

      write(stemp,*) 'jam_get_clim_date: using ', sfile
      __diag(kfile_diag,TRIM(stemp))

      call DATE_AND_TIME(TIME=ttemp)
      write(stemp,*) 'jam_get_clim_date: file read start: ', TRIM(ttemp)
      __diag(kfile_diag,TRIM(stemp))

!     * get first date from climate data

      kdpm(:) = 0
      k_first_dyear = zhead(3) / 10000
      zyear = zhead(3) / 10000               ! actual year
      zmonth = MOD(zhead(3) / 100, 100)      ! actual month
      kdpm(zmonth) = MOD(zhead(3), 100)      ! count number of days per month
      if (kdpm(zmonth) .ne. 1) then
        if (kdpm(zmonth) .ne. 0) then
          write(*,*) 'ERROR: climate file data start not with the first day of month'
          stop
        endif
        kdpm(zmonth) = 1                     ! using monmeans
      endif

!     * count the date from climate data

      ktspd  = 0
      zlyear = zyear                         ! init last leap year
      zndays = 0
      do while (ioerr .eq. 0)               ! run up to the end of file
        zts = 0
        do while ((kdpm(zmonth) .eq. MOD(zhead(3), 100)) .and.         &
     &            (zmonth .eq. MOD(zhead(3) / 100, 100)) .and.         &
     &            (ioerr .eq. 0))
          zts = zts + 1
          __globe_read_srvnc(tmp_fileID,zhead,ztemp)
          __globe_read_srvnc(tmp_fileID,zhead,iostat=ioerr)
!         if (ioerr .gt. 0) __error_read(ioerr,sfile)
        enddo
        if (ktspd .eq. 0) ktspd = zts
        if (zts .ne. ktspd) then
          write(*,*) 'ERROR: climate file data are not continuous in timestep'
          stop
        endif
        if (ioerr .eq. 0) then
!         * at the first day of the new year
          if (zhead(3) / 10000 .ne. zyear) then
!           * check for continuous years
            if (zhead(3) / 10000 .ne. zyear + 1) then
              write(*,*) 'ERROR: climate file data are not continuous in year'
              stop
            endif
!           * check for leap year
            if (((zndays .eq. 28) .or. ((zndays .eq. 0) .and. (kdpm(1) .eq. 31))) .and. (kdpm(2) .eq. 29)) then
              if (kdpm(13) .eq. 0) then
                kdpm(13) = zyear             ! first data leap year
              endif
              if (zyear - zlyear .gt. 4) kdpm(14) = zyear - 4  ! data year without leap day
              zlyear = zyear                 ! last leap year
            endif
!           * count and reset
            __diag_num(kfile_diag,'finished year ',zyear)
            zyear = zyear + 1
            n_data_years = n_data_years + 1
            zndays = kdpm(2)
            kdpm(1:12) = 0
          endif
!         * count the days per month
          zmonth = MOD(zhead(3) / 100, 100)
          if (MOD(zhead(3), 100) .ne. kdpm(zmonth) + 1) then
            if (MOD(zhead(3), 100) .ne. 0) then
              write(*,*) 'ERROR: climate file data are not continuous in day'
              stop
            endif
            kdpm(zmonth) = 1                 ! using monmeans
          else
            kdpm(zmonth) = MOD(zhead(3),100) ! number of days per month
          endif
        endif
      enddo
      call globe_rewind(tmp_fileID)

      call DATE_AND_TIME(TIME=ttemp)
      write(stemp,*) 'jam_get_clim_date: file read end: ', TRIM(ttemp)
      __diag(kfile_diag,TRIM(stemp))

!     * check for leap year in day numbers per month
      if ((kdpm(2) .eq. 29) .and. (kdpm(13) .ne. 0)) kdpm(2) = 28

!     * check the number of monthes
      if (kdpm(12) .eq. 0) then
        write(*,*) 'ERROR: Climate file include less than 12 months'
        stop
      endif

      deallocate(ztemp)

!     * output the calendar

      __diag_num(kfile_diag,'First year in the data file: ',k_first_dyear)
      if (kdpm(1) == 1) then
        __diag(kfile_diag,'Data file includes monthly means.')
        __diag(kfile_diag,'Filling every month with 30 days.')
      else
        __diag_num(kfile_diag,'Number of time steps per day: ',ktspd)
        __diag_num(kfile_diag,'Number of days per year: ',SUM(kdpm(1:12)))
        if (kdpm(13) .ne. 0) then
          __diag_num(kfile_diag,'First leap day in year: ',kdpm(13))
          if (kdpm(14) .ne. 0) then
            __diag_num(kfile_diag,'Year with skipped leap day: ',kdpm(14))
          endif
        else
          __diag(kfile_diag,'No leap days.')
        endif
      endif
      __diag_num(kfile_diag,'Number of years in the data file: ',n_data_years)

!     * end

      __diag(kfile_diag,'jam_get_clim_date: done')

      return
      end subroutine jam_get_clim_date
