#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam_fio.F90
!> \brief I/O related routines of JAM

!> \file jam_fio.F90
!> This file includes the subroutines which handle the file input
!! , like the namelist and the input of climate forcing.

!     ******************************************************************
!     JAM_READ_NAMELIST
!     ******************************************************************

!> \brief Reads the JAM namelist file

!> \detail This routine reads the namelist parameter file of JAM and
!! broadcasts the parameters to all processors.\n

!> \b Variables \b set: kCO2_dyn, kRandYear,
!! pDeltaT, pDeltaP, pTEMP_init, fRADs_as_init, fRADl_as_init,
!! fH2Ol_as_init, cH2Og_a_init, xCO2g_a_init

!     ******************************************************************

      subroutine jam_read_namelist ()
      use jam_mod
      implicit none

      integer :: ioerr

!     * Definition variable parameters
      namelist /jampar/ kCO2_dyn,       &
     &                  kRandYear,      &
     &                  pDeltaT,        &
     &                  pDeltaP,        &
     &                  pTEMP_init,     &
     &                  fRADs_as_init,  &
     &                  fRADl_as_init,  &
     &                  fH2Ol_as_init,  &
     &                  cH2Og_a_init,   &
     &                  xCO2g_a_init

      __diag(kfile_diag,'jam_read_namelist')

      if (kmypid == kroot) then

!       * input of variable parameters from the parameter file

        open(kfile_namelist, FILE=sfile_namelist, STATUS='old',        &
       &                     FORM='formatted', IOSTAT=ioerr)
        if (ioerr .eq. 0) then
          read(kfile_namelist, jampar, iostat=ioerr)
!          read(kfile_namelist, jampar)

          if (kdiag .eq. 1) then
            write(kfile_diag, '(" *************************************************")')
            write(kfile_diag, '(" * Namelist JAMPAR from ",A)') sfile_namelist
            write(kfile_diag, '(" *************************************************")')
          endif
          close(kfile_namelist)
        else
          if (kdiag .eq. 1) then
            write(kfile_diag, '(" ******************************************")')
            write(kfile_diag, '(" * ERROR reading file ",A)') sfile_namelist
            write(kfile_diag, '(" * Using default values")')
            write(kfile_diag, '(" ******************************************")')
          endif
        endif
        if (kdiag .eq. 1) then
          write(kfile_diag, jampar)
          flush(kfile_diag)
        endif

      endif

      __diag(kfile_diag,'jam_read_namelist: done')

!     * broadcast namelist parameters

      __globe_mpbc(kCO2_dyn)
      __globe_mpbc(kRandYear)
      __globe_mpbc(pDeltaT)
      __globe_mpbc(pDeltaP)

      __globe_mpbc(pTEMP_init)
      __globe_mpbc(fRADs_as_init)
      __globe_mpbc(fRADl_as_init)
      __globe_mpbc(fH2Ol_as_init)
      __globe_mpbc(cH2Og_a_init)
      __globe_mpbc(xCO2g_a_init)

      return
      end subroutine jam_read_namelist

!     ******************************************************************
!     JAM_INPUT_OPEN
!     ******************************************************************

!> \brief Opens a climate forcing file

!> \detail This routine tries to open a climate input file. At success the
!! input field count is incremented.\n

!> \b Variables \b modified: ndatafields, kfieldIDs

!> \param kfile File ID of the climate input file
!> \param sfile File name of the climate input file
!> \param input Status flag if file exist

!     ******************************************************************

      subroutine jam_input_open (kfile, sfile, input)
      use jam_mod
      implicit none

      integer,          INTENT(in)  :: kfile
      character(len=*), INTENT(in)  :: sfile
      logical,          INTENT(out) :: input

      integer :: ioerr

      if (kmypid == kroot) then

        inquire(FILE=sfile, EXIST=input)
        if (input) then
          call globe_open_input(TRIM(sfile), kfile, kfile_diag)
          ndatafields = ndatafields + 1
          kfieldIDs(ndatafields) = kfile
          call globe_rewind(kfile)
        endif

      endif

      return
      end subroutine jam_input_open

!     ******************************************************************
!     JAM_INPUT_STOP
!     ******************************************************************

!> \brief Closes all climate forcing files

!> \detail This routine closes all the open climate forcing files.

!     ******************************************************************

      subroutine jam_input_stop ()
      use jam_mod
      implicit none

      __diag(kfile_diag,'jam_input_stop')

      if (kmypid == kroot) then

        if (input_pTEMP)     call globe_close_input(kfile_temp2m, kfile_diag)
        if (input_fRADs_as)  call globe_close_input(kfile_solrad, kfile_diag)
        if (input_fRADl_as)  call globe_close_input(kfile_terrad, kfile_diag)
        if (input_fH2Ol_as)  call globe_close_input(kfile_precip, kfile_diag)
        if (input_cH2Og_a)   call globe_close_input(kfile_relhum, kfile_diag)
        if (input_xCO2g_a)   call globe_close_input(kfile_co2a, kfile_diag)

        if (input_pco2list)  close(kfile_pco2list)

      endif

      __diag(kfile_diag,'jam_input_stop: done')

      return
      end subroutine jam_input_stop

!     ******************************************************************
!     JAM_SPT_PREREAD_CLIMATE
!     ******************************************************************

!> \brief Pre-read climate forcing files

!> \detail This routine is only used by the single grid point version.
!! It reads all climate forcing data of this grid point from all
!! existing files to memory.\n

!> \b Variables \b set: preclim, nmaxpoints

!     ******************************************************************

      subroutine jam_spt_preread_climate ()
      use jam_mod
      implicit none

      integer :: i, np, kfile, zhead(8), ioerr
      real,allocatable :: ztmp(:)

      __diag(kfile_diag,'jam_spt_preread_climate')

      if (kmypid == kroot) then

        allocate(ztmp(nlon*nlat))

!       * fill field preclim(:,:) and set nmaxpoints

        nmaxpoints = 0
        do i = 1, ndatafields
          kfile = kfieldIDs(i)
          np = 0
          ioerr = 0
          do while (ioerr .eq. 0)
            __globe_read_srvnc(kfile,zhead,iostat=ioerr)
            if (ioerr .eq. 0) then
              __globe_read_srvnc(kfile,zhead,ztmp)
              np = np + 1
              preclim(i,np) = ztmp(kLPID(1))
            endif
          enddo
          if (nmaxpoints .eq. 0) nmaxpoints = np
          if (nmaxpoints .ne. np) then
            write(*,*) 'ERROR: The climate data files differ in size'
            stop
          endif
        enddo

        deallocate(ztmp)

      endif

      __diag(kfile_diag,'jam_spt_preread_climate: done')

      return
      end subroutine jam_spt_preread_climate

!     ******************************************************************
!     JAM_READ_CLIMATE_FILE
!     ******************************************************************

!> \brief Reads data from climate forcing files

!> \detail This routine reads one data field from a climate input file,
!! or in the single point mode gets the data point from the pre-read
!! variable in memory.\n

!> \b Variables \b modified: nfieldcount, nspointer

!> \param climdata  Climate data field
!> \param dsize     Field size of \c climdata
!> \param kfile     The file ID of the climate input file
!> \param sfile     The file name of the climate input file
!> \param zdatayear The actual year to check file data consistency

!     ******************************************************************

      subroutine jam_read_climate_file (climdata,dsize,kfile,sfile,zdatayear)
      use jam_mod
      implicit none

      integer,          INTENT(in)  :: dsize, kfile, zdatayear
      real,             INTENT(out) :: climdata(dsize)
      character(len=*), INTENT(in)  :: sfile

      integer :: zhead(8), ioerr, field, i
      real,allocatable :: ztmp(:)

      if (kmypid == kroot) then

!       * get pre-read data for single point version

        if (ksingle .ne. 0 .and. kRandYear .eq. 0) then
          do i = 1, size(kfieldIDs)
            if (kfieldIDs(i) .eq. kfile) field = i
          enddo

          climdata(1) = preclim(field,nspointer)

          nfieldcount = nfieldcount + 1
          if (nfieldcount .gt. ndatafields) then
            nfieldcount = 1
            nspointer = nspointer + 1
            if (nspointer .gt. nmaxpoints) nspointer = 1
          endif

          return
        endif

        allocate(ztmp(nlon*nlat))

!       * check the header

        __globe_read_srvnc(kfile,zhead,iostat=ioerr)
        if (ioerr .ne. 0) then
          call globe_rewind(kfile)
        endif

!       * read the data

        __globe_read_srvnc(kfile,zhead,ztmp)
        climdata(:) = ztmp(kLPID(1:nLPts))

!       * check the data year

        if (zdatayear .ne. zhead(3)/10000 .and. kRandYear .eq. 0) then
          write(*,*) 'ERROR: Computed year different to data set'
          stop
        endif

        deallocate(ztmp)

      endif

      return
      end subroutine jam_read_climate_file

!     ******************************************************************
!     JAM_SEEK_CLIM_FILE
!     ******************************************************************

!> \brief Seeks climate forcing files

!> \detail This routine seeks all existing climate input files for the
!! actual date, or sets the pointer of the pre-read memory variable.\n

!> \b Variables \b set: nspointer

!     ******************************************************************

      subroutine jam_seek_clim_file ()
      use jam_mod
      implicit none

      integer :: i
      integer :: zlyear = 0        ! first counted leap year
      integer :: zldays = 0        ! number of leap days to seek forward
      integer :: zryears           ! rest of data years by mod with actual year
      integer :: zdpy
      integer :: zhead(8)
      real    :: random
      real,allocatable :: ztmp(:)
      character(len=135) :: msg

      __diag(kfile_diag,'jam_seek_clim_file')

!     * Rest of the years that are not a multiple of the number of years in the climate data set
!     * (kyear is the next year to be computed -> skip kyear-1 years)

      if (kRandYear .eq. 0) then
        zryears = MOD(kyear - 1, n_data_years)
      else
        call RANDOM_NUMBER(random)
        zryears = INT(random * real(n_data_years))
        if (zryears .ge. n_data_years) zryears = n_data_years - 1
        write(msg,*) "JAM seek ",zryears," years"
        __diag(kfile_diag,TRIM(msg))
      endif

!     * Compute the number of leap days to seek forward

      if (kdpm(13) .ne. 0) then
        zlyear = kdpm(13) - k_first_dyear + 1
        if (zlyear .le. zryears) zldays = (zryears - zlyear + 4) / 4
        if (kdpm(14) .ne. 0) then
          zlyear = kdpm(14) - k_first_dyear + 1
          if (zlyear .le. zryears) zldays = zldays - 1
        endif
      endif

!     * "normal" number of days per year

      zdpy = SUM(kdpm(1:12))

!     * set pointer for pre-read climate (single point)

      if (ksingle .ne. 0 .and. kRandYear .eq. 0) then
        nspointer = (zdpy * zryears + zldays) * ktspd + 1
        return
      endif

!     * Seek the file pointers to the right position

      allocate(ztmp(nlon*nlat))

      do i = 1, (zdpy * zryears + zldays) * ktspd
        if (input_pTEMP)    __globe_read_srvnc(kfile_temp2m,zhead,ztmp)
        if (input_fRADs_as) __globe_read_srvnc(kfile_solrad,zhead,ztmp)
        if (input_fRADl_as) __globe_read_srvnc(kfile_terrad,zhead,ztmp)
        if (input_fH2Ol_as) __globe_read_srvnc(kfile_precip,zhead,ztmp)
        if (input_cH2Og_a)  __globe_read_srvnc(kfile_relhum,zhead,ztmp)
        if (input_xCO2g_a)  __globe_read_srvnc(kfile_co2a,zhead,ztmp)
      enddo

      deallocate(ztmp)

      __diag(kfile_diag,'jam_seek_clim_file: done')

      return
      end subroutine jam_seek_clim_file

!     ******************************************************************
!     JAM_READ_CO2_LIST
!     ******************************************************************

!> \brief Reads the annual CO<sub>2</sub> list file

!> \detail This routine reads one CO<sub>2</sub> value from a text file
!! with an annual CO<sub>2</sub> list.

!> \param zoutyear Requested data year
!> \param zparam   CO<sub>2</sub> value from list

!     ******************************************************************
!     * read data from zoutyear in formatted file to parameter zparam
!     * example file format
!     * 1859 280.0000
!     * 1860 280.1000
!     * 1861 282.3000

      subroutine jam_read_co2_list (zparam, zoutyear)
      use jam_mod
      implicit none

      integer, INTENT(in)  :: zoutyear
      real,    INTENT(out) :: zparam

      integer :: ioerr, zlistyear
      real    :: zlistparam

      character(len=80) :: format
      character(len=135) :: msg

      if (kmypid == kroot) then

        rewind(kfile_pco2list, IOSTAT=ioerr)

        do while(ioerr .eq. 0)
          read(kfile_pco2list, *, IOSTAT=ioerr) zlistyear, zlistparam
          if (ioerr .ne. 0) then
            format = '("Year ",I4," not found in ",A,". Parameter not changed: ",A)'
            write(msg,format) zoutyear, sfile_pco2list, globe_string(zparam)
            __diag(kfile_diag, TRIM(msg))
          else
            if (zlistyear .eq. zoutyear) then
              zparam = zlistparam
              format = '("New parameter value ",F8.3," found for year ",I4)'
              write(msg,format) zparam, zoutyear
              __diag(kfile_diag, TRIM(msg))
              exit
            endif
          endif
        enddo

      endif

      return
      end subroutine jam_read_co2_list
