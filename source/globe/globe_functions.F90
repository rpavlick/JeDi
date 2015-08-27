#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_functions.F90
!> \brief Collection of global routines

!> \file globe_functions.F90
!> This file includes a collection of useful subroutines frequently
!! used by all models of GLOBE. These include file error handling,
!! output of results, diagnostic output, data file reading, and the
!! calculation of the day length.

!     ******************************************************************
!     ERROR HANDLNG
!     ******************************************************************

!> \brief Error handling at file opening

!> \detail This routine shows a message in case of an error during the
!! opening of a new file and stops the program.

!> \param ioerr Status message of the command \c open
!> \param sfile  File name of the file to open

!     ******************************************************************

      subroutine error_open_new (ioerr, sfile)
      implicit none

      integer,          INTENT(in) :: ioerr
      character(len=*), INTENT(in) :: sfile

      if (ioerr .ne. 0) then
        write(*,*) 'ERROR: Creating file ', TRIM(sfile)
        stop
      endif

      return
      end subroutine error_open_new

!     ******************************************************************

!> \brief Error handling at file opening

!> \detail This routine shows a message in case of an error during the
!! reopening of an existing file and stops the program.

!> \param ioerr Status message of the command \c open
!> \param sfile  File name of the file to open

!     ******************************************************************

      subroutine error_open_old (ioerr, sfile)
      implicit none

      integer,          INTENT(in) :: ioerr
      character(len=*), INTENT(in) :: sfile

      if (ioerr .ne. 0) then
        write(*,*) 'ERROR: Opening file ', TRIM(sfile)
        stop
      endif

      return
      end subroutine error_open_old

!     ******************************************************************

!> \brief Error handling at file reading

!> \detail This routine shows a message in case of an error during the
!! reading of a file and stops the program.

!> \param ioerr Status message of the command \c read
!> \param sfile  File name of the file to read
!> \param f90    File name of the Fortran file the \c read was called
!> \param line   Line number in the Fortran file of the \c read command

!     ******************************************************************

      subroutine error_read (ioerr, sfile, f90, line)
      use globe_mod, error_read_local => error_read
      implicit none

      integer,                    INTENT(in) :: ioerr
      character(len=*),           INTENT(in) :: sfile
      character(len=*), OPTIONAL, INTENT(in) :: f90
      integer,          OPTIONAL, INTENT(in) :: line

      if (ioerr .lt. 0) then
        if (present(f90) .and. present(line)) then
          write(*,*) f90,': ', TRIM(globe_string(line)),               &
     &               ': ERROR: End of file in ', TRIM(sfile)
        else
          write(*,*) 'ERROR: End of file in ', TRIM(sfile)
        endif
        write(*,'("ERROR at model date (yymmddtt): ",I2,I2,I2,I4)')    &
     &           kglobe_year, kglobe_month, kglobe_day, kglobe_ts
        stop
      endif
      if (ioerr .gt. 0) then
        if (present(f90) .and. present(line)) then
          write(*,*) f90,': ', TRIM(globe_string(line)), ': ERROR: Reading file ', TRIM(sfile)
        else
          write(*,*) 'ERROR: Reading file ', TRIM(sfile)
        endif
        write(*,'("ERROR at model date (yymmddtttt): ",I2,I2,I2,I4)')    &
     &           kglobe_year, kglobe_month, kglobe_day, kglobe_ts
        stop
      endif

      return
      end subroutine error_read

#ifdef __NETCDF
!     ******************************************************************

!> \brief Error handling of NetCDF

!> \detail This routine shows a message in case of an error during the
!! different file operations of NetCDF and stops the program.

!> \param ioerr The NetCDF error status.

!     ******************************************************************

      subroutine nc_check(ioerr)
      use netcdf
      implicit none

      integer, INTENT(in) :: ioerr

      if(ioerr /= NF90_NOERR) then
        write(*,*) TRIM(nf90_strerror(ioerr))
        stop
      endif

      return
      end subroutine nc_check
#endif /* __NETCDF */

!     ******************************************************************
!     GLOBE_LANDAVG
!     ******************************************************************

!> \brief Global average of fields with land area

!> \detail This routine calculates the global average of data fields,
!! scaled by the land area at the data fields. Because of the map
!! projection used, the area of the map fields depends on its latitude.

!> \param pField   The data field to average
!> \param pAverage The global averaged of the data field

!     ******************************************************************

      subroutine globe_landavg (pField, pAverage)
      use globe_mod
      implicit none

      real, INTENT(in)    :: pField(nglobe_cgpt)
      real, INTENT(out)   :: pAverage

      integer             :: k
      real                :: pTotalArea

      __globe_mpga_to_from(globe_tmp_nGPts2,pField)

      if (kglobe_mypid == kglobe_kroot) then
        pTotalArea = SUM(pglobe_gLP_area(:))
        pAverage   = SUM(pglobe_gLP_area(:) * globe_tmp_nGPts2(1:nglobe_landpts))
        pAverage   = pAverage / pTotalArea
      endif
      __globe_mpbc(pAverage)

      return
      end subroutine globe_landavg

!     ******************************************************************
!     GLOBE_DAYLENGTH
!     ******************************************************************

!> \brief Calculates the day length for all map fields

!> \detail This routine calculates the length of the actual day for all
!! given map fields.

!> \param dpy     Number of days per year
!> \param diy     Actual day of the year continuous counted
!> \param p_lat   Field with the latitudes
!> \param cgpt    Number of field elements
!> \param zDayLen Resulting field with the day length of field elements

!     ******************************************************************

      subroutine globe_daylength (dpy, diy, p_lat, zDayLen, cgpt)
      use globe_mod
      implicit none

      integer, INTENT(in)  :: dpy, diy, cgpt
      real,    INTENT(in)  :: p_lat(cgpt)
      real,    INTENT(out) :: zDayLen(cgpt)

      integer :: i
      real, parameter :: PI2 = 6.28318530719586

      real,allocatable    :: zADay(:)
      allocate(zADay(nglobe_cgpt))

      zADay(:)  = TAN(PI2 / 360.0 * p_lat(:)) * TAN(PI2 / 360.0 * 23.4 * COS(PI2 / dpy * (diy + 10.0)))
      do i = 1 , nglobe_cgpt
        if (zADAY(i) .ge. 1.0) then
          zDayLen(i) = 0.0
        else
          if (zADAY(i) .le. -1.0) then
            zDayLen(i) = 24.0
          else
            zDayLen(i) = 2.0 * 24.0 / PI2 * ACOS(zADAY(i))
          endif
        endif
      enddo

      deallocate(zADay)

      return
      end subroutine globe_daylength

!     ******************************************************************
!     GLOBE_PLOT_MAP
!     ******************************************************************

!> \brief Output of a data field as an ASCII map

!> \detail This routine creates an ASCII map of a land points field and
!! writes it to a text file.

!> \param pField   Land points field to output as ASCII map
!> \param kFile_ID File ID of an open text file

!     ******************************************************************

      subroutine globe_plot_map (kFile_ID, pField)
      use globe_mod
      implicit none

      integer, INTENT(in) :: kFile_ID
      real,    INTENT(in) :: pField(nglobe_gpts2)

      real      :: zMaxVal, zMinVal
      integer   :: iLon, jLat, k

      character,allocatable :: iField(:), outiField(:)

      if (kglobe_mypid == kglobe_kroot) then

!       * allocate temporary fields
        allocate(iField(nglobe_gpts2))
        allocate(outiField(nglobe_nlon*nglobe_nlat))

        zMaxVal    = MAXVAL(pField(:))
        zMinVal    = MINVAL(pField(:))

        iField(:)  = CHAR(INT(10.0 * (pField(:) - zMinVal) / MAX(1.E-20, (zMaxVal - zMinVal))) + 48)
        where(iField .eq. ":") iField = "*"
        outiField(:) = "."

        write(kFile_ID,*) 'min = ', zMinVal, ' max = ', zMaxVal

        do k = 1, nglobe_landpts
          do iLon = 1, nglobe_nlon
            do jLat = 1, nglobe_nlat
              if ((jLat-1) * nglobe_nlon + iLon .eq. kglobe_gLPIDs(k)) then
                if (pField(k) .eq. zMaxVal) then
                  outiField((jLat-1) * nglobe_nlon + iLon) = "M"
                else if (pField(k) .eq. zMinVal) then
                  outiField((jLat-1) * nglobe_nlon + iLon) = "m"
                else
                  outiField((jLat-1) * nglobe_nlon + iLon) = iField(k)
                endif
              endif
            enddo
          enddo
        enddo

        do jLat = 1, nglobe_nlat
          do iLON = 1, nglobe_nlon
            write(kFile_ID,'(A1)',ADVANCE="NO") outiField((jLat-1) * nglobe_nlon + iLon)
          enddo
          write(kFile_ID,*) ""
        enddo

!       * deallocate temporary fields
        deallocate(iField)
        deallocate(outiField)

      endif

      return
      end subroutine globe_plot_map

!     ******************************************************************
!     GLOBE_OPEN_INPUT
!     ******************************************************************

!> \brief Opens the input files

!> \detail This routine is called by the different models to open and
!! initialize the input files containing the climate forcing in the SRV
!! format.

!> \param sfile File name of the input file
!> \param kfile The file ID of the input file
!> \param kdiag The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_open_input (sfile, kfile, kdiag)
      use globe_mod
      implicit none

      character(len=*), INTENT(in)    :: sfile
      integer,          INTENT(inout) :: kfile, kdiag

      character(len=80) :: filename
      integer           :: ioerr
      logical           :: exist

      if (kglobe_mypid == kglobe_kroot) then

        __diag(kdiag,'globe_open_input')

!       write(filename,'(A,A)') TRIM(sfile), ".srv"
        write(filename,'(A)') TRIM(sfile)
        inquire(FILE=TRIM(filename), EXIST=exist)
        if (.not. exist) then
          write(*,*) 'ERROR: ',TRIM(filename),': File not found'
          stop
        endif

!       * open srv-input file

        open(kfile, FILE=TRIM(filename), STATUS='old',                 &
     &              FORM='unformatted', IOSTAT=ioerr)
        call error_open_old(ioerr, filename)

        __diag(kdiag,'globe_open_input: done')

      endif  ! (kglobe_mypid == kglobe_kroot)

      return
      end subroutine globe_open_input

!     ******************************************************************
!     GLOBE_CLOSE_INPUT
!     ******************************************************************

!> \brief Closes the input files

!> \detail This routine is called by the different models to close the
!! input files.

!> \param kfile The file ID of the input file
!> \param kdiag The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_close_input (kfile, kdiag)
      use globe_mod
      implicit none

      integer, INTENT(in) :: kfile, kdiag

      if (kglobe_mypid == kglobe_kroot) then

        __diag(kdiag,'globe_close_input')

        close(kfile)

        __diag(kdiag,'globe_close_input: done')

      endif

      return
      end subroutine globe_close_input

!     ******************************************************************
!     GLOBE_REWIND
!     ******************************************************************

!> \brief Set the file pointer to the start of file

!> \detail This routine sets the file pointer to the start of the file.
!! SRV input files are read in direct access mode by using its own file
!! pointer.

!     ******************************************************************

      subroutine globe_rewind (kfile)
      use globe_mod
      implicit none

      integer :: kfile

      if (kglobe_mypid == kglobe_kroot) then
        rewind(kfile)
      endif

      return
      end subroutine globe_rewind

!     ******************************************************************
!     GLOBE_OPEN_OUTPUT
!     ******************************************************************

!> \brief Opens the output files

!> \detail This routine is called by the different models to open and
!! initialize the output files containing the results in the SRV or
!! NetCDF format.
!> In case of SRV format and output of land fields only
!! (\c kglobe_mapoutput switch in the GLOBE namelist), the land index
!! list is written as the first element  of the file.
!> In case of NetCDF output (\c kglobe_netcdf switch in the GLOBE
!! namelist), the NetCDF file (version 3 or 4) is initialized.

!> \param sfile File name of the results output file
!> \param kfile The file ID of the results output file
!> \param kdiag The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_open_output (sfile, kfile, kdiag)
      use globe_mod
#ifdef __NETCDF
      use netcdf
#endif /* __NETCDF */
      implicit none

      character(len=*), INTENT(in)    :: sfile
      integer,          INTENT(inout) :: kfile, kdiag

      character(len=80) :: filename
      integer           :: ioerr
      logical           :: exist

      if (kglobe_mypid == kglobe_kroot) then

        __diag(kdiag,'globe_open_output')

        if (kglobe_netcdf .eq. 0) then
          write(filename,'(A,A)') TRIM(sfile), ".srv"
        else
          write(filename,'(A,A)') TRIM(sfile), ".nc"
        endif

        inquire(FILE=TRIM(filename), EXIST=exist)

        if (kglobe_netcdf .eq. 0) then

!         * open srv-output file

          if (kglobe_restart .gt. 0 .and. exist) then
            open(kfile, FILE=TRIM(filename), IOSTAT=ioerr,                &
     &                  STATUS='old', FORM='unformatted', POSITION='append')
            call error_open_old(ioerr, filename)
          else
            open(kfile, FILE=TRIM(filename), IOSTAT=ioerr,              &
   &                    STATUS='replace', FORM='unformatted')
            call error_open_new(ioerr, filename)
            if (kglobe_mapoutput .eq. 0) then
              globe_tmp_nLPs(:) = real(kglobe_gLPIDs(:))
              __globe_writeoutput(kfile,globe_tmp_nLPs,kcode_ls)
            endif
          endif

#ifdef __NETCDF
        else

!         * open NetCDF-output file

          if (kglobe_restart .gt. 0 .and. exist) then
            call nc_check(nf90_open(TRIM(filename), NF90_WRITE, kfile))
          else
            if (kglobe_netcdf .eq. 1) then
              call nc_check(nf90_create(TRIM(filename), NF90_CLOBBER, kfile))
            else
              call nc_check(nf90_create(TRIM(filename), NF90_NETCDF4+NF90_CLASSIC_MODEL, kfile))
            endif
            call nc_check(nf90_enddef(kfile))
          endif
#endif /* __NETCDF */

        endif  ! (kglobe_netcdf .eq. 0)

        __diag(kdiag,'globe_open_output: done')

      endif  ! (kglobe_mypid == kglobe_kroot)

      return
      end subroutine globe_open_output

!     ******************************************************************
!     GLOBE_CLOSE_OUTPUT
!     ******************************************************************

!> \brief Closes the output files

!> \detail This routine is called by the different models to close the
!! output files.

!> \param kfile The file ID of the results output file
!> \param kdiag The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_close_output (kfile, kdiag)
      use globe_mod
#ifdef __NETCDF
      use netcdf
#endif /* __NETCDF */
      implicit none

      integer, INTENT(in) :: kfile, kdiag

      if (kglobe_mypid == kglobe_kroot) then

        __diag(kdiag,'globe_close_output')

        if (kglobe_netcdf .eq. 0) then
          close(kfile)
#ifdef __NETCDF
        else
          call nc_check(nf90_close(kfile))
#endif /* __NETCDF */
        endif

        __diag(kdiag,'globe_close_output: done')

      endif

      return
      end subroutine globe_close_output

!     ******************************************************************
!     GLOBE_OPEN_DIAG
!     ******************************************************************

!> \brief Opens the diagnostic output files

!> \detail This routine is called by the different models to open their
!! diagnostic output text files.

!> \param sfile File name of the diagnostic output file
!> \param kfile The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_open_diag (sfile, kfile)
      use globe_mod
      implicit none

      character(len=*), INTENT(in) :: sfile
      integer,          INTENT(in) :: kfile

      integer :: ioerr
      logical :: exist

      if (kglobe_mypid == kglobe_kroot .and. kglobe_diag .eq. 1) then

        inquire(FILE=TRIM(sfile), EXIST=exist)
        if (exist) then
          open(kfile, FILE=TRIM(sfile), STATUS='old',                  &
     &                IOSTAT=ioerr, FORM='formatted', POSITION='append')
          call error_open_old(ioerr,TRIM(sfile))
        else
          open(kfile, FILE=TRIM(sfile), STATUS='replace',              &
     &                IOSTAT=ioerr, FORM='formatted')
          call error_open_new(ioerr,TRIM(sfile))
        endif

      endif

      return
      end subroutine globe_open_diag

!     ******************************************************************
!     GLOBE_WRITE_DIAG
!     ******************************************************************

!> \brief Writes the diagnostic output to the file

!> \detail This routine is usually called by the different models using
!! a macro (defined in globe_macros.f90) to write the diagnostic output
!! to a text file.

!> \param kfile The file ID of the diagnostic output file
!> \param msg   The message to write to the diagnostic output file
!> \param num   A second message, appended to the first one, usually used
!! to append numeric values to the message

!     ******************************************************************

      subroutine globe_write_diag (kfile, msg, num)
      use globe_mod
      implicit none

      integer,          INTENT(in) :: kfile
      character(len=*), INTENT(in) :: msg, num

      if (kglobe_mypid == kglobe_kroot .and. kglobe_diag .eq. 1) then

        write(kfile,*) msg, TRIM(num)
        flush(kfile)

      endif

      return
      end subroutine globe_write_diag

!     ******************************************************************
!     GLOBE_CLOSE_DIAG
!     ******************************************************************

!> \brief Closes the diagnostic output files

!> \detail This routine is called by the different models to close their
!! diagnostic output text files.

!> \param kfile The file ID of the diagnostic output file

!     ******************************************************************

      subroutine globe_close_diag (kfile)
      use globe_mod
      implicit none

      integer, INTENT(in) :: kfile

      if (kglobe_mypid == kglobe_kroot .and. kglobe_diag .eq. 1) then

        close(kfile)

      endif

      return
      end subroutine globe_close_diag

!     ******************************************************************
!     GLOBE_READ_SRVNC
!     ******************************************************************

!> \brief Reads data fields from SRV files

!> \detail This routine reads data from SRV files, independent of the
!! SRV file format and converts the data to the requested variable type.
!! This routine is usually called by the interface routine
!! globe_read_srvnc_ (in globe_mod_func.F90). The output is a byte
!! field, which is converted to the output variable type during transfer
!! at the subroutine call.\n

!> The first part of the routine reads the data header and determines the
!! data format (real*4 or real*8). If the optional parameter iostat is
!! set, this part of the routine passes the file read status to it and returns.
!> The second part of the routine determines the size of the output field
!! according to the given variable type \c dtype.
!> The third part reads the data section from the file.
!> The last part converts the data to the requested variable type, and
!! transfers it to the output data field.

!> \param kfile  File ID of the input file
!> \param f90    Name of the Fortran file the routine is called from
!> \param line   Line number in the Fortran file of the call
!> \param head   The read data header
!> \param data   The read and converted data field
!> \param nbyte  Byte size of the data field
!> \param dtype  Data type of the corresponding data field variable
!> \param iostat Status message of the read file operation of the header

!     ******************************************************************

      subroutine globe_read_srvnc (kfile, f90, line, head, data, nbyte, dtype, iostat)
      use globe_mod, globe_read_srvnc_local => globe_read_srvnc
      implicit none

      !     * function declarations
#ifndef __GCC
      integer   :: fseek
#endif /* __GCC */
      integer*8 :: ftell8

      integer,             INTENT(in)  :: kfile, line
      character(len=*),    INTENT(in)  :: f90
      integer,             INTENT(out) :: head(8)
      integer,             INTENT(in)  :: nbyte
      integer,   OPTIONAL, INTENT(in)  :: dtype  ! = range(data)*kind(data)
      integer*1, OPTIONAL, INTENT(out) :: data(nbyte)
      integer,   OPTIONAL, INTENT(out) :: iostat

!     * data types for the transfer of the data to its requested type
      integer*1               :: b
      integer*1,  allocatable :: i_1(:)
      integer*2,  allocatable :: i_2(:)
      integer*4,  allocatable :: i_4(:)
      integer*8,  allocatable :: i_8(:)
      real*4,     allocatable :: r_4(:)
      real*8,     allocatable :: r_8(:)
      complex*8,  allocatable :: c_8(:)
      complex*16, allocatable :: c_16(:)

!     * type definitions for the SRV data and header field
      integer*4           :: zhead_4(8)
      integer*8           :: zhead_8(8)
      real*4, allocatable :: zdata_4(:)
      real*8, allocatable :: zdata_8(:)

      integer      :: dsize, ioerr
      integer*8    :: fpos1_8, fpos2_8
      character*80 :: sfile

      if (kglobe_mypid .ne. kglobe_kroot) return

      inquire(kfile, NAME=sfile)

!     * read the SRV header (first try integer*8 and than integer*4)

      fpos1_8 = ftell8(kfile)
      read(kfile, IOSTAT=ioerr) zhead_8
      if (ioerr .ne. 0) then
        fpos2_8 = ftell8(kfile)
#ifdef __GCC
        call fseek(kfile, fpos1_8-fpos2_8, 1, ioerr)
#else /* __GCC */
        ioerr = fseek(kfile, fpos1_8-fpos2_8, 1)
#endif /* __GCC */
        if (present(iostat) .and. (ioerr .ne. 0)) then
          iostat = ioerr
          return
        endif
        call error_read(ioerr, sfile, f90, line)
        read(kfile, IOSTAT=ioerr) zhead_4
        if (present(iostat)) then
          if (ioerr .ne. 0) then
            iostat = ioerr
            return
          endif
          head(:) = zhead_4(:)
          fpos2_8 = ftell8(kfile)
#ifdef __GCC
          call fseek(kfile, fpos1_8-fpos2_8, 1, ioerr)
#else /* __GCC */
          ioerr = fseek(kfile, fpos1_8-fpos2_8, 1)
#endif /* __GCC */
          if (ioerr .ne. 0) then
            iostat = ioerr
            return
          endif
          iostat = 0
          return
        endif
        call error_read(ioerr, sfile, f90, line)
        head(:) = zhead_4(:)
      else
        head(:) = zhead_8(:)
        if (present(iostat)) then
          fpos2_8 = ftell8(kfile)
#ifdef __GCC
          call fseek(kfile, fpos1_8-fpos2_8, 1, ioerr)
#else /* __GCC */
          ioerr = fseek(kfile, fpos1_8-fpos2_8, 1)
#endif /* __GCC */
          if (ioerr .ne. 0) then
            iostat = ioerr
            return
          endif
          iostat = 0
          return
        endif
      endif
      if (.not. present(data)) then
        fpos2_8 = ftell8(kfile)
#ifdef __GCC
        call fseek(kfile, fpos1_8-fpos2_8, 1, ioerr)
#else /* __GCC */
        ioerr = fseek(kfile, fpos1_8-fpos2_8, 1)
#endif /* __GCC */
        call error_read(ioerr, sfile, f90, line)
        return
      endif
!     * get the number of elements of the data field

      select case(dtype)
      case (range(i_1)*1)
        dsize = nbyte / 1
      case (range(i_2)*2)
        dsize = nbyte / 2
      case (range(i_4)*4)
        dsize = nbyte / 4
      case (range(i_8)*8)
        dsize = nbyte / 8
      case (range(r_4)*4)
        dsize = nbyte / 4
      case (range(r_8)*8)
        dsize = nbyte / 8
      case (range(c_8)*8)
        dsize = nbyte / 8 * 2  ! read real and imaginary part separate
      case (range(c_16)*16)
        dsize = nbyte / 16 * 2  ! read real and imaginary part separate
      end select

!     * check the correct size of data variable and file

      if (dsize .ne. head(5)*head(6)) then
        write(*,*) f90, ": ", TRIM(globe_string(line)), ": ERROR: Variable size differ to the field size at the file"
        stop
      endif

!     * read the SRV data (first try real*8 and than real*4)

      allocate(zdata_4(dsize))
      allocate(zdata_8(dsize))

      fpos1_8 = ftell8(kfile)
      read(kfile, IOSTAT=ioerr) zdata_8
      if (ioerr .ne. 0) then
        fpos2_8 = ftell8(kfile)
#ifdef __GCC
        call fseek(kfile, fpos1_8-fpos2_8, 1, ioerr)
#else /* __GCC */
        ioerr = fseek(kfile, fpos1_8-fpos2_8, 1)
#endif /* __GCC */
        call error_read(ioerr, sfile, f90, line)
        read(kfile, IOSTAT=ioerr) zdata_4
        call error_read(ioerr, sfile, f90, line)
        zdata_8(:) = zdata_4(:)
      endif

!     * convert the data to the requested data type

      select case(dtype)
      case (range(i_1)*1)
        allocate(i_1(dsize))
!        i_1(:) = IIDNNT(zdata_8(:))
        i_1(:) = NINT(zdata_8(:),1)
        data = transfer(i_1,b,nbyte)
        deallocate(i_1)
      case (range(i_2)*2)
        allocate(i_2(dsize))
!        i_2(:) = IIDNNT(zdata_8(:))
        i_2(:) = NINT(zdata_8(:),2)
        data = transfer(i_2,b,nbyte)
        deallocate(i_2)
      case (range(i_4)*4)
        allocate(i_4(dsize))
!        i_4(:) = JIDNNT(zdata_8(:))
        i_4(:) = NINT(zdata_8(:),4)
        data = transfer(i_4,b,nbyte)
        deallocate(i_4)
      case (range(i_8)*8)
        allocate(i_8(dsize))
!        i_8(:) = KIDNNT(zdata_8(:))
        i_8(:) = NINT(zdata_8(:),8)
        data = transfer(i_8,b,nbyte)
        deallocate(i_8)
      case (range(r_4)*4)
        allocate(r_4(dsize))
        r_4(:) = zdata_8(:)
        data = transfer(r_4,b,nbyte)
        deallocate(r_4)
      case (range(r_8)*8)
        allocate(r_8(dsize))
        r_8(:) = zdata_8(:)
        data = transfer(r_8,b,nbyte)
        deallocate(r_8)
      case (range(c_8)*8)
        allocate(c_8(dsize/2))
        c_8 = transfer(zdata_8,c_16,dsize/2)
        data = transfer(c_8,b,nbyte)
        deallocate(c_8)
      case (range(c_16)*16)
        allocate(c_16(dsize/2))
        c_16 = transfer(zdata_8,c_16,dsize/2)
        data = transfer(c_16,b,nbyte)
        deallocate(c_16)
      end select

      deallocate(zdata_4)
      deallocate(zdata_8)

      return
      end subroutine globe_read_srvnc

!     ******************************************************************
!     FTELL8
!     ******************************************************************

!> \brief Implementation of a 64 bit version of ftell()

!> \detail This function implements ftell8() as a wrapper to the 64 bit
!! version of ftell(), which is defined as ftelli8() in the Intel Fortran
!! Compiler, but as ftell() in the Portland Group Fortran 90/95 Compiler.

!> \param kfile File unit of the file to return the file pointer

!> \return returns the position of the file pointer

!     ******************************************************************

      integer*8 function ftell8 (kfile)
      implicit none

      integer, INTENT(in) :: kfile

#ifdef __IFORT
      integer*8 :: ftelli8

      ftell8 = ftelli8(kfile)

#else /* __IFORT */
      integer*8 :: ftell

      ftell8 = ftell(kfile)

#endif /* __IFORT */
      return
      end function ftell8
