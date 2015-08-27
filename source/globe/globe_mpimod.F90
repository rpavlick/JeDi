#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_mpimod.F90
!> \brief MPI routines for parallelization

!> \file globe_mpimod.F90
!> This file includes all routines which deal with MPI to perform
!! the inter-process data exchange, the distribution of read data, and
!! the collection of data for output.

!     ******************************************************************
!     GLOBE_MPIMOD
!     ******************************************************************


      module globe_mpimod
      use globe_mod
#ifdef __MPI
      use mpi

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                  Basic data types of MPI
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer :: mpi_itype = MPI_INTEGER4
      integer :: mpi_rtype = MPI_REAL4
      integer :: mpi_btype = MPI_BYTE

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
#endif /* __MPI */
__SEC()                        MPI handler
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer :: mpinfo  = 0

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      end module globe_mpimod

!     ******************************************************************
!     GLOBE_MPSTART
!     ******************************************************************

!> \brief Initialization of MPI

!> \detail This routine initializes MPI by setting the MPI handler and
!! writing out the processor configuration to the standard output.\n

!> \b Variables \b set: kglobe_myworld, nglobe_npro, kglobe_mypid

!     ******************************************************************

      subroutine globe_mpstart ()
      use globe_mpimod
      implicit none

      integer :: ilen,i,j
      integer :: itest
      real    :: rtest
      integer, allocatable :: node(:,:)
      character(len=80) :: yname

#ifdef __MPI
      if (kind(itest) == 8) mpi_itype = MPI_INTEGER8
      if (kind(rtest) == 8) mpi_rtype = MPI_REAL8

      call mpi_init(mpinfo)
      kglobe_myworld = MPI_COMM_WORLD

      call mpi_comm_size(kglobe_myworld, nglobe_npro, mpinfo)
      call mpi_comm_rank(kglobe_myworld, kglobe_mypid, mpinfo)

      call mpi_get_processor_name(yname, ilen, mpinfo)
      if (ilen > 40) ilen = 40

      allocate(node(40,nglobe_npro))
      node(:,:) = ICHAR(' ')
      do i = 1, ilen
        node(i,1) = ICHAR(yname(i:i))
      enddo

      if (nglobe_npro > 1) then
        call mpi_gather(node(1,1), 40, mpi_itype,                      &
     &                  node(1,1), 40, mpi_itype,                      &
     &                  kglobe_kroot, kglobe_myworld, mpinfo)
      endif

      if (kglobe_mypid == kglobe_kroot) then
        write (*,'(61("*"))')
        do j = 1, nglobe_npro
          write (*,'("Process nr.",i3," runs on ",40a1)') j-1, (char(node(i,j)),i=1,40)
        enddo
        write (*,'(61("*"))')
      endif

      deallocate (node)
#endif /* __MPI */

      return
      end subroutine globe_mpstart

!     ******************************************************************
!     GLOBE_MPSTOP
!     ******************************************************************

!> \brief Shutdown of MPI

!> \detail This routine shuts down the MPI environment.

!     ******************************************************************

      subroutine globe_mpstop ()
      use globe_mpimod
      implicit none

#ifdef __MPI
      call mpi_barrier(kglobe_myworld, mpinfo)
      call mpi_finalize(mpinfo)
#endif /* __MPI */

      return
      end subroutine globe_mpstop

!     ******************************************************************
!     GLOBE_MPCPUSYNC
!     ******************************************************************

!> \brief Waits until the execution of all processes has reached the same stage (debug)

!> \detail This routine waits until the program execution of all
!! processes has reached the same stage. (for debuging only)

!     ******************************************************************

      subroutine globe_mpcpusync ()
      use globe_mpimod
      implicit none

#ifdef __MPI
      call mpi_barrier(kglobe_myworld, mpinfo)
#endif /* __MPI */

      return
      end subroutine globe_mpcpusync

!     ******************************************************************
!     GLOBE_MPBC
!     ******************************************************************

!> \brief Broadcasts a scalar or array of any type

!> \detail This routine is used to send a scalar or array variable to all
!! processes of the program. The root process sends the data variable
!! and all other processes receive it.\n
!! The data variable is of type \c byte to make the routine variable
!! type independent, because byte is the base of all variable types.
!! The conversion of the variable is done automatically by hand over at
!! the procedure call.

!> \param pp Variable to broadcast
!> \param nb Size of the variable field in byte

!     ******************************************************************

      subroutine globe_mpbc (pp, nb)
      use globe_mpimod
      implicit none

      integer,   INTENT(in)    :: nb
      integer*1, INTENT(inout) :: pp(nb)

#ifdef __MPI
      if (nglobe_npro > 1) then
        call mpi_bcast(pp, nb, mpi_btype, kglobe_kroot, kglobe_myworld, mpinfo)
      endif
#endif /* __MPI */

      return
      end subroutine globe_mpbc

!     ******************************************************************
!     GLOBE_MPSC
!     ******************************************************************

!> \brief Scatters an array of any type

!> \detail This routine is used to split up a global variable field
!! \c pf to equal parts and distribute the parts to the variable \c pp
!! of all subprocesses of the program.\n
!! A check of the variable sizes is done during the split up of \c pf.
!! The input variable size \c nb1pf has to be a multiple of the number
!! of subprocesses and the output variable size \c nb1pp.\n
!! The data variables are of type \c byte to make the routine variable
!! type independent, because byte is the base of all variable types.
!! The conversion of the variables is done automatically by hand over at
!! the procedure call.

!> \param sfile Name of the Fortran file with the procedure call
!> \param nr    Line number of the procedure call in the Fortran file
!> \param pf    Global variable field to scatter to the subprocesses
!> \param nb1pf Size of the first dimension of the input variable
!> \param d2pf  Size of the second dimension of the input variable
!> \param pp    This variable include a part of the global variable \c pf
!> \param nb1pp Size of the first dimension of the output variable
!> \param d2pp  Size of the second dimension of the output variable

!     ******************************************************************

      subroutine globe_mpsc (sfile, nr, pf, nb1pf, d2pf, pp, nb1pp, d2pp)
      use globe_mpimod
      implicit none

      integer,          INTENT(in)  :: nr, nb1pf, d2pf, nb1pp, d2pp
      integer*1,        INTENT(in)  :: pf(nb1pf, d2pf)
      integer*1,        INTENT(out) :: pp(nb1pp, d2pp)
      character(len=*), INTENT(in)  :: sfile

      character(80) :: line, snr
      integer :: ii

!     * error checking

      if (kglobe_mypid == kglobe_kroot) then
        if (d2pf .ne. d2pp) then
          write(snr,*) nr
          write(line,*) sfile, "(", TRIM(snr), ")"
          write(*,*) "ERROR: globe_mpsc: ", TRIM(line), ":",  &
     &      "incompartible second dimension of source and destination"
          stop
        endif
        if (nb1pf .ne. nb1pp * nglobe_npro) then
          write(snr,*) nr
          write(line,*) sfile, "(", TRIM(snr), ")"
          write(*,*) "ERROR: globe_mpsc: ", TRIM(line), ":",  &
     &      "incompatible first dimension or kind of source and destination"
          stop
        endif
      endif

!     * scattering

      if (nglobe_npro == 1) then
        pp(:,:) = pf(:,:)
#ifdef __MPI
      else
        do ii = 1, d2pf
          call mpi_scatter(pf(1,ii), nb1pp, mpi_btype,                 &
     &                     pp(1,ii), nb1pp, mpi_btype,                 &
     &                     kglobe_kroot, kglobe_myworld, mpinfo)
        enddo
#endif /* __MPI */
      endif

      return
      end subroutine globe_mpsc

!     ******************************************************************
!     GLOBE_MPGA
!     ******************************************************************

!> \brief Gathers an array of any type

!> \detail This routine is used to collect the variable \c pp from all
!! subprocesses of the program and merge it to the global variable field
!! \c pf.\n
!! A check of the variable sizes is done during the collection of \c pp.
!! The output variable size \c nb1pf has to be a multiple of the number
!! of subprocesses and the input variable size \c nb1pp.\n
!! The data variables are of type \c byte to make the routine variable
!! type independent, because byte is the base of all variable types.
!! The conversion of the variables is done automatically by hand over at
!! the procedure call.

!> \param sfile Name of the Fortran file with the procedure call
!> \param nr    Line number of the procedure call in the Fortran file
!> \param pf    Global variable field to gather from the subprocesses
!> \param nb1pf Size of the first dimension of the output variable
!> \param d2pf  Size of the second dimension of the output variable
!> \param pp    This variable include a part of the global variable \c pf
!> \param nb1pp Size of the first dimension of the input variable
!> \param d2pp  Size of the second dimension of the input variable

!     ******************************************************************

      subroutine globe_mpga (sfile, nr, pf, nb1pf, d2pf, pp, nb1pp, d2pp)
      use globe_mpimod
      implicit none

      integer,          INTENT(in)  :: nr, nb1pf, d2pf, nb1pp, d2pp
      integer*1,        INTENT(out) :: pf(nb1pf, d2pf)
      integer*1,        INTENT(in)  :: pp(nb1pp, d2pp)
      character(len=*), INTENT(in)  :: sfile

      character(80) :: line, snr
      integer :: ii

!     * error checking

      if (kglobe_mypid == kglobe_kroot) then
        if (d2pf .ne. d2pp) then
          write(snr,*) nr
          write(line,*) sfile, "(", TRIM(snr), ")"
          write(*,*) "ERROR: globe_mpga: ", TRIM(line), ":",           &
     &      "incompartible second dimension of source and destination"
          stop
        endif
        if (nb1pf .ne. nb1pp*nglobe_npro) then
          write(snr,*) nr
          write(line,*) sfile, "(", TRIM(snr), ")"
          write(*,*) "ERROR: globe_mpga: ", TRIM(line), ":",           &
     &      "incompatible first dimension or kind of source and destination"
          stop
        endif
      endif

!     * gathering

      if (nglobe_npro == 1) then
        pf(:,:) = pp(:,:)
#ifdef __MPI
      else
        do ii = 1 , d2pf
          call mpi_gather(pp(1,ii), nb1pp, mpi_btype,                  &
     &                    pf(1,ii), nb1pp, mpi_btype,                  &
     &                    kglobe_kroot, kglobe_myworld, mpinfo)
        enddo
#endif /* __MPI */
      endif

      return
      end subroutine globe_mpga

!     ******************************************************************
!     GLOBE_MPREADGP
!     ******************************************************************

!> \brief Reads and scatters restart data

!> \detail This routine is used to read a variable field from the
!! restart file and scatter it to all subprocesses of the program.
!! Restart files are used to save all important data for the restart of
!! the program.\n
!! All saved data are of type \c real. For this reason the output
!! variable type is real.

!> \param kfile  File ID of the restart file
!> \param pp     Subprocess variable to get the restart data
!> \param klev   Size of the second dimension of the output variable
!> \param nCPts  Subprocess variable size
!> \param nGPts  Global variable size (size of variable in the file)
!> \param nGPts2 Global variable size, extended as multiple of \c nCPts

!     ******************************************************************

      subroutine globe_mpreadgp (kfile, pp, klev, nCPts, nGPts, nGPts2)
      use globe_mpimod
      implicit none

      integer, INTENT(in)  :: kfile, klev, nCPts, nGPts, nGPts2
      real,    INTENT(out) :: pp(nCPts,klev)

      integer :: j, ioerr

      real,allocatable :: zz(:,:)
      real,allocatable :: zp(:,:)

      allocate(zz(nGPts2, klev))
      allocate(zp(nGPts, klev))

      if (kglobe_mypid == kglobe_kroot) then
        read (kfile, IOSTAT=ioerr) zp
        __error_read(ioerr,'subroutine mpreadgp')
        zz(1:nGPts,:) = zp(:,:)
        do j= nGPts+1, nGPts2
          zz(j,:) = zz(nGPts,:)
        enddo
      endif
      if (nCPts .eq. nGPts2) then
        pp(:,:) = zz(:,:)
#ifdef __MPI
      else
        __globe_mpsc_from_to(zz,pp)
#endif /* __MPI */
      endif

      deallocate(zz)
      deallocate(zp)

      return
      end subroutine globe_mpreadgp

!     ******************************************************************
!     GLOBE_MPWRITEGP
!     ******************************************************************

!> \brief Gathers and writes restart data

!> \detail This routine is used to gather a variable from all
!! subprocesses of the program and write it to the restart file. Restart
!! files are used to save all important data for the restart of the
!! program.\n
!! All saved data are of type \c real. For this reason the input
!! variable type is real.

!> \param kfile  File ID of the restart file
!> \param pp     Subprocess variable which hold the restart data
!> \param klev   Size of the second dimension of the input variable
!> \param nCPts  Subprocess variable size
!> \param nGPts  Global variable size (size of variable at file)
!> \param nGPts2 Global variable size, extended as multiple of \c nCPts

!     ******************************************************************

! -> see function TRANSFER: Change type but maintain bit representation
! -> see function RESHAPE: Change the shape of an array
!     or function UNPACK: Unpack a rank-one array into an array of multiple dimensions

      subroutine globe_mpwritegp (kfile, pp, klev, nCPts, nGPts, nGPts2)
      use globe_mpimod
      implicit none

      integer, INTENT(in) :: kfile, klev, nCPts, nGPts, nGPts2
      real,    INTENT(in) :: pp(nCPts,klev)

      real,allocatable :: zz(:,:)
      real,allocatable :: zp(:,:)

      allocate(zz(nGPts2,klev))
      allocate(zp(nGPts,klev))

      if (nCPts .eq. nGPts2) then
        zz(:,:) = pp(:,:)
#ifdef __MPI
      else
        __globe_mpga_to_from(zz,pp)
#endif /* __MPI */
      endif
      if (kglobe_mypid == kglobe_kroot) then
        zp(:,:) = zz(1:nGPts,:)
        write(kfile) zp
      endif

      deallocate(zz)
      deallocate(zp)

      return
      end subroutine globe_mpwritegp

!     ******************************************************************
!     GLOBE_MPREADARRAY
!     ******************************************************************

!> \brief Reads input data and scatter it (not used)

!> \detail This routine reads data array maps and scatters them to all
!! processes.

!     ******************************************************************

      subroutine globe_mpreadarray (kfile, pp, nCPts, nGPts, nGPts2, code)
      use globe_mpimod
      implicit none

      integer, INTENT(in)  :: kfile, nCPts, nGPts, nGPts2, code
      real,    INTENT(out) :: pp(nCPts)

      integer :: i, j, k, ioerr
      integer :: numX, numY
      integer(kind=4) :: zhead_4(8)

      real,allocatable :: zz(:)

      allocate(zz(nGPts2))

      if (kglobe_mypid == kglobe_kroot) then
        read (kfile) zhead_4
        numX = zhead_4(5)
        numY = zhead_4(6)

!       * error checking
        if (numX .ne. nglobe_nlon .or. numY .ne. nglobe_nlat) then
          write (*,*) 'mpreadarray: *** fatal error ***'
          write (*,*) 'mpreadarray: mismatch in file dimension'
          write (*,*) 'mpreadarray: numX = ', numX, ' numY = ', numY
          write (*,*) 'mpreadarray: nLon = ', nglobe_nlon, ' nLat = ', nglobe_nlat
          write (*,*) 'mpreadarray: end'
          stop
        endif
        read (kfile, IOSTAT=ioerr) globe_tmp_lonlat_4
        __error_read(ioerr,'subroutine mpreadarray')

! -> see function PACK: Pack array into rank-one array

        do k = 1, nGPts
          i = MOD(kglobe_gLPIDs(k), nglobe_nlon)
          j = (kglobe_gLPIDs(k) - i) / nglobe_nlon + 1
          zz(k) = globe_tmp_lonlat_4((j-1) * nglobe_nlon + i)
        enddo

        do k = nGPts + 1, nGPts2
          zz(k) = zz(nGPts)
        enddo
      endif

      if (nCPts .eq. nGPts2) then
        pp(:) = zz(:)
#ifdef __MPI
      else
        __globe_mpsc_from_to(zz,pp)
#endif /* __MPI */
      endif

      deallocate(zz)

      return
      end subroutine globe_mpreadarray

!     ******************************************************************
!     GLOBE_MPWRITEARRAY
!     ******************************************************************

!> \brief Gathers data fields and writes it to the output file

!> \detail This routine writes one and two dimensional data arrays to the
!! output file in the SRV or NetCDF format, optional as a two
!! dimensional map.\n
!! Every second dimension is written as a separate data layer, labeled
!! with numbers starting from one.\n
!! After the setup of the map size, and the writing of the SRV-header in
!! case of SRV-file output, distributed data fields are gathered from
!! all subprocesses. As next step, the land only data are distributed
!! over a map, where sea points are filled with zeros.\n
!! The last step is the writing to the output file in SRV or NetCDF
!! format. Only the SRV format allows the output of land only data. The
!! NetCDF format only supports map output.

!> \param kfile File ID of the output file
!> \param data  Data field for output
!> \param ndata First dimension of the data field
!> \param nall  Overall size of the data field (over all dimensions)
!> \param code  Numerical variable identifier (1 - 9999)

!     ******************************************************************

      subroutine globe_mpwritearray (kfile, data, ndata, nall, code)
      use globe_mpimod
      implicit none

      integer, INTENT(in) :: kfile, ndata, nall, code
      real,    INTENT(in) :: data(ndata,nall/ndata)

      integer   :: nlon, nlat, level, nlevels, onelevel
      integer*8 :: rec_8
      real*4    :: r_4

!     * number of levels in the data field
      nlevels = nall / ndata

!     * set level (used in the header) to zero for output of only one level field
      onelevel = 0
      if (nlevels == 1) onelevel = 1

!     * loop through all levels
      do level = 1, nlevels

!       * write the SRV header to the output file in case of no NetCDF output

        if (kglobe_mypid == kglobe_kroot .and. kglobe_netcdf == 0) then
!         * set correct field size to header in case of land points output only
          if (kglobe_mapoutput == 0) then
            nlon = nglobe_landpts
            nlat = 1
          else
            nlon = nglobe_nlon
            nlat = nglobe_nlat
          endif
!         * create and write the header
          call globe_mpwritesrvheader(kfile, level - onelevel, code, nlon, nlat)
        endif

!       * fill global field with different kind of input data

        if (ndata .eq. nglobe_gpts2) then
          if (kglobe_mypid == kglobe_kroot) globe_tmp_nGPts2(:) = data(:,level)
        elseif (ndata .eq. nglobe_landpts) then
          if (kglobe_mypid == kglobe_kroot) globe_tmp_nGPts2(1:nglobe_landpts) = data(:,level)
#ifdef __MPI
        elseif (ndata .eq. nglobe_cgpt) then
          __globe_mpga_to_from(globe_tmp_nGPts2,data(:,level))
#endif /* __MPI */
        elseif (ndata .eq. nglobe_nlon*nglobe_nlat) then
          if (kglobe_mypid == kglobe_kroot) globe_tmp_lonlat(:) = data(:,level)
        else
          if (kglobe_mypid == kglobe_kroot) then
            write(*,*) "ERROR: globe_mpwritearray: code ",code,": Wrong array size"
            stop
          endif
        endif

!       * fill up the whole map

        if (kglobe_mypid == kglobe_kroot) then
          if (ndata .ne. nglobe_nlon*nglobe_nlat) then
            globe_tmp_lonlat(:) = 0.0
            globe_tmp_lonlat(kglobe_gLPIDs(:)) = globe_tmp_nGPts2(1:nglobe_landpts)
          endif
        endif

!       * write the data to the output file

        if (kglobe_mypid == kglobe_kroot) then
          if (kglobe_netcdf == 0) then
            if (kglobe_mapoutput == 0) then
              globe_tmp_nLPs_4(:) = globe_tmp_lonlat(kglobe_gLPIDs(:))
              write(kfile) globe_tmp_nLPs_4
            else
              globe_tmp_lonlat_4(:) = globe_tmp_lonlat(:)
              write(kfile) globe_tmp_lonlat_4
            endif
#ifdef __NETCDF
          else
            call globe_mpwritearray_nc(kfile, level, nlevels, code)
#endif /* __NETCDF */
          endif
        endif

      enddo  ! level

      return
      end subroutine globe_mpwritearray

!     ******************************************************************
!     GLOBE_MPWRITESRVHEADER
!     ******************************************************************

!> \brief Routine to create and write the SRV-header for output

!> \detail This routine creates the SRV-header and writes it to the output
!! file. The header elements 7 and 8 are not predefined by the SRV
!! format and therefore used internally to save the actual time step and
!! the map size (used in case of land only output).

!> \param kfile File ID of the output file
!> \param level Value of the second dimension of the data field
!> \param code  Numerical variable identifier (1 - 9999)
!> \param lon   Number of elements in longitudinal direction
!> \param lat   Number of elements in latitudinal direction

!     ******************************************************************

      subroutine globe_mpwritesrvheader (kfile, level, code, lon, lat)
      use globe_mod
      implicit none

      integer, INTENT(in) :: kfile, level, code, lon, lat

      integer(kind=4) :: zhead_4(8)
      integer         :: date

      date = __date_string()

      zhead_4(1) = code
      zhead_4(2) = level
      zhead_4(3) = date
      zhead_4(4) = 0                                   ! time
      zhead_4(5) = lon
      zhead_4(6) = lat
      zhead_4(7) = kglobe_timestep                     ! disp1 (free)
      zhead_4(8) = (nglobe_nlon * 10000) + nglobe_nlat ! disp2 (free)
      write(kfile) zhead_4

      return
      end subroutine globe_mpwritesrvheader

#ifdef __NETCDF
!     ******************************************************************
!     GLOBE_MPWRITEARRAY_NC
!     ******************************************************************

!> \brief Writes the output for the NetCDF format

!> \detail This routine writes the output to a file in NetCDF format.\n
!! The first step is to check the definition of all axes of the map
!! output, get the actual time step of the file, and create one if no
!! valid definition could be found. The check of the map axes is done
!! every time, because the definitions are created if needed. This is
!! necessary, because the number of levels changes for different
!! variables, and every number of levels needs its own definition.\n
!! The second step is the determination of the variable ID, and its
!! creation if the variable is written the first time.\n
!! At the last step the data is written to the last time position in
!! the NetCDF file.\n
!! The data field is transferred by the global variable \c globe_tmp_lonlat.

!> \param kfile  File ID of the output file
!> \param level  The actual level
!> \param levels The whole number of levels (second dimension of data array)
!> \param code   Numerical variable identifier

!     ******************************************************************

      subroutine globe_mpwritearray_nc (kfile, level, levels, code)
      use globe_mpimod
      use netcdf
      implicit none

      integer, INTENT(in) :: kfile, level, levels, code

      integer :: stat_var, var_id, time_pos
      integer :: start0(3), ncount0(3), dim_ids0(3)
      integer :: start(4), ncount(4), dim_ids(4)
      character(len=10) :: var_name

!     * check the correct axes definition of the variable
!       (returns actual time position in the file and the IDs of the axes)

      if (levels == 1) then
        call globe_nc_check_var_axes(kfile, time_pos, dim_ids0, 3, 0)
      else
        call globe_nc_check_var_axes(kfile, time_pos, dim_ids, 4, levels)
      endif

!     * get or define the variable ID

      write(var_name,'("code",i6.6)') code
      stat_var = nf90_inq_varid(kfile, var_name, var_id)

      if (stat_var /= NF90_NOERR) then
        call nc_check(nf90_redef(kfile))
        if (levels == 1) then
          call nc_check(nf90_def_var(kfile, var_name, NF90_REAL, dim_ids0, var_id))
        else
          call nc_check(nf90_def_var(kfile, var_name, NF90_REAL, dim_ids, var_id))
        endif
        if (kglobe_netcdf .gt. 2) then
          call nc_check(nf90_def_var_deflate(kfile, var_id, 1, 1, 1))
        endif
        call nc_check(nf90_enddef(kfile))
      endif

!     * write the data at last time position

      if (levels == 1) then
        start0 = (/ 1, 1, time_pos /)
        ncount0 = (/ nglobe_nlon, nglobe_nlat, 1 /)
      else
        start = (/ 1, 1, level, time_pos /)
        ncount = (/ nglobe_nlon, nglobe_nlat, 1, 1 /)
      endif

      if (levels == 1) then
        call nc_check(nf90_put_var(kfile, var_id, globe_tmp_lonlat, start0, ncount0))
      else
        call nc_check(nf90_put_var(kfile, var_id, globe_tmp_lonlat, start, ncount))
      endif

      return
      end subroutine globe_mpwritearray_nc

!     ******************************************************************
!     GLOBE_NC_CHECK_VAR_AXES
!     ******************************************************************

!> \brief Checks or creates axes of the arrays for the output to NetCDF

!> \detail This routine checks and sets up definitions for the axes of
!! the maps for the NetCDF output.\n
!! The first step is to get or define the dimension IDs of the axes
!! variables (longitude, latitude, time, level), and the creation of the
!! returned dimension IDs.\n
!! The second step is to check or define the longitude, latitude, time,
!! and level axes variables itself.\n
!! The last step is to check and set the actual time at the last time
!! position in the NetCDF file.

!> \param kfile    File ID of the output file
!> \param time_pos Actual time variable position in the NetCDF file
!> \param dim_ids  Dimension IDs of the axes
!> \param dim1     Number of dimensions (3 or 4, with levels)
!> \param levels   Number of levels (second dimension of data field)

!     ******************************************************************

      subroutine globe_nc_check_var_axes (kfile, time_pos, dim_ids, dim1, levels)
      use globe_mpimod
      use netcdf
      implicit none

      integer, INTENT(in)  :: kfile, dim1, levels
      integer, INTENT(out) :: time_pos
      integer, INTENT(out) :: dim_ids(dim1)

      integer :: stat_lon, stat_lat, stat_lev, stat_time
      integer :: lon_id, lat_id, lev_id, time_id, var_id, tmp_id
      integer :: ytime(1), vtime(1)
      character*7 :: lname

!     * determine the label for the level

      write(lname,'("lev_",I03)') levels

!     * get or define the dimension IDs of the variables (axes)

      stat_lon  = nf90_inq_dimid(kfile, "lon", lon_id)
      stat_lat  = nf90_inq_dimid(kfile, "lat", lat_id)
      stat_time = nf90_inq_dimid(kfile, "time", time_id)

      if (stat_lon /= NF90_NOERR .and. stat_lat /= NF90_NOERR .and. stat_time /= NF90_NOERR) then
        call nc_check(nf90_redef(kfile))
        call nc_check(nf90_def_dim(kfile, "lon", nglobe_nlon, lon_id))
        call nc_check(nf90_def_dim(kfile, "lat", nglobe_nlat, lat_id))
        call nc_check(nf90_def_dim(kfile, "time", NF90_UNLIMITED, time_id))
        call nc_check(nf90_enddef(kfile))
      else
        call nc_check(stat_lon)
        call nc_check(stat_lat)
        call nc_check(stat_time)
      endif

      if (levels .gt. 0) then
        stat_lev = nf90_inq_dimid(kfile, lname, lev_id)
        if (stat_lev /= NF90_NOERR) then
          call nc_check(nf90_redef(kfile))
          call nc_check(nf90_def_dim(kfile, lname, levels, lev_id))
          call nc_check(nf90_enddef(kfile))
        else
          call nc_check(stat_lev)
        endif
      endif

      if (levels == 0) then
        dim_ids = (/ lon_id, lat_id, time_id /)
      else
        dim_ids = (/ lon_id, lat_id, lev_id, time_id /)
      endif

!     * check or define the longitude, latitude, and time variables (axes)

      stat_lon = nf90_inq_varid(kfile, "lon", tmp_id)
      stat_lat = nf90_inq_varid(kfile, "lat", tmp_id)
      stat_time = nf90_inq_varid(kfile, "time", var_id)
      if (stat_lon /= NF90_NOERR .and. stat_lat /= NF90_NOERR .and. stat_time /= NF90_NOERR) then
        call globe_nc_set_axes_values(kfile, dim_ids, dim1, 0)
      else
        call nc_check(stat_lon)
        call nc_check(stat_lat)
        call nc_check(stat_time)
      endif

      if (levels .gt. 0) then
        stat_lev = nf90_inq_varid(kfile, lname, tmp_id)
        if (stat_lev /= NF90_NOERR) then
          call globe_nc_set_axes_values(kfile, dim_ids, dim1, levels)
        else
          call nc_check(stat_lev)
        endif
      endif

!     * check and set actual time at last time position

      ytime(1) = __date_string()
      call nc_check(nf90_inquire_dimension(kfile, time_id, LEN=time_pos))
      call nc_check(nf90_inq_varid(kfile, "time", var_id))
      call nc_check(nf90_get_var(kfile, var_id, vtime, start=(/time_pos/), count=(/1/)))
      if (ytime(1) .ne. vtime(1)) then
        call nc_check(nf90_put_var(kfile, var_id, ytime, start=(/time_pos+1/), count=(/1/)))
      endif

      call nc_check(nf90_inquire_dimension(kfile, time_id, LEN=time_pos))

      return
      end subroutine globe_nc_check_var_axes

!     ******************************************************************
!     GLOBE_NC_SET_AXES_VALUES
!     ******************************************************************

!> \brief Creates and sets values of newly created array axes for the NetCDF output

!> \detail This routine defines the values of newly created axes for the
!! NetCDF file. The definition is done separately for the levels.\n
!! In the first part the dimension IDs are assigned.\n
!! The second part defines the variables (longitude, latitude, time, and
!! level. This is done in the redefinition mode of NetCDF to modify the
!! file header.\n
!! The last part fills the axes variables with values. The longitude and
!! latitude axes are equally filled with nglobe_nlon and nglobe_nlat
!! values. The time axes is an open axes and only initialized with the
!! actual time. This axes is growing with the addition of data during
!! following time steps.

!> \param kfile   File ID of the output file
!> \param dim_ids Dimension IDs of the axes
!> \param dim1    Number of dimensions (3 or 4, with levels)
!> \param levels  Number of levels (second dimension of data field)

!     ******************************************************************

      subroutine globe_nc_set_axes_values (kfile, dim_ids, dim1, levels)
      use globe_mpimod
      use netcdf
      implicit none

      integer, INTENT(in) :: kfile, dim1, levels
      integer, INTENT(in) :: dim_ids(dim1)

      integer :: lon_id, lat_id, lev_id, time_id
      integer :: var_id, var1_id, var2_id, var3_id
      integer :: ytime(1), ii
      character*7 :: lname
      real,allocatable :: lonlats(:)

      allocate(lonlats(nglobe_nlon))

!     * determine the label for the level

      write(lname,'("lev_",I03)') levels

!     * assign IDs

      lon_id  = dim_ids(1)
      lat_id  = dim_ids(2)
      if (dim1 == 3) then
        time_id = dim_ids(3)
      else
        lev_id  = dim_ids(3)
        time_id = dim_ids(4)
      endif

      call nc_check(nf90_redef(kfile))

      if (levels == 0) then

!       * define the longitude variable

        call nc_check(nf90_def_var(kfile, "lon", NF90_DOUBLE, lon_id, var1_id))
        call nc_check(nf90_put_att(kfile, var1_id, "long_name", "longitude"))
        call nc_check(nf90_put_att(kfile, var1_id, "units", "degrees_east"))
        call nc_check(nf90_put_att(kfile, var1_id, "standard_name", "longitude"))

!       * define the latitude variable

        call nc_check(nf90_def_var(kfile, "lat", NF90_DOUBLE, lat_id, var2_id))
        call nc_check(nf90_put_att(kfile, var2_id, "long_name", "latitude"))
        call nc_check(nf90_put_att(kfile, var2_id, "units", "degrees_north"))
        call nc_check(nf90_put_att(kfile, var2_id, "standard_name", "latitude"))

!       * define the time variable

        call nc_check(nf90_def_var(kfile, "time", NF90_DOUBLE, time_id, var_id))
        call nc_check(nf90_put_att(kfile, var_id, "units", "day as %Y%m%d.%f"))
        call nc_check(nf90_put_att(kfile, var_id, "calendar", "proleptic_gregorian"))

      else

!       * define the level variable

        call nc_check(nf90_def_var(kfile, lname, NF90_DOUBLE, lev_id, var3_id))
        call nc_check(nf90_put_att(kfile, var3_id, "long_name", "level"))
        call nc_check(nf90_put_att(kfile, var3_id, "units", "count"))
        call nc_check(nf90_put_att(kfile, var3_id, "standard_name", "level"))

      endif

      call nc_check(nf90_enddef(kfile))

      if (levels == 0) then

!       * fill the longitude variable

        do ii = 1, nglobe_nlon
          lonlats(ii) = 360.0 / nglobe_nlon * (real(ii) - 0.5)
        enddo
        call nc_check(nf90_put_var(kfile, var1_id, lonlats(1:nglobe_nlon)))

!       * fill the latitude variable

        do ii = 1, nglobe_nlat
          lonlats(ii) = 180.0 / nglobe_nlat * (real(nglobe_nlat - ii) + 0.5) - 90.0
        enddo
        call nc_check(nf90_put_var(kfile, var2_id, lonlats(1:nglobe_nlat)))

!       * set the time variable

        ytime(1) = __date_string()
        call nc_check(nf90_put_var(kfile, var_id, ytime))

      else

!       * fill the level variable

        do ii = 1, levels
          lonlats(ii) = real(ii)
        enddo
        call nc_check(nf90_put_var(kfile, var3_id, lonlats(1:levels)))

      endif

      deallocate(lonlats)

      return
      end subroutine globe_nc_set_axes_values
#endif /* __NETCDF */
