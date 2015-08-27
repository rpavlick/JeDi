#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_fio.F90
!> \brief I/O-routines of GLOBE

!> \file globe_fio.F90
!> This file includes the input/output routines of GLOBE which read the
!! model configuration input, read and create initialization data, and
!! write out configuration, read and create initialization data, and
!! write out the process status.

!     ******************************************************************
!     GLOBE_READ_NAMELIST
!     ******************************************************************

!> \brief Reads the GLOBE namelist file

!> \detail This routine reads the namelist parameter file of GLOBE and
!! broadcasts the parameters to all processors.\n

!> \b Variables \b set: kglobe_firstyear, kglobe_lastyear,
!! kglobe_yrsskip, kglobe_yrsout, kglobe_status, kglobe_diag,
!! kglobe_restart, kglobe_mapoutput, kglobe_netcdf

!     ******************************************************************

      subroutine globe_read_namelist ()
      use globe_mod
      implicit none

      integer :: ioerr

!     * definition of variable parameters
      namelist /globepar/ kglobe_firstyear,    &
     &                    kglobe_lastyear,     &
     &                    kglobe_yrsskip,      &
     &                    kglobe_yrsout,       &
     &                    kglobe_status,       &
     &                    kglobe_diag,         &
     &                    kglobe_restart,      &
     &                    kglobe_mapoutput,    &
     &                    kglobe_netcdf,       &
     &                    kglobe_spt,          &
     &                    kglobe_spt_gpid,     &
     &                    pglobe_spt_lon,      &
     &                    pglobe_spt_lat,      &
     &                    pglobe_one_deltalon, &
     &                    pglobe_one_deltalat


!     * input of variable parameters from the parameter file

      if (kglobe_mypid == kglobe_kroot) then
        open(kglobe_fnamelist, FILE=sglobe_fnamelist, FORM='formatted', &
     &                         STATUS='old', IOSTAT=ioerr)
        if (ioerr .eq. 0) then
          read(kglobe_fnamelist, globepar)
          close(kglobe_fnamelist)
#ifndef __NETCDF
          kglobe_netcdf = 0
#endif /* __NETCDF */
        endif

!       * open diagnostic and status output files
        call globe_open_diag(sglobe_fdiag, kglobe_fdiag)
        __diag(kglobe_fdiag,'globe_init')
        __diag(kglobe_fdiag,'globe_read_namelist')

        if (kglobe_diag .eq. 1) then
          if (ioerr .eq. 0) then
            write(kglobe_fdiag,'(" *******************************************")')
            write(kglobe_fdiag,'(" * Namelist GLOBEPAR from ",A)') sglobe_fnamelist
            write(kglobe_fdiag,'(" *******************************************")')
          else
            write(kglobe_fdiag,'(" *******************************************")')
            write(kglobe_fdiag,'(" * ERROR reading file ",A)') sglobe_fnamelist
            write(kglobe_fdiag,'(" * Using default values")')
            write(kglobe_fdiag,'(" *******************************************")')
          endif
          write(kglobe_fdiag,globepar)
          flush(kglobe_fdiag)
        endif
      endif

!     * single point version works only with one processor

      if (kglobe_spt .ne. 0 .and. nglobe_npro .gt. 1) then
        if (kglobe_mypid == kglobe_kroot)                              &
     &    write(*,*) 'ERROR: The single point version works only with one processor'
        stop
      endif

!     * broadcast parameters

      __globe_mpbc(kglobe_firstyear) ! last year to simulate
      __globe_mpbc(kglobe_lastyear)  ! last year to simulate
      __globe_mpbc(kglobe_yrsskip)   ! number of years to skip output
      __globe_mpbc(kglobe_yrsout)    ! number years to output

      __globe_mpbc(kglobe_status)    ! flag for status output
      __globe_mpbc(kglobe_diag)      ! diagnostic flag
      __globe_mpbc(kglobe_restart)   ! restart flag

!     * end

      __diag(kglobe_fdiag,'globe_read_namelist: done')

      return
      end subroutine globe_read_namelist

!     ******************************************************************
!     GLOBE_READ_SURFACE_FILE
!     ******************************************************************

!> \brief Reads the surface parameter file

!> \detail This routine reads the surface parameter file and initializes
!! all included data fields. It also reads and checks the land mask
!! file, and also sets up the requested land field of the single point
!! version of the model.\n

!> \b Variables \b set: nglobe_nlon, nglobe_nlat, nglobe_landpts,
!! kglobe_gLPIDs, pglobe_gLP_elevation,
!> pglobe_gLP_paw, pglobe_gLP_numnn,
!! pglobe_gLP_idnn, pglobe_gLP_distnn, kglobe_lpt_nn,
!> pglobe_longitude, pglobe_latitude, pglobe_gLP_lon, pglobe_gLP_lat,
!! pglobe_gLP_len, pglobe_gLP_wid, pglobe_gLP_area, pglobe_gLP_topograd\n

!> \b Variables \b modified: kglobe_spt_gpid

!     ******************************************************************

    subroutine globe_read_surface_file ()
    use globe_mod
    implicit none

    integer :: ihead(8), ioerr, done, i, j, k
    real    :: reslon, reslat
    logical :: fexist
    logical :: is_elevation
    logical :: is_paw

    integer(kind=4)           :: mhead_4(8)
    real(kind=4), allocatable :: zmask_4(:)

    __diag(kglobe_fdiag,'globe_read_surface_file')

    if (kglobe_mypid == kglobe_kroot) then

!       * open surface file

      call globe_open_input(sglobe_fsurface, kglobe_fsurface, kglobe_fdiag)

!       * read header to get the array size

      __globe_read_srvnc(kglobe_fsurface,ihead,IOSTAT=ioerr)
      nglobe_nlon    = ihead(5)
      nglobe_nlat    = ihead(6)
      nglobe_landpts = ihead(8)
      call globe_rewind(kglobe_fsurface)

      allocate(zmask_4(nglobe_nlon*nglobe_nlat))
      zmask_4 = 0.0

      if (kglobe_spt .gt. 0) then
        nglobe_landpts = 1
      else

!         * if the a land mask file exists -> get the number of land points

        inquire(FILE=sglobe_flandmask, EXIST=fexist)
        if (fexist) then
          __diag(kglobe_fdiag,'globe_read_surface_file: land mask found')
          call globe_open_input(sglobe_flandmask, kglobe_flandmask, kglobe_fdiag)
          __globe_read_srvnc(kglobe_flandmask,mhead_4,zmask_4)
          call globe_close_input(kglobe_flandmask, kglobe_fdiag)
          nglobe_landpts = INT(SUM(MIN(ABS(REAL(zmask_4(:))), 1.0)) + 0.5)
        endif

      endif  ! (kglobe_spt .gt. 0) else

    endif  ! (kglobe_mypid == kglobe_kroot)

!     * broadcast parameters

    __globe_mpbc(nglobe_landpts)
    __globe_mpbc(nglobe_nlon)
    __globe_mpbc(nglobe_nlat)

!     * allocate land point fields

    call globe_landsea_alloc

    if (kglobe_mypid == kglobe_kroot) then
      is_elevation    = .FALSE.
      is_paw       = .FALSE.

!       * read the longitude and latitude fields from the surface file

      __diag(kglobe_fdiag,'globe_read_surface_file: read long lat fields')

      ioerr = 0
      do while (ioerr .eq. 0)
        __globe_read_srvnc(kglobe_fsurface,ihead,IOSTAT=ioerr)
        if (ioerr .eq. 0) then
          __globe_read_srvnc(kglobe_fsurface,ihead,globe_tmp_lonlat)

          SELECT CASE (ihead(1))
          CASE (kcode_longitude)
            pglobe_longitude(:) = globe_tmp_lonlat(:)
          CASE (kcode_latitude)
            pglobe_latitude(:) = globe_tmp_lonlat(:)
          END SELECT
        endif
      enddo
      call globe_rewind(kglobe_fsurface)

!       * read the index mask and check the single point

      __diag(kglobe_fdiag,'globe_read_surface_file: read indexes')

      ioerr = 0
      done = 0
      do while (ioerr == 0 .and. done == 0)
        __globe_read_srvnc(kglobe_fsurface,ihead,IOSTAT=ioerr)
        if (ioerr == 0) then
          __globe_read_srvnc(kglobe_fsurface,ihead,globe_tmp_lonlat)

          if (ihead(1) == kcode_index) then
            if (kglobe_spt .gt. 0) then

!               * check if there is land at the given index

              if (globe_tmp_lonlat(kglobe_spt_gpid) .eq. 0.0) then
                write(0,*) 'ERROR: globe_read_surface_file: There is no land at the given coordinates'
                stop
              endif

              kglobe_gLPIDs(1) = kglobe_spt_gpid
              j = 1
            else  ! (kglobe_spt .gt. 0)

!               * get the land point index field

              if (SUM(zmask_4(:)) .gt. 0.0) then
                globe_tmp_lonlat(:) = globe_tmp_lonlat(:) * zmask_4(:)
              endif
              j = 0
              do i = 1, nglobe_nlon*nglobe_nlat
                if (globe_tmp_lonlat(i) .ne. 0.0) then
                  j = j + 1
                  kglobe_gLPIDs(j) = i
                endif
              enddo

            endif  ! (kglobe_spt .gt. 0) else
            done = 1
          endif  ! (ihead(1) == kcode_index)

        endif  ! (ioerr .eq. 0)
      enddo  ! while (ioerr .eq. 0)
      call globe_rewind(kglobe_fsurface)

!       * stop if no index field is found or something is wrong with the indexes

      if (done == 0) then
        write(0,*) "ERROR: globe_read_surface_file: Index mask not", &
   &               " found in the created file ",TRIM(sglobe_fsurface)
        stop
      else
        if (j .ne. nglobe_landpts) then
          write(0,*) "ERROR: globe_read_surface_file: Wrong land masking."
          write(0,*) "Maybe the land masking file cover ocean fields."
          stop
        endif
      endif

!       * read all other fields from the surface file

      __diag(kglobe_fdiag,'globe_read_surface_file: read surface fields')

      ioerr = 0
      do while (ioerr .eq. 0)
        __globe_read_srvnc(kglobe_fsurface,ihead,IOSTAT=ioerr)
        if (ioerr .eq. 0) then
          __globe_read_srvnc(kglobe_fsurface,ihead,globe_tmp_lonlat)

          SELECT CASE (ihead(1))
          CASE (kcode_elevation)
            is_elevation = .TRUE.
            pglobe_gLP_elevation(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_paw)
           is_paw = .TRUE.
           pglobe_gLP_paw(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_longitude)
            pglobe_gLP_lon(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_latitude)
            pglobe_gLP_lat(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_length)
            pglobe_gLP_len(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_width)
            pglobe_gLP_wid(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_area)
            pglobe_gLP_area(:) = globe_tmp_lonlat(kglobe_gLPIDs)
          CASE (kcode_topography)
            pglobe_gLP_topograd(:) = globe_tmp_lonlat(kglobe_gLPIDs)
         CASE (kcode_numnn)
           pglobe_gLP_numnn(:) = globe_tmp_lonlat(kglobe_gLPIDs)
         CASE (kcode_idnn)
           pglobe_gLP_idnn(:,ihead(2)) = globe_tmp_lonlat(kglobe_gLPIDs)
         CASE (kcode_distnn)
           pglobe_gLP_distnn(:,ihead(2)) = globe_tmp_lonlat(kglobe_gLPIDs)
          END SELECT
        endif
      enddo

      call globe_close_input(kglobe_fsurface, kglobe_fdiag)

!       * set not read fields to its default value

      if (.not. is_elevation)    pglobe_gLP_elevation(:) = pglobe_gLP_elevation0
      if (.not. is_paw)       pglobe_gLP_paw(:)          = pglobe_gLP_paw0

!       * localize nearest neighbors as grid points in land points numbering

      kglobe_lpt_nn(:,:) = 0
      do i = 1, nglobe_landpts
        do j = 1, 8
          do k = 1, nglobe_landpts
            if (pglobe_gLP_idnn(i,j) .eq. kglobe_gLPIDs(k)) kglobe_lpt_nn(i,j) = k
          enddo
        enddo
      enddo

    endif  ! (kglobe_mypid == kglobe_kroot)

    __deallocate(zmask_4)

    __diag(kglobe_fdiag,'globe_read_surface_file: done')

    return
    end subroutine globe_read_surface_file

!     ******************************************************************
!     GLOBE_READ_LANDPOINT
!     ******************************************************************

!> \brief Reads the land points file

!> \detail This routine reads the land points file, and initializes the
!! corresponding data fields. It also sets up the requested land field
!! of the single point version of the model.\n

!> \b Variables \b set: nglobe_nlon, nglobe_nlat, nglobe_landpts,
!! kglobe_gLPIDs, pglobe_gLP_elevation,
!> pglobe_gLP_paw, pglobe_gLP_numnn,
!! pglobe_gLP_idnn, pglobe_gLP_distnn, kglobe_lpt_nn,
!> pglobe_gLP_lon, pglobe_gLP_lat, pglobe_gLP_len, pglobe_gLP_wid,
!! pglobe_gLP_area, pglobe_gLP_topograd\n

!> \b Variables \b modified: kglobe_spt_gpid

!     ******************************************************************

      subroutine globe_read_landpoint ()
      use globe_mod
      implicit none

      integer :: i, j, k, ioerr
      real    :: reslon, reslat
      integer :: zn_landpts

      character(80) :: GPfileformat

      __diag(kglobe_fdiag,'globe_read_landpoint')

!     * read the land points file

      if (kglobe_mypid == kglobe_kroot) then
        open(kglobe_fgridpoint, FILE=sglobe_fgridpoint, STATUS='old',  &
     &                          FORM='formatted', IOSTAT=ioerr)
        call error_open_old(ioerr, sglobe_fgridpoint)
        read(kglobe_fgridpoint,'(3I6)', IOSTAT=ioerr)                 &
     &       nglobe_landpts, nglobe_nlon, nglobe_nlat
        __error_read(ioerr, sglobe_fgridpoint)
        zn_landpts = nglobe_landpts

!       * single point version -> compute the gridpoint-ID from coordinates

        if (kglobe_spt .gt. 0) then
          if (kglobe_spt_gpid .eq. 0) then
            reslon = 360.0 / nglobe_nlon
            reslat = 180.0 / nglobe_nlat
            if (pglobe_spt_lon < 0) then
              i = INT((pglobe_spt_lon + 360.0) / reslon + 0.5 + 0.5)
            else
              i = INT(pglobe_spt_lon / reslon + 0.5 + 0.5)
            endif
            j = INT(ABS(pglobe_spt_lat -  90.0) / reslat + 0.5)
            kglobe_spt_gpid = j * nglobe_nlon + i
          endif
          nglobe_landpts = 1
        endif

      endif

!     * broadcast parameters

      __globe_mpbc(nglobe_landpts)
      __globe_mpbc(nglobe_nlon)
      __globe_mpbc(nglobe_nlat)

!     * allocate land point fields

      call globe_landsea_alloc

      if (kglobe_mypid == kglobe_kroot) then

        i = 1
        do j = 1, zn_landpts
!         * single point version -> index i always 1
          if (kglobe_spt == 0) i = j
!         * single point version -> read up to the right line (i than set to 0)
          if (i .gt. 0) then
            GPfileformat = '(I6,f10.4,f9.4,2f10.1,f10.1,f8.1,f8.5,I6,8(I6,f8.1))'
            read(kglobe_fgridpoint,GPfileformat, IOSTAT=ioerr)         &
     &        kglobe_gLPIDs(i), pglobe_gLP_lon(i), pglobe_gLP_lat(i),  &
     &        pglobe_gLP_len(i), pglobe_gLP_wid(i),                    &
     &        pglobe_gLP_area(i), pglobe_gLP_elevation(i),             &
     &        pglobe_gLP_paw(i),                                       &
     &        pglobe_gLP_numnn(i),                                     &
     &        pglobe_gLP_idnn(i,1), pglobe_gLP_distnn(i,1),            &
     &        pglobe_gLP_idnn(i,2), pglobe_gLP_distnn(i,2),            &
     &        pglobe_gLP_idnn(i,3), pglobe_gLP_distnn(i,3),            &
     &        pglobe_gLP_idnn(i,4), pglobe_gLP_distnn(i,4),            &
     &        pglobe_gLP_idnn(i,5), pglobe_gLP_distnn(i,5),            &
     &        pglobe_gLP_idnn(i,6), pglobe_gLP_distnn(i,6),            &
     &        pglobe_gLP_idnn(i,7), pglobe_gLP_distnn(i,7),            &
     &        pglobe_gLP_idnn(i,8), pglobe_gLP_distnn(i,8)
            __error_read(ioerr, sglobe_fgridpoint)
!           * single point version -> stop reading, if LPID found
            if (kglobe_spt .gt. 0 .and. kglobe_gLPIDs(i) .eq. kglobe_spt_gpid) i = 0
          endif
        enddo

        close(kglobe_fgridpoint)

!       * single point version -> was there a land point at the given coordinates ?

        if (kglobe_spt .gt. 0 .and. i .eq. 1) then
          write(*,*) 'ERROR: There is no land at the given coordinates'
          stop
        endif

#ifndef __JEDI
!       * localize nearest neighbors as grid points in land points numbering

        kglobe_lpt_nn(:,:) = 0
        do i=1,nglobe_landpts
          do j=1,8
            do k=1,nglobe_landpts
              if (pglobe_gLP_idnn(i,j) .eq. kglobe_gLPIDs(k)) kglobe_lpt_nn(i,j) = k
            enddo
          enddo
        enddo
#endif /* __JEDI */

      endif

!     * end

      __diag(kglobe_fdiag,'globe_read_landpoint: done')

      return
      end subroutine globe_read_landpoint

!     ******************************************************************
!     GLOBE_READ_RESTART
!     ******************************************************************

!> \brief Reads the restart data of GLOBE

!> \detail This routine calls the subroutines to read the restart data
!! fields in case the model run is continued.\n

!> \b Variables \b set:
!! - \b global \b variables: nglobe_nlon, nglobe_nlat, nglobe_landpts,
!! kglobe_year, kglobe_timestep, kglobe_gLPIDs, pglobe_gLP_lon,
!! pglobe_gLP_lat, pglobe_gLP_elevation, pglobe_gLP_len,
!> pglobe_gLP_paw, pglobe_gLP_numnn, pglobe_gLP_idnn,
!! pglobe_gLP_distnn, kglobe_lpt_nn,
!> pglobe_gLP_wid, pglobe_gLP_area, pglobe_gLP_topograd, globe_FVEG
!>
!! - \b carbon \b variables: xglobe_CO2g_a, fglobe_CO2d_as,
!! fglobe_CO2g_sa, fglobe_CO2g_av, fglobe_CO2g_va, fglobe_CO2g_vs,
!! fglobe_CORG_vs, fglobe_CO2d_sr, fglobe_CO2d_sv
!>
!! - \b water \b variables: xglobe_H2Ol_v, xglobe_H2Og_a, fglobe_H2Ol_sv,
!! fglobe_H2Ol_sr, fglobe_H2Ol_sb, fglobe_H2Ol_br, fglobe_H2Og_sa,
!! fglobe_H2Og_sa_pot, fglobe_H2Og_va

!     ******************************************************************

      subroutine globe_read_restart ()
      use globe_mod
      implicit none

      integer :: ioerr

      __diag(kglobe_fdiag,'globe_read_restart')

!     * read GLOBE restart file

      if (kglobe_mypid == kglobe_kroot) then
        open(kglobe_frestart, FILE=sglobe_frestart, STATUS='old',      &
     &                        FORM='unformatted', IOSTAT=ioerr)
        call error_open_old(ioerr, sglobe_frestart)

!       * read scalars

        read(kglobe_frestart, IOSTAT=ioerr) nglobe_nlon
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) nglobe_nlat
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) nglobe_landpts
        __error_read(ioerr, sglobe_frestart)

        read(kglobe_frestart, IOSTAT=ioerr) kglobe_year
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) kglobe_timestep
        __error_read(ioerr, sglobe_frestart)

        if (kglobe_restart .eq. 2) kglobe_year = 1
      endif

!     * broadcast parameters

      __globe_mpbc(nglobe_nlon)
      __globe_mpbc(nglobe_nlat)
      __globe_mpbc(nglobe_landpts)
      __globe_mpbc(kglobe_year)
      __globe_mpbc(kglobe_timestep)

!     * allocate the fields

      call globe_landsea_alloc

!     * read the global fields

      if (kglobe_mypid == kglobe_kroot) then
        read(kglobe_frestart, IOSTAT=ioerr) kglobe_gLPIDs
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_longitude
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_latitude
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_lon
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_lat
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_elevation
        __error_read(ioerr, sglobe_frestart)
     read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_paw
     __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_len
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_wid
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_area
        __error_read(ioerr, sglobe_frestart)
        read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_topograd
        __error_read(ioerr, sglobe_frestart)

     read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_numnn
     __error_read(ioerr, sglobe_frestart)
     read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_idnn
     __error_read(ioerr, sglobe_frestart)
     read(kglobe_frestart, IOSTAT=ioerr) pglobe_gLP_distnn
     __error_read(ioerr, sglobe_frestart)
     read(kglobe_frestart, IOSTAT=ioerr) kglobe_lpt_nn
     __error_read(ioerr, sglobe_frestart)
      endif

      call globe_npro_fields
      call globe_alloc

!     * read the distributed fields (except climate forcing data -> jam)

!     * surface state

      call globe_mpreadgp(kglobe_frestart,globe_FVEG,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

!     * water balance

      call globe_mpreadgp(kglobe_frestart,xglobe_H2Ol_v     ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,xglobe_H2Og_a     ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpreadgp(kglobe_frestart,fglobe_H2Ol_sv    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Ol_sr    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Ol_sb    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Ol_br    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Og_sa    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Og_sa_pot,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_H2Og_va    ,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

!     * carbon balance

      call globe_mpreadgp(kglobe_frestart,xglobe_CO2g_a, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpreadgp(kglobe_frestart,fglobe_CO2d_as,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_CO2g_sa,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpreadgp(kglobe_frestart,fglobe_CO2g_av,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_CO2g_va,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_CO2g_vs,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_CORG_vs,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpreadgp(kglobe_frestart,fglobe_CO2d_sr,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpreadgp(kglobe_frestart,fglobe_CO2d_sv,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      if (kglobe_mypid == kglobe_kroot) close(kglobe_frestart,STATUS='DELETE')

!     * end
      __diag(kglobe_fdiag,'globe_read_restart: done')

      return
      end subroutine globe_read_restart

!     ******************************************************************
!     GLOBE_WRITE_RESTART
!     ******************************************************************

!> \brief Writes the restart data of GLOBE to a file

!> \detail This routine calls the subroutines to write restart data to
!! a file. The restart data is needed to continue a model run after
!! restart.

!     ******************************************************************

      subroutine globe_write_restart ()
      use globe_mod
      implicit none

      integer :: ioerr

      __diag(kglobe_fdiag,'globe_write_restart')

!     * write GLOBE restart file

      if (kglobe_mypid == kglobe_kroot) then
        open(kglobe_frestart, FILE=sglobe_frestart, STATUS='replace',  &
     &                        FORM='unformatted', IOSTAT=ioerr)
        call error_open_new(ioerr,sglobe_frestart)

!       * write scalars

        write(kglobe_frestart) nglobe_nlon
        write(kglobe_frestart) nglobe_nlat
        write(kglobe_frestart) nglobe_landpts

        write(kglobe_frestart) kglobe_year
        write(kglobe_frestart) kglobe_timestep

!       * write global fields

        write(kglobe_frestart) kglobe_gLPIDs
        write(kglobe_frestart) pglobe_longitude
        write(kglobe_frestart) pglobe_latitude
        write(kglobe_frestart) pglobe_gLP_lon
        write(kglobe_frestart) pglobe_gLP_lat
        write(kglobe_frestart) pglobe_gLP_elevation
        write(kglobe_frestart) pglobe_gLP_paw
        write(kglobe_frestart) pglobe_gLP_len
        write(kglobe_frestart) pglobe_gLP_wid
        write(kglobe_frestart) pglobe_gLP_area
        write(kglobe_frestart) pglobe_gLP_topograd

       write(kglobe_frestart) pglobe_gLP_numnn
       write(kglobe_frestart) pglobe_gLP_idnn
       write(kglobe_frestart) pglobe_gLP_distnn
       write(kglobe_frestart) kglobe_lpt_nn
      endif

!     * write distributed fields (except climate forcing data -> jam)

!     * surface state

      call globe_mpwritegp(kglobe_frestart,globe_FVEG,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

!     * water balance

      call globe_mpwritegp(kglobe_frestart,xglobe_H2Ol_v,      1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,xglobe_H2Og_a,      1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpwritegp(kglobe_frestart, fglobe_H2Ol_sv,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Ol_sr,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Ol_sb,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Ol_br,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Og_sa,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Og_sa_pot,1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart, fglobe_H2Og_va,    1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

!     * carbon balance

      call globe_mpwritegp(kglobe_frestart,xglobe_CO2g_a,  1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpwritegp(kglobe_frestart,fglobe_CO2d_as, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,fglobe_CO2g_sa, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpwritegp(kglobe_frestart,fglobe_CO2g_av, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,fglobe_CO2g_va, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,fglobe_CO2g_vs, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,fglobe_CORG_vs, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      call globe_mpwritegp(kglobe_frestart,fglobe_CO2d_sr, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)
      call globe_mpwritegp(kglobe_frestart,fglobe_CO2d_sv, 1,nglobe_cgpt,nglobe_landpts,nglobe_gpts2)

      if (kglobe_mypid == kglobe_kroot) close(kglobe_frestart)

!     * end
      __diag(kglobe_fdiag,'globe_write_restart: done')

      return
      end subroutine globe_write_restart

!     ******************************************************************
!     GLOBE_STATUS
!     ******************************************************************

!> \brief Writes the time step to a file

!> \detail This routine writes at a monthly basis the time step to a
!! status file.

!     ******************************************************************

      subroutine globe_status ()
      use globe_mod
      implicit none

      integer :: ioerr

      if (kglobe_mypid == kglobe_kroot .and. kglobe_status == 1) then
        open(kglobe_fstatus, FILE=sglobe_fstatus, STATUS='replace',    &
     &                       FORM='formatted', IOSTAT=ioerr)
        call error_open_new(ioerr,sglobe_fstatus)
        write(kglobe_fstatus,'(A,I7,A,I7,A,I7,A,I7,A)')                &
     &    'globe_status: day in year = ', kglobe_diy,                  &
     &    ' month = ', kglobe_month,                                   &
     &    ' year = ', kglobe_year + kglobe_firstyear - 1,              &
     &    ' of ', kglobe_lastyear, ' years'
        close(kglobe_fstatus)
      endif

      return
      end subroutine globe_status
