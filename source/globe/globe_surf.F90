#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_surf.F90
!> \brief Routines used to create the surface state file

!> \file globe_surf.F90
!> This file includes all routines, which are only used to create the
!! surface state file, which includes all constant fields describing the
!! state of the global model, originating from the global climate
!! forcing data.\n
!! Most important is the index mask, which includes the field index
!! of the land points, because only land points are processed.
!! Other important fields are the longitude, latitude, length, width,
!! and area of the fields.
!> Fields only used by the JEDI module are index, count, and
!! distance of nearest land neighbors.


!     ******************************************************************
!     GLOBE_CREATE_SURFACE_FILE
!     ******************************************************************

!> \brief Main routine to create the surface file

!> \detail This routine calls the different subroutines, which create and
!! write the fields to the surface file.

!     ******************************************************************

      subroutine globe_create_surface_file ()
      use globe_mod
      implicit none

      integer          :: ioerr

      logical          :: zlon_exist
      logical          :: zlat_exist
      real,allocatable :: zindexmask(:)
      real,allocatable :: zelevation(:)
      real,allocatable :: zlatitude(:)
      real,allocatable :: zlength(:)
      real,allocatable :: zwidth(:)
      real,allocatable :: zlonres(:)
      real,allocatable :: zlatres(:)
      real,allocatable :: zidnn(:,:)
      real,allocatable :: zdistnn(:,:)

      __diag(kglobe_fdiag,'globe_create_surface_file')

      if (kglobe_mypid == kglobe_kroot) then

!       * open surface file

        open(kglobe_fsurface, FILE=sglobe_fsurface, STATUS='replace',  &
     &                        FORM='unformatted', IOSTAT=ioerr)
        call error_open_new(ioerr, sglobe_fsurface)

!       * create and insert index mask to surface file

        call globe_insert_indexmask

!       * allocate local fields

        allocate(zindexmask(nglobe_nlon*nglobe_nlat))
        allocate(zelevation(nglobe_nlon*nglobe_nlat))
        allocate(zlatitude(nglobe_nlon*nglobe_nlat))
        allocate(zlength(nglobe_nlon*nglobe_nlat))
        allocate(zwidth(nglobe_nlon*nglobe_nlat))
        allocate(zlonres(nglobe_nlon*nglobe_nlat))
        allocate(zlatres(nglobe_nlon*nglobe_nlat))
        allocate(zidnn(nglobe_nlon*nglobe_nlat,8))
        allocate(zdistnn(nglobe_nlon*nglobe_nlat,8))

!       * globe_tmp_lonlat is read by routine globe_insert_indexmask
        zindexmask(:) = globe_tmp_lonlat(:)

!       * add the elevation to the surface file

        call globe_insert_map(kglobe_felevation, sglobe_felevation, kcode_elevation)
!       * globe_tmp_lonlat is read by routine globe_insert_map
!         (if no elevation is given, globe_tmp_lonlat is filled with zeros)
        zelevation(:) = globe_tmp_lonlat(:)
!       * zero out the elevation for unused map points
        where (zindexmask(:) .eq. 0.0) zelevation(:) = 0.0

   !    * add plant available water to the surface file

        call globe_insert_map(kglobe_fpaw, sglobe_fpaw, kcode_paw)

!       * add longitude, latitude, area, and extension of fields to surface file
!         (output parameters: zlatitude, zlength, zwidth)

        call globe_landsea_area(zindexmask, zlatitude, zlength, zwidth, &
     &                          zlon_exist, zlat_exist, zlonres, zlatres)

!       * calculate immediate neighbors indexes and count
   !      and add it to the surface file
!         (putput parameter: zidnn)

        call globe_landsea_idnn(zindexmask, zlon_exist, zidnn)

!       * calculate distance to immediate neighbors
   !      and add it to the surface file
!         (output parameter: zdistnn)

        call globe_landsea_distance(zlatitude, zlength, zwidth,        &
     &                              zlonres, zlatres, zdistnn)

!       * calculate and add topographic gradient to surface file

        call globe_topographic_gradient(zelevation, zidnn, zdistnn,       &
     &                                  zlon_exist, zlat_exist)

!       * clean up

        close(kglobe_fsurface)

!       * global variable globe_tmp_lonlat was only temporary allocated
        deallocate(globe_tmp_lonlat)

        deallocate(zindexmask)
        deallocate(zelevation)
        deallocate(zlatitude)
        deallocate(zlength)
        deallocate(zwidth)
        deallocate(zlonres)
        deallocate(zlatres)
        deallocate(zidnn)
        deallocate(zdistnn)

      endif  ! (kglobe_mypid == kglobe_kroot)

      __diag(kglobe_fdiag,'globe_create_surface_file: done')

      return
      end subroutine globe_create_surface_file

!     ******************************************************************
!     GLOBE_INSERT_INDEXMASK
!     ******************************************************************

!> \brief Creates and writes the index mask to the surface file

!> \detail This routine creates the land index mask from a given
!! land-sea-mask-file, and writes it to the surface file. The
!! land-sea-mask is a map where land points are marked with one and sea
!! points with zero.\n
!! The first part reads the land-sea-mask file and gets the map size.\n
!! The second part changes the land-sea-mask to contain only 0 and 1
!! by splitting the values at 0.5.\n
!! The third part fills the mask with a continuous index, starting with
!! the index of one from the upper left corner of the map.\n
!! The last part writes the index mask to the surface file.

!     ******************************************************************

      subroutine globe_insert_indexmask ()
      use globe_mod
      implicit none

      integer :: ihead(8), ioerr, i

      if (kglobe_mypid == kglobe_kroot) then

!       * read the land-sea mask to globe_tmp_lonlat

        call globe_open_input(sglobe_flandsea, kglobe_flandsea, kglobe_fdiag)
        __globe_read_srvnc(kglobe_flandsea,ihead,iostat=ioerr)

!       * get the map size
        nglobe_nlon = ihead(5)
        nglobe_nlat = ihead(6)
!       * temporary allocate globe_tmp_lonlat, real allocation is done later
        allocate(globe_tmp_lonlat(nglobe_nlon*nglobe_nlat))

        __globe_read_srvnc(kglobe_flandsea,ihead,globe_tmp_lonlat)

        call globe_close_input(kglobe_flandsea, kglobe_fdiag)

!       * set the land-sea mask to 0 and 1 only

        where (globe_tmp_lonlat(:) .lt. 0.5)
          globe_tmp_lonlat(:) = 0.0
        elsewhere
          globe_tmp_lonlat(:) = 1.0
        end where

!       * get the nunber of land points

        nglobe_landpts = INT(SUM(globe_tmp_lonlat(:)) + 0.5)

!       * fill the land-sea mask with the field indexes

        do i = 1, nglobe_nlon * nglobe_nlat
          if (globe_tmp_lonlat(i) .ne. 0.0) globe_tmp_lonlat(i) = real(i)
        enddo

!       * adapt and write the field header to the surface file

        ihead(1) = kcode_index
        ihead(2) = 0
        ihead(3) = 0
        ihead(4) = 0
        ihead(8) = nglobe_landpts
        write(kglobe_fsurface) ihead

!       * write the data to the surface file

        write(kglobe_fsurface) globe_tmp_lonlat

!       * do not deallocate globe_tmp_lonlat, because it is also used by
!         the caller routine

      endif

      return
      end subroutine globe_insert_indexmask

!     ******************************************************************
!     GLOBE_INSERT_MAP
!     ******************************************************************

!> \brief Inserts data from a given file to the surface file

!> \detail This routine adds the variable field from the given file to
!! the surface file. If the given file does not exist, nothing is
!! written.\n
!! The first part reads and checks the header from the given file.\n
!! The last part adopts the header for the surface file, inserts the
!! given \c code number, and writes the data to the surface file.

!> \param kfile File ID of the input file
!> \param sfile File name of the input file
!> \param code  Code number of the data field in the surface file

!     ******************************************************************

      subroutine globe_insert_map (kfile, sfile, code)
      use globe_mod
      implicit none

      integer,          INTENT(in) :: kfile, code
      character(len=*), INTENT(in) :: sfile

      integer :: ihead(8), ioerr
      logical :: fexist

      if (kglobe_mypid == kglobe_kroot) then

        inquire(FILE=sfile, EXIST=fexist)
        if (fexist) then

!         * open the input file
          call globe_open_input(sfile, kfile, kglobe_fdiag)

!         * read and check the header
          __globe_read_srvnc(kfile,ihead,iostat=ioerr)
          if (nglobe_nlon .ne. ihead(5) .or. nglobe_nlat .ne. ihead(6)) then
            write(0,*) 'ERROR: globe_insert_map: Array-size of land-sea-mask and ', sfile, ' differ'
            stop
          endif

!         * read the data
          __globe_read_srvnc(kfile,ihead,globe_tmp_lonlat)

!         * adapt and write the header to the surface file
          ihead(1) = code
          ihead(2) = 0
          ihead(3) = 0
          ihead(4) = 0
          ihead(8) = nglobe_landpts
          write(kglobe_fsurface) ihead

!         * write the data to the surface file
          write(kglobe_fsurface) globe_tmp_lonlat
          call globe_close_input(kfile, kglobe_fdiag)

        else
          globe_tmp_lonlat(:) = 0.0
        endif

      endif

      return
      end subroutine globe_insert_map

!     ******************************************************************
!     GLOBE_LANDSEA_AREA
!     ******************************************************************

!> \brief Calculates and writes field positions and sizes to the surface file

!> \detail This routine calculates for every grid cell of the map, defined
!! by the index mask, its central longitude, and latitude, plus its
!! \c area [km<sup>2</sup>] and extension [km] as length and width, and
!! writes the data to the surface file.

!> \param zindexmask Index mask calculated by globe_insert_indexmask
!> \param zlatitude  Field with latitudes
!> \param zlength    Field with the extension in longitude of grid cells
!> \param zwidth     Field with the extension in latitude of grid cells
!> \param zlon_exist Flag for longitudes read from file
!> \param zlat_exist Flag for latitudes read from file
!> \param zlonres    Field with the resolution of the grid cells in longitudes
!> \param zlatres    Field with the resolution of the grid cells in latitudes

!     ******************************************************************

      subroutine globe_landsea_area (zindexmask, zlatitude, zlength, zwidth, &
     &                               zlon_exist, zlat_exist, zlonres, zlatres)
      use globe_mod
      implicit none

      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(in)  :: zindexmask
      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(out) :: zlatitude
      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(out) :: zlength
      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(out) :: zwidth
      logical,                                 INTENT(out) :: zlon_exist
      logical,                                 INTENT(out) :: zlat_exist
      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(out) :: zlonres
      real,dimension(nglobe_nlon*nglobe_nlat), INTENT(out) :: zlatres

      real,parameter   :: radius = 6378.137  ! GRS-80 spherical approximation
      double precision :: pi, deg2rad
      integer          :: ihead(8), i, j, k, x1, x2, x3, x4, ioerr

      real,allocatable :: elevation(:)
      real,allocatable :: longitude(:)
      real,allocatable :: area(:)
      real,allocatable :: deltalon(:)
      real,allocatable :: deltalat(:)

      if (kglobe_mypid == kglobe_kroot) then

!       * initialization

        allocate(elevation(nglobe_nlon*nglobe_nlat))
        allocate(longitude(nglobe_nlon*nglobe_nlat))
        allocate(area(nglobe_nlon*nglobe_nlat))
        allocate(deltalon(nglobe_nlon*nglobe_nlat))
        allocate(deltalat(nglobe_nlon*nglobe_nlat))

!       * get the land-sea-mask from the index mask
        globe_tmp_lonlat(:) = MIN(1.0, zindexmask(:))

        elevation(:) = 0.0
        longitude(:) = 0.0
        zlatitude(:) = 0.0

        pi = 4.0 * ATAN(1.0)
        deg2rad = pi / 180.0

!       * check if file sglobe_longitude exist and read it

        inquire(FILE=sglobe_longitude, EXIST=zlon_exist)
        if (zlon_exist) then
!         * open the input file
          call globe_open_input(sglobe_longitude, kglobe_longitude, kglobe_fdiag)
!         * read and check the header
          __globe_read_srvnc(kglobe_longitude,ihead,iostat=ioerr)
          if (ihead(5) .ne. nglobe_nlon .or. (ihead(6) .ne. 1 .and. ihead(6) .ne. nglobe_nlat)) then
            write(0,*) 'ERROR: globe_landsea_area: Number of longitudes of land-sea-mask and ', sglobe_longitude, ' differ'
            stop
          endif
!         * read the data
          if (ihead(6) .ne. 1) then
            __globe_read_srvnc(kglobe_longitude,ihead,longitude(:))
          else
            __globe_read_srvnc(kglobe_longitude,ihead,longitude(1:nglobe_nlon))
!           * fill the whole longitude field
            do i = 1, nglobe_nlat-1
              x1 = nglobe_nlon * i + 1
              x2 = nglobe_nlon * i + nglobe_nlon
              longitude(x1:x2) = longitude(1:nglobe_nlon)
            enddo
          endif
          call globe_close_input(kglobe_longitude, kglobe_fdiag)
        endif

!       * check if file sglobe_latitude exist and read it

        inquire(FILE=sglobe_latitude, EXIST=zlat_exist)
        if (zlat_exist) then
!         * open the input file
          call globe_open_input(sglobe_latitude, kglobe_latitude, kglobe_fdiag)
!         * read and check the header
          __globe_read_srvnc(kglobe_latitude,ihead,iostat=ioerr)
          if ((ihead(5) .ne. nglobe_nlat .and. ihead(5) .ne. nglobe_nlon) .or. (ihead(6) .ne. 1 .and. ihead(6) .ne. nglobe_nlat)) then
            write(0,*) 'ERROR: globe_landsea_area: Number of latitudes of land-sea-mask and ', sglobe_latitude, ' differ'
            stop
          endif
!         * read the data
          if (ihead(6) .ne. 1) then
            __globe_read_srvnc(kglobe_latitude,ihead,zlatitude(:))
          else
            __globe_read_srvnc(kglobe_latitude,ihead,zlatitude(1:nglobe_nlat))
!           * fill the whole latitude field
            do i = nglobe_nlat, 1, -1
              zlatitude((i-1)*nglobe_nlon+1:i*nglobe_nlon) = zlatitude(i)
            enddo
          endif
          call globe_close_input(kglobe_latitude, kglobe_fdiag)
!         * check orientation of latitude (north has to be up)
          if (zlatitude(1) .lt. zlatitude(nglobe_nlon*nglobe_nlat)) then
            write(*,*) "ERROR: ",TRIM(sglobe_latitude), ": Wrong map orientation: north has to be on top"
            stop
          endif
        endif

!       * calculate the longitude and latitude for the fields

        if (.not. zlon_exist .or. .not. zlat_exist) then
          zlonres(:) = 360.0 / nglobe_nlon
          zlatres(:) = 180.0 / nglobe_nlat
          k = 0
          do j = 1 , nglobe_nlat
            do i = 1 , nglobe_nlon
              k = k + 1
              if (.not. zlon_exist) then
                longitude(k) = (i-1) * zlonres(1) + zlonres(1) / 2.0
                if (longitude(k) .gt. 180.0) longitude(k) = longitude(k) - 360.0
              endif
              if (.not. zlat_exist) then
                zlatitude(k)  = 90.0 - (j-1) * zlatres(1) - zlatres(1) / 2.0
              endif
            enddo
          enddo
        endif

!       * calculate the grid resolution in longitude and latitude

        if (zlon_exist) then
          if (nglobe_nlon == 1) then
            zlonres(:) = pglobe_one_deltalon
          else
            do i = 1, nglobe_nlat
              x1 = (i-1) * nglobe_nlon + 1
              x2 = (i-1) * nglobe_nlon + nglobe_nlon
              zlonres(x1:x2-1) = longitude(x1+1:x2) - longitude(x1:x2-1)
              zlonres(x2) = zlonres(x2-1)
              if (nglobe_nlon .gt. 2) then
                zlonres(x1+1:x2-1) = (zlonres(x1:x2-2) + zlonres(x1+1:x2-1)) / 2.0
              endif
            enddo
            where (zlonres(:) == 0.0) zlonres(:) = pglobe_one_deltalon
          endif
        endif

        if (zlat_exist) then
          if (nglobe_nlat == 1) then
            zlatres(:) = pglobe_one_deltalat
          else
            do i = 1, nglobe_nlat-1
              x1 = (i-1) * nglobe_nlon + 1
              x2 = (i-1) * nglobe_nlon + nglobe_nlon
              x3 = i * nglobe_nlon + 1
              x4 = i * nglobe_nlon + nglobe_nlon
              zlatres(x1:x2) = zlatitude(x1:x2) - zlatitude(x3:x4)
            enddo
            zlatres(x3:x4) = zlatres(x1:x2)
            if (nglobe_nlat .gt. 2) then
              do i = 2, nglobe_nlat-1
                x1 = (i-2) * nglobe_nlon + 1
                x2 = (i-2) * nglobe_nlon + nglobe_nlon
                x3 = (i-1) * nglobe_nlon + 1
                x4 = (i-1) * nglobe_nlon + nglobe_nlon
                zlatres(x3:x4) = (zlatres(x1:x2) + zlatres(x3:x4)) / 2.0
              enddo
            endif
            where (zlatres(:) == 0.0) zlatres(:) = pglobe_one_deltalat
          endif
        endif

!       * calculate the area [km^2] and extension [km] of the fields

        elevation(:) = deg2rad * zlatitude(:)
        deltalat(:)  = deg2rad * zlatres(:)
        deltalon(:)  = deg2rad * zlonres(:)

        zlength(:)   = (deg2rad * COS(elevation(:)) * radius)          &
     &               * zlonres(:) * globe_tmp_lonlat(:)
        zwidth(:)    = deg2rad * radius * zlatres(:) * globe_tmp_lonlat(:)

        area(:)      = 2.0 * radius**2 * deltalon(:) * COS(elevation(:)) &
     &               * SIN(deltalat(:) / 2.0) * globe_tmp_lonlat(:)

!       * write the fields to the surface file

        ihead(:) = 0
        ihead(5) = nglobe_nlon
        ihead(6) = nglobe_nlat
        ihead(8) = nglobe_landpts

        ihead(1) = kcode_longitude
        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) longitude

        ihead(1) = kcode_latitude
        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) zlatitude

        ihead(1) = kcode_length
        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) zlength

        ihead(1) = kcode_width
        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) zwidth

        ihead(1) = kcode_area
        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) area

        deallocate(elevation)
        deallocate(longitude)
        deallocate(area)
        deallocate(deltalon)
        deallocate(deltalat)

        zlatitude(:) = zlatitude(:) * globe_tmp_lonlat(:)

      endif

      return
      end subroutine globe_landsea_area

!     ******************************************************************
!     GLOBE_LANDSEA_IDNN
!     ******************************************************************

!> \brief Calculates the index and count of the nearest land neighbors

!> \detail This routine calculates a layered field of all 8 nearest
!! neighbors of a grid cell from the index mask
   !>, and writes it to the surface file
!>. Since sea points are zero in the index mask, this indexes are also
!! zero in the layered field.
   !> The count of nearest neighbors not equal
   !! zero is also written to the surface file.

!> \param zindexmask Index mask calculated by globe_insert_indexmask
!> \param zlon_exist Flag for longitudes read from file
!> \param zidnn      Index field with all nearest land neighbours

!     ******************************************************************

      subroutine globe_landsea_idnn (zindexmask, zlon_exist, zidnn)
      use globe_mod
      implicit none

      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zindexmask
      logical,                                   INTENT(in)  :: zlon_exist
      real,dimension(nglobe_nlon*nglobe_nlat,8), INTENT(out) :: zidnn

      integer :: ihead(8), i, nlon, nlat

      integer, allocatable :: idx(:,:)
      real,    allocatable :: znumnn(:)

      if (kglobe_mypid == kglobe_kroot) then

        allocate(idx(nglobe_nlat,2))
        allocate(znumnn(nglobe_nlon*nglobe_nlat))

        nlon = nglobe_nlon
        nlat = nglobe_nlat

!       * shift all immediate neighbours over the center -> neighbour-index

        zidnn(:,:) = 0.0
        zidnn(1+nlon+1:nlon*nlat       ,1) = zindexmask(1       :nlon*nlat-nlon-1)
        zidnn(1+nlon  :nlon*nlat       ,2) = zindexmask(1       :nlon*nlat-nlon  )
        zidnn(1+nlon-1:nlon*nlat       ,3) = zindexmask(1       :nlon*nlat-nlon+1)
        zidnn(1       :nlon*nlat     -1,4) = zindexmask(1     +1:nlon*nlat       )
        zidnn(1       :nlon*nlat-nlon-1,5) = zindexmask(1+nlon+1:nlon*nlat       )
        zidnn(1       :nlon*nlat-nlon  ,6) = zindexmask(1+nlon  :nlon*nlat       )
        zidnn(1       :nlon*nlat-nlon+1,7) = zindexmask(1+nlon-1:nlon*nlat       )
        zidnn(1     +1:nlon*nlat       ,8) = zindexmask(1       :nlon*nlat     -1)

!       * create an index-vector for the left and right edge

        do i = 1, nlat
          idx(i,1) = 1    + (i-1) * nlon
          idx(i,2) = nlon + (i-1) * nlon
        enddo

!       * correct the left and right edge of the world map

        if (zlon_exist) then
          zidnn(idx(1  :nlat-1,1),1) = 0.0
          zidnn(idx(1+1:nlat  ,2),3) = 0.0
          zidnn(idx(1+1:nlat  ,2),4) = 0.0
          zidnn(idx(1+1:nlat  ,2),5) = 0.0
          zidnn(idx(1  :nlat-1,1),7) = 0.0
          zidnn(idx(1  :nlat-1,1),8) = 0.0
        else
          zidnn(idx(1  :nlat-1,1),1) = zidnn(idx(1+1:nlat  ,1),1)
          zidnn(idx(1+1:nlat  ,2),3) = zidnn(idx(1  :nlat-1,2),3)
          zidnn(idx(1+1:nlat  ,2),4) = zidnn(idx(1  :nlat-1,2),4)
          zidnn(idx(1+1:nlat  ,2),5) = zidnn(idx(1  :nlat-1,2),5)
          zidnn(idx(1  :nlat-1,1),7) = zidnn(idx(1+1:nlat  ,1),7)
          zidnn(idx(1  :nlat-1,1),8) = zidnn(idx(1+1:nlat  ,1),8)
        endif

!       * write the field in levels to the surface file

        ihead(:) = 0
        ihead(1) = kcode_idnn
        ihead(5) = nglobe_nlon
        ihead(6) = nglobe_nlat
        ihead(8) = nglobe_landpts

        do i = 1, 8
          ihead(2) = i
          write(kglobe_fsurface) ihead
          write(kglobe_fsurface) zidnn(:,i)
        enddo

!       * count immediate land neighbors -> neighbor-count

        znumnn(:) = 0.0
        where (zidnn(:,1) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,2) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,3) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,4) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,5) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,6) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,7) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0
        where (zidnn(:,8) .ne. 0.0) znumnn(:) = znumnn(:) + 1.0

!       * write the field to the surface file

        ihead(:) = 0.0
        ihead(1) = kcode_numnn
        ihead(5) = nglobe_nlon
        ihead(6) = nglobe_nlat
        ihead(8) = nglobe_landpts

        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) znumnn

        deallocate(idx)
        deallocate(znumnn)

      endif

      return
      end subroutine globe_landsea_idnn

!     ******************************************************************
!     GLOBE_LANDSEA_DISTANCE
!     ******************************************************************

!> \brief Calculates the distance to immediate neighbors

!> \detail This routine calculates a layered map containing the distance
!! to all nearest neighbors
!> and writes it to the surface file
!> .

!> \param zlatitude Field with the latitude calculated by globe_landsea_area
!> \param zlength   Field with the length calculated by globe_landsea_area
!> \param zwidth    Field with the width calculated by globe_landsea_area
!> \param zlonres    Field with the resolution of the grid cells in longitudes
!> \param zlatres    Field with the resolution of the grid cells in latitudes
!> \param zdistnn   Field with the distance to all nearest neighbours

!     ******************************************************************

      subroutine globe_landsea_distance (zlatitude, zlength, zwidth, zlonres, zlatres, zdistnn)
      use globe_mod
      implicit none

      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zlatitude
      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zlength
      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zwidth
      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zlonres
      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in)  :: zlatres
      real,dimension(nglobe_nlon*nglobe_nlat,8), INTENT(out) :: zdistnn

      integer :: ihead(8), i
      real,allocatable :: distance(:)

      real,             allocatable :: T1(:)
      double precision, allocatable :: T2_8(:)
      double precision, allocatable :: T3_8(:)

      double precision, parameter :: PI_8 = 3.14159265358979d0

      if (kglobe_mypid == kglobe_kroot) then

        allocate(distance(nglobe_nlon*nglobe_nlat))
        allocate(T1(nglobe_nlon*nglobe_nlat))
        allocate(T2_8(nglobe_nlon*nglobe_nlat))
        allocate(T3_8(nglobe_nlon*nglobe_nlat))

        T1(:)   = 111.0 * zlatres(:)
        T2_8(:) = COS((zlatitude(:) + zlatres(:) / 2.0) / 180.0 * PI_8)
        T3_8(:) = zlonres(:) / zlatres(:)

        distance(:) = T1 / COS(ATAN(T2_8(:) * T3_8))

        zdistnn(:,1) = distance(:)
        zdistnn(:,2) = zwidth(:)
        zdistnn(:,3) = distance(:)
        zdistnn(:,4) = zlength(:)
        zdistnn(:,5) = distance(:)
        zdistnn(:,6) = zwidth(:)
        zdistnn(:,7) = distance(:)
        zdistnn(:,8) = zlength(:)

!       * write the field in levels to the surface file

        ihead(:) = 0
        ihead(1) = kcode_distnn
        ihead(5) = nglobe_nlon
        ihead(6) = nglobe_nlat
        ihead(8) = nglobe_landpts

        do i = 1, 8
          ihead(2) = i
          write(kglobe_fsurface) ihead
          write(kglobe_fsurface) zdistnn(:,i)
        enddo

        deallocate(distance)
        deallocate(T1)
        deallocate(T2_8)
        deallocate(T3_8)

      endif

      return
      end subroutine globe_landsea_distance

!     ******************************************************************
!     GLOBE_TOPOGRAPHIC_GRADIENT
!     ******************************************************************

!> \brief Adds the mean topographic gradient to immediate neighbors to
!! the surface file

!> \detail This routine calculates the mean topographic gradient from
!! the difference of the elevation to the nearest neighbors, and writes it
!! to the surface file.

!> \param elevation  elevation map as read from file \c sglobe_felevation
!> \param zidnn   Index map with all nearest land neighbors
!> \param zdistnn Field with the distance to all nearest neighbors
!> \param zlon_exist Flag for longitudes read from file
!> \param zlat_exist Flag for latitudes read from file

!     ******************************************************************

      subroutine globe_topographic_gradient (elevation, zidnn, zdistnn, zlon_exist, zlat_exist)
      use globe_mod
      implicit none

      real,dimension(nglobe_nlon*nglobe_nlat),   INTENT(in) :: elevation
      real,dimension(nglobe_nlon*nglobe_nlat,8), INTENT(in) :: zidnn
      real,dimension(nglobe_nlon*nglobe_nlat,8), INTENT(in) :: zdistnn
      logical,                                   INTENT(in) :: zlon_exist
      logical,                                   INTENT(in) :: zlat_exist

      integer :: ihead(8), i, j, n, x1, x2, x3, x4
      real    :: z_topograd, zelevation

      real,allocatable :: ztopograd(:)

      if (kglobe_mypid == kglobe_kroot) then

        allocate(ztopograd(nglobe_nlon*nglobe_nlat))

        ztopograd(:) = 0.0

        do j = 1, nglobe_nlon * nglobe_nlat
          if (elevation(j) .ne. 0.0) then
            n = 0
            do i = 1, 8
              zelevation = 0.0
              if (zidnn(j,i) > 0) zelevation = elevation(INT(zidnn(j,i) + 0.5))
              if (zelevation > elevation(j)) then
                z_topograd = 0.0
              else
                z_topograd = ( elevation(j) - zelevation ) / zdistnn(j,i)
                n          = n + 1
              endif
              ztopograd(j) = ztopograd(j) + z_topograd
            enddo
            if (n > 0) then
              ztopograd(j) = ztopograd(j) / n
            else
              ztopograd(j) = 0.0
            endif
          endif
        enddo

!       * set ztopograd for land points at edges of site maps to the value of the nearest cell

        if (zlon_exist) then
          do i = 1, nglobe_nlat
            x1 = (i-1) * nglobe_nlon + 1
            x2 = (i-1) * nglobe_nlon + nglobe_nlon
            ztopograd(x1) = ztopograd(x1+1)
            ztopograd(x2) = ztopograd(x2-1)
          enddo
        endif
        if (zlat_exist) then
          x1 = 1
          x2 = nglobe_nlon
          x3 = nglobe_nlon + 1
          x4 = nglobe_nlon + nglobe_nlon
          ztopograd(x1:x2) = ztopograd(x3:x4)
          x1 = (nglobe_nlat-1) * nglobe_nlon + 1
          x2 = (nglobe_nlat-1) * nglobe_nlon + nglobe_nlon
          x3 = (nglobe_nlat-2) * nglobe_nlon + 1
          x4 = (nglobe_nlat-2) * nglobe_nlon + nglobe_nlon
          ztopograd(x1:x2) = ztopograd(x3:x4)
        endif
        if (zlon_exist .or. zlat_exist) then
          where (elevation(:) .eq. 0.0) ztopograd(:) = 0.0
        endif

!       * write the field to the surface file

        ihead(:) = 0
        ihead(1) = kcode_topography
        ihead(5) = nglobe_nlon
        ihead(6) = nglobe_nlat
        ihead(8) = nglobe_landpts

        write(kglobe_fsurface) ihead
        write(kglobe_fsurface) ztopograd

        deallocate(ztopograd)

      endif

      return
      end subroutine globe_topographic_gradient
