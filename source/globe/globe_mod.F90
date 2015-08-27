#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file globe_mod.F90
!> \brief GLOBE variable definitions

!> \file globe_mod.F90
!> This file contains the \c module \c globe_mod which defines all
!! global variables used by the GLOBE routines.

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      module globe_mod

      use globe_functions
      use globe_flux_mod
      use globe_stat_mod

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                         Grid size
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer             :: nglobe_nlon     = 64 !< number of longitudes
      integer             :: nglobe_nlat     = 32 !< number of latitudes
      integer             :: nglobe_landpts  = 0  !< number of land points
      integer,allocatable :: kglobe_gLPIDs(:)     !< land point index field

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                 Surface parameter fields
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real, allocatable   :: pglobe_latitude(:)          !< latitudes of the grid
      real, allocatable   :: pglobe_longitude(:)         !< longitudes of the grid
      real, allocatable   :: pglobe_gLP_lat(:)           !< latitude of the land points
      real, allocatable   :: pglobe_gLP_lon(:)           !< longitude of the land points
      real, allocatable   :: pglobe_gLP_elevation(:)     !< elevation
      real                :: pglobe_gLP_elevation0 = 0.0 !< elevation: default value
      real, allocatable   :: pglobe_gLP_paw(:)           !< plant available water (mm water per mm soil depth)
      real                :: pglobe_gLP_paw0 = 0.10      !< plant available water: default value
      real, allocatable   :: pglobe_gLP_topograd(:)      !< mean topographic gradient
      real, allocatable   :: pglobe_gLP_len(:)           !< length of grid cell
      real, allocatable   :: pglobe_gLP_wid(:)           !< width of grid cell
      real, allocatable   :: pglobe_gLP_area(:)          !< area of grid cell

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()         Field size as multiple of number of CPUs
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer,allocatable :: kglobe_gLPIDs2(:)    !< land point index field
      real, allocatable   :: pglobe_gLP_lat2(:)   !< latitude of the land points
      real, allocatable   :: pglobe_gLP_lon2(:)   !< longitude of the land points
      real, allocatable   :: pglobe_gLP_hght2(:)  !< elevation
      real, allocatable   :: pglobe_gLP_paw2(:)   !< plant available water
      real, allocatable   :: pglobe_gLP_togr2(:)  !< mean topographic gradient
      real, allocatable   :: pglobe_gLP_len2(:)   !< length of grid cell
      real, allocatable   :: pglobe_gLP_wid2(:)   !< width of grid cell
      real, allocatable   :: pglobe_gLP_area2(:)  !< area of grid cell

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()            Nearest neighbors and their indices
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer,allocatable :: pglobe_gLP_numnn(:)    !< number of nearest neighbors
      integer,allocatable :: pglobe_gLP_idnn(:,:)   !< map-index of nearest neighbors
      real,   allocatable :: pglobe_gLP_distnn(:,:) !< distance to nearest neighbors
      integer,allocatable :: kglobe_lpt_nn(:,:)     !< land-index of nearest neighbors

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                      Parallelization
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer, parameter  :: kglobe_kroot   = 0   !< master node number
      integer             :: kglobe_myworld = 0   !< MPI identification
      integer             :: kglobe_mypid   = 0   !< process ID
      integer             :: nglobe_npro    = 1   !< number of processes
      integer             :: nglobe_cgpt    = 0   !< land points per CPU
      integer             :: nglobe_gpts2   = 0   !< nglobe_cgpt * nglobe_npro

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     GLOBE run control
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer             :: kglobe_status = 1    !< switch (0/1) to create status file
      integer             :: kglobe_diag = 1      !< switch (0/1) to create diagnostic output
      integer             :: kglobe_restart = 0   !< switch for restart: 0 = false, 1 = normal, 2 = reset kglobe_year (works only in special cases), -1 = no creation of restart files
      integer             :: kglobe_firstyear = 1 !< first year
      integer             :: kglobe_lastyear = 10 !< last year
      integer             :: kglobe_yrsskip = 0   !< number of years to skip output (repeating)
      integer             :: kglobe_yrsout = 1    !< number of years to output
      integer             :: kglobe_mapoutput = 1 !< switch for srv-output: 1 = as map, 0 = only land points
      integer             :: kglobe_netcdf = 0    !< switch to output netCDF format: 1 = version3, 2 = version4, 3 = version4 compressed

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                   Single point version
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer             :: kglobe_spt = 0       !< switch (0/1) for for single point version
      integer             :: kglobe_spt_gpid = 0  !< spt: grid point id for single point
      real                :: pglobe_spt_lon = 0   !< spt: give longitude to determine grid point id
      real                :: pglobe_spt_lat = 0   !< spt: give latitude to determine grid point id
      real                :: pglobe_one_deltalon = 0 !< extension of the grid points in longitude [degree] (overwrites calculated values if .ne. 0, but required for fields with only one longitude)
      real                :: pglobe_one_deltalat = 0 !< extension of the grid points in latitude [degree] (overwrites calculated values if .ne. 0, but required for fields with only one latitude)

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                      Loop variables
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer             :: kglobe_year = 1      !< current year
      integer             :: kglobe_dpy = 0       !< number of days in the year
      integer             :: kglobe_diy           !< day within year

      integer             :: kglobe_month         !< current month of year
      integer             :: kglobe_dpm(14) = 0   !< days per month
      integer             :: kglobe_dim           !< max. number of dpm = field size

      integer             :: kglobe_day           !< current day of month
      integer             :: kglobe_leapday = 0   !< flag of leap day in month

      integer             :: kglobe_ts            !< current time step in day
      integer             :: kglobe_tspd = 1      !< number of time steps per day
      integer             :: kglobe_timestep = 0  !< overall time step count

      real                :: kglobe_spts = 86400. !< seconds per time step (initialized as seconds per day)

      integer             :: kglobe_file_pointer_8(9999) = -1 !< to save the position after last file read for SRV data files

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()           Predefined indicators for step handling
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer, parameter  :: GLOBE_STEP_PREYEAR  = 1 !< execute before the first month of year
      integer, parameter  :: GLOBE_STEP_PREMONTH = 2 !< execute before the first day of month
      integer, parameter  :: GLOBE_STEP_PREDAY   = 3 !< execute before the first time step per day
      integer, parameter  :: GLOBE_STEP_TIMESTEP = 4 !< execute at every time step
      integer, parameter  :: GLOBE_STEP_DAY      = 5 !< execute after the last time step of day
      integer, parameter  :: GLOBE_STEP_MONTH    = 6 !< execute after the last day of month
      integer, parameter  :: GLOBE_STEP_YEAR     = 7 !< execute after the last month of year

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     GLOBE file units
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer             :: kglobe_fnamelist    = 90  !< namelist file
      integer             :: kglobe_frestart     = 91  !< restart file
      integer             :: kglobe_fsurface     = 92  !< input of surface parameters
      integer             :: kglobe_fgridpoint   = 921 !< grid points file
      integer             :: kglobe_flandsea     = 922 !< land sea mask
      integer             :: kglobe_felevation   = 923 !< elevation for land sea mask
      integer             :: kglobe_fpaw         = 925 !< plant available water
      integer             :: kglobe_longitude    = 926 !< longitude values in case of site or irregular grid
      integer             :: kglobe_latitude     = 927 !< latitude values in case of site or irregular grid
      integer             :: kglobe_flandmask    = 93  !< masking
      integer             :: kglobe_fdiag        = 94  !< diagnostic output
      integer             :: kglobe_fstatus      = 96  !< GLOBE status file
      integer             :: kglobe_foutfldpos   = 981 !< data field positions in the output files

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     GLOBE file names
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      character(len=*),parameter :: sglobe_fnamelist   = "globe_namelist"   !< namelist file
      character(len=*),parameter :: sglobe_frestart    = "globe_restart"    !< restart file
      character(len=*),parameter :: sglobe_fsurface    = "globe_surf.srv"   !< input of surface parameters
      character(len=*),parameter :: sglobe_fgridpoint  = "gridpoints.txt"   !< grid points file
      character(len=*),parameter :: sglobe_flandsea    = "landsea.srv"      !< land sea mask
      character(len=*),parameter :: sglobe_felevation  = "elevation.srv"    !< elevation for land sea mask
      character(len=*),parameter :: sglobe_fpaw        = "paw.srv"          !< plant available water file
      character(len=*),parameter :: sglobe_longitude   = "longitude.srv"    !< longitude values in case of site or irregular grid
      character(len=*),parameter :: sglobe_latitude    = "latitude.srv"     !< latitude values in case of site or irregular grid
      character(len=*),parameter :: sglobe_flandmask   = "landmask.srv"     !< land mask
      character(len=*),parameter :: sglobe_fdiag       = "globe_diag.txt"   !< diagnostic output
      character(len=*),parameter :: sglobe_fstatus     = "globe_status.txt" !< GLOBE status file
      character(len=*),parameter :: sglobe_foutfldpos  = "globe_fldpos.srv" !< data field positions in the output files

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     GLOBE field codes
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer, parameter  :: kcode_ls = 172       !< land-sea mask

!     * field codes of the surface file

      integer, parameter  :: kcode_index      = 901 !< index at landsea mask
      integer, parameter  :: kcode_elevation  = 902 !< elevation map
      integer, parameter  :: kcode_paw        = 904 !< plant available water
      integer, parameter  :: kcode_longitude  = 905 !< longitude of field
      integer, parameter  :: kcode_latitude   = 906 !< latitude of field
      integer, parameter  :: kcode_length     = 907 !< extension in longitude of field [ km ]
      integer, parameter  :: kcode_width      = 908 !< extension in latitude of field [ km ]
      integer, parameter  :: kcode_area       = 909 !< area size of field [ km<sup>2</sup> ]
      integer, parameter  :: kcode_idnn       = 910 !< index of all 8 nearest neighbors (8 levels)
      integer, parameter  :: kcode_distnn     = 911 !< distance to nearest neighbours (8 levels)
      integer, parameter  :: kcode_topography = 912 !< mean topographic gradient of field
      integer, parameter  :: kcode_numnn      = 913 !< number of nearest neighbors with land

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()           Temporary variables for sub-procedures
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()
!     * (prevents heavy allocation and deallocation - may save time)

      real(kind=4), allocatable :: globe_tmp_lonlat_4(:) !< size: nglobe_nlon * nglobe_nlat
      real(kind=4), allocatable :: globe_tmp_nLPs_4(:)   !< size: nglobe_landpts
      real,         allocatable :: globe_tmp_lonlat(:)   !< size: nglobe_nlon * nglobe_nlat
      real,         allocatable :: globe_tmp_nLPs(:)     !< size: nglobe_landpts
      real,         allocatable :: globe_tmp_nGPts2(:)   !< size: nglobe_gpts2

#ifndef __PLASIM
__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                   Variables for PLASIM
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

!     * indicator for the first run of landmod intern in PLASIM

      integer           :: globe_plasim_init_done = 0

!     * variables of PLASIM

      real              :: globe_plasim_ALS
      real              :: globe_plasim_ALV
      integer           :: globe_plasim_NHOR
      integer           :: globe_plasim_NLEP
      integer           :: globe_plasim_NLEV
      integer           :: globe_plasim_NLPP
      integer           :: globe_plasim_NLSOIL
      real              :: globe_plasim_TMELT
      real              :: globe_plasim_co2
      real, allocatable :: globe_plasim_dalb(:)
      real, allocatable :: globe_plasim_dcsoil(:)
      real, allocatable :: globe_plasim_dcveg(:)
      real              :: globe_plasim_delt
      real, allocatable :: globe_plasim_devap(:)
      real, allocatable :: globe_plasim_dpet(:)
      real, allocatable :: globe_plasim_dflux(:,:)
      real, allocatable :: globe_plasim_dforest(:)
      real, allocatable :: globe_plasim_dglac(:)
      real, allocatable :: globe_plasim_dgpp(:)
      real, allocatable :: globe_plasim_dgppl(:)
      real, allocatable :: globe_plasim_dgppw(:)
      real, allocatable :: globe_plasim_dlai(:)
      real, allocatable :: globe_plasim_dlhfl(:)
      real, allocatable :: globe_plasim_dlitter(:)
      real, allocatable :: globe_plasim_dls(:)
      real, allocatable :: globe_plasim_dnogrow(:)
      real, allocatable :: globe_plasim_dnpp(:)
      real, allocatable :: globe_plasim_dp(:)
      real, allocatable :: globe_plasim_dprc(:)
      real, allocatable :: globe_plasim_dprl(:)
      real, allocatable :: globe_plasim_dprs(:)
      real, allocatable :: globe_plasim_dq(:,:)
      real, allocatable :: globe_plasim_dres(:)
      real, allocatable :: globe_plasim_drhs(:)
      real, allocatable :: globe_plasim_drunoff(:)
      real, allocatable :: globe_plasim_dshfl(:)
      real, allocatable :: globe_plasim_dsmelt(:)
      real, allocatable :: globe_plasim_dsndch(:)
      real, allocatable :: globe_plasim_dsnow(:)
      real, allocatable :: globe_plasim_dswfl(:,:)
      real, allocatable :: globe_plasim_dt(:,:)
      real, allocatable :: globe_plasim_dtd2(:)
      real, allocatable :: globe_plasim_dtd3(:)
      real, allocatable :: globe_plasim_dtd4(:)
      real, allocatable :: globe_plasim_dtd5(:)
      real, allocatable :: globe_plasim_dtsoil(:)
      real, allocatable :: globe_plasim_dveg(:)
      real, allocatable :: globe_plasim_dwatc(:)
      real, allocatable :: globe_plasim_dwmax(:)
      real, allocatable :: globe_plasim_dz0(:)
      integer           :: globe_plasim_nperpetual
      integer           :: globe_plasim_nprhor
      integer           :: globe_plasim_nprint
      real              :: globe_plasim_psurf
      real              :: globe_plasim_ra1
      real              :: globe_plasim_ra2
      real              :: globe_plasim_ra4
      real              :: globe_plasim_rdbrv
      real              :: globe_plasim_ww

!     * variables of landmod

      integer           :: globe_plasim_nlandt
      integer           :: globe_plasim_nlandw
      integer           :: globe_plasim_nbiome
      integer           :: globe_plasim_nmaxevap
      integer           :: globe_plasim_newsurf
      real              :: globe_plasim_albland
      real              :: globe_plasim_albsmin
      real              :: globe_plasim_albsmax
      real              :: globe_plasim_albsminf
      real              :: globe_plasim_albsmaxf
      real              :: globe_plasim_albgmin
      real              :: globe_plasim_albgmax
      real              :: globe_plasim_dz0land
      real              :: globe_plasim_drhsland
      real              :: globe_plasim_drhsfull
      real              :: globe_plasim_dzglac
      real              :: globe_plasim_dztop
      real              :: globe_plasim_dsmax
      real              :: globe_plasim_wsmax
      integer           :: globe_plasim_ncveg
      real              :: globe_plasim_rlue
      real              :: globe_plasim_co2conv
      real              :: globe_plasim_tau_veg
      real              :: globe_plasim_tau_soil
      real              :: globe_plasim_riniveg
      real              :: globe_plasim_rinisoil
      real              :: globe_plasim_rnbiocats
      real              :: globe_plasim_forgrow
      real              :: globe_plasim_rlaigrow
      real              :: globe_plasim_gs
      real              :: globe_plasim_z0_max
      real, allocatable :: globe_plasim_doro(:)
      real, allocatable :: globe_plasim_dts(:)
      real, allocatable :: globe_plasim_dtsm(:)
      real, allocatable :: globe_plasim_dqs(:)
      real, allocatable :: globe_plasim_driver(:)
      real, allocatable :: globe_plasim_duroff(:)
      real, allocatable :: globe_plasim_dvroff(:)
      real, allocatable :: globe_plasim_darea(:)
      real, allocatable :: globe_plasim_dsoilz(:)
      real, allocatable :: globe_plasim_dsoilt(:,:)
      real, allocatable :: globe_plasim_dsnowt(:)
      real, allocatable :: globe_plasim_dtclsoil(:)
      real, allocatable :: globe_plasim_dsnowz(:)
      real, allocatable :: globe_plasim_dwater(:)
      real, allocatable :: globe_plasim_pgrow(:)
      real, allocatable :: globe_plasim_plai(:)
      real, allocatable :: globe_plasim_pgs(:)
      real, allocatable :: globe_plasim_pz0_max(:)
      real, allocatable :: globe_plasim_dtcl(:,:)
      real, allocatable :: globe_plasim_dwcl(:,:)
      real, allocatable :: globe_plasim_dtclim(:)
      real, allocatable :: globe_plasim_dwclim(:)
      real, allocatable :: globe_plasim_dz0clim(:)
      real, allocatable :: globe_plasim_dz0climo(:)
      real, allocatable :: globe_plasim_dalbclim(:)
      real, allocatable :: globe_plasim_dh2ol_added(:)
      real, allocatable :: globe_plasim_geomask(:)
#endif /* __PLASIM */

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      end module globe_mod
