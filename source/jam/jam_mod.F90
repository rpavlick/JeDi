#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jam_mod.F90
!> \brief Variable definitions of JAM

!> \file jam_mod.F90
!> This file contains the \c module \c jam_mod, which defines all the
!! global variables used by the JAM routines only.

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()

      module jam_mod
      use globe_functions

!_vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     Constants for JAM
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real, parameter   :: cTMelt  = 273.15       !< melting temperature
      real, parameter   :: cRg     = 8.3145       !< gas constant (J/mol/K)
      real, parameter   :: cRH2Og  = 461.5        !< gas constant for water vapor
      real, parameter   :: cRhoH2O = 1000.0       !< density of water (kg/m<sup>3</sup>)

!     ==================================================================
!     ------------------------------------------------------------------
!     * parameters provided by GLOBE
!     ------------------------------------------------------------------
!     ==================================================================

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                   GLOBE parallel stuff
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer           :: kroot  = 0             !< master node ID
      integer           :: kmypid = 0             !< process ID

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                    GLOBE grid handling
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer, allocatable :: kLPID(:)            !< land point index field
      integer              :: nLPts  = 0          !< number of actual land points
      integer              :: nCPts  = 0          !< number of land points per CPU
      integer              :: nLPts2 = 0          !< npro * nCPts
      integer              :: nlon                !< resolution longitude
      integer              :: nlat                !< resolution latitude

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                    GLOBE date handling
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer           :: ktsim                  !< time step in the actual month
      integer           :: kts                    !< actual time step
      integer           :: kday                   !< actual day
      integer           :: kmonth                 !< actual month
      integer           :: kyear                  !< actual year

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                     GLOBE run control
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer           :: krestart = 0           !< same as kglobe_restart
      integer           :: kdiag    = 1           !< same as kglobe_diag
      integer           :: ksingle  = 0           !< same as kglobe_spt

!     ==================================================================

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                 GLOBE exchange variables
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real              :: pdt = 86400.0          !< time step in seconds (initialized as seconds per day)
      integer           :: ktspd                  !< time steps per day
      integer           :: kdpm(14) = 0           !< days of month from the climate files
      integer           :: kdpy                   !< days per year
      integer           :: nleapday               !< has the month a leap day

      integer           :: mdim                   !< max. number of days per month

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                          Energy
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real, allocatable :: pTEMP(:,:)             !< surface temperature
      real, allocatable :: fRADs_as(:,:)          !< short wave radiation
      real, allocatable :: fRADl_as(:,:)          !< long wave radiation

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                           Water
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real, allocatable :: fH2Ol_as(:,:)          !< precipitation
      real, allocatable :: cH2Og_a(:,:)           !< H<sub>2</sub>O in boundary layer (rel. hum.)

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                          Carbon
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real, allocatable :: xCO2g_a(:,:)           !< CO<sub>2</sub> in atmosphere (ppm)

!     ==================================================================
!     ------------------------------------------------------------------
!     * JAM parameters
!     ------------------------------------------------------------------
!     ==================================================================

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                  JAM namelist parameters
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real              :: pDeltaT   = 0.0        !< add temperature anomaly
      real              :: pDeltaP   = 0.0        !< multiply precipitation anomaly

      integer           :: kCO2_dyn  = 0          !< switch (0/2) to calculate dynamic CO<sub>2</sub> (2 = read yearly data from text file)
      integer           :: kRandYear = 0          !< switch (0/1) to use random year of climate forcing

      real              :: pTEMP_init     = 300.0 !< initializes surface temperature
      real              :: fRADs_as_init  = 380.0 !< initializes short wave radiation
      real              :: fRADl_as_init  = 40.0  !< initializes long wave radiation
      real              :: fH2Ol_as_init  = 20.0e-07 !< initializes precipitation
      real              :: cH2Og_a_init   = 0.7   !< initializes H<sub>2</sub>O in boundary layer (rel. hum.)
      real              :: xCO2g_a_init   = 360.0 !< initializes CO<sub>2</sub> in atmosphere (ppm)

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                   JAM status parameters
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer           :: kmonmeans = 0          !< data is monthly means ? 0/1
      integer           :: n_data_years = 1       !< number of years in the data
      integer           :: k_first_dyear          !< first year in data set

      integer           :: nspointer = 1          !< data set pointer for pre-read data (single point)
      integer           :: nmaxpoints = 0         !< number of data sets in climate files
      integer           :: ndatafields = 0        !< number of climate variables to read
      integer           :: nfieldcount = 1        !< count the assigned fields (single point)
      integer           :: kfieldIDs(12) = 0      !< reference of fields to data file IDs (single point)

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()            Atmospheric forcing data (monthly)
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()
!     The field size of kfieldIDs must be the number of those data fields

      real,allocatable  :: pTEMP2(:,:)            !< surface temperature
      real,allocatable  :: fRADs_as2(:,:)         !< short wave radiation
      real,allocatable  :: fRADl_as2(:,:)         !< long wave radiation
      real,allocatable  :: fH2Ol_as2(:,:)         !< precipitation
      real,allocatable  :: cH2Og_a2(:,:)          !< H<sub>2</sub>O in boundary layer (rel. hum.)
      real,allocatable  :: xCO2g_a2(:,:)          !< CO<sub>2</sub> in atmosphere (ppm)

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()      Atmospheric forcing data (pre-read single point)
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      real,allocatable  :: preclim(:,:)           !< variable to store whole atmospheric forcing

!     ==================================================================
!     ------------------------------------------------------------------
!     * file handling
!     ------------------------------------------------------------------
!     ==================================================================

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                Switch if input file exist
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      logical           :: input_pTEMP
      logical           :: input_fRADs_as
      logical           :: input_fRADl_as
      logical           :: input_fH2Ol_as
      logical           :: input_cH2Og_a
      logical           :: input_xCO2g_a

      logical           :: input_pco2list

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                        File units
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      integer :: kfile_namelist = 800             !< namelist file
      integer :: kfile_diag     = 810             !< diagnostic output
      integer :: kfile_temp2m   = 841             !< temperature
      integer :: kfile_solrad   = 842             !< solar radiation
      integer :: kfile_terrad   = 843             !< terrestrial radiation
      integer :: kfile_precip   = 844             !< precipitation
      integer :: kfile_relhum   = 845             !< relative humidity
      integer :: kfile_co2a     = 855             !< CO<sub>2</sub> in atmosphere
      integer :: kfile_pco2list = 860             !< CO<sub>2</sub> concentration

__vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv()
__SEC()                        File names
__oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo()

      character(len=*),parameter :: sfile_namelist = 'jam_namelist' !< namelist file
      character(len=*),parameter :: sfile_diag     = 'jam_diag.txt' !< diagnostic output
      character(len=*),parameter :: sfile_temp2m   = 'tas.srv'     !< temperature
      character(len=*),parameter :: sfile_solrad   = 'rsds.srv'   !< downward solar radiation
      character(len=*),parameter :: sfile_terrad   = 'rlns.srv'   !< net longwave radiation
      character(len=*),parameter :: sfile_precip   = 'pr.srv'     !< precipitation
      character(len=*),parameter :: sfile_relhum   = 'hurs.srv'   !< relative humidity
      character(len=*),parameter :: sfile_co2a     = 'co2a.srv'     !< CO<sub>2</sub> in atmosphere
      character(len=*),parameter :: sfile_pco2list = 'pco2.txt'     !< CO<sub>2</sub> concentration

      end module jam_mod
