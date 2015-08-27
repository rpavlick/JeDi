#define __ACTIVATE
#include "../globe/globe_macros.f90"

      module jedi_mod
      use globe_functions

!     ------------------------------------------------------------------
!     * global parameters
!     ------------------------------------------------------------------

      parameter ( PI2     = 6.283185307179586 )
      parameter ( cGamma  = 65.0              )
      parameter ( cLambda = 2.5E6             )
      parameter ( cTiny   = 1.E-6             )
      real(8) :: r8_tiny = 0.0123456
      real(8) :: r8_epsilon = 0.0123456

!     ------------------------------------------------------------------
!     * parameters provided by GLOBE
!     ------------------------------------------------------------------

!     * parallel stuff

      integer :: kroot    = 0                   ! Master node
      integer :: kmypid   = 0                   ! Process ID

      integer :: gNHOR    = 2048                ! global field number
      integer :: NHOR     = 2048                ! local field number
      integer :: kNumGPts = 2048                ! number of land points

!     * grid handling

      integer :: nlat     = 32                  ! Number of latitudes
      integer :: nlon     = 64                  ! Number of longitudes

!     * flags

      integer :: krestart = 0                   ! 1 for true, 0 for false
      integer :: kdiag    = 1
      integer :: kspinup  = 1                   ! 1 for true, 0 for false
      integer :: ktopspecies = 0
      integer :: kdynamic = 0
      integer :: kspec_yrout = 0
      integer :: koutput = 1

!     * date handling

      integer             :: YEAR
      integer             :: OUTYEAR
      integer             :: MONTH
      integer             :: kdiy
      integer             :: kdpy

!     * climate variables and units:
!       - temperature:     in centigrade
!       - precipitation:   in mm/d
!       - shortwave down:  in W/m2
!       - longwave net:    in W/m2
      real,allocatable    :: pjedi_rain   (:)   ! rainfall
      real,allocatable    :: pjedi_snow   (:)   ! snowfall
      real,allocatable    :: pjedi_temp   (:)   ! temperature
      real,allocatable    :: pjedi_swdown (:)   ! shortwave down
      real,allocatable    :: pjedi_lwnet  (:)   ! longwave net
      real,allocatable    :: cCO2g_a      (:)   ! atmospheric CO2 concentration

      integer,allocatable :: gkGPID       (:)   ! land point index field
      real,allocatable    :: pjedi_lat    (:)   ! latitude for daylen

!     * area handling
      real,allocatable    :: pjedi_area   (:)   ! total area of landpoints
      integer,allocatable :: klpt_nn(:,:)       ! land-index of nearest neighbors

      real,allocatable    :: pjedi_paw    (:)   ! plant available water (mm water per mm soil)

!     * Seed for random numbers

      integer             :: seed = 123456      ! base of random seed
      integer             :: rand_seed          ! seed to generate random numbers

!     ------------------------------------------------------------------
!     * file units
!     ------------------------------------------------------------------

      integer :: kFile_Grid       = 5000        ! grid-based fields
      integer :: kFile_SPP        = 5100        ! species-specific fields
      integer :: kFile_SPPParm    = 5200        ! species parameter files
      integer :: kFile_Namelist   = 5300        ! namelist file
      integer :: kFile_Restart    = 5400        ! restart file
      integer :: kFile_Diag       = 5500        ! diagnostic text output
      integer :: kFile_SPPSucc    = 5600        ! species-specific output
      integer :: kFile_SPPGrid    = 5610        ! species-specific grid output
      integer :: kFile_OptiSPP    = 5700        ! optimum species parameters
      integer :: kFile_Opti       = 5800        ! optimization log
      integer :: kFile_OptiOut    = 5900        ! optimization output
      integer :: kFile_SPPExtinct = 5650        ! extinct species files
      integer :: kFile_GASPP      = 5210        ! community-aggregated species parameter file
      integer :: kFile_SPPTopGrid = 5620        ! species-specific grid input

!     ------------------------------------------------------------------
!     * file names
!     ------------------------------------------------------------------

      character(len=*),parameter :: sfile_Grid       = 'jedi_output'
      character(len=*),parameter :: sfile_SPP        = 'jedi_species.txt'
      character(len=*),parameter :: sfile_SPPParm    = 'jedi_specparms.txt'
      character(len=*),parameter :: sfile_SPPGrid    = 'jedi_species'
      character(len=*),parameter :: sfile_Namelist   = 'jedi_namelist'
      character(len=*),parameter :: sfile_Restart    = 'jedi_restart'
      character(len=*),parameter :: sfile_Diag       = 'jedi_diag.txt'
      character(len=*),parameter :: sfile_SPPSucc    = 'jedi_success.txt'
      character(len=*),parameter :: sfile_OptiSPP    = 'jedi_opti_spp.txt'
      character(len=*),parameter :: sfile_Opti       = 'jedi_opti_log.txt'
      character(len=*),parameter :: sfile_OptiOut    = 'jedi_opti'
      character(len=*),parameter :: sfile_SPPExtinct = 'jedi_exspecies.txt'
      character(len=*),parameter :: sfile_GASPP      = 'jedi_gaspp.txt'
      character(len=*),parameter :: sfile_SPPTopGrid = 'jedi_topspecies.srv'

!     ------------------------------------------------------------------
!     * field codes
!     ------------------------------------------------------------------

      parameter (kcode_CORG_v  = 5000)           ! vegetation biomass
      parameter (kcode_CO2g_av = 5100)           ! gross primary production

      parameter (kcode_goalfct = 5900)           ! optimization goal function
      parameter (kcode_optispp = 5901)           ! optimum species

      parameter (kcode_5300    = 5300)
      parameter (kcode_5301    = 5301)
      parameter (kcode_5302    = 5302)
      parameter (kcode_5304    = 5304)
      parameter (kcode_5370    = 5370)
      parameter (kcode_5310    = 5310)
      parameter (kcode_5311    = 5311)
      parameter (kcode_5308    = 5308)
      parameter (kcode_5182    = 5182)
      parameter (kcode_5174    = 5174)
      parameter (kcode_5200    = 5200)
      parameter (kcode_5199    = 5199)

      parameter (kcode_5229    = 5229)
      parameter (kcode_5378    = 5378)
      parameter (kcode_5379    = 5379)
      parameter (kcode_5600    = 5600)
      parameter (kcode_5611    = 5611)
      parameter (kcode_5612    = 5612)
      parameter (kcode_5613    = 5613)
      parameter (kcode_5614    = 5614)
      parameter (kcode_5621    = 5621)
      parameter (kcode_5622    = 5622)

      parameter (kcode_5623    = 5623)
      parameter (kcode_5624    = 5624)
      parameter (kcode_5625    = 5625)

      parameter (kcode_5410    = 5410)
      parameter (kcode_5411    = 5411)
      parameter (kcode_5412    = 5412)
      parameter (kcode_5413    = 5413)
      parameter (kcode_5414    = 5414)
      parameter (kcode_5415    = 5415)

      parameter (kcode_5416    = 5416)
      parameter (kcode_5417    = 5417)
      parameter (kcode_5402    = 5402)
      parameter (kcode_5403    = 5403)
      parameter (kcode_5404    = 5404)


!     ------------------------------------------------------------------
!     * namelist parameters
!     ------------------------------------------------------------------

!     * JEDI version and options

      integer :: kPop_dyn    = 1       ! 0 for individuals, 1 for population dynamics
      integer :: kSpec_dyn   = 0       ! 1 for Speciation on, 0 for off
      integer :: kMig_dyn    = 0       ! 1 for Migration on, 0 off
      integer :: kDist_dyn   = 0
      integer :: kOpti       = 0       ! 1 for maximization of productivity
      integer :: nonepoint   = 0       ! 1 run only 1 point, 0 normal output
      integer :: kFrPl       = 1       ! first plant for simulation
      integer :: kToPl       = 200     ! last plant for simulation
      integer :: kMaxSPP     = 200     ! maximal number of modeled species (important for MES)
      integer :: BareSPP     = 1       ! empty SPP with no plant for bare soil
      integer :: FirstSPP    = 2       ! first plant SPP

!     * ecophysiological parameters

      real    :: pq10R       = 2.0     ! q10 value for resp
      real    :: pQ10H       = 2.0     ! value of Q10 for heterotrophic respiration
      real    :: pC_GPP      = 1.2E-1  ! light use efficiency c_GPP in gC/W d
      real    :: pC_RES      = 2.0E-3  ! specific respiration rate in gC/d gC
      real    :: pC_LAI      = 0.030   ! specific leaf area in m^2 gC^-1
      real    :: pC_TRS      = 0.5     ! specific supply rate in mm d^-1 gC^-1
      real    :: pC_WSMAX    = 40.0    ! specific root depth in (mm gC^-1)^-1/2
      real    :: pC_WLMAX    = 0.2     ! specific canopy capacity in mm m^2
      real    :: pC_FFOR     = 0.002   ! forest ratio conversion in gC
      real    :: pPAW0       = 50.0    ! minimum soil moisture capacity in mm
      real    :: pAlb0       = 0.2     ! unvegetated albedo
      real    :: PAlb1       = 0.12    ! vegetated albedo
      real    :: pT_CRIT     = 10.0    ! critical temperature in K
      real    :: pMortTau    = 150.0   ! mortality par (years)
      real    :: pDistTau    = 3.0     ! disturbances par (years)
      real    :: pRAbdTau    = 1.0     ! relative abundance timescale (years)
      real    :: pT_GPP1     = 10.0    ! critical temperature
      real    :: pT_GPP2     = 8.0     ! critical temperature
      real    :: cPETCorr    = 1.0

!     * parameters related to litter and soil pools

      real   :: pLitterTau      = 2.0   ! litter residence time (years)
      real   :: pCWDTau         = 25.0  ! woody debris residence time (years)
      real   :: pSoilTau        = 100.0 ! soil carbon residence time (years)
      real   :: pLitterFrac2Atm = 0.77  ! fraction of litter decay going to atmosphere
      real   :: pCWDFrac2Atm    = 0.20  ! fraction of woody decay going to atmosphere

      real    :: pN_RES      = 0.218   ! specific respiration in gC/gN/d
      real    :: pCN_Wood    = 0.00303 ! 1/C:N Ratio for Woody Tissue
      real    :: pCN_Root    = 0.02    ! 1/C:N Ratio for Fine Roots
      real    :: pCN_Leaf    = 0.08    ! 1/C:N Ratio for Leaf Roots
      real    :: cT_MIN      = -5.0

      real    :: pA0         = 1.0     ! initial seed size A0
      real    :: pSeedTau    = 10.0    ! seed decomposition time scale (years)
      real    :: pSeedFix    = 0.01    ! seed fix constant (unitless)

      integer :: kGermFix    = 1       ! flag for germination fix
      integer :: kCbal       = 1       ! flag for germination fix
      integer :: kCWT        = 0       ! flag for community-weighted traits

!     * parameters for speciation and migration and Disturbances

      real    :: pMig        = 0.01    ! fraction of seeds that disperse to neigbors
      real    :: pSpec       = 0.1     ! variance for variation of species

!     * CO2 fertilization factor defined as GPP ~ (1 + beta ln(C/C0))
!     beta = 0.3-0.6, C0 = 360ppm - Harvey 1989

      real     :: pCO2sens   = 0.0
      real     :: pPCO2      = 364.0
      real, parameter :: co2_ref = 364.0 ! CO2 reference concentration in ppm
      real, parameter :: pjedi_defaultpaw = 0.25 ! default soil texture if paw.srv not present

      real, parameter :: cTmin    = -2.5
      real, parameter :: cTref    = 20.0
      real, parameter :: cTmax350 = 28.0
      real, parameter :: cTmax700 = 33.0
      real     :: Tmax

!     * initialization

      integer :: initGP      = 0       ! from which GP to start initially
      real    :: pdaylen     = 0.0     ! sets fixed day length
      real    :: pfGerm      = 1.0     ! set competition for bare space (seeds)
      real    :: pfComp      = 1.0     ! set competition for occ space (dominance)

!     ------------------------------------------------------------------
!     * tags for handling migration, speciation
!     ------------------------------------------------------------------

      integer :: jedi_migruns  = 0
      integer :: jedi_specruns = 0

!     ------------------------------------------------------------------
!     * constants
!     ------------------------------------------------------------------

      real, parameter :: c2PI         = 6.28318512
      real, parameter :: cK           = 0.5
      real, parameter :: cWSUB_MAX    = 1000.0
      real, parameter :: cDieDeltaC   = 1.0
      real, parameter :: cGResLeaf    = 1.0 / (1.0 + 0.33)
      real, parameter :: cGResWood    = 1.0 / (1.0 + 0.25)
      real, parameter :: cGResRoot    = 1.0 / (1.0 + 0.49)
      real, parameter :: cGResSeed    = 1.0 / (1.0 + 0.66)
      real, parameter :: cWRFrac      = 0.05
      real, parameter :: cDaysPerYear = 365.0

      real, parameter :: cT0          = 273.15
      real, parameter :: cPAR         = 0.50

!     ------------------------------------------------------------------
!     * plant species parameters
!     ------------------------------------------------------------------

      real, allocatable :: p01 (:)
      real, allocatable :: p02 (:)
      real, allocatable :: p05 (:)
      real, allocatable :: p06 (:)
      real, allocatable :: p07 (:)
      real, allocatable :: p08 (:)
      real, allocatable :: p09 (:)
      real, allocatable :: p10 (:)
      real, allocatable :: p13 (:)
      real, allocatable :: p16 (:)
      real, allocatable :: p14 (:)
      real, allocatable :: p15 (:)
      real, allocatable :: p04 (:)
      real, allocatable :: p17 (:)
      real, allocatable :: p12 (:)
      real, allocatable :: p11 (:)
      real, allocatable :: p03 (:)
      real, allocatable :: p18 (:)
      real, allocatable :: p19 (:)
      real, allocatable :: p20 (:)

!     ------------------------------------------------------------------
!     * plant species handling
!     ------------------------------------------------------------------

      integer, allocatable :: kSPPID(:) ! pointer to species
      integer :: kNumSPPID              ! last set SPPID Number
      integer :: kNumSPP                ! last set species

!     ------------------------------------------------------------------
!     * state variables
!     ------------------------------------------------------------------

!     * state variables -- pools

      real, allocatable :: rW        (:,:)
      real, allocatable :: rWS       (:,:)
      real, allocatable :: rWL       (:,:)
      real, allocatable :: rWSUB     (:,:)
      real, allocatable :: rCA       (:,:)
      real, allocatable :: rCL       (:,:)
      real, allocatable :: rCR       (:,:)
      real, allocatable :: rCWL      (:,:)
      real, allocatable :: rCWR      (:,:)
      real, allocatable :: rCS       (:,:)
      real, allocatable :: rCtot     (:,:)

      real, allocatable :: rLIT_CL   (:)
      real, allocatable :: rLIT_CR   (:)
      real, allocatable :: rLIT_CWL  (:)
      real, allocatable :: rLIT_CWR  (:)
      real, allocatable :: rCSLOW    (:)

      real, allocatable :: rArea     (:,:)
      real, allocatable :: rAreaBare (:)

!     * state variables -- life time

      real, allocatable :: rLIVE     (:,:)

!     * state variables -- memory

      real, allocatable :: rGrowW    (:,:)
      real, allocatable :: rGrowG    (:,:)
      real, allocatable :: rGrowT    (:,:)
      real, allocatable :: rDie      (:,:)

!     * variables for gather

      real, allocatable :: grArea     (:,:) ! Area
      real, allocatable :: grAreaBare (:)   ! Bare Area
      real, allocatable :: grCA       (:,:)
      real, allocatable :: grCS       (:,:)
      real, allocatable :: grCWL      (:,:)
      real, allocatable :: grCWR      (:,:)
      real, allocatable :: grCL       (:,:)
      real, allocatable :: grCR       (:,:)

      real, allocatable :: grLIT_CL   (:)
      real, allocatable :: grLIT_CR   (:)
      real, allocatable :: grLIT_CWL  (:)
      real, allocatable :: grLIT_CWR  (:)
      real, allocatable :: grCSLOW    (:)

      real, allocatable :: grLIVE     (:,:)

!     ------------------------------------------------------------------
!     * output variables
!     ------------------------------------------------------------------

      real, allocatable, target :: BM_total  (:)
      real, allocatable, target :: Richness  (:)
      real, allocatable         :: dRAbd     (:,:)
      real, allocatable         :: dGGermCnt (:)
      real, allocatable         :: dGDeadCnt (:)

      real, allocatable :: dGAGPP    (:)
      real, allocatable :: dGANPP    (:)
      real, allocatable :: dGARES    (:)
      real, allocatable :: dGARESH   (:)
      real, allocatable :: dGALIT    (:)

      real, allocatable :: dGASOLRAD (:)
      real, allocatable :: dGALH     (:)
      real, allocatable :: dGASH     (:)

      real, allocatable :: dGACALOSS (:)
      real, allocatable :: dGACLLOSS (:)
      real, allocatable :: dGACRLOSS (:)
      real, allocatable :: dGACWLLOSS(:)
      real, allocatable :: dGACWRLOSS(:)

      real, allocatable :: dGACLALLOC (:)
      real, allocatable :: dGACRALLOC (:)
      real, allocatable :: dGACWLALLOC(:)
      real, allocatable :: dGACWRALLOC(:)
      real, allocatable :: dGACSALLOC (:)

      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CS
      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CL
      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CR
      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CWL
      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CWR
      REAL,allocatable,dimension(:) :: dGALOSS_LIT_CSLOW

      real, allocatable :: dGAFH2O   (:)
      real, allocatable :: dGAFT     (:)

      real, allocatable :: dGAET     (:)
      real, allocatable :: dGALEVAP     (:)
      real, allocatable :: dGATRANS     (:)
      real, allocatable :: dGAESOIL     (:)
      real, allocatable :: dGAQTOT     (:)
      real, allocatable :: dGAQSURF     (:)

      real, allocatable :: dGAALB    (:)
      real, allocatable :: dGALAI    (:)
      real, allocatable :: dGAFVEG   (:)
      real, allocatable :: dGAFFOR   (:)
      real, allocatable :: dGAWMAX   (:)

      real, allocatable :: dGAtau    (:)
      real, allocatable :: dGAtauM   (:)
      real, allocatable :: dGAGerm   (:)
      real, allocatable :: dGACol    (:)
      real, allocatable :: dGAExcl   (:)
      real, allocatable :: dGAMort   (:)
      real, allocatable :: dGADist   (:)
      real, allocatable :: dGASpec   (:)
      real, allocatable :: dGAMig    (:)
      real, allocatable :: dGAExt    (:)

      real              :: nAccuCount = 0

!     * variables for gather

      real, allocatable,target :: gdGAGPP  (:)
      real, allocatable,target :: gdGANPP  (:)
      real, allocatable,target :: gdGANEE  (:)
      real, allocatable,target :: gdGARES  (:)
      real, allocatable,target :: gdGARESH (:)
      real, allocatable,target :: gdGALIT  (:)
      real, allocatable        :: gdGARESE (:)

      real, allocatable,target :: gdGACALOSS  (:)
      real, allocatable,target :: gdGACLLOSS  (:)
      real, allocatable,target :: gdGACRLOSS  (:)
      real, allocatable,target :: gdGACWLLOSS (:)
      real, allocatable,target :: gdGACWRLOSS (:)

      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CS
      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CL
      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CR
      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CWL
      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CWR
      REAL,allocatable,dimension(:) :: gdGALOSS_LIT_CSLOW

      real, allocatable        :: gdGASOLRAD  (:)
      real, allocatable        :: gdGALH      (:)
      real, allocatable        :: gdGASH      (:)

      real, allocatable,target :: gdGACLALLOC  (:)
      real, allocatable,target :: gdGACRALLOC  (:)
      real, allocatable,target :: gdGACWLALLOC (:)
      real, allocatable,target :: gdGACWRALLOC (:)
      real, allocatable,target :: gdGACSALLOC  (:)

      real, allocatable,target :: gdGAFT       (:)
      real, allocatable,target :: gdGAFH2O     (:)

      real, allocatable,target :: gdGAET       (:)
      real, allocatable,target :: gdGAESOIL       (:)
      real, allocatable,target :: gdGALEVAP       (:)
      real, allocatable,target :: gdGATRANS       (:)
      real, allocatable,target :: gdGAQSURF       (:)
      real, allocatable        :: gdGAQTOT       (:)

      real, allocatable :: gdGAALB     (:)
      real, allocatable :: gdGALAI     (:)
      real, allocatable :: gdGAFVEG    (:)
      real, allocatable :: gdGAFFOR    (:)
      real, allocatable :: gdGAWMAX    (:)

      real, allocatable,target :: gdGAGerm    (:)
      real, allocatable,target :: gdGACol     (:)
      real, allocatable,target :: gdGAExcl    (:)
      real, allocatable,target :: gdGAMort    (:)
      real, allocatable,target :: gdGADist    (:)
      real, allocatable,target :: gdGASpec    (:)
      real, allocatable        :: gdGAtau     (:)
      real, allocatable        :: gdGAtauM    (:)

      real, allocatable        :: gdGAMig     (:)
      real, allocatable,target :: gdGAExt     (:)

!     * species-specific fields

      real, allocatable :: dSARAbd     (:,:)
      real, allocatable :: dSAGPP      (:,:)
      real, allocatable :: dSANPP      (:,:)
      real, allocatable :: dSACS       (:,:)
      real, allocatable :: dSACA       (:,:)
      real, allocatable :: dSACL       (:,:)
      real, allocatable :: dSACR       (:,:)
      real, allocatable :: dSACWL      (:,:)
      real, allocatable :: dSACWR      (:,:)

      real, allocatable :: dSAtau      (:,:)
      real, allocatable :: dSAtauM     (:,:)

!     * variables for gather

      real, allocatable :: gdSARAbd    (:,:)
      real, allocatable :: gdSAGPP     (:,:)
      real, allocatable :: gdSANPP     (:,:)
      real, allocatable :: gdSACS      (:,:)
      real, allocatable :: gdSACA      (:,:)
      real, allocatable :: gdSACL      (:,:)
      real, allocatable :: gdSACR      (:,:)
      real, allocatable :: gdSACWL     (:,:)
      real, allocatable :: gdSACWR     (:,:)

      real, allocatable :: gdSAtau     (:,:)
      real, allocatable :: gdSAtauM    (:,:)
      real, allocatable :: gdSAMort    (:,:)
      real, allocatable :: gdSAGerm    (:,:)
      real, allocatable :: gdSACol     (:,:)
      real, allocatable :: gdSAExcl    (:,:)

      real              :: nSACount = 1

!     * variables for single point configuration

      real, allocatable :: dSALAI      (:,:)
      real, allocatable :: dSAWMAX     (:,:)
      real, allocatable :: dSARES      (:,:)
      real, allocatable :: dSALIT      (:,:)
      real, allocatable :: dSAFVEG     (:,:)
      real, allocatable :: dSAFH2O     (:,:)
      real, allocatable :: dSAFT       (:,:)

      real, allocatable :: dSAET       (:,:)
      real, allocatable :: dSASEVAP    (:,:)
      real, allocatable :: dSALEVAP    (:,:)
      real, allocatable :: dSABEVAP    (:,:)
      real, allocatable :: dSATRANS    (:,:)

      real, allocatable :: dSAMort     (:,:)
      real, allocatable :: dSAGerm     (:,:)
      real, allocatable :: dSACol      (:,:)
      real, allocatable :: dSAExcl     (:,:)

!     ------------------------------------------------------------------
!     * optimizing version
!     ------------------------------------------------------------------

      integer              :: kOptiYears   = 50
      integer              :: nOptiCnt     = 0
      integer              :: kOptiSeed    = 1604
      real                 :: pOptiFrac    = 0.5
      real                 :: pOptiVar     = 0.05 ! variation of existing species
      real                 :: pOptiFracRan = 0.1  ! propability for random new species

      real, allocatable    :: dOptiGoalAcc      (:,:)
      real, allocatable    :: gdOptiGoalAcc     (:,:)
      real, allocatable    :: gdOptiGoalMax     (:)

      integer, allocatable :: kOptiSorted       (:,:)
      integer, allocatable :: kOptiSortedGlobal (:)
      integer, allocatable :: kOptiBest         (:)

      real, allocatable    :: gOptiBestMap      (:)

      end module jedi_mod

!     ******************************************************************
!     JEDI_ALLOC
!     ******************************************************************

      subroutine jedi_alloc ()
      use jedi_mod
      implicit none

      if (kmypid == kroot) write(kFile_Diag,*) 'jedi_init: jedi_alloc'

      __allocate(klpt_nn,(kNumGPts,8))

!     * allocate climate data fields

      __allocate(pjedi_rain,(NHOR))
      __allocate(pjedi_snow,(NHOR))
      __allocate(pjedi_temp,(NHOR))
      __allocate(pjedi_swdown,(NHOR))
      __allocate(pjedi_lwnet,(NHOR))
      __allocate(cCO2g_a,(NHOR))

!     * allocate plants
      __allocate(kSPPID,(kMaxSPP))

      __allocate(p01,(kMaxSPP))
      __allocate(p02,(kMaxSPP))
      __allocate(p05,(kMaxSPP))
      __allocate(p06,(kMaxSPP))
      __allocate(p07,(kMaxSPP))
      __allocate(p08,(kMaxSPP))
      __allocate(p09,(kMaxSPP))
      __allocate(p10,(kMaxSPP))
      __allocate(p13,(kMaxSPP))
      __allocate(p16,(kMaxSPP))
      __allocate(p14,(kMaxSPP))
      __allocate(p15,(kMaxSPP))
      __allocate(p04,(kMaxSPP))
      __allocate(p17,(kMaxSPP))
      __allocate(p12,(kMaxSPP))
      __allocate(p11,(kMaxSPP))
      __allocate(p03,(kMaxSPP))
      __allocate(p18,(kMaxSPP))
      __allocate(p19,(kMaxSPP))
      __allocate(p20,(kMaxSPP))

!     * allocate process fields

      __allocate(grCA,(gNHOR,kMaxSPP))
      __allocate(grCS,(gNHOR,kMaxSPP))
      __allocate(grCWL,(gNHOR,kMaxSPP))
      __allocate(grCWR,(gNHOR,kMaxSPP))
      __allocate(grCL,(gNHOR,kMaxSPP))
      __allocate(grCR,(gNHOR,kMaxSPP))
      __allocate(grLIVE,(gNHOR,kMaxSPP))

      __allocate(rLIT_CL,(NHOR))
      __allocate(rLIT_CR,(NHOR))
      __allocate(rLIT_CWL,(NHOR))
      __allocate(rLIT_CWR,(NHOR))
      __allocate(rCSLOW,(NHOR))

      __allocate(grLIT_CL,(gNHOR))
      __allocate(grLIT_CR,(gNHOR))
      __allocate(grLIT_CWL,(gNHOR))
      __allocate(grLIT_CWR,(gNHOR))
      __allocate(grCSLOW,(gNHOR))

      __allocate(rW,(NHOR,kMaxSPP))
      __allocate(rWS,(NHOR,kMaxSPP))
      __allocate(rWL,(NHOR,kMaxSPP))
      __allocate(rWSUB,(NHOR,kMaxSPP))
      __allocate(rCA,(NHOR,kMaxSPP))
      __allocate(rCL,(NHOR,kMaxSPP))
      __allocate(rCR,(NHOR,kMaxSPP))
      __allocate(rCWL,(NHOR,kMaxSPP))
      __allocate(rCWR,(NHOR,kMaxSPP))
      __allocate(rCS,(NHOR,kMaxSPP))
      __allocate(rCtot,(NHOR,kMaxSPP))
      __allocate(rLIVE,(NHOR,kMaxSPP))
      __allocate(rGrowW,(NHOR,kMaxSPP))
      __allocate(rGrowG,(NHOR,kMaxSPP))
      __allocate(rGrowT,(NHOR,kMaxSPP))
      __allocate(rDie,(NHOR,kMaxSPP))

      if (kPop_dyn .ge. 1) then
        __allocate(grArea,(gNHOR,kMaxSPP))
        __allocate(grAreaBare,(gNHOR))
        __allocate(rArea,(NHOR,kMaxSPP))
        __allocate(rAreaBare,(NHOR))
      endif

      if (kOpti == 1) then
        __allocate(dOptiGoalAcc,(NHOR,kMaxSPP))
        __allocate(gdOptiGoalAcc,(gNHOR,kMaxSPP))
        __allocate(gdOptiGoalMax,(gNHOR))
        __allocate(kOptiSorted,(gNHOR,kMaxSPP))
        __allocate(kOptiSortedGlobal,(kMaxSPP))
        __allocate(kOptiBest,(kMaxSPP))
        __allocate(gOptiBestMap,(gNHOR))
      endif

      if (kmypid == kroot) write(kFile_Diag,*) 'jedi_alloc: done'

      return
      end subroutine jedi_alloc

!     ******************************************************************
!     JEDI_DEALLOC
!     ******************************************************************

      subroutine jedi_dealloc ()
      use jedi_mod
      implicit none

      if (kmypid == kroot) write(kFile_Diag,*) 'jedi_stop: jedi_dealloc'

      __deallocate(pjedi_rain)
      __deallocate(pjedi_snow)
      __deallocate(pjedi_temp)
      __deallocate(pjedi_swdown)
      __deallocate(pjedi_lwnet)
      __deallocate(cCO2g_a)

      __deallocate(gkGPID)
      __deallocate(pjedi_lat)

      __deallocate(pjedi_area)
      __deallocate(klpt_nn)

      __deallocate(pjedi_paw)

      __deallocate(p01)
      __deallocate(p02)
      __deallocate(p05)
      __deallocate(p06)
      __deallocate(p07)
      __deallocate(p08)
      __deallocate(p09)
      __deallocate(p10)
      __deallocate(p13)
      __deallocate(p16)
      __deallocate(p14)
      __deallocate(p15)
      __deallocate(p04)
      __deallocate(p17)
      __deallocate(p12)
      __deallocate(p11)
      __deallocate(p03)
      __deallocate(p18)
      __deallocate(p19)
      __deallocate(p20)

      __deallocate(kSPPID)

      __deallocate(rW)
      __deallocate(rWS)
      __deallocate(rWL)
      __deallocate(rWSUB)
      __deallocate(rCA)
      __deallocate(rCL)
      __deallocate(rCR)
      __deallocate(rCWL)
      __deallocate(rCWR)
      __deallocate(rCS)
      __deallocate(rCtot)
      __deallocate(rArea)
      __deallocate(rAreaBare)

      __deallocate(rLIVE)

      __deallocate(rGrowW)
      __deallocate(rGrowG)
      __deallocate(rGrowT)
      __deallocate(rDie)

      __deallocate(grArea)
      __deallocate(grAreaBare)
      __deallocate(grCA)
      __deallocate(grCS)
      __deallocate(grCWL)
      __deallocate(grCWR)
      __deallocate(grCL)
      __deallocate(grCR)

      __deallocate(rLIT_CL)
      __deallocate(rLIT_CR)
      __deallocate(rLIT_CWL)
      __deallocate(rLIT_CWR)
      __deallocate(rCSLOW)

      __deallocate(grLIT_CL)
      __deallocate(grLIT_CR)
      __deallocate(grLIT_CWL)
      __deallocate(grLIT_CWR)
      __deallocate(grCSLOW)

      __deallocate(grLIVE)

      __deallocate(dRAbd)
      __deallocate(dGGermCnt)
      __deallocate(dGDeadCnt)

      __deallocate(dGAGPP)
      __deallocate(dGANPP)
      __deallocate(dGARES)
      __deallocate(dGARESH)
      __deallocate(dGALIT)

      __deallocate(dGASOLRAD)
      __deallocate(dGALH)
      __deallocate(dGASH)

      __deallocate(dGACALOSS)
      __deallocate(dGACLLOSS)
      __deallocate(dGACRLOSS)
      __deallocate(dGACWLLOSS)
      __deallocate(dGACWRLOSS)

      __deallocate(dGALOSS_LIT_CS)
      __deallocate(dGALOSS_LIT_CL)
      __deallocate(dGALOSS_LIT_CR)
      __deallocate(dGALOSS_LIT_CWL)
      __deallocate(dGALOSS_LIT_CWR)
      __deallocate(dGALOSS_LIT_CSLOW)

      __deallocate(dGACLALLOC)
      __deallocate(dGACRALLOC)
      __deallocate(dGACWLALLOC)
      __deallocate(dGACWRALLOC)
      __deallocate(dGACSALLOC)

      __deallocate(dGAFH2O)
      __deallocate(dGAFT)

      __deallocate(dGAET)
      __deallocate(dGAQTOT)
      __deallocate(dGAQSURF)
      __deallocate(dGALEVAP)
      __deallocate(dGAESOIL)
      __deallocate(dGATRANS)

      __deallocate(dGAALB)
      __deallocate(dGALAI)
      __deallocate(dGAFVEG)
      __deallocate(dGAFFOR)
      __deallocate(dGAWMAX)

      __deallocate(dGAtau)
      __deallocate(dGAtauM)
      __deallocate(dGAGerm)
      __deallocate(dGACol)
      __deallocate(dGAExcl)
      __deallocate(dGAMort)
      __deallocate(dGADist)
      __deallocate(dGASpec)
      __deallocate(dGAMig)
      __deallocate(dGAExt)

      __deallocate(gdGAGPP)
      __deallocate(gdGANPP)
      __deallocate(gdGANEE)
      __deallocate(gdGARES)
      __deallocate(gdGARESH)
      __deallocate(gdGALIT)
      __deallocate(gdGARESE)

      __deallocate(gdGASOLRAD)
      __deallocate(gdGALH)
      __deallocate(gdGASH)

      __deallocate(gdGACALOSS)
      __deallocate(gdGACLLOSS)
      __deallocate(gdGACRLOSS)
      __deallocate(gdGACWLLOSS)
      __deallocate(gdGACWRLOSS)

      __deallocate(gdGALOSS_LIT_CS)
      __deallocate(gdGALOSS_LIT_CL)
      __deallocate(gdGALOSS_LIT_CR)
      __deallocate(gdGALOSS_LIT_CWL)
      __deallocate(gdGALOSS_LIT_CWR)
      __deallocate(gdGALOSS_LIT_CSLOW)

      __deallocate(gdGACLALLOC)
      __deallocate(gdGACRALLOC)
      __deallocate(gdGACWLALLOC)
      __deallocate(gdGACWRALLOC)
      __deallocate(gdGACSALLOC)

      __deallocate(gdGAFT)
      __deallocate(gdGAFH2O)

      __deallocate(gdGAET)
      __deallocate(gdGAQTOT)
      __deallocate(gdGAQSURF)
      __deallocate(gdGALEVAP)
      __deallocate(gdGAESOIL)
      __deallocate(gdGATRANS)

      __deallocate(gdGAALB)
      __deallocate(gdGALAI)
      __deallocate(gdGAFVEG)
      __deallocate(gdGAFFOR)
      __deallocate(gdGAWMAX)

      __deallocate(gdGAtau)
      __deallocate(gdGAtauM)

      __deallocate(gdGAGerm)
      __deallocate(gdGACol)
      __deallocate(gdGAExcl)
      __deallocate(gdGAMort)
      __deallocate(gdGADist)
      __deallocate(gdGASpec)
      __deallocate(gdGAMig)
      __deallocate(gdGAExt)

      __deallocate(dSARAbd)
      __deallocate(dSAGPP)
      __deallocate(dSANPP)
      __deallocate(dSACS)
      __deallocate(dSACA)
      __deallocate(dSACL)
      __deallocate(dSACR)
      __deallocate(dSACWL)
      __deallocate(dSACWR)

      __deallocate(dSAtau)
      __deallocate(dSAtauM)
      __deallocate(dSAMort)
      __deallocate(dSAGerm)
      __deallocate(dSACol)
      __deallocate(dSAExcl)

      __deallocate(gdSARAbd)
      __deallocate(gdSAGPP)
      __deallocate(gdSANPP)
      __deallocate(gdSACS)
      __deallocate(gdSACA)
      __deallocate(gdSACL)
      __deallocate(gdSACR)
      __deallocate(gdSACWL)
      __deallocate(gdSACWR)

      __deallocate(gdSAtau)
      __deallocate(gdSAtauM)
      __deallocate(gdSAMort)
      __deallocate(gdSAGerm)
      __deallocate(gdSACol)
      __deallocate(gdSAExcl)

      __deallocate(dSALAI)
      __deallocate(dSAWMAX)
      __deallocate(dSARES)
      __deallocate(dSALIT)
      __deallocate(dSAFVEG)
      __deallocate(dSAFH2O)
      __deallocate(dSAFT)

      __deallocate(dSAET)
      __deallocate(dSASEVAP)
      __deallocate(dSALEVAP)
      __deallocate(dSABEVAP)
      __deallocate(dSATRANS)

      __deallocate(dOptiGoalAcc)
      __deallocate(gdOptiGoalAcc)
      __deallocate(gdOptiGoalMax)

      __deallocate(kOptiSorted)
      __deallocate(kOptiSortedGlobal)
      __deallocate(kOptiBest)

      __deallocate(gOptiBestMap)

      __deallocate(BM_total)
      __deallocate(Richness)

      return
      end subroutine jedi_dealloc
