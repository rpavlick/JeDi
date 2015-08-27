#define __ACTIVATE
#include "../globe/globe_macros.f90"

      module jedi_dyn_mod

      REAL,allocatable,dimension(:) :: zWMAX
      REAL,allocatable,dimension(:) :: zALB
      REAL,allocatable,dimension(:) :: zFVEG
      REAL,allocatable,dimension(:) :: zSLA
      REAL,allocatable,dimension(:) :: zLAI
      REAL,allocatable,dimension(:) :: zFFOR
      REAL,allocatable,dimension(:) :: zFSNOW
      REAL,allocatable,dimension(:) :: zWLMAX
      REAL,allocatable,dimension(:) :: zFW
      REAL,allocatable,dimension(:) :: BareFW
      REAL,allocatable,dimension(:) :: zFT2
      REAL,allocatable,dimension(:) :: zFT3
      REAL,allocatable,dimension(:) :: zFT
      REAL,allocatable,dimension(:) :: zA
      REAL,allocatable,dimension(:) :: zB
      REAL,allocatable,dimension(:) :: zThFall
      REAL,allocatable,dimension(:) :: zTRANS
      REAL,allocatable,dimension(:) :: zRunoff
      REAL,allocatable,dimension(:) :: zInfil
      REAL,allocatable,dimension(:) :: zMelt
      REAL,allocatable,dimension(:) :: zDrain
      REAL,allocatable,dimension(:) :: zSET
      REAL,allocatable,dimension(:) :: zFH2O
      REAL,allocatable,dimension(:) :: zPAR
      REAL,allocatable,dimension(:) :: zGPP
      REAL,allocatable,dimension(:) :: zAMAX
      REAL,allocatable,dimension(:) :: zNPP
      REAL,allocatable,dimension(:,:) :: zRES
      REAL,allocatable,dimension(:) :: zRESH
      REAL,allocatable,dimension(:) :: zBETA

!     * phenology

      REAL,allocatable,dimension(:) :: ff_GrowL
      REAL,allocatable,dimension(:) :: ff_GrowR
      REAL,allocatable,dimension(:) :: ff_GrowT
      REAL,allocatable,dimension(:) :: ff_GrowW
      REAL,allocatable,dimension(:) :: ff_GrowG
      REAL,allocatable,dimension(:) :: ff_Grow
      REAL,allocatable,dimension(:) :: ff_Die
      REAL,allocatable,dimension(:) :: ff_Seed
      REAL,allocatable,dimension(:,:) :: fGerm
      REAL,allocatable,dimension(:,:) :: fGermFix
      REAL,allocatable,dimension(:,:) :: fLIT

!     * allocation

      REAL,allocatable,dimension(:) :: ff_A
      REAL,allocatable,dimension(:) :: ff_AL
      REAL,allocatable,dimension(:) :: ff_AS
      REAL,allocatable,dimension(:) :: ff_AR
      REAL,allocatable,dimension(:) :: ff_LW
      REAL,allocatable,dimension(:) :: ff_RW
      REAL,allocatable,dimension(:) :: zALLOC
      REAL,allocatable,dimension(:) :: zALLOC_CL
      REAL,allocatable,dimension(:) :: zALLOC_CR
      REAL,allocatable,dimension(:) :: zALLOC_CWL
      REAL,allocatable,dimension(:) :: zALLOC_CWR
      REAL,allocatable,dimension(:) :: zALLOC_CS

      REAL,allocatable,dimension(:) :: zLOSS_CA
      REAL,allocatable,dimension(:) :: zLOSS_CS
      REAL,allocatable,dimension(:) :: zLOSS_CL
      REAL,allocatable,dimension(:) :: zLOSS_CR
      REAL,allocatable,dimension(:) :: zLOSS_CWL
      REAL,allocatable,dimension(:) :: zLOSS_CWR

      REAL,allocatable,dimension(:) :: zGALOSS_CA
      REAL,allocatable,dimension(:) :: zGALOSS_CL
      REAL,allocatable,dimension(:) :: zGALOSS_CR
      REAL,allocatable,dimension(:) :: zGALOSS_CWL
      REAL,allocatable,dimension(:) :: zGALOSS_CWR

      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CS
      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CL
      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CR
      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CWL
      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CWR
      REAL,allocatable,dimension(:) :: zGALOSS_LIT_CSLOW

      REAL,allocatable,dimension(:) :: fff_
      REAL,allocatable,dimension(:) :: ffx_
      REAL,allocatable,dimension(:) :: zDW

      REAL,allocatable,dimension(:) :: zPET
      REAL,allocatable,dimension(:) :: zLEVAP
      REAL,allocatable,dimension(:) :: zBEVAP
      REAL,allocatable,dimension(:) :: zSSG
      REAL,allocatable,dimension(:) :: zzSSG
      REAL,allocatable,dimension(:) :: zSALB
      REAL,allocatable,dimension(:) :: zSALBMAX
      REAL,allocatable,dimension(:) :: zSALBMIN

      REAL,allocatable,dimension(:) :: zSNOWALB
      REAL,allocatable,dimension(:) :: zRNET
      REAL,allocatable,dimension(:) :: zSRNET
      REAL,allocatable,dimension(:) :: zSEVAP
      REAL,allocatable,dimension(:) :: zWL
      REAL,allocatable,dimension(:) :: zSubDrain

      REAL,allocatable,dimension(:) :: zdlen         ! day length

!     * abundance dynamics

      REAL,allocatable,dimension(:,:) :: zDeltaRAbd
      REAL,allocatable,dimension(:,:) :: zDeltaCtot
      REAL,allocatable,dimension(:)   :: zgCtot

!     * area dynamics

      real :: InitArea
      real :: pCompTau

      REAL,allocatable,dimension(:)   :: dfAreaGerm
      REAL,allocatable,dimension(:)   :: dfAreaComp
      REAL,allocatable,dimension(:)   :: dfAreaExcl
      REAL,allocatable,dimension(:)   :: dfAreaM
      REAL,allocatable,dimension(:)   :: sumArea

      REAL,allocatable,dimension(:,:) :: dominance
      REAL,allocatable,dimension(:,:) :: comp
      REAL,allocatable,dimension(:)   :: sumCab
      REAL,allocatable,dimension(:)   :: maxCab

      REAL,allocatable,dimension(:)   :: NPP0
      REAL,allocatable,dimension(:)   :: tau
      REAL,allocatable,dimension(:,:) :: mTau
      REAL,allocatable,dimension(:)   :: fG
      REAL,allocatable,dimension(:)   :: sumfG

      REAL :: dfAreaBare
      REAL :: sumdfAreaG

!     * migration and speciation

      real,allocatable,dimension(:,:) :: fseed       ! Seedflux
      real,allocatable,dimension(:,:) :: gfseed      ! Seedflux
      real,allocatable,dimension(:,:) :: gAfseed     ! Seedflux
      real,allocatable,dimension(:,:) :: Afseed      ! Accumulated Seedflux for speciation

      real,allocatable,dimension(:,:) :: fseedout    ! Seed out flow via migration
      real,allocatable,dimension(:,:) :: gfseedout
      real,allocatable,dimension(:,:) :: gfseedin

      integer,allocatable,dimension(:):: inFreePosSPP
      integer                         :: firstFree=0, lastFree=0, count=0

      real,allocatable,dimension(:)   :: countSpec, countExt
      real,allocatable,dimension(:,:) :: Ext

      end module

!     ******************************************************************
!     JEDI_DYN_ALLOC
!     ******************************************************************

      subroutine jedi_dyn_alloc ()
      use jedi_dyn_mod
      use jedi_mod
      implicit none

!     * allocate temporary variables

      __allocate(zWMAX,(NHOR))
      __allocate(zALB,(NHOR))
      __allocate(zSLA,(NHOR))
      __allocate(zFVEG,(NHOR))
      __allocate(zLAI,(NHOR))
      __allocate(zFFOR,(NHOR))
      __allocate(zFSNOW,(NHOR))
      __allocate(zWLMAX,(NHOR))
      __allocate(zFW,(NHOR))
      __allocate(BareFW,(NHOR))
      __allocate(zFT2,(NHOR))
      __allocate(zFT3,(NHOR))
      __allocate(zFT,(NHOR))
      __allocate(zA,(NHOR))
      __allocate(zB,(NHOR))
      __allocate(zThFall,(NHOR))
      __allocate(zTRANS,(NHOR))
      __allocate(zRunoff,(NHOR))
      __allocate(zInfil,(NHOR))
      __allocate(zMelt,(NHOR))
      __allocate(zDrain,(NHOR))
      __allocate(zSET,(NHOR))
      __allocate(zFH2O,(NHOR))
      __allocate(zPAR,(NHOR))
      __allocate(zGPP,(NHOR))
      __allocate(zAMAX,(NHOR))

      __allocate(zNPP,(NHOR))
      __allocate(zRES,(NHOR,kMaxSPP))
      __allocate(zRESH,(NHOR))
      __allocate(zBETA,(NHOR))

      __allocate(ff_GrowT,(NHOR))
      __allocate(ff_GrowW,(NHOR))
      __allocate(ff_GrowG,(NHOR))
      __allocate(ff_Grow,(NHOR))
      __allocate(ff_Die,(NHOR))
      __allocate(ff_Seed,(NHOR))
      __allocate(ff_A,(NHOR))
      __allocate(ff_AL,(NHOR))
      __allocate(ff_AS,(NHOR))
      __allocate(ff_AR,(NHOR))
      __allocate(ff_GrowL,(NHOR))
      __allocate(ff_GrowR,(NHOR))
      __allocate(ff_LW,(NHOR))
      __allocate(ff_RW,(NHOR))
      __allocate(zALLOC,(NHOR))
      __allocate(zALLOC_CL,(NHOR))
      __allocate(zALLOC_CR,(NHOR))
      __allocate(zALLOC_CWL,(NHOR))
      __allocate(zALLOC_CWR,(NHOR))
      __allocate(zALLOC_CS,(NHOR))

      __allocate(zLOSS_CA,(NHOR))
      __allocate(zLOSS_CS,(NHOR))
      __allocate(zLOSS_CL,(NHOR))
      __allocate(zLOSS_CR,(NHOR))
      __allocate(zLOSS_CWL,(NHOR))
      __allocate(zLOSS_CWR,(NHOR))

      __allocate(zGALOSS_CA,(NHOR))
      __allocate(zGALOSS_CL,(NHOR))
      __allocate(zGALOSS_CR,(NHOR))
      __allocate(zGALOSS_CWL,(NHOR))
      __allocate(zGALOSS_CWR,(NHOR))

      __allocate(zGALOSS_LIT_CS,(NHOR))
      __allocate(zGALOSS_LIT_CL,(NHOR))
      __allocate(zGALOSS_LIT_CR,(NHOR))
      __allocate(zGALOSS_LIT_CWL,(NHOR))
      __allocate(zGALOSS_LIT_CWR,(NHOR))
      __allocate(zGALOSS_LIT_CSLOW,(NHOR))

      __allocate(fff_,(NHOR))
      __allocate(ffx_,(NHOR))
      __allocate(zDW,(NHOR))

      __allocate(zPET,(NHOR))
      __allocate(zLEVAP,(NHOR))
      __allocate(zBEVAP,(NHOR))
      __allocate(zSSG,(NHOR))
      __allocate(zzSSG,(NHOR))
      __allocate(zSALB,(NHOR))
      __allocate(zSALBMAX,(NHOR))
      __allocate(zSALBMIN,(NHOR))

      __allocate(zSNOWALB,(NHOR))
      __allocate(zRNET,(NHOR))
      __allocate(zSRNET,(NHOR))

      __allocate(zSEVAP,(NHOR))
      __allocate(zWL,(NHOR))
      __allocate(zSubDrain,(NHOR))

      __allocate(zdlen,(NHOR))
      __allocate(fLIT,(NHOR,kMaxSPP))
      __allocate(fGerm,(NHOR,kMaxSPP))
      __allocate(fGermFix,(NHOR,kMaxSPP))

!     * abundance dynamics

      if (kPop_dyn .eq. 1) then
        __allocate(zDeltaRAbd,(NHOR,kMaxSPP))
        __allocate(zDeltaCtot,(NHOR,kMaxSPP))
        __allocate(zgCtot,(NHOR))
      endif

!     * area dynamics

      if (kPop_dyn .eq. 2) then
        __allocate(Comp,(kMaxSPP,kMaxSPP))
        __allocate(NPP0,(kMaxSPP))
        __allocate(tau,(kMaxSPP))
        __allocate(mTau,(NHOR,kMaxSPP))
        __allocate(dfAreaM,(kMaxSPP))
        __allocate(dfAreaGerm,(kMaxSPP))
        __allocate(dfAreaComp,(kMaxSPP))
        __allocate(dfAreaExcl,(kMaxSPP))
        __allocate(fG,(kMaxSPP))
        __allocate(sumfG,(NHOR))
        __allocate(sumArea,(NHOR))
        __allocate(sumCab,(NHOR))
        __allocate(maxCab,(NHOR))
        __allocate(dominance,(NHOR,kMaxSPP))
      endif

      if (kMig_dyn .gt. 0) then
        __allocate(gfseedout,(gNHOR,kMaxSPP))
        __allocate(gfseedin,(gNHOR,kMaxSPP))
        __allocate(fseedout,(NHOR,kMaxSPP))
        fseedout(:,:) = 0.0
      endif

      if (kMig_dyn .gt. 0 .or. kSpec_dyn .gt. 0) then
        __allocate(fseed,(NHOR,kMaxSPP))
        __allocate(gfseed,(NHOR,kMaxSPP))
      endif

      if (kSpec_dyn .gt. 0) then
        __allocate(gAfseed,(gNHOR,kMaxSPP))
        __allocate(Afseed,(NHOR,kMaxSPP))

        __allocate(inFreePosSPP,(kMaxSPP))
        __allocate(Ext,(NHOR,kMaxSPP))
        Ext(:,:) = 0.0
        __allocate(countExt,(NHOR))
        countExt(:) = 0.0
        __allocate(countSpec,(NHOR))
        countSpec(:) = 0.0
      endif

      return
      end subroutine jedi_dyn_alloc

!     ******************************************************************
!     JEDI_DYN_DEALLOC
!     ******************************************************************

      subroutine jedi_dyn_dealloc ()
      use jedi_dyn_mod
      implicit none

!     * deallocate temporary variables

      __deallocate(zWMAX)
      __deallocate(zALB)
      __deallocate(zFVEG)
      __deallocate(zSLA)
      __deallocate(zLAI)
      __deallocate(zFFOR)
      __deallocate(zFSNOW)
      __deallocate(zWLMAX)
      __deallocate(zFW)
      __deallocate(BareFW)
      __deallocate(zFT2)
      __deallocate(zFT3)
      __deallocate(zFT)
      __deallocate(zA)
      __deallocate(zB)
      __deallocate(zThFall)
      __deallocate(zTRANS)
      __deallocate(zRunoff)
      __deallocate(zInfil)
      __deallocate(zMelt)
      __deallocate(zDrain)
      __deallocate(zSET)
      __deallocate(zFH2O)
      __deallocate(zPAR)
      __deallocate(zGPP)
      __deallocate(zAMAX)
      __deallocate(zNPP)
      __deallocate(zRES)
      __deallocate(zRESH)
      __deallocate(zBETA)

      __deallocate(ff_GrowL)
      __deallocate(ff_GrowR)
      __deallocate(ff_GrowT)
      __deallocate(ff_GrowW)
      __deallocate(ff_GrowG)
      __deallocate(ff_Grow)
      __deallocate(ff_Die)
      __deallocate(ff_Seed)
      __deallocate(fGerm)
      __deallocate(fGermFix)
      __deallocate(fLIT)

      __deallocate(ff_A)
      __deallocate(ff_AL)
      __deallocate(ff_AS)
      __deallocate(ff_AR)
      __deallocate(ff_LW)
      __deallocate(ff_RW)
      __deallocate(zALLOC)
      __deallocate(zALLOC_CL)
      __deallocate(zALLOC_CR)
      __deallocate(zALLOC_CWL)
      __deallocate(zALLOC_CWR)
      __deallocate(zALLOC_CS)

      __deallocate(zLOSS_CL)
      __deallocate(zLOSS_CR)
      __deallocate(zLOSS_CWL)
      __deallocate(zLOSS_CWR)
      __deallocate(zLOSS_CA)
      __deallocate(zLOSS_CS)

      __deallocate(zGALOSS_CA)
      __deallocate(zGALOSS_CL)
      __deallocate(zGALOSS_CR)
      __deallocate(zGALOSS_CWL)
      __deallocate(zGALOSS_CWR)

      __deallocate(zGALOSS_LIT_CS)
      __deallocate(zGALOSS_LIT_CL)
      __deallocate(zGALOSS_LIT_CR)
      __deallocate(zGALOSS_LIT_CWL)
      __deallocate(zGALOSS_LIT_CWR)
      __deallocate(zGALOSS_LIT_CSLOW)

      __deallocate(fff_)
      __deallocate(ffx_)
      __deallocate(zDW)

      __deallocate(zPET)
      __deallocate(zLEVAP)
      __deallocate(zBEVAP)
      __deallocate(zSSG)
      __deallocate(zzSSG)
      __deallocate(zSALB)
      __deallocate(zSALBMAX)
      __deallocate(zSALBMIN)

      __deallocate(zSNOWALB)
      __deallocate(zRNET)
      __deallocate(zSRNET)

      __deallocate(zSEVAP)
      __deallocate(zWL)
      __deallocate(zSubDrain)

      __deallocate(zdlen)

!     * abundance dynamics

      __deallocate(zDeltaRAbd)
      __deallocate(zDeltaCtot)
      __deallocate(zgCtot)

!     * area dynamics

      __deallocate(Comp)
      __deallocate(dfAreaM)
      __deallocate(dfAreaGerm)
      __deallocate(dfAreaComp)
      __deallocate(dfAreaExcl)
      __deallocate(NPP0)
      __deallocate(tau)
      __deallocate(mTau)
      __deallocate(fG)
      __deallocate(sumfG)
      __deallocate(sumArea)
      __deallocate(sumCab)
      __deallocate(maxCab)
      __deallocate(dominance)

!     * Migration & Speciation

      __deallocate(fseed)
      __deallocate(gfseed)
      __deallocate(gAfseed)
      __deallocate(Afseed)

      __deallocate(fseedout)
      __deallocate(gfseedout)
      __deallocate(gfseedin)

      __deallocate(inFreePosSPP)

      __deallocate(countSpec)
      __deallocate(countExt)
      __deallocate(Ext)

      return
      end subroutine jedi_dyn_dealloc
