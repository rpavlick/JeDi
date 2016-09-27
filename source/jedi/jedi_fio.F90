#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jedi_fio.F90
!> \brief JEDI

!     ******************************************************************
!     JEDI_READ_NAMELIST
!     ******************************************************************

      subroutine jedi_read_namelist ()
      use jedi_mod
      implicit none

      integer :: ioerr

!     * define namelist

      namelist /jedipar/ &
     &  kPop_dyn,        &
     &  kSpec_dyn,       &
     &  kMig_dyn,        &
     &  kDist_dyn,       &
     &  kOpti,           &
     &  kOptiYears,      &
     &  pOptiFrac,       &
     &  pOptiVar,        &
     &  pOptiFracRan,    &
     &  pA0,             &
     &  pSeedFix,        &
     &  pC_GPP,          &
     &  pC_RES,          &
     &  pC_LAI,          &
     &  pC_TRS,          &
     &  pC_WSMAX,        &
     &  pC_WLMAX,        &
     &  pC_FFOR,         &
     &  kFrPl,           &
     &  kToPl,           &
     &  kMaxSPP,         &
     &  pPAW0,           &
     &  pAlb0,           &
     &  pAlb1,           &
     &  pT_CRIT,         &
     &  pQ10H,           &
     &  pQ10R,           &
     &  pSeedTau,        &
     &  pMortTau,        &
     &  pDistTau,        &
     &  pSoilTau,        &
     &  pRAbdTau,        &
     &  kGermFix,        &
     &  kCbal,           &
     &  kCWT,            &
     &  pT_GPP1,         &
     &  pT_GPP2,         &
     &  cPETCorr,        &
     &  pN_RES,          &
     &  pCN_Wood,        &
     &  pCN_Leaf,        &
     &  pCN_Root,        &
     &  cT_MIN,          &
     &  pmig,            &
     &  pSpec,           &
     &  initGP,          &
     &  pdaylen,         &
     &  pfGerm,          &
     &  pfComp,          &
     &  pCO2sens,        &
     &  kspinup,         &
     &  seed,            &
     &  ktopspecies,     &
     &  kdynamic,        &
     &  kspec_yrout,     &
     &  koutput


      __diag(kFile_Diag,'jedi_init: read namelist parameters')

!     * read namelist parameters

      if (kmypid == kroot) then
        open(kFile_Namelist, FILE=sfile_Namelist, STATUS='old',        &
     &                       FORM='formatted', IOSTAT=ioerr)
        if (ioerr .eq. 0) then
          read(kFile_Namelist, jedipar)
          if (kdiag .eq. 1) then
            write(kFile_Diag,'(" *****************************************")')
            write(kFile_Diag,'(" * Namelist JEDIPAR from ",A)') sfile_Namelist
            write(kFile_Diag,'(" *****************************************")')
          endif
          close(kFile_Namelist)
        else
          if (kdiag .eq. 1) then
            write(kFile_Diag,'(" **************************************")')
            write(kFile_Diag,'(" * ERROR reading file ",A)') sfile_Namelist
            write(kFile_Diag,'(" * Using default values")')
            write(kFile_Diag,'(" **************************************")')
          endif
        endif
        if (kdiag .eq. 1) then
          write(kFile_Diag,jedipar)
          flush(kFile_Diag)
        endif

        if (kCWT == 1) then
          kFrPl = 1
          kToPl = gNHOR
        endif

        kNumSPP = (kToPl - kFrPl) + 1
        if (kSpec_dyn .gt. 0) then
          kMaxSPP = MAX(kMaxSPP, kNumSPP)
        else
          kMaxSPP = kNumSPP
        endif

!       * add one for bare soil
        kMaxSPP = kMaxSPP + 1

!       * convert namelist parameter timescales
        pseedtau = 1.0 / (pseedtau * cDaysPerYear)
        pmorttau = 1.0 / (pmorttau * cDaysPerYear)
        pdisttau = 1.0 / (pdisttau * cDaysPerYear)
        pRAbdTau = 1.0 / (pRAbdTau * cDaysPerYear)
        pSoilTau = 1.0 / (psoiltau * cDaysPerYear)

      endif

!     * broadcast namelist parameters

      __globe_mpbc(kMaxSPP)
      __globe_mpbc(kNumSPP)

      __globe_mpbc(kPop_dyn)
      __globe_mpbc(kSpec_dyn)
      __globe_mpbc(kMig_dyn)
      __globe_mpbc(kDist_dyn)

      __globe_mpbc(kOpti)
      __globe_mpbc(kOptiYears)
      __globe_mpbc(pOptiFrac)
      __globe_mpbc(pOptiVar)
      __globe_mpbc(pOptiFracRan)

      __globe_mpbc(pq10R)
      __globe_mpbc(kFrPl)
      __globe_mpbc(kToPl)
      __globe_mpbc(kMaxSPP)

      __globe_mpbc(pA0)
      __globe_mpbc(pSeedFix)
      __globe_mpbc(kGermFix)
      __globe_mpbc(kCbal)
      __globe_mpbc(kCWT)

      __globe_mpbc(pC_GPP)
      __globe_mpbc(pC_RES)
      __globe_mpbc(pC_LAI)
      __globe_mpbc(pC_TRS)
      __globe_mpbc(pC_WSMAX)
      __globe_mpbc(pC_WLMAX)
      __globe_mpbc(pC_FFOR)
      __globe_mpbc(pPAW0)
      __globe_mpbc(pAlb0)
      __globe_mpbc(PAlb1)
      __globe_mpbc(pT_CRIT)
      __globe_mpbc(pQ10H)

      __globe_mpbc(pSeedTau)
      __globe_mpbc(pMortTau)
      __globe_mpbc(pDistTau)
      __globe_mpbc(pSoilTau)
      __globe_mpbc(pRAbdTau)

      __globe_mpbc(pT_GPP1)
      __globe_mpbc(pT_GPP2)
      __globe_mpbc(cPETCorr)

      __globe_mpbc(pN_RES)
      __globe_mpbc(pCN_Leaf)
      __globe_mpbc(pCN_Root)
      __globe_mpbc(pCN_Wood)
      __globe_mpbc(cT_MIN)

      __globe_mpbc(pmig)
      __globe_mpbc(pSpec)
      __globe_mpbc(initGP)
      __globe_mpbc(pdaylen)
      __globe_mpbc(pfGerm)
      __globe_mpbc(pfComp)
      __globe_mpbc(pCO2sens)

      __globe_mpbc(kspinup)
      __globe_mpbc(seed)
      __globe_mpbc(ktopspecies)
      __globe_mpbc(kdynamic)
      __globe_mpbc(kspec_yrout)
      __globe_mpbc(koutput)

      return
      end subroutine jedi_read_namelist

!     ******************************************************************
!     JEDI_READ_SPECIES
!     ******************************************************************

      subroutine jedi_read_species ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, i
      integer :: ioerr
      character :: dummy

      __diag(kFile_Diag,'jedi_init: read specparms file')

      if (kmypid == kroot) then
        open (kFile_SPPParm, FILE=sfile_SPPParm, STATUS='old',         &
     &                       FORM='FORMATTED', IOSTAT=ioerr)
        call error_open_old(ioerr,sfile_SPPParm)
        do i = 1, kFrPl - 1
          read(kFile_SPPParm,'(A1)', IOSTAT=ioerr) dummy
          __error_read(ioerr,sfile_SPPParm)
        enddo
        do i = kFrPl, kToPl
          iSPP = i - kFrPl + FirstSPP
          read(kFile_SPPParm,'(20f6.3)', IOSTAT=ioerr)                 &
     &          p01(iSPP), p02(iSPP), p03(iSPP), p04(iSPP), p05(iSPP), &
     &          p06(iSPP), p07(iSPP), p08(iSPP), p09(iSPP), p10(iSPP), &
     &          p11(iSPP), p12(iSPP), p13(iSPP), p14(iSPP), p15(iSPP), &
     &          p16(iSPP), p17(iSPP), p18(iSPP), p19(iSPP), p20(iSPP)
          __error_read(ioerr,sfile_SPPParm)
          kSPPID(iSPP) = iSPP                !* Fill Species ID
        enddo
        close(kFile_SPPParm)
        if (FirstSPP .gt. 1) then
          do iSPP = 1, FirstSPP - 1
            p01(iSPP) = 0.0
            p02(iSPP) = 0.0
            p05(iSPP) = 0.0
            p06(iSPP) = 0.0
            p07(iSPP) = 0.0
            p08(iSPP) = 0.0
            p09(iSPP) = 0.0
            p10(iSPP) = 0.0
            p13(iSPP) = 0.0
            p16(iSPP) = 0.0
            p14(iSPP) = 0.0
            p15(iSPP) = 0.0
            p04(iSPP) = 0.0
            p17(iSPP) = 0.0
            p12(iSPP) = 0.0
            p11(iSPP) = 0.0
            p03(iSPP) = 0.0
            p18(iSPP) = 0.0
            p19(iSPP) = 0.0
            p20(iSPP) = 0.0
          enddo
        endif
      endif

!     * broadcast species

      __globe_mpbc(p01)
      __globe_mpbc(p02)
      __globe_mpbc(p05)
      __globe_mpbc(p06)
      __globe_mpbc(p07)
      __globe_mpbc(p08)
      __globe_mpbc(p09)
      __globe_mpbc(p10)
      __globe_mpbc(p13)
      __globe_mpbc(p16)
      __globe_mpbc(p14)
      __globe_mpbc(p15)
      __globe_mpbc(p04)
      __globe_mpbc(p17)
      __globe_mpbc(p12)
      __globe_mpbc(p11)
      __globe_mpbc(p03)
      __globe_mpbc(p18)
      __globe_mpbc(p19)
      __globe_mpbc(p20)

      kNumSPPID = kNumSPP

      do i = kFrPl, kToPl
        iSPP = i - kFrPl + FirstSPP
        call jedi_transform_parms(iSPP)
      enddo

!     * fill up list with free Positions ini species arrays

      if ( (kSpec_dyn .ge. 1) .and. (kmypid == kroot) ) then
        firstFree = 1 ! index for first free position in SPP arrays (p01,p02,...)
        do i = kNumSPP+1, kMaxSPP
          if (lastFree .eq. kMaxSPP) then
            lastFree = 1
          else
            lastFree = lastFree + 1
          endif
          inFreePosSPP(lastFree) = i
          count = count + 1
        enddo
      endif

      return
      end subroutine jedi_read_species

!     ******************************************************************
!     JEDI_READ_RESTART
!     ******************************************************************

      subroutine jedi_read_restart ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, ioerr

      __diag(kFile_Diag,'jedi_init: read restart')

      if (kmypid == kroot) then
        open(kFile_Restart, FILE=sfile_Restart, STATUS='old',          &
     &                      FORM='unformatted', IOSTAT=ioerr)
        call error_open_old(ioerr,sfile_Restart)

        read(kFile_Restart) kMaxSPP
        read(kFile_Restart) kNumSPP
        read(kFile_Restart) kNumSPPID
      endif

      __globe_mpbc(kMaxSPP)
      __globe_mpbc(kNumSPP)
      __globe_mpbc(kNumSPPID)

!     * allocate fields

      call jedi_alloc
      if (kPop_dyn .ge. 1) call jedi_dyn_alloc
      call jedi_output_init

      if (kmypid == kroot) then
        read(kFile_Restart) kSPPID

        read(kFile_Restart) p01
        read(kFile_Restart) p02
        read(kFile_Restart) p03
        read(kFile_Restart) p04
        read(kFile_Restart) p05
        read(kFile_Restart) p06
        read(kFile_Restart) p07
        read(kFile_Restart) p08
        read(kFile_Restart) p09
        read(kFile_Restart) p10
        read(kFile_Restart) p11
        read(kFile_Restart) p12
        read(kFile_Restart) p13
        read(kFile_Restart) p14
        read(kFile_Restart) p15
        read(kFile_Restart) p16
        read(kFile_Restart) p17
        read(kFile_Restart) p18
        read(kFile_Restart) p19
        read(kFile_Restart) p20

        do iSPP = FirstSPP, kMaxSPP
          call jedi_transform_parms(iSPP)
        enddo

        read(kFile_Restart) nOptiCnt

        if (kSpec_dyn .gt. 0) read(kFile_Restart) count
        if (kSpec_dyn .gt. 0) read(kFile_Restart) firstFree
        if (kSpec_dyn .gt. 0) read(kFile_Restart) lastFree
        if (kSpec_dyn .gt. 0) read(kFile_Restart) inFreePosSPP

        read(kFile_Restart) nSACount
      endif

      __globe_mpbc(kSPPID)

      __globe_mpbc(p01)
      __globe_mpbc(p02)
      __globe_mpbc(p05)
      __globe_mpbc(p06)
      __globe_mpbc(p07)
      __globe_mpbc(p08)
      __globe_mpbc(p09)
      __globe_mpbc(p10)
      __globe_mpbc(p13)
      __globe_mpbc(p16)
      __globe_mpbc(p14)
      __globe_mpbc(p15)
      __globe_mpbc(p04)
      __globe_mpbc(p17)
      __globe_mpbc(p12)
      __globe_mpbc(p11)
      __globe_mpbc(p03)
      __globe_mpbc(p18)
      __globe_mpbc(p19)
      __globe_mpbc(p20)

      __globe_mpbc(nOptiCnt)

      __globe_mpbc(nSACount)

      call globe_mpreadgp(kFile_Restart,rW,    kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rWS,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rWL,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rWSUB, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,rCA,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCL,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCR,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCWL,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCWR,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCS,   kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,rLIT_CL,   1,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rLIT_CR,   1,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rLIT_CWL,  1,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rLIT_CWR,  1,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rCSLOW,    1,NHOR,kNumGPts,gNHOR)

      if (kPop_dyn .ge. 1) call globe_mpreadgp(kFile_Restart,rArea,kMaxSPP,NHOR,kNumGPts,gNHOR)
      if (kPop_dyn .ge. 1) call globe_mpreadgp(kFile_Restart,rAreaBare,1,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,rLIVE, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,rGrowW,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rGrowG,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rGrowT,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,rDie,  kMaxSPP,NHOR,kNumGPts,gNHOR)

      if (kMig_dyn .gt. 0 .or. kSpec_dyn .gt. 0)                       &
     &  call globe_mpreadgp(kFile_Restart,fseed,kMaxSPP,NHOR,kNumGPts,gNHOR)
      if (kMig_dyn .gt. 0) call globe_mpreadgp(kFile_Restart,fseedout,kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,dSARAbd,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSAGPP, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSANPP, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACS,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACA,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACL,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACR,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACWL, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACWR, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpreadgp(kFile_Restart,dSACTOT, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpreadgp(kFile_Restart,dSACVEG, kMaxSPP,NHOR,kNumGPts,gNHOR)

      if (kmypid == kroot) close(kFile_Restart, STATUS='DELETE')

      return
      end subroutine jedi_read_restart

!     ******************************************************************
!     JEDI_WRITE_RESTART
!     ******************************************************************

      subroutine jedi_write_restart ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, ioerr

      __diag(kFile_Diag,'jedi_stop: write restart')

      if (kmypid == kroot) then
        open(kFile_Restart, FILE=sfile_Restart, STATUS='replace',      &
     &                      FORM='unformatted', IOSTAT=ioerr)
        call error_open_new(ioerr,sfile_Restart)

        write(kFile_Restart) kMaxSPP
        write(kFile_Restart) kNumSPP
        write(kFile_Restart) kNumSPPID

        write(kFile_Restart) kSPPID

        do iSPP = FirstSPP, kMaxSPP
          call jedi_transform_parms_back(iSPP)
        enddo

        write(kFile_Restart) p01
        write(kFile_Restart) p02
        write(kFile_Restart) p03
        write(kFile_Restart) p04
        write(kFile_Restart) p05
        write(kFile_Restart) p06
        write(kFile_Restart) p07
        write(kFile_Restart) p08
        write(kFile_Restart) p09
        write(kFile_Restart) p10
        write(kFile_Restart) p11
        write(kFile_Restart) p12
        write(kFile_Restart) p13
        write(kFile_Restart) p14
        write(kFile_Restart) p15
        write(kFile_Restart) p16
        write(kFile_Restart) p17
        write(kFile_Restart) p18
        write(kFile_Restart) p19
        write(kFile_Restart) p20

        write(kFile_Restart) nOptiCnt

        if (kSpec_dyn .gt. 0) write(kFile_Restart) count
        if (kSpec_dyn .gt. 0) write(kFile_Restart) firstFree
        if (kSpec_dyn .gt. 0) write(kFile_Restart) lastFree
        if (kSpec_dyn .gt. 0) write(kFile_Restart) inFreePosSPP

        write(kFile_Restart) nSACount
      endif

      call globe_mpwritegp(kFile_Restart,rW,    kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rWS,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rWL,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rWSUB, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,rCA,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCL,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCR,   kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCWL,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCWR,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCS,   kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,rLIT_CL,   1,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rLIT_CR,   1,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rLIT_CWL,  1,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rLIT_CWR,  1,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rCSLOW,    1,NHOR,kNumGPts,gNHOR)

      if (kPop_dyn .ge. 1) call globe_mpwritegp(kFile_Restart,rArea,kMaxSPP,NHOR,kNumGPts,gNHOR)
      if (kPop_dyn .ge. 1) call globe_mpwritegp(kFile_Restart,rAreaBare,1,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,rLIVE, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,rGrowW,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rGrowG,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rGrowT,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,rDie,  kMaxSPP,NHOR,kNumGPts,gNHOR)

      if (kPop_dyn .gt. 2) call globe_mpwritegp(kFile_Restart,fseed,kMaxSPP,NHOR,kNumGPts,gNHOR)
      if (kMig_dyn .gt. 0) call globe_mpwritegp(kFile_Restart,fseedout,kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,dSARAbd,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSAGPP, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSANPP, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACS,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACA,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACL,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACR,  kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACWL, kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACWR, kMaxSPP,NHOR,kNumGPts,gNHOR)

      call globe_mpwritegp(kFile_Restart,dSACTOT,kMaxSPP,NHOR,kNumGPts,gNHOR)
      call globe_mpwritegp(kFile_Restart,dSACVEG,kMaxSPP,NHOR,kNumGPts,gNHOR)

      if (kmypid == kroot) close(kFile_Restart)

      return
      end subroutine jedi_write_restart

!     ******************************************************************
!     JEDI_OUTPUT_SPECIES
!     ******************************************************************
!
!     This subroutine writes out the 6 carbon pools characterizing a species
!     to a file for those species that fulfill the success criterion
!
      subroutine jedi_output_species ()
      use jedi_mod
      implicit none

      integer  iLon, jLat, jGrid, iSPP, i, j, k, ioerr

      __diag(kFile_Diag,'jedi_output_species: start')

!     * gather fields

      __globe_mpga_to_from(gdSARAbd,dSARAbd)
      __globe_mpga_to_from(gdSAGPP,dSAGPP)
      __globe_mpga_to_from(gdSANPP,dSANPP)
      __globe_mpga_to_from(gdSACS,dSACS)
      __globe_mpga_to_from(gdSACA,dSACA)
      __globe_mpga_to_from(gdSACL,dSACL)
      __globe_mpga_to_from(gdSACR,dSACR)
      __globe_mpga_to_from(gdSACWL,dSACWL)
      __globe_mpga_to_from(gdSACWR,dSACWR)

      __globe_mpga_to_from(gdSACTOT,dSACTOT)
      __globe_mpga_to_from(gdSACVEG,dSACVEG)

      if (kPop_dyn .eq. 2) then
        __globe_mpga_to_from(gdSAtau,dSAtau)
      endif

!     * averaging

      if ((kmypid == kroot) .and. (nSACount .gt. 0)) then
        gdSARAbd(:,:)   = gdSARAbd(:,:) / REAL(nSACount)
        gdSAGPP(:,:)    = gdSAGPP(:,:)  / REAL(nSACount)
        gdSANPP(:,:)    = gdSANPP(:,:)  / REAL(nSACount)
        gdSACS(:,:)     = gdSACS(:,:)   / REAL(nSACount)
        gdSACA(:,:)     = gdSACA(:,:)   / REAL(nSACount)
        gdSACL(:,:)     = gdSACL(:,:)   / REAL(nSACount)
        gdSACR(:,:)     = gdSACR(:,:)   / REAL(nSACount)
        gdSACWL(:,:)    = gdSACWL(:,:)  / REAL(nSACount)
        gdSACWR(:,:)    = gdSACWR(:,:)  / REAL(nSACount)
        gdSAtau(:,:)    = gdSAtau(:,:)  / REAL(nSACount)
        gdSACTOT(:,:)   = gdSACTOT(:,:) / REAL(nSACount)
        gdSACVEG(:,:)   = gdSACVEG(:,:) / REAL(nSACount)
      endif

!     * write out successful species
!       format:  lon, lat, sppnum, relabd, gpp, npp, cs, ca, cl, cr, cwl, cwr, (ra)

      if (kmypid == kroot) then
        open(kFile_SPP, FILE=sfile_SPP, STATUS='replace',              &
     &                  FORM='formatted', IOSTAT=ioerr)
        call error_open_new(ioerr,sfile_SPP)
        do k = 1, kNumGPts
          do i = 1, nlon
            do j = 1, nlat
              if (((j - 1) * nlon + i) .eq. gkGPID(k)) then
                iLon = i
                jLat = j
              endif
            enddo
          enddo
          select case(kPop_dyn)
    !     case(0)
    !         do iSPP = FirstSPP, kMaxSPP
    !           if (grCS(k,iSPP) .ge. p14(iSPP) * pA0) then
    !             write(kFile_SPP,'(2I5,I10,9E12.4)') iLon, jLat,        &
    !  &            kFrPl + iSPP - 1, gdSARAbd(k,iSPP), gdSAGPP(k,iSPP), &
    !  &            gdSANPP(k,iSPP), gdSACS(k,iSPP), gdSACA(k,iSPP),     &
    !  &            gdSACL(k,iSPP), gdSACR(k,iSPP), gdSACWL(k,iSPP),     &
    !  &            gdSACWR(k,iSPP)
    !           endif
    !         enddo
            case(1)
            do iSPP = FirstSPP, kMaxSPP
              if (grCA(k,iSPP) + grCL(k,iSPP) + grCR(k,iSPP) + grCWL(k,iSPP) + grCWR(k,iSPP) .gt. pA0) then
                write(kFile_SPP,'(2I5,I10,9E12.4)') iLon, jLat,        &
     &            kFrPl + iSPP - 1, gdSARAbd(k,iSPP), gdSAGPP(k,iSPP), &
     &            gdSANPP(k,iSPP), gdSACS(k,iSPP), gdSACA(k,iSPP),     &
     &            gdSACL(k,iSPP), gdSACR(k,iSPP), gdSACWL(k,iSPP),     &
     &            gdSACWR(k,iSPP)
              endif
            enddo
          case(2)
            do iSPP = FirstSPP, kMaxSPP
              if ( (grArea(k,iSPP) .gt. 0.0) .and. (grCA(k,iSPP)       &
     &            + grCL(k,iSPP) + grCR(k,iSPP) + grCWL(k,iSPP)        &
     &            + grCWR(k,iSPP) .gt. cTiny)) then
                write(kFile_SPP,'(2I5,I10,10E12.4)') iLon, jLat,       &
     &            kFrPl + iSPP - 1, gdSARAbd(k,iSPP), gdSAGPP(k,iSPP), &
     &            gdSANPP(k,iSPP), gdSACS(k,iSPP), gdSACA(k,iSPP),     &
     &            gdSACL(k,iSPP), gdSACR(k,iSPP), gdSACWL(k,iSPP),     &
     &            gdSACWR(k,iSPP), gdSAtau(k,iSPP)
              endif
            enddo
          end select
        enddo
        close(kFile_SPP)
      endif

      __diag(kFile_Diag,'jedi_output_species: end')

      return
      end subroutine jedi_output_species

!     ******************************************************************
!     JEDI_SUCCESS
!     ******************************************************************

!     * writes out parameters of successful species into the "success.data" file and the end of the simulation
!     * writes out extinct species after simulation

      subroutine jedi_success ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, k, iLon, jLat, ioerr
      integer :: kSuccess
      real    :: zp01, zp02, zp13, zp14, zp04, zp12, zp11

!     * gather fields

      __globe_mpga_to_from(grCA,rCA)
      __globe_mpga_to_from(grCS,rCS)
      __globe_mpga_to_from(grCWL,rCWL)
      __globe_mpga_to_from(grCWR,rCWR)
      __globe_mpga_to_from(grCL,rCL)
      __globe_mpga_to_from(grCR,rCR)

      __globe_mpga_to_from(grLive,rLive)

      if ((kPop_dyn .eq. 2) .or. (kPop_dyn .eq. 3)) then
        __globe_mpga_to_from(grArea,rArea)
      endif

!     * write out extinct species

      if (kSpec_dyn .gt. 0) call jedi_extinction

!     * output successful species

      if (kmypid == kroot) then
        open(kFile_SPPSucc, FILE=sfile_SPPSucc, STATUS='replace',      &
     &                      FORM='formatted', IOSTAT=ioerr)
        call error_open_new(ioerr,sfile_SPPSucc)
        __diag(kFile_Diag,'jedi_success: write success file')
        do iSPP = FirstSPP, kMaxSPP
          kSuccess = 0

!         * check if species iSPP is successful anywhere on the grid

          select case(kPop_dyn)
          ! case(0)
          !   do k = 1, kNumGPts
          !     if (grCS(k,iSPP) .ge. p14(iSPP) * pA0) then
          !       kSuccess = 1
          !     endif
          !   enddo

          case(1)
            do k = 1, kNumGPts
              if (grCA(k,iSPP) + grCL(k,iSPP) + grCR(k,iSPP) + grCWL(k,iSPP) + grCWR(k,iSPP) .gt. cTiny) then
                kSuccess = 1
              endif
            enddo


          case(2)
            do k = 1, kNumGPts
              if (grArea(k,iSPP) .gt. cTiny) then
                kSuccess = 1
              endif
            enddo
          end select
!         * if successful then write parameters to file

          if (kSuccess .gt. 0) then
            call jedi_transform_parms_back(iSPP)
            write(kFile_SPPSucc,'(i8,20f6.3)') kSPPID(iSPP),           &
     &        p01(iSPP), p02(iSPP), p03(iSPP), p04(iSPP), p05(iSPP),   &
     &        p06(iSPP), p07(iSPP), p08(iSPP), p09(iSPP), p10(iSPP),   &
     &        p11(iSPP), p12(iSPP), p13(iSPP), p14(iSPP), p15(iSPP),   &
     &        p16(iSPP), p17(iSPP), p18(iSPP), p19(iSPP), p20(iSPP)
            call jedi_transform_parms(iSPP)
          endif
        enddo
        close(kFile_SPPSucc)
      endif

      return
      end subroutine jedi_success

!     ******************************************************************
!     JEDI_OUTPUT_SPT
!     ******************************************************************

      subroutine jedi_output_spt ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, kFile_onepoint

      __diag(kFile_Diag,'jedi_diag: jedi_output_onepoint')

      if (kmypid == kroot) then
        rCtot(:,:) = dSACA(:,:) + dSACL(:,:) + dSACR(:,:) + dSACWL(:,:) + dSACWR(:,:)

        if (nAccuCount .gt. 0) then
          dSARAbd(:,:)  = dSARAbd(:,:)  / REAL(nAccuCount)
          dSAGPP(:,:)   = dSAGPP(:,:)   / REAL(nAccuCount)
          dSACS(:,:)    = dSACS(:,:)    / REAL(nAccuCount)
          dSACA(:,:)    = dSACA(:,:)    / REAL(nAccuCount)
          dSACL(:,:)    = dSACL(:,:)    / REAL(nAccuCount)
          dSACR(:,:)    = dSACR(:,:)    / REAL(nAccuCount)
          dSACWL(:,:)   = dSACWL(:,:)   / REAL(nAccuCount)
          dSACWR(:,:)   = dSACWR(:,:)   / REAL(nAccuCount)
          rCtot(:,:)    = rCtot(:,:)    / REAL(nAccuCount)

          dSAtau(:,:)   = dSAtau(:,:)   / REAL(nAccuCount)

          dSALAI(:,:)   = dSALAI(:,:)   / REAL(nAccuCount)
          dSAWMAX(:,:)  = dSAWMAX(:,:)  / REAL(nAccuCount)
          dSANPP(:,:)   = dSANPP(:,:)   / REAL(nAccuCount)
          dSARES(:,:)   = dSARES(:,:)   / REAL(nAccuCount)
          dSALIT(:,:)   = dSALIT(:,:)   / REAL(nAccuCount)
          dSAFVEG(:,:)  = dSAFVEG(:,:)  / REAL(nAccuCount)
          dSAFH2O(:,:)  = dSAFH2O(:,:)  / REAL(nAccuCount)
          dSAFT(:,:)    = dSAFT(:,:)    / REAL(nAccuCount)

          dSAET(:,:)    = dSAET(:,:)    / REAL(nAccuCount)
          dSASEVAP(:,:) = dSASEVAP(:,:) / REAL(nAccuCount)
          dSALEVAP(:,:) = dSALEVAP(:,:) / REAL(nAccuCount)
          dSABEVAP(:,:) = dSABEVAP(:,:) / REAL(nAccuCount)
          dSATRANS(:,:) = dSATRANS(:,:) / REAL(nAccuCount)
        endif

        do iSPP = FirstSPP, kMaxSPP
          kFile_onepoint = 100 + iSPP
          write(kFile_onepoint,'(2(i5," "),14(e14.7," "))') OUTYEAR, kdiy, &
     &      dSARAbd(:,iSPP), rCtot(:,iSPP), dSAGPP(:,iSPP),             &
     &      dSANPP(:,iSPP), dSALAI(:,iSPP), dSAWMAX(:,iSPP), dSAFVEG(:,iSPP), &
     &      dSACS(:,iSPP), dSACA(:,iSPP), dSACL(:,iSPP), dSACR(:,iSPP), &
     &      dSACWL(:,iSPP), dSACWR(:,iSPP), dSAtau(:,iSPP)

!    &      dSAGPP(:,iSPP), dSANPP(:,iSPP), dSARES(:,iSPP), dSALIT(:,iSPP),      &
!    &      dSAET(:,iSPP), dSASEVAP(:,iSPP), dSALEVAP(:,iSPP), dSABEVAP(:,iSPP), &
!    &      dSATRANS(:,iSPP), dSALAI(:,iSPP),  dSAWMAX(:,iSPP), dSAFVEG(:,iSPP), &
!    &      dSAFH2O(:,iSPP), dSAFT(:,iSPP), dSARAbd(:,iSPP),  dSACS(:,iSPP),     &
!    &      dSACA(:,iSPP), dSACL(:,iSPP), dSACR(:,iSPP), dSACWL(:,iSPP),         &
!    &      dSACWR(:,iSPP), rCtot(:,iSPP)
        enddo
      endif

      return
      end subroutine jedi_output_spt

!     ******************************************************************
!     JEDI_WRITE_EXTINCT
!     ******************************************************************

!     * writes out parameters of EXTINCT species into the "jedi_exspecies.txt" file

      subroutine jedi_write_extinct (iSPP)
      use jedi_mod
      implicit none

      integer :: iSPP

      if (kmypid == kroot) then
        call jedi_transform_parms_back(iSPP)

        write(kFile_SPPExtinct,'(i8,20f6.3)') kSPPID(iSPP),            &
     &    p01(iSPP), p02(iSPP), p05(iSPP), p06(iSPP), p07(iSPP),       &
     &    p08(iSPP), p09(iSPP), p10(iSPP), p13(iSPP), p16(iSPP),       &
     &    p14(iSPP), p15(iSPP), p04(iSPP), p17(iSPP), p12(iSPP),       &
     &    p11(iSPP), p03(iSPP), p18(iSPP), p19(iSPP), p20(iSPP)
      endif

      return
      end subroutine jedi_write_extinct

!     ******************************************************************
!     JEDI_OUTPUT_SPPGRID
!     ******************************************************************
!
!     This subroutine writes out the 6 carbon pools characterizing a species
!     to a file for those species that fulfill the success criterion
!
      subroutine jedi_output_sppgrid ()
      use jedi_mod
      implicit none

      integer iSPP, ioerr
      character(len=80) :: sfile

      __diag(kFile_Diag,'jedi_output_sppgrid: start')

!     * gather fields

      __globe_mpga_to_from(gdSARAbd,dSARAbd)
!      __globe_mpga_to_from(gdSARES,dSARES)
      __globe_mpga_to_from(gdSAGPP,dSAGPP)
      __globe_mpga_to_from(gdSANPP,dSANPP)
!      __globe_mpga_to_from(gdSACS,dSACS)
!      __globe_mpga_to_from(gdSACA,dSACA)
!      __globe_mpga_to_from(gdSACL,dSACL)
!      __globe_mpga_to_from(gdSACR,dSACR)
!      __globe_mpga_to_from(gdSACWL,dSACWL)
!      __globe_mpga_to_from(gdSACWR,dSACWR)
      __globe_mpga_to_from(gdSACTOT,dSACTOT)
      __globe_mpga_to_from(gdSACVEG,dSACVEG)

!     * averaging
      if ((kmypid == kroot) .and. (nSACount .gt. 0)) then
        gdSARAbd(:,:)   = gdSARAbd(:,:) / REAL(nSACount)
!        gdSARES(:,:)    = dSARES(:,:)   / REAL(nAccuCount)
        gdSAGPP(:,:)    = gdSAGPP(:,:)  / REAL(nSACount)
        gdSANPP(:,:)    = gdSANPP(:,:)  / REAL(nSACount)
!        gdSACS(:,:)     = gdSACS(:,:)   / REAL(nSACount)
!        gdSACA(:,:)     = gdSACA(:,:)   / REAL(nSACount)
!        gdSACL(:,:)     = gdSACL(:,:)   / REAL(nSACount)
!        gdSACR(:,:)     = gdSACR(:,:)   / REAL(nSACount)
!        gdSACWL(:,:)    = gdSACWL(:,:)  / REAL(nSACount)
!        gdSACWR(:,:)    = gdSACWR(:,:)  / REAL(nSACount)
        gdSACTOT(:,:)   = gdSACTOT(:,:)  / REAL(nSACount)
        gdSACVEG(:,:)   = gdSACVEG(:,:)  / REAL(nSACount)
!        gdSAtau(:,:)    = gdSAtau(:,:)  / REAL(nSACount)
      endif

!     * write out successful species
!       format:  lon, lat, sppnum, relabd, gpp, npp, cs, ca, cl, cr, cwl, cwr, (ra)

      if (kmypid == kroot) then
        write(sfile,'(A)') sfile_SPPGrid
        if (kspec_yrout == 1) then
          write(sfile,'(A,"_",I4)') sfile_SPPGrid, YEAR
          if (YEAR .lt. 1000) write(sfile,'(A,"_0",I3)') sfile_SPPGrid, YEAR
          if (YEAR .lt. 100) write(sfile,'(A,"_00",I2)') sfile_SPPGrid, YEAR
          if (YEAR .lt. 10) write(sfile,'(A,"_000",I1)') sfile_SPPGrid, YEAR
        endif

        call globe_open_output(TRIM(sfile), kFile_SPPGrid, kFile_Diag)

        __globe_writeoutput(kFile_SPPGrid,gdSARAbd(1:kNumGPts,:),5400)
        __globe_writeoutput(kFile_SPPGrid,gdSAGPP(1:kNumGPts,:),5300)
        __globe_writeoutput(kFile_SPPGrid,gdSANPP(1:kNumGPts,:),5301)
!        __globe_writeoutput(kFile_SPPGrid,gdSACS(1:kNumGPts,:),5372)
!        __globe_writeoutput(kFile_SPPGrid,gdSACA(1:kNumGPts,:),5373)
!        __globe_writeoutput(kFile_SPPGrid,gdSACL(1:kNumGPts,:),5374)
!        __globe_writeoutput(kFile_SPPGrid,gdSACR(1:kNumGPts,:),5375)
!        __globe_writeoutput(kFile_SPPGrid,gdSACWL(1:kNumGPts,:),5376)
!        __globe_writeoutput(kFile_SPPGrid,gdSACWR(1:kNumGPts,:),5377)
        __globe_writeoutput(kFile_SPPGrid,gdSACVEG(1:kNumGPts,:),5304)
        __globe_writeoutput(kFile_SPPGrid,gdSACTOT(1:kNumGPts,:),5399)

        call globe_close_output(kFile_SPPGrid, kFile_Diag)
      endif

      if (kspec_yrout == 1) call jedi_output_species_reset

      __diag(kFile_Diag,'jedi_output_sppgrid: end')

      return
      end subroutine jedi_output_sppgrid

!     ******************************************************************
!     JEDI_OUTPUT_GASPP
!     ******************************************************************

!     * writes out community-weighted trait parameters

      subroutine jedi_output_gaspp ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, i, ioerr
      real    :: zp01,zp02,zp05,zp06,zp07,zp08,zp09,zp10,zp13,zp16
      real    :: zp14,zp15,zp04,zp17,zp12,zp11,zp03,zp18,zp19,zp20
      real    :: zAreaSum, zAL

!     * gather relative abundances

      __globe_mpga_to_from(gdSARAbd,dRAbd)

!     * output community-aggregated trait parameters

      if (kmypid == kroot) then

!     * averaging

      gdSARAbd(:,:) = gdSARAbd(:,:) / REAL(nSACount)

      open(kFile_GASPP, FILE=sfile_GASPP, STATUS='replace',            &
     &                  FORM='formatted', IOSTAT=ioerr)
      call error_open_new(ioerr,sfile_GASPP)
      __diag(kFile_Diag,'jedi_gaspp: write community-aggregated species parameters')

!     * loop over land grid cells

      do i = 1, gNHOR
        zAreaSum = SUM(gdSARAbd(i,FirstSPP:kMaxSPP))  ! sum relative abundances within each cell
        if (zAreaSum .gt. cTiny) then

!         * aggregate traits based on relative abundances

          zp01 = SUM(p01(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp02 = SUM(p02(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp05 = SUM(p05(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp06 = SUM(p06(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp07 = SUM(p07(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp08 = SUM(p08(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp09 = SUM(p09(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp10 = SUM(p10(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp13 = SUM(p13(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp16 = SUM(p16(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp14 = SUM(p14(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp15 = SUM(p15(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp04 = SUM(p04(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp17 = SUM(p17(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp12 = SUM(p12(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp11 = SUM(p11(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp03 = SUM(p03(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp18 = SUM(p18(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp19 = SUM(p19(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum
          zp20 = SUM(p20(FirstSPP:kMaxSPP) * gdSARAbd(i,FirstSPP:kMaxSPP)) / zAreaSum

!         * transform to between zero to one

          zAL  = (zp05 + zp06 + zp07 + zp08)
          zp01 = (LOG10(zp01) + 2.0) / 4.0
          zp02 = (LOG10(zp02) + 2.0) / 4.0
          zp05 = zp05 * zAL
          zp06 = zp06 * zAL
          zp07 = zp07 * zAL
          zp08 = zp08 * zAL
          zp13 = (LOG10(zp13) + 1.0) / 4.0
          zp15 = (LOG10(zp15) + 2.0)
          zp04 = (LOG10(zp04) + 3.0) / 3.0
          zp12 = (LOG10(1.0 / (zp12* (365/12.0)))) / 2.0
          zp11 = 1.0 / (80.0 * 365.0 * zp11)
          zp03 = zp03 / 15.00
        else  ! if no living species, set to default value 0.5
          zp01 = 0.5
          zp02 = 0.5
          zp05 = 0.5
          zp06 = 0.5
          zp07 = 0.5
          zp08 = 0.5
          zp09 = 0.5
          zp10 = 0.5
          zp13 = 0.5
          zp16 = 0.5
          zp14 = 0.5
          zp15 = 0.5
          zp04 = 0.5
          zp17 = 0.5
          zp12 = 0.5
          zp11 = 0.5
          zp03 = 0.5
          zp18 = 0.5
          zp19 = 0.5
          zp20 = 0.5
        endif

!       * write to community-aggregate trait parameter file

        write(kFile_GASPP,'(20f6.3)')                                  &
     &    zp01, zp02, zp03, zp04, zp05, zp06, zp07, zp08, zp09, zp10,  &
     &    zp11, zp12, zp13, zp14, zp15, zp16, zp17, zp18, zp19, zp20
        enddo

        close(kFile_GASPP)
      endif

      return
      end subroutine jedi_output_gaspp
