#define __ACTIVATE
#include "../globe/globe_macros.f90"

!> \file jedi_dyn.F90
!> \brief JEDI 

!     ******************************************************************
!     JEDI_STEP_MODEL_DYN
!     ******************************************************************

      subroutine jedi_step_model_dyn ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP , i

!     ------------------------------------------------------------------
!     * environmental variables
!     ------------------------------------------------------------------

      fLIT(:,:)  = 0.0

!-    * initialize seed pool

      if ((YEAR .lt. 10) .and. (kspinup .eq. 1) .and. (kgermfix .eq. 0)) rCS(:,:) = rCS(:,:) + (pA0 / 365.0)
      rCS(:,:)   = MAX(rCS(:,:), 0.0)

      fGerm(:,:) = 0.0

      zGALOSS_CA(:)     = 0.0
      zGALOSS_CL(:)     = 0.0
      zGALOSS_CR(:)     = 0.0
      zGALOSS_CWL(:)    = 0.0
      zGALOSS_CWR(:)    = 0.0
      zGALOSS_LIT_CS(:) = 0.0

      zLOSS_CA(:)   = 0.0
      zLOSS_CL(:)   = 0.0
      zLOSS_CR(:)   = 0.0
      zLOSS_CWL(:)  = 0.0
      zLOSS_CWR(:)  = 0.0

      zALLOC_CL(:)  = 0.0
      zALLOC_CR(:)  = 0.0
      zALLOC_CWL(:) = 0.0
      zALLOC_CWR(:) = 0.0
      zALLOC_CS(:)  = 0.0

      ! TODO: document this
      zSSG(:)   = cLambda * EXP(17.269 * pjedi_temp(:) / (237.3 + pjedi_temp(:)) ) / (237.3 + pjedi_temp(:))**2
      zSSG(:)   = zSSG(:) / (zSSG(:) + cGamma) / cLambda

      ! TODO: document this and make permanent fix for non-global forcing datasets
      call globe_daylength(kdpy, kdiy, pjedi_lat, zdlen, NHOR)

      ! TODO: check cPAR against
      zPAR(:)   = cPAR * pjedi_swdown(:)

      ! TODO: provide reference to papers for these equations
      zBETA(:)  = 3.822E-6 * (1.0 - (100.0 / pPCO2))                   &
     &          * (1.0 - 0.015 * (MAX(pjedi_temp(:),15.0) - 15.0) * (350.0 / pPCO2))

      Tmax = cTmax350 + (cTmax700 - cTmax350) * ((pPCO2 - 350.0) / (700.0 - 350.0))
      Tmax = 0.5 * ((3.0 * Tmax) - cTmin)

      ! TODO: clean up variable names, double check constant values correspond to Pavlick et al 2013
      zFT(:) = (pjedi_temp(:) - cTmin)**2.0 * (Tmax - pjedi_temp(:))
      zFT(:) = zFT(:) / ((cTref - cTmin)**2.0 * (Tmax - cTref))
      zFT(:) = zFT(:) * (3.0 / (1.0 + (700.0 / pPCO2)))
      where(pjedi_temp(:) .le. cTmin) zFT(:) = 0.0
      zFT(:) = MAX(0.0, zFT(:))

      ! TODO: provide paper and equation number
      zFT2(:)   = pQ10R**(0.10 * (pjedi_temp(:) - 20.0))

!     ------------------------------------------------------------------

!     * update relative abundances

      ! TODO: document kcwt feature
      if (kPop_dyn .eq. 1 .and. kdynamic .eq. 1) then
        if (kcwt == 1) then
          call jedi_relabd_cwt
        else
          call jedi_relabd_dyn
        endif
      endif

      do iSPP = BareSPP, kMaxSPP
        call jedi_dyn_environment(iSPP)
        call jedi_dyn_plants(iSPP)
        call jedi_dyn_accum(iSPP)  ! update accumulated variables
      enddo

!     * soil carbon
      call jedi_dyn_soil

!     TODO: ask kristin about removing kpop_dyn=2 to simplify code
!     * competition
      if (kPop_dyn .eq. 2 .and. (YEAR > 69 .or. kspinup .eq. 0) ) call jedi_dyn_areastep

!     TODO: ask kristin about removing migration to simplify code
!     * migration
      if (jedi_migruns .eq. 1) call jedi_collectSeeds

      ! TODO: harmonize variable names
      dGASOLRAD    = dGASOLRAD(:) + pjedi_swdown(:) 
      dGAFH2O(:)   = dGAFH2O(:)   + pjedi_rain(:) + pjedi_snow(:)

      ! TODO: add better comment describing these counters and various loops
      nAccuCount     = nAccuCount + 1
      nSACount       = nSACount   + 1

      return
      end subroutine jedi_step_model_dyn

!     ******************************************************************
!     JEDI_DYN_environment
!     ******************************************************************

      subroutine jedi_dyn_environment (iSPP)
      use jedi_dyn_mod
      use jedi_mod
      implicit none

      ! TODO: explain or improve these variable names
      integer :: iSPP, i
      integer :: StartHOR, EndHOR

      StartHOR = 1
      EndHOR   = NHOR

      ! TODO: document or remove kcwt feature
      if (kcwt == 1 .and. iSPP .ge. FirstSPP) then
        StartHOR = iSPP - 1
        EndHOR   = iSPP - 1
      endif

!     ----------------------------------------------------------------
!     * calculate land surface parameters
!     ----------------------------------------------------------------
      ! TODO: provide paper/equation number where possible
      do i = StartHOR, EndHOR
        zWMAX(i) = SQRT(pC_WSMAX * rCWR(i,iSPP)) * pjedi_paw(i)
        zWMAX(i) = MIN(cWSUB_MAX - 1.0, MAX(pPAW0, zWMAX(i)))

        if (iSPP .gt. BareSPP) then
          zSLA(i)   = pC_LAI * (p12(iSPP) * 365.0)**(0.46)
          zLAI(i)   = zSLA(i) * rCL(i,iSPP)
          zFVEG(i)  = 1.0 - EXP(-cK * zLAI(i))
          zFFOR(i)  = 1.0 - EXP(- rCWL(i,iSPP) * pC_FFOR)
        else
          zSLA(i)   = 0.0
          zLAI(i)   = 0.0
          zFVEG(i)  = 0.0
          zFFOR(i)  = 0.0
        endif
        pAlb1     = 3.216 * p15(iSPP) + 0.02         ! Hollinger 2010 GCB
        zALB(i)   = pAlb0 * (1.0 - zFVEG(i)) + pAlb1 * zFVEG(i)
        zWLMAX(i) = pC_WLMAX * zLAI(i)
        zFSNOW(i) = rWS(i,iSPP) * 0.066666
        zFSNOW(i) = MIN(1.0, zFSNOW(i))

!       ----------------------------------------------------------------
!       water balance
!       ----------------------------------------------------------------
        ! TODO: document this kludge fix related to jedi_relabd_dyn scaling
        zFW(i)     = MAX(0.0, rW(i,iSPP) - zWMAX(i))
        rW(i,iSPP) = rW(i,iSPP) - zFW(i)
        zFW(i)     = MIN(0.999, rW(i,iSPP) / zWMAX(i)) !fixme

!       * water flux calculations

        zMelt(i)   = MIN(3.22 * MAX(0.0, pjedi_temp(i)), pjedi_snow(i) + rWS(i,iSPP))

!       * interception by canopy

        zThFall(i) = (rWL(i,iSPP) + pjedi_rain(i) + zMelt(i)) - zWLMAX(i)
        zThFall(i) = MAX(0.0, zThFall(i))

!       * calculate surface runoff
        ! TODO: implement slope heterogeneity from ECHAM4 here
        zA(i) = MAX(0.0, (1.0 - zFW(i))**(1.0 / (1.0 + 1.0)))
        zB(i) = zThFall(i) / (1.0 + (1.0 * zWMAX(i)))
        zA(i) = MAX(0.0, (zA(i) - zB(i))**(1.0+1.0))

        zRunoff(i) = MIN(zThFall(i), zThFall(i) - (zWMAX(i) - rW(i,iSPP)) + zWMAX(i) * zA(i))
        zRunoff(i) = MAX(0.0, zRunOff(i))

        zInfil(i)  = zThFall(i) - zRunoff(i)

!       * no infiltration with frozen soil

        if (pjedi_temp(i) .lt. 0.0) then
          zRunoff(i) = zRunoff(i) + zInfil(i)
          zInfil(i)  = 0.0
        endif

!       * soil drainage
        ! TODO: implement drainage based on soil texture from LPJ here

        zDrain(i)    = 0.0
        if (zFW(i) .gt. 0.05) zDrain(i) = MAX(0.0, 0.0242 * zFW(i))

        zA(i)        = 2.42 * ((MAX(0.0, (zFW(i) - 0.90)) * 10.0)**1.5)

        zDrain(i)    = zDrain(i) + zA(i)
        zB(i)        = rW(i,iSPP) + (zInfil(i) - zDrain(i)) - zWMAX(i)
        zDrain(i)    = zDrain(i) + MAX(0.0, zB(i))
        zDrain(i)    = MIN(zDrain(i), rW(i,iSPP) + zInfil(i))
        zDrain(i)    = MAX(0.0, zDrain(i))
        zSubDrain(i) = rWSUB(i,iSPP) - (cWSUB_MAX - zWMAX(i))
        zSubDrain(i) = MAX(0.0, zSubDrain(i))

!       * update state variables

        rWS(i,iSPP)   = rWS(i,iSPP)   + pjedi_snow(i)  - zMelt(i)
        rWL(i,iSPP)   = rWL(i,iSPP)   + pjedi_rain(i)  + zMelt(i) - zThFall(i)
        rW(i,iSPP)    = rW(i,iSPP)    + zInfil(i) - zDrain(i)
        rWSUB(i,iSPP) = rWSUB(i,iSPP) + zDrain(i) - zSubDrain(i)

!       * calculate snow albedo

        zSALBMIN(i) = 0.3 * zFFOR(i) + 0.4 * (1.0 - zFFOR(i))
        zSALBMAX(i) = 0.4 * zFFOR(i) + 0.8 * (1.0 - zFFOR(i))
        zA(i)       = (pjedi_temp(i) + 10.0) * 0.10
        zA(i)       = MIN(1.0, MAX(0.0, zA(i)))
        zSALB(i)    = zSALBMAX(i) - (zSALBMAX(i) - zSALBMIN(i)) * zA(i)
        zSNOWALB(i) = zALB(i) + (zSALB(i) - zALB(i)) * rWS(i, iSPP)    &
     &              / (rWS(i, iSPP) + 10.0)
        zSRNET(i)   = pjedi_swdown(i) * (1.0 - zSNOWALB(i)) + zdlen(i) &
     &              / 24.0 * pjedi_lwnet(i)
        zSRNET(i)   = MAX(0.0, zSRNET(i)) * zFSNOW(i)

!       * snow evaporation

        zPET(i)    = cPETCorr * zSSG(i) * zSRNET(i) * 86400.0
        zPET(i)    = MAX(0.0, zPET(i))
        zSEVAP(i)  = MIN(rWS(i,iSPP), zPET(i))

!       * evapotranspiration from snow free areas

        zRNET(i)   = pjedi_swdown(i) * (1.0 - zALB(i)) + zdlen(i)      &
     &             / 24.0 * pjedi_lwnet(i)
        zRNET(i)   = MAX(0.0, zRNET(i)) * (1.0 - zFSNOW(i))
        zPET(i)    = cPETCorr * zSSG(i) * zRNET(i) * 86400.0
        zPET(i)    = MAX(1.e-6, zPET(i))

!       * column average net radiation

        zRNET(i)   = zRNET(i) + zSRNET(i)

!       * evaporation from skin reservoir of leaves

        zLEVAP(i)  = MIN(rWL(i, iSPP), zPET(i))

!       * evaporation from bare soil

        zWL(i)     = MIN(pPAW0, zWMAX(i))
        zA(i)      = rW(i, iSPP) - (zWMAX(i) - zWL(i))
        if (zA(i) .gt. 0.0) then
          zB(i)    = 0.5 * (1.0 - COS(c2PI / 2 * zA(i) / zWL(i)))
        else
          zB(i)    = 0.0
        endif
        zBEVAP(i)  = zB(i) * (zPET(i) - zLEVAP(i)) * (1.0 - zFVEG(i))
        zBEVAP(i)  = MIN(zBEVAP(i), rW(i, iSPP))

!       ----------------------------------------------------------------
!       plant fluxes: transpiration and carbon uptake
!       ----------------------------------------------------------------

!       * transpiration

        if (rCA(i,iSPP) .gt. 0.0) then

!         * different water stress formulation to ECHAM

          zSET(i)   = pC_TRS * rCR(i,iSPP) * zFW(i)
          zFH2O(i)  = 1.0 - EXP(-zSET(i) / zPET(i))
          zTRANS(i) = zFH2O(i) * (zPET(i) - zLEVAP(i)) * zFVEG(i)
          zTRANS(i) = MIN(zTRANS(i), rW(i,iSPP) - zBEVAP(i))

!         * Canopy N-Amax relationship (Ollinger et al. 2008, PNAS)
!         * Non-rectangular hyperbola from
!         * Rabinowitch (1951), Johnson and Thornley (1984), Franklin (2007)

!         TODO: test Sands 1995 scaling as implemented in https://github.com/mdekauwe/GDAY

          zAMAX(i)  = (0.0059196 * p15(iSPP) + 0.000113844) * zFT(i)
          zGPP(i)   = zBETA(i) * zPAR(i) * zFVEG(i)
          zGPP(i)   = (zGPP(i) + zAMAX(i)) - SQRT( (zGPP(i)+ zAMAX(i))**2.0 - (4.0 * zGPP(i) * zAMAX(i) * pC_GPP))
          zGPP(i)   = ((3600.0 * zdlen(i)) / (2.0 * pC_GPP)) * zGPP(i) * zFH2O(i)

          zRES(i,iSPP)   = cWRFrac * pCN_Wood * (rCWR(i,iSPP) + rCWL(i,iSPP))
          zRES(i,iSPP)   = zRES(i,iSPP) + p15(iSPP) * (rCL(i,iSPP) + rCR(i,iSPP))
          zRES(i,iSPP)   = zRES(i,iSPP) * pN_RES * zFT2(i)
          zRES(i,iSPP)   = zRES(i,iSPP) + pSeedTau  * rCA(i,iSPP)

          zNPP(i)   = zGPP(i) - zRES(i,iSPP)

        else
          zFH2O(i)     = 0.0
          zTRANS(i)    = 0.0
          zGPP(i)      = 0.0
          zRES(i,iSPP) = 0.0
          zNPP(i)      = 0.0
        endif

!       * update soil water budget

        rWS(i, iSPP) = rWS(i, iSPP) - zSEVAP(i)
        rWL(i, iSPP) = rWL(i, iSPP) - zLEVAP(i)
        rW(i, iSPP)  = rW(i, iSPP)  - zBEVAP(i) - zTRANS(i)

!       ----------------------------------------------------------------
!       plant phenology
!       ----------------------------------------------------------------

!       * determine growth conditions

        BareFW(i) = rW(i, BareSPP) / pPAW0

        rGrowW(i, iSPP) = (zFW(i)    + p01(iSPP) * rGrowW(i,iSPP)) / (p01(iSPP) + 1.0)
        rGrowG(i, iSPP) = (BareFW(i) + p01(iSPP) * rGrowG(i,iSPP)) / (p01(iSPP) + 1.0)
        rGrowT(i, iSPP) = (pjedi_temp(i)    + p02(iSPP) * rGrowT(i,iSPP)) / (p02(iSPP) + 1.0)

        ff_GrowW(i) = 0.0
        ff_GrowG(i) = 0.0
        ff_GrowT(i) = 0.0

        if (rGrowW(i,iSPP) .gt. 0.5)       ff_GrowW(i)  = 1.0
        if (rGrowG(i,iSPP) .gt. 0.5)       ff_GrowG(i)  = 1.0
        if (rGrowT(i,iSPP) .gt. p03(iSPP)) ff_GrowT(i) = 1.0

        ff_Grow(i)  = ff_GrowT(i) * ff_GrowW(i)
        ff_GrowG(i) = ff_GrowG(i) * ff_GrowT(i)

!       ----------------------------------------------------------------
!-      germination
!       ----------------------------------------------------------------


!-      * add virtual seeds to seed pool and GPP flux

        if (OUTYEAR .gt. 2 .and. kGermFix .eq. 1) then
!         * set min. value for rCS, because the automated setting to zero
!           for small values by compiler depends on processor, compiler,
!           and the use of the real8 switch
          if (rCA(i,iSPP) .le. r8_tiny .and. rCS(i,iSPP) .le. 1.0E-8) then
            rCS(i,iSPP) = pA0 * ff_GrowG(i)
!           dGAGPP(i) = dGAGPP(i) + pA0 * ff_GrowG(i)
          endif
        endif

!-      * calculate germination flux
        fGerm(i,iSPP) = rCS(i,iSPP) * p04(iSPP) * ff_GrowG(i)

!-      * add germination flux to assimilate pool
!-      * seed decay goes to litter pool
        if (dRAbd(i,iSPP) .gt. pSeedFix) then
          rCA(i,iSPP) = rCA(i,iSPP) + (fGerm(i,iSPP) / dRAbd(i,iSPP))
          if (YEAR .gt. 150 .or. kspinup .eq. 0)                       &
     &      zGALOSS_LIT_CS(i) = zGALOSS_LIT_CS(i) + rCS(i,iSPP) * pSeedTau
        else
          rCA(i,iSPP) = rCA(i,iSPP) + (fGerm(i,iSPP) / pSeedFix)
          if (YEAR .gt. 150 .or. kspinup .eq. 0)                       &
     &      zGALOSS_LIT_CS(i) = zGALOSS_LIT_CS(i) + (rCS(i,iSPP) * pSeedTau) &
     &                        + (fGerm(i,iSPP) * (1.0 - (dRAbd(i,iSPP) / pSeedFix)))
        endif

!-      * subtract germination and seed decay from seed pool
        rCS(i,iSPP) = (rCS(i,iSPP) * (1.0 - pSeedTau)) - fGerm(i,iSPP)

      enddo

      return
      end subroutine jedi_dyn_environment

!     ******************************************************************
!     JEDI_DYN_PLANTS
!     ******************************************************************

      subroutine jedi_dyn_plants (iSPP)
      use jedi_dyn_mod
      use jedi_mod
      implicit none

      integer :: iSPP, i
      integer :: StartHOR, EndHOR

      StartHOR = 1
      EndHOR   = NHOR

      if (kcwt == 1 .and. iSPP .ge. FirstSPP) then
        StartHOR = iSPP - 1
        EndHOR   = iSPP - 1
      endif

      do i = StartHOR, EndHOR

!       ----------------------------------------------------------------
!       loop over all living plants
!       ----------------------------------------------------------------

        if (rCA(i,iSPP) .gt. 0.0) then
          rLive(i,iSPP) = 1.0      ! Species population is alive

!         * seed allocation

          ff_Seed(i) = 0.0
          if (zNPP(i) .gt. 0.0)  ff_Seed(i) = 1.0

!         * marginal growth/plasticity -- not implemented yet

          ff_GrowL(i) = zFVEG(i)
          ff_GrowR(i) = 1.0

!-        * allocation parameters

          ff_AS(i) = p05(iSPP) * ff_Seed(i) * ff_Grow(i)
          ff_AL(i) = p06(iSPP) * ff_Grow(i)
          ff_AR(i) = p07(iSPP) * ff_Grow(i)
          ff_LW(i) = p09(iSPP) * ff_GrowL(i)
          ff_RW(i) = p10(iSPP) * ff_GrowL(i)
          ff_A(i)  = ff_AS(i) + ff_AR(i) + ff_AL(i)

!-        * allocation from assimilate pool

          rCA(i,iSPP) = rCA(i,iSPP) + zNPP(i)
          zALLOC(i) = MAX(0.0, rCA(i,iSPP))

!-        * allocation to plant tissue pools

          zALLOC_CL(i)   = ff_AL(i) * zALLOC(i) * (1.0 - ff_LW(i)) * cGResLeaf
          zALLOC_CR(i)   = ff_AR(i) * zALLOC(i) * (1.0 - ff_RW(i)) * cGResLeaf
          zALLOC_CWL(i)  = ff_AL(i) * zALLOC(i) * ff_LW(i) * cGResWood
          zALLOC_CWR(i)  = ff_AR(i) * zALLOC(i) * ff_RW(i) * cGResRoot
          zALLOC_CS(i)   = ff_AS(i) * zALLOC(i) * cGResSeed

          rCA(i,iSPP)  = rCA(i, iSPP) - (zALLOC(i) * ff_A(i))

!-        * add growth respiration to autotrophic respiration

          zRES(i,iSPP) = zRES(i,iSPP)                                  &
     &            + (zALLOC_CL(i)  / cGResLeaf) * (1.0 - cGResLeaf)    &
     &            + (zALLOC_CR(i)  / cGResLeaf) * (1.0 - cGResLeaf)    &
     &            + (zALLOC_CWL(i) / cGResWood) * (1.0 - cGResWood)    &
     &            + (zALLOC_CWR(i) / cGResRoot) * (1.0 - cGResRoot)    &
     &            + (zALLOC_CS(i)  / cGResSeed) * (1.0 - cGResSeed)

          zNPP(i) = zGPP(i) - zRES(i,iSPP)

!-        * senescence

          ff_Die(i)    = 0.0
          rDie(i,iSPP) = (zNPP(i) + p13(iSPP) * rDie(i, iSPP)) / (p13(iSPP) + 1.0)
          if (rDie(i,iSPP) .lt. 0.0 .and. zNPP(i) .lt. 0.0) ff_Die(i) = 1.0 - ff_Grow(i)

!-        * turnover

          zLOSS_CL(i)  = rCL(i,iSPP) * (p12(iSPP) + 0.10 * ff_Die(i) * p14(iSPP))
          zLOSS_CR(i)  = rCR(i,iSPP) * (p12(iSPP) + 0.10 * ff_Die(i) * (1.0 - p14(iSPP)))
          zLOSS_CWL(i) = p11(iSPP) * rCWL(i,iSPP)
          zLOSS_CWR(i) = p11(iSPP) * rCWR(i,iSPP)

!-        * update plant carbon pools

          rCL(i,iSPP)  = MAX(0.0, rCL(i,iSPP) + zALLOC_CL(i) - zLOSS_CL(i))
          rCR(i,iSPP)  = MAX(0.0, rCR(i,iSPP)  + zALLOC_CR(i)  - zLOSS_CR(i))
          rCWL(i,iSPP) = MAX(0.0, rCWL(i,iSPP) + zALLOC_CWL(i) - zLOSS_CWL(i))
          rCWR(i,iSPP) = MAX(0.0, rCWR(i,iSPP) + zALLOC_CWR(i) - zLOSS_CWR(i))
          rCS(i,iSPP)  = rCS(i,iSPP)  + (zALLOC_CS(i) * dRAbd(i,iSPP))

!-        * death

          if (rCA(i,iSPP) .le. cTiny) then
            zLOSS_CL(i)    = zLOSS_CL(i)  + rCL(i,iSPP)
            zLOSS_CR(i)    = zLOSS_CR(i)  + rCR(i,iSPP)
            zLOSS_CWL(i)   = zLOSS_CWL(i) + rCWL(i,iSPP)
            zLOSS_CWR(i)   = zLOSS_CWR(i) + rCWR(i,iSPP)
            rCA(i,iSPP)    = 0.0
            rCL(i,iSPP)    = 0.0
            rCR(i,iSPP)    = 0.0
            rCWL(i,iSPP)   = 0.0
            rCWR(i,iSPP)   = 0.0
            rDie(i,iSPP)   = 0.0
            rLive(i,iSPP)  = 0.0
            rArea(i,iSPP)  = 0.0
            if (kPop_dyn .eq. 2) then
              rAreaBare(i)  = rAreaBare(i) + rArea(i,iSPP)
              rArea(i,iSPP) = 0.0
              rLive(i,iSPP) = 0.0                              ! tag for dead population
              if (rCS(i,iSPP) .le. 0.0)  rLive(i,iSPP) = -1.0  ! tag for locally extinct population
              if (kSpec_dyn .gt. 0) then
                if ((rLive(i,iSPP) .eq. -1.0) .and. (fLIT(i,iSPP) .gt. 0.0)) Ext(i,iSPP) = 1.0
              endif
            endif
          endif

!         * transfer water from sublayer to rooting zone layer because of root growth
!         * or vice versa because of root loss

          zDW(i)      = SQRT(pC_WSMAX * rCWR(i,iSPP)) * pjedi_paw(i)
          zDW(i)      = MIN(cWSUB_MAX - 1.0, MAX(pPAW0, zDW(i))) - zWMAX(i)

          zA(i)       = MAX(0.0, zDW(i)) * rWSUB(i,iSPP) / (cWSUB_MAX - zWMAX(i))
          zB(i)       = MIN(0.0, zDW(i)) * rW(i,iSPP) / zWMAX(i)

          rW(i,iSPP)    = rW(i,iSPP)    + zA(i) + zB(i)
          rWSUB(i,iSPP) = rWSUB(i,iSPP) - zA(i) - zB(i)
        endif
      enddo

      return
      end subroutine jedi_dyn_plants

!     ******************************************************************
!     JEDI_RELABD
!     ******************************************************************
      ! TODO: need better comments

      subroutine jedi_relabd_dyn ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP , i
      integer :: StartHOR, EndHOR

      StartHOR = 1
      EndHOR   = NHOR

      zgCtot(:)       = 0.0
      dRAbd(:,:)      = 0.0
      zDeltaRAbd(:,:) = 0.0

      rCtot(:,:) = rCWL(:,:) + rCL(:,:) + rCR(:,:) + rCWR(:,:) + rCA(:,:)

!     * calculate equilibrium relative abundances

      zgCtot(:)  = SUM(rCtot(:,:),dim=2)
      do iSPP = FirstSPP, kMaxSPP
        where (zgCtot(:) .gt. cTiny) !fixme cTiny
          dRAbd(:,iSPP) = rCtot(:,iSPP) / zgCtot(:)
        end where
      enddo

!-    * calculate timestep change in relative abundances

      zDeltaRAbd(:,:) = (dRAbd(:,:) - rArea(:,:)) * pRAbdTau
      dGACOL(:) = dGACOL(:) + SUM(ABS(zDeltaRAbd(:,:)), DIM=2)

      if (kCbal .eq. 1) then


        where ((rArea(:,:) + zDeltaRAbd(:,:)) .gt. rArea(:,:)) !fixme cTiny
          zDeltaCtot(:,:) = rArea(:,:) / (rArea(:,:) + zDeltaRAbd(:,:))
        end where

        zDeltaCtot(:,:) = MAX(0.0,zDeltaCtot(:,:))
        zDeltaCtot(:,:) = MIN(1.0,zDeltaCtot(:,:))
       ! rCA(:,:)        = rCA(:,:)  * zDeltaCtot(:,:)
        rCL(:,:)        = MAX(0.0, rCL(:,:)  * zDeltaCtot(:,:))
        rCR(:,:)        = MAX(0.0, rCR(:,:)  * zDeltaCtot(:,:))
        rCWL(:,:)       = MAX(0.0, rCWL(:,:) * zDeltaCtot(:,:))
        rCWR(:,:)       = MAX(0.0, rCWR(:,:) * zDeltaCtot(:,:))

!-      * virtual carbon goes to litter

        zDeltaCtot(:,:) = MAX(0.0, -zDeltaRAbd(:,:))

        zGALOSS_CA(:)  = zGALOSS_CA(:)  + SUM(rCA(:,:)  * zDeltaCtot(:,:), dim=2)
        zGALOSS_CL(:)  = zGALOSS_CL(:)  + SUM(rCL(:,:)  * zDeltaCtot(:,:), dim=2)
        zGALOSS_CR(:)  = zGALOSS_CR(:)  + SUM(rCR(:,:)  * zDeltaCtot(:,:), dim=2)
        zGALOSS_CWL(:) = zGALOSS_CWL(:) + SUM(rCWL(:,:) * zDeltaCtot(:,:), dim=2)
        zGALOSS_CWR(:) = zGALOSS_CWR(:) + SUM(rCWR(:,:) * zDeltaCtot(:,:), dim=2)

      endif

!-    * update area
      rArea(:,:) = rArea(:,:) + zDeltaRAbd(:,:)
      rArea(:,BareSPP) = 1.0 - SUM(rArea(:,FirstSPP:kMaxSPP),dim=2)
      rArea(:,:) = MAX(0.0,rArea(:,:))
      rArea(:,:) = MIN(1.0,rArea(:,:))
      dRAbd(:,:) = rArea(:,:)

      return
      end subroutine jedi_relabd_dyn

!     ******************************************************************
!     JEDI_RELABD
!     ******************************************************************

      subroutine jedi_relabd_cwt ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP

      rCtot(:,:) = rCWL(:,:) + rCL(:,:) + rCR(:,:) + rCWR(:,:) + rCA(:,:)

      if (kmypid == kroot) then
        grArea(:,:) = 0.0
        iSPP = FirstSPP
        do iSPP = FirstSPP, kMaxSPP
          grArea(iSPP-1,iSPP) = 1.0
        enddo
      endif

      __globe_mpsc_from_to(grArea,rArea)
      dRAbd(:,:) = rArea(:,:)

      return
      end subroutine jedi_relabd_cwt

!     ******************************************************************
!     JEDI_DYN_SOIL
!     ******************************************************************

      subroutine jedi_dyn_soil ()
      use jedi_dyn_mod
      use jedi_mod
      implicit none

      dGACALOSS(:)   = dGACALOSS(:)   +  zGALOSS_CA(:)
      dGACLLOSS(:)   = dGACLLOSS(:)   +  zGALOSS_CL(:)
      dGACRLOSS(:)   = dGACRLOSS(:)   +  zGALOSS_CR(:)
      dGACWLLOSS(:)  = dGACWLLOSS(:)  +  zGALOSS_CWL(:)
      dGACWRLOSS(:)  = dGACWRLOSS(:)  +  zGALOSS_CWR(:)

      zFT3(:)    = pQ10H**(0.1 * pjedi_temp(:))

      pLitterTau = 2.00
      pCWDTau    = 30.0
      pSoilTau   = 125.0

      zGALOSS_LIT_CL(:)    = rLIT_CL(:)  * zFT3(:) * (1.0 / (pLitterTau * cDaysPerYear))
      zGALOSS_LIT_CR(:)    = rLIT_CR(:)  * zFT3(:) * (1.0 / (pLitterTau * cDaysPerYear))
      zGALOSS_LIT_CWL(:)   = rLIT_CWL(:) * zFT3(:) * (1.0 / (pCWDTau    * cDaysPerYear))
      zGALOSS_LIT_CWR(:)   = rLIT_CWR(:) * zFT3(:) * (1.0 / (pCWDTau    * cDaysPerYear))
      zGALOSS_LIT_CSLOW(:) = rCSLOW(:)   * zFT3(:) * (1.0 / (pSoilTau   * cDaysPerYear))

      dGALOSS_LIT_CS(:)    = dGALOSS_LIT_CS(:)    + zGALOSS_LIT_CS(:)
      dGALOSS_LIT_CL(:)    = dGALOSS_LIT_CL(:)    + zGALOSS_LIT_CL(:)
      dGALOSS_LIT_CR(:)    = dGALOSS_LIT_CR(:)    + zGALOSS_LIT_CR(:)
      dGALOSS_LIT_CWL(:)   = dGALOSS_LIT_CWL(:)   + zGALOSS_LIT_CWL(:)
      dGALOSS_LIT_CWR(:)   = dGALOSS_LIT_CWR(:)   + zGALOSS_LIT_CWR(:)
      dGALOSS_LIT_CSLOW(:) = dGALOSS_LIT_CSLOW(:) + zGALOSS_LIT_CSLOW(:)

!-    * calculate heterotrophic flux to atmosphere

      dGARESH(:) = dGARESH(:) + pLitterFrac2Atm * (zGALOSS_LIT_CL(:)  + zGALOSS_LIT_CR(:))
      dGARESH(:) = dGARESH(:) + pCWDFrac2Atm    * (zGALOSS_LIT_CWL(:) + zGALOSS_LIT_CWR(:))
      dGARESH(:) = dGARESH(:) + pLitterFrac2Atm * zGALOSS_LIT_CS(:)
      dGARESH(:) = dGARESH(:) + zGALOSS_LIT_CSLOW(:)

!     * update litter and soil pools

      rLIT_CL(:)   = rLIT_CL(:)  + zGALOSS_CL(:)  + zGALOSS_CA(:) - zGALOSS_LIT_CL(:)
      rLIT_CR(:)   = rLIT_CR(:)  + zGALOSS_CR(:)  - zGALOSS_LIT_CR(:)
      rLIT_CWL(:)  = rLIT_CWL(:) + zGALOSS_CWL(:) - zGALOSS_LIT_CWL(:)
      rLIT_CWR(:)  = rLIT_CWR(:) + zGALOSS_CWR(:) - zGALOSS_LIT_CWR(:)

      rCSLOW(:) = rCSLOW(:) + (1.0 - pLitterFrac2Atm) * (zGALOSS_LIT_CL(:)  + zGALOSS_LIT_CR(:))
      rCSLOW(:) = rCSLOW(:) + (1.0 - pCWDFrac2Atm)    * (zGALOSS_LIT_CWL(:) + zGALOSS_LIT_CWR(:))
      rCSLOW(:) = rCSLOW(:) + (1.0 - pLitterFrac2Atm) * zGALOSS_LIT_CS(:)
      rCSLOW(:) = rCSLOW(:) - zGALOSS_LIT_CSLOW(:)

      return
      end subroutine jedi_dyn_soil

!     ******************************************************************
!     JEDI_DYN_ACCUM
!     ******************************************************************

      subroutine jedi_dyn_accum (iSPP)
      use jedi_dyn_mod
      use jedi_mod
      implicit none

      integer :: iSPP, i
      integer :: StartHOR, EndHOR

      StartHOR = 1
      EndHOR   = NHOR

      if (kcwt == 1 .and. iSPP .ge. FirstSPP) then
        StartHOR = iSPP - 1
        EndHOR   = iSPP - 1
      endif

      do i = StartHOR, EndHOR

!       * species-based accumulation for one point (distributed)

        if (nonepoint .eq. 1) then
          dSALAI(i,iSPP)   = dSALAI(i,iSPP)   + zLAI(i)
          dSAWMAX(i,iSPP)  = dSAWMAX(i,iSPP)  + zWMAX(i)
          dSARES(i,iSPP)   = dSARES(i,iSPP)   + zRES(i,iSPP)
          dSALIT(i,iSPP)   = dSALIT(i,iSPP)   + fLIT(i,iSPP)
          dSAFVEG(i,iSPP)  = dSAFVEG(i,iSPP)  + zFVEG(i)
          dSAFH2O(i,iSPP)  = dSAFH2O(i,iSPP)  + zFH2O(i)
          dSAFT(i,iSPP)    = dSAFT(i,iSPP)    + zFT(i)
          dSAET(i,iSPP)    = dSAET(i,iSPP)    + zSEVAP(i) + zLEVAP(i) + zBEVAP(i) + zTRANS(i)
          dSASEVAP(i,iSPP) = dSASEVAP(i,iSPP) + zSEVAP(i)
          dSALEVAP(i,iSPP) = dSALEVAP(i,iSPP) + zLEVAP(i)
          dSABEVAP(i,iSPP) = dSABEVAP(i,iSPP) + zBEVAP(i)
          dSATRANS(i,iSPP) = dSATRANS(i,iSPP) + zTRANS(i)
        endif

!       * species-based accumulation (distributed) (rCtot is *not* weighted by RA so as similar
!         as possible to the value used by the richness in jedi.F90)

        dSARAbd(i,iSPP) = dSARAbd(i,iSPP) + dRAbd(i,iSPP)
        dSAGPP(i,iSPP)  = dSAGPP(i,iSPP)  + dRAbd(i,iSPP) * zGPP(i)
        dSANPP(i,iSPP)  = dSANPP(i,iSPP)  + dRAbd(i,iSPP) * zNPP(i)
!        dSACS(i,iSPP)   = dSACS(i,iSPP)   + dRAbd(i,iSPP) * rCS(i,iSPP)
!        dSACA(i,iSPP)   = dSACA(i,iSPP)   + dRAbd(i,iSPP) * rCA(i,iSPP)
!        dSACL(i,iSPP)   = dSACL(i,iSPP)   + dRAbd(i,iSPP) * rCL(i,iSPP)
!        dSACR(i,iSPP)   = dSACR(i,iSPP)   + dRAbd(i,iSPP) * rCR(i,iSPP)
!        dSACWL(i,iSPP)  = dSACWL(i,iSPP)  + dRAbd(i,iSPP) * rCWL(i,iSPP)
!        dSACWR(i,iSPP)  = dSACWR(i,iSPP)  + dRAbd(i,iSPP) * rCWR(i,iSPP)
        dSACTOT(i,iSPP) = dSACTOT(i,iSPP) + rCtot(i,iSPP)
        dSACVEG(i,iSPP) = dSACVEG(i,iSPP) + dRAbd(i,iSPP) * (rCA(i,iSPP) + rCWL(i,iSPP) + rCWR(i,iSPP) + rCL(i,iSPP) + rCR(i,iSPP))

!       * grid-based accumulation
!       using relative abundances of previous time step
!       these are reset monthly in jedi_output_reset

        if(iSPP /= 1) then ! don't accumulate bare earth
        dGARAbd(i)     = dGARAbd(i)     + dRAbd(i,iSPP)
        endif
        dGAGPP(i)      = dGAGPP(i)      + dRAbd(i,iSPP) * zGPP(i)
        dGANPP(i)      = dGANPP(i)      + dRAbd(i,iSPP) * zNPP(i)
        dGARES(i)      = dGARES(i)      + dRAbd(i,iSPP) * zRES(i,iSPP)

        dGACA(i)   = dGACA(i)   + rCA(i,iSPP)
        dGACL(i)   = dGACL(i)   + rCL(i,iSPP)
        dGACR(i)   = dGACR(i)   + rCR(i,iSPP)
        dGACWL(i)  = dGACWL(i)  + rCWL(i,iSPP)
        dGACWR(i)  = dGACWR(i)  + rCWR(i,iSPP)

        dGACVEG(i)   = dGACVEG(i) + dRAbd(i,iSPP) * (rCA(i,iSPP) + rCWL(i,iSPP) + rCWR(i,iSPP) + rCL(i,iSPP) + rCR(i,iSPP))

        dGACLALLOC(i)  = dGACLALLOC(i)  + dRAbd(i,iSPP) * zALLOC_CL(i)
        dGACRALLOC(i)  = dGACRALLOC(i)  + dRAbd(i,iSPP) * zALLOC_CR(i)
        dGACWLALLOC(i) = dGACWLALLOC(i) + dRAbd(i,iSPP) * zALLOC_CWL(i)
        dGACWRALLOC(i) = dGACWRALLOC(i) + dRAbd(i,iSPP) * zALLOC_CWR(i)
        dGACSALLOC(i)  = dGACSALLOC(i)  + dRAbd(i,iSPP) * zALLOC_CS(i)

!       zGALOSS_CA(i)    = zGALOSS_CA(i)   + dRAbd(i,iSPP) * zLOSS_CA(i)
        zGALOSS_CL(i)    = zGALOSS_CL(i)   + dRAbd(i,iSPP) * zLOSS_CL(i)
        zGALOSS_CR(i)    = zGALOSS_CR(i)   + dRAbd(i,iSPP) * zLOSS_CR(i)
        zGALOSS_CWL(i)   = zGALOSS_CWL(i)  + dRAbd(i,iSPP) * zLOSS_CWL(i)
        zGALOSS_CWR(i)   = zGALOSS_CWR(i)  + dRAbd(i,iSPP) * zLOSS_CWR(i)

        dGAALB(i)      = dGAALB(i)      + zSLA(i) * dRAbd(i,iSPP)

        dGALAI(i)      = dGALAI(i)      + zLAI(i)  * dRAbd(i,iSPP)
        dGAFVEG(i)     = dGAFVEG(i)     + p15(iSPP) * dRAbd(i,iSPP)
        dGAWMAX(i)     = dGAWMAX(i)     + zWMAX(i) * dRAbd(i,iSPP)
        dGAQTOT(i)     = dGAQTOT(i)     + (zRunoff(i) + zSubDrain(i)) * dRAbd(i,iSPP)
        dGAFT(i)       = dGAFT(i)       + dRAbd(i,iSPP) * zRNET(i)

        dGAET(i)       = dGAET(i)       + dRAbd(i,iSPP) * (zSEVAP(i)   &
     &                 + zLEVAP(i) + zBEVAP(i) + zTRANS(i))


        dGATRANS(i)     = dGATRANS(i)     + dRAbd(i,iSPP) * zTRANS(i)
        dGAQSURF(i)     = dGAQSURF(i)     + dRAbd(i,iSPP) * zRunoff(i)
        dGAESOIL(i)     = dGAESOIL(i)     + dRAbd(i,iSPP) * (zSEVAP(i) + zBEVAP(i))
        dGALEVAP(i)     = dGALEVAP(i)     + dRAbd(i,iSPP) * zLEVAP(i)


!       * optimizing version

        if (kOpti == 1) dOptiGoalAcc(i,iSPP) = dOptiGoalAcc(i,iSPP) + zGPP(i)

      enddo

      return
      end subroutine jedi_dyn_accum
