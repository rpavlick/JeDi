! #define __ACTIVATE
! #include "../globe/globe_macros.f90"
!
! !     ******************************************************************
! !     JEDI_STEP_MODEL_IND
! !     ******************************************************************
! 
!       subroutine jedi_step_model_ind ()
!       use jedi_mod
!       implicit none
!
!       integer :: iSPP , i, j, k, n, zS, zT, iT
!
!       REAL,allocatable,dimension(:) :: zWMAX
!       REAL,allocatable,dimension(:) :: zALB
!       REAL,allocatable,dimension(:) :: zFVEG
!       REAL,allocatable,dimension(:) :: zLAI
!       REAL,allocatable,dimension(:) :: zFFOR
!       REAL,allocatable,dimension(:) :: zFSNOW
!       REAL,allocatable,dimension(:) :: zWLMAX
!       REAL,allocatable,dimension(:) :: zFW
!       REAL,allocatable,dimension(:) :: zFT2
!       REAL,allocatable,dimension(:) :: zFT3
!       REAL,allocatable,dimension(:) :: zFT
!       REAL,allocatable,dimension(:) :: zA
!       REAL,allocatable,dimension(:) :: zB
!       REAL,allocatable,dimension(:) :: zThFall
!       REAL,allocatable,dimension(:) :: zTRANS
!       REAL,allocatable,dimension(:) :: zRunoff
!       REAL,allocatable,dimension(:) :: zInfil
!       REAL,allocatable,dimension(:) :: zMelt
!       REAL,allocatable,dimension(:) :: zDrain
!       REAL,allocatable,dimension(:) :: zSET
!       REAL,allocatable,dimension(:) :: zFH2O
!       REAL,allocatable,dimension(:) :: zPAR
!       REAL,allocatable,dimension(:) :: zGPP
!       REAL,allocatable,dimension(:) :: zNPP
!       REAL,allocatable,dimension(:) :: zRES
!       REAL,allocatable,dimension(:) :: zLIT
!       REAL,allocatable,dimension(:) :: ff_GrowT, ff_GrowW
!       REAL,allocatable,dimension(:) :: ff_Grow
!       REAL,allocatable,dimension(:) :: ff_Die
!       REAL,allocatable,dimension(:) :: ff_Seed
!       REAL,allocatable,dimension(:) :: ff_A
!       REAL,allocatable,dimension(:) :: ff_AL
!       REAL,allocatable,dimension(:) :: ff_AS
!       REAL,allocatable,dimension(:) :: ff_AR
!       REAL,allocatable,dimension(:) :: ff_GrowL
!       REAL,allocatable,dimension(:) :: ff_GrowR
!       REAL,allocatable,dimension(:) :: ff_LW
!       REAL,allocatable,dimension(:) :: ff_RW
!       REAL,allocatable,dimension(:) :: ff_LD
!       REAL,allocatable,dimension(:) :: ff_RD
!       REAL,allocatable,dimension(:) :: ZALLOC
!       REAL,allocatable,dimension(:) :: fff_
!       REAL,allocatable,dimension(:) :: ffx_
!       REAL,allocatable,dimension(:) :: zDW
!
!       REAL,allocatable,dimension(:) :: zPET
!       REAL,allocatable,dimension(:) :: zLEVAP
!       REAL,allocatable,dimension(:) :: zRAIN
!       REAL,allocatable,dimension(:) :: zBEVAP
!       REAL,allocatable,dimension(:) :: zSSG, zzSSG
!       REAL,allocatable,dimension(:) :: zSALB
!       REAL,allocatable,dimension(:) :: zSALBMAX
!       REAL,allocatable,dimension(:) :: zSALBMIN
!
!       REAL,allocatable,dimension(:) :: zSNOW
!       REAL,allocatable,dimension(:) :: zSNOWALB
!       REAL,allocatable,dimension(:) :: zRNET
!       REAL,allocatable,dimension(:) :: zSEVAP
!       REAL,allocatable,dimension(:) :: zWL
!       REAL,allocatable,dimension(:) :: zSubDrain
!
!       REAL,allocatable,dimension(:) :: zdt           ! correction : dt - 273.15 kelvin
!       REAL,allocatable,dimension(:) :: zdlen         ! daylen
!
!       REAL,allocatable,dimension(:) :: zLifePlants, zSuccessPlants
!
! !     * allocate temporary fields
!
!       __allocate(zWMAX,(NHOR))
!       __allocate(zALB,(NHOR))
!       __allocate(zFVEG,(NHOR))
!       __allocate(zLAI,(NHOR))
!       __allocate(zFFOR,(NHOR))
!       __allocate(zFSNOW,(NHOR))
!       __allocate(zWLMAX,(NHOR))
!       __allocate(zFW,(NHOR))
!       __allocate(zFT2,(NHOR))
!       __allocate(zFT3,(NHOR))
!       __allocate(zFT,(NHOR))
!       __allocate(zA,(NHOR))
!       __allocate(zB,(NHOR))
!       __allocate(zThFall,(NHOR))
!       __allocate(zTRANS,(NHOR))
!       __allocate(zRunoff,(NHOR))
!       __allocate(zInfil,(NHOR))
!       __allocate(zMelt,(NHOR))
!       __allocate(zDrain,(NHOR))
!       __allocate(zSET,(NHOR))
!       __allocate(zFH2O,(NHOR))
!       __allocate(zPAR,(NHOR))
!       __allocate(zGPP,(NHOR))
!       __allocate(zNPP,(NHOR))
!       __allocate(zRES,(NHOR))
!       __allocate(zLIT,(NHOR))
!       __allocate(ff_GrowT,(NHOR))
!       __allocate(ff_GrowW,(NHOR))
!       __allocate(ff_Grow,(NHOR))
!       __allocate(ff_Die,(NHOR))
!       __allocate(ff_Seed,(NHOR))
!       __allocate(ff_A,(NHOR))
!       __allocate(ff_AL,(NHOR))
!       __allocate(ff_AS,(NHOR))
!       __allocate(ff_AR,(NHOR))
!       __allocate(ff_GrowL,(NHOR))
!       __allocate(ff_GrowR,(NHOR))
!       __allocate(ff_LW,(NHOR))
!       __allocate(ff_RW,(NHOR))
!       __allocate(ff_LD,(NHOR))
!       __allocate(ff_RD,(NHOR))
!       __allocate(ZALLOC,(NHOR))
!       __allocate(fff_,(NHOR))
!       __allocate(ffx_,(NHOR))
!       __allocate(zDW,(NHOR))
!
!       __allocate(zPET,(NHOR))
!       __allocate(zLEVAP,(NHOR))
!       __allocate(zRAIN,(NHOR))
!       __allocate(zBEVAP,(NHOR))
!       __allocate(zSSG,(NHOR))
!       __allocate(zzSSG,(NHOR))
!       __allocate(zSALB,(NHOR))
!       __allocate(zSALBMAX,(NHOR))
!       __allocate(zSALBMIN,(NHOR))
!
!       __allocate(zSNOW,(NHOR))
!       __allocate(zSNOWALB,(NHOR))
!       __allocate(zRNET,(NHOR))
!       __allocate(zSEVAP,(NHOR))
!       __allocate(zWL,(NHOR))
!       __allocate(zSubDrain,(NHOR))
!
!       __allocate(zdt,(NHOR))
!       __allocate(zdlen,(NHOR))
!
!       __allocate(zLifePlants,(NHOR))
!       __allocate(zSuccessPlants,(NHOR))
!
!       zLIT(:) = 0.0
!
! !     ------------------------------------------------------------------
! !     * environmental variables
! !     ------------------------------------------------------------------
!
!       zdt(:)    = pjedi_temp(:) - cT0
!       zFT(:)    = 0.5 * (TANH((zdt(:) - 15.0) / 8.0) + 1.0)          ! KM2000
!       zFT2(:)   = pQ10R**((zdt(:) - 10.0)/10.0)
!       zFT3(:)   = MAX(0.0, zdt(:) - cT_MIN) /( MAX(0.0, zdt(:) - cT_MIN) + (pT_CRIT - cT_MIN))
!
!       if (pdaylen .eq. 0.0) then
!         call globe_daylength(kdpy,kdiy,pjedi_lat,zdlen,NHOR)
!       else
!         zdlen   = pdaylen
!       endif
!
!       zSSG(:)   = cLambda * EXP(17.269 * zdt(:)/(237.3 + zdt(:)) ) / (237.3 + zdt(:))**2
!       zSSG(:)   = zSSG(:) / (zSSG(:) + cGamma) / cLambda
!
! !     ------------------------------------------------------------------
!
!       do iSPP = FirstSPP, kMaxSPP
!
! !       ----------------------------------------------------------------
! !       * calculate land surface parameters
! !       ----------------------------------------------------------------
!
! !       zWMAX(:)  = MAX(0.5 * pjedi_paw(:), pC_WSMAX * SQRT(MAX(0.,rCWR(:,iSPP))) * pjedi_paw(:)/100.)
!         zWMAX(:)  = pC_WSMAX * SQRT(MAX(0.0, rCWR(:,iSPP))) * pjedi_paw(:)
!         where (zWMAX(:) .lt. pPAW0)
!           zWMAX(:) = pPAW0
!         end where
!         zLAI(:)   = pC_LAI * rCL(:,iSPP)
!         zFVEG(:)  = 1.0 - EXP(-cK * zLAI(:))
!         zFFOR(:)  = 1.0 - EXP(-rCWL(:,iSPP) / pC_FFOR)
!         zALB(:)   = pAlb0 * (1. - zFVEG(:)) + pAlb1 * zFVEG(:)
!         zWLMAX(:) = pC_WLMAX * ((1. - zFVEG(:)) + zFVEG(:) * zLAI(:))
!         zFSNOW(:) = MIN(1.0, rWS(:,iSPP) / 15.0)
!         zFW(:)    = rW(:,iSPP) / zWMAX(:)
!
! !       ----------------------------------------------------------------
! !       water balance
! !       ----------------------------------------------------------------
!
! !       *** may not conserve water balance by this:
!
!         rW(:, iSPP) = MIN(rW(:,iSPP),  zWMAX(:))
!         rWL(:, iSPP)= MIN(rWL(:,iSPP), zWLMAX(:))
!
! !       * water flux calculations
! !       * let it rain/snow
!
!         zA(:)     = 3.3 - zdt(:)
!         zA(:)     = MIN(1.0, MAX(0.0, zA(:) / 4.4))
!
!         zSnow(:)  = pjedi_prec(:) * zA(:)
!         zRain(:)  = pjedi_prec(:) - zSnow(:)
!
!         zMelt(:)  = MIN(3.22 * MAX(0.0, zdt(:)), zSnow(:) + rWS(:,iSPP))
!
! !       * interception by canopy
!
!         zThFall(:)= MAX(0.0, rWL(:,iSPP)- zWLMAX(:)+ zRain(:) + zMelt(:) )
!
! !       * calculate surface runoff
!
!         zRunoff(:)= MAX(0.0, rW(:,iSPP) - zWMAX(:))
!         zInfil(:) = zThFall(:) - zRunoff(:)
!
! !       * no infiltration with frozen soil
!
!         where (zdt(:) .lt. 0.0)
!           zRunoff(:) = zRunoff(:) + zInfil(:)
!           zInfil(:)  = 0.0
!         end where
!
! !       * soil drainage
!
!         zDrain(:)    = 0.0
!         zSubDrain(:) = MAX(0.0, rWSUB(:,iSPP) - cWSUB_MAX)
!
! !       * update state variables
!
!         rWS(:,iSPP)   = rWS(:,iSPP)   + zSnow(:)  - zMelt(:)
!         rWL(:,iSPP)   = rWL(:,iSPP)   + zRain(:)  + zMelt(:)  - zThFall(:)
!         rW(:,iSPP)    = rW(:,iSPP)    + zInfil(:) - zDrain(:)
!         rWSUB(:,iSPP) = rWSUB(:,iSPP) + zDrain(:) - zSubDrain(:)
!
! !       * calculate snow albedo
!
!         zSALBMIN(:) = 0.3 * zFFOR(:) + 0.4 * (1.0 - zFFOR(:))
!         zSALBMAX(:) = 0.4 * zFFOR(:) + 0.8 * (1.0 - zFFOR(:))
!         zA(:)       = MIN(1.0, MAX(0.0,(zdt(:) - (-10.0)) / (0.0 - (-10.0))))
!         zSALB(:)    = zSALBMAX(:) - (zSALBMAX(:) - zSALBMIN(:)) * zA(:)
!         zSNOWALB(:) = zALB(:) + (zSALB(:) - zALB(:)) * rWS(:, iSPP) / (rWS(:, iSPP) + 10.0)
!         zRNET(:)    = pjedi_swdown(:) * (1.0 - zSNOWALB(:)) + zdlen(:) / 24.0 * pjedi_lwnet(:)
!         zRNET(:)    = MAX(0.0, zRNET(:)) * 86400.0 * zFSNOW(:)
!
! !       * snow evaporation
!
!         zPET(:)    = MAX(0.0, cPETCorr * zSSG(:) * zRNET(:))
!         zSEVAP(:)  = MIN(rWS(:,iSPP), zPET(:))
!
! !       * evapotranspiration from snow free areas
!
!         zRNET(:)   = pjedi_swdown(:) * (1.0 - zALB(:)) + zdlen(:) / 24.0 * pjedi_lwnet(:)
!         zRNET(:)   = MAX(0.0, zRNET(:)) * 86400.0 * (1.0 - zFSNOW(:))
!         zPAR(:)    = 0.55 * pjedi_swdown(:)
!         zPET(:)    = MAX(1.0E-6, cPETCorr * zSSG(:) * zRNET(:))
!
! !       * evaporation from skin reservoir of leaves
!
!         zLEVAP(:)  = MIN(rWL(:,iSPP), zPET(:))
!
! !       * evaporation from bare soil
!
! !       zWL(:)     = MIN(0.5 * pjedi_paw(:), zWMAX(:))
!         zWL(:)     = MIN(pPAW0, zWMAX(:))
!         zA(:)      = rW(:, iSPP) - (zWMAX(:) - zWL(:))
!         zB(:)      = 0.5 * (1.0 - COS(c2PI / 2 * zA(:) / zWL(:)))
!         zBEVAP(:)  = zB(:) * (zPET(:) - zLEVAP(:)) * (1.0 - zFVEG(:))
!         zBEVAP(:)  = MIN(zBEVAP(:), rW(:,iSPP))
!
! !       ----------------------------------------------------------------
! !       plant fluxes: transpiration and carbon uptake
! !       ----------------------------------------------------------------
!
! !       * transpiration
!
!         where (rLIVE(:,iSPP) .ge. 0.)
!
! !         * different water stress formulation to ECHAM
!
!           zSET(:)   = pC_TRS * rCR(:,iSPP) * zFW(:)
!           zFH2O(:)  = 1.0 - EXP(-zSET(:) / zPET(:))
!           zTRANS(:) = zFH2O(:) * (zPET(:) - zLEVAP(:)) * zFVEG(:)
!           zTRANS(:) = MIN(zTRANS(:), rW(:,iSPP) - zBEVAP(:))
!
! !         * calculation of net primary production
!
!           zGPP(:)   = pC_GPP * p15(iSPP) * zFT(:)  * zFH2O(:)   * zFVEG(:)  * zPAR(:)
!           zRES(:)   = pC_RES * p15(iSPP) * zFT2(:) * (rCL(:,iSPP) + rCR(:,iSPP) + cWRFrac * (rCWL(:,iSPP) + rCWR(:,iSPP)))
!           zNPP(:)   = zGPP(:) - zRES(:)
!
!         elsewhere
!
!           zFH2O(:)  = 0.0
!           zTRANS(:) = 0.0
!           zGPP(:)   = 0.0
!           zRES(:)   = 0.0
!           zNPP(:)   = 0.0
!
!         end where
!
! !       * update soil water budget
!
!         rWS(:,iSPP) = rWS(:,iSPP) - zSEVAP(:)
!         rWL(:,iSPP) = rWL(:,iSPP) - zLEVAP(:)
!         rW(:,iSPP)  = rW(:,iSPP)  - zBEVAP(:) - zTRANS(:)
!
! !       ----------------------------------------------------------------
! !       plant phenology
! !       ----------------------------------------------------------------
!
! !       * determine growth conditions
!
!         rGrowW(:,iSPP) = (p01(iSPP) * zFW(:) + rGrowW(:,iSPP)) / (p01(iSPP) + 1.0)
!         rGrowT(:,iSPP) = (p02(iSPP) * zFT3(:) + rGrowT(:,iSPP)) / (p02(iSPP) + 1.0)
!
!         where (rGrowW(:,iSPP) .gt. 0.5)
!           ff_GrowW(:) = 1.0
!         elsewhere
!           ff_GrowW(:) = 0.0
!         end where
!
!         where (rGrowT(:,iSPP) .gt. 0.5)
!           ff_GrowT(:) = 1.0
!         elsewhere
!           ff_GrowT(:) = 0.0
!         end where
!
!         ff_Grow(:)   = ff_GrowT(:) * ff_GrowW(:)
!         rDie(:,iSPP) = (p13(iSPP) * zNPP(:) + rDie(:,iSPP)) / (p13(iSPP) + 1.0)
!         ff_Die(:)    = rDie(:,iSPP)
!
! !       ----------------------------------------------------------------
! !       germination
! !       ----------------------------------------------------------------
!
! !       * start a new life?
!
!         where(ff_Grow(:) .gt. 0.5 .and. rLIVE(:,iSPP) .lt. 0.0)
!           rLIVE(:,iSPP) = 0.0
!           rCA(:,iSPP)   = p14(iSPP) * pA0
!           dGGermCnt(:)  = dGGermCnt(:) + 1.0
!         end where
!
! !       ----------------------------------------------------------------
! !       loop over all living plants
! !       ----------------------------------------------------------------
!
!         where(rLIVE(:,iSPP) .ge. 0.0)
!
! !         * seed allocation
!
!           where (rCA(:,iSPP) .gt. p14(iSPP) * pA0)
!             ff_Seed(:) = 1.0
!           end where
!
! !         * marginal growth/plasticity -- not implemented yet
!
!           ff_GrowL(:) = 1.0
!           ff_GrowR(:) = 1.0
!
! !         * allocation parameters
!
!           ff_A(:)  = p05(iSPP) + p06(iSPP) + p07(iSPP) + p08(iSPP)
!           ff_AS(:) = p05(iSPP) / ff_A(:) * ff_Seed(:)
!           ff_AL(:) = p06(iSPP) / ff_A(:) * ff_Grow(:) * ff_GrowL(:)
!           ff_AR(:) = p07(iSPP) / ff_A(:) * ff_Grow(:) * ff_GrowR(:)
!           ff_LW(:) = p09(iSPP)
!           ff_RW(:) = p10(iSPP)
!           ff_LD(:) = 0.0
!           ff_RD(:) = 0.0
!           ff_A(:)  = ff_AS(:) + ff_AR(:) + ff_AL(:)
!
! !         * senescence
!
!           where((ff_Die(:) .lt. 0.0) .and. (zNPP(:) .lt. 0.0))
!             ff_LD(:) = p16(iSPP)
!             ff_RD(:) = (1.0 - p16(iSPP))
!           end where
!
! !         * solve equations of state
!
!           where(ff_A(:) .gt. 0.0)
!             fff_(:)  = ff_A(:)
!             ffx_(:)  = zNPP(:) / fff_(:)
!
!             zALLOC(:)    = zNPP(:) + (rCA(:,iSPP) - zNPP(:) / fff_(:)) * (1.0 - EXP(-fff_(:)))
!             rCA(:,iSPP) = ffx_(:) + (rCA(:,iSPP) - ffx_(:)) * EXP(-fff_(:))
!
! !           * this should not occur !!!
!             zALLOC(:) = MIN(rCA(:,iSPP) + zNPP(:), zALLOC(:))
!             zALLOC(:) = MAX(0.0, zALLOC(:))
!
!             rCL(:,iSPP)  = rCL(:,iSPP)  + ff_AL(:) / ff_A(:) * zALLOC(:) * (1.0 - ff_LW(:)) * cGResLeaf - ff_LD(:) * cDieDeltaC
!             rCR(:,iSPP)  = rCR(:,iSPP)  + ff_AR(:) / ff_A(:) * zALLOC(:) * (1.0 - ff_RW(:)) * cGResRoot - ff_RD(:) * cDieDeltaC
!             rCWL(:,iSPP) = rCWL(:,iSPP) + ff_AL(:) / ff_A(:) * zALLOC(:) * ff_LW(:) * cGResWood
!             rCWR(:,iSPP) = rCWR(:,iSPP) + ff_AR(:) / ff_A(:) * zALLOC(:) * ff_RW(:) * cGResWood
!             rCS(:,iSPP)  = rCS(:,iSPP)  + ff_AS(:) / ff_A(:) * zALLOC(:)
!           end where
!
!           where(ff_A(:) .le. 0.)
!             rCA(:,iSPP) = rCA(:,iSPP) + zNPP(:)
!             zALLOC(:)   = 0.0
!             rCL(:,iSPP) = rCL(:,iSPP) - ff_LD(:) * cDieDeltaC
!             rCR(:,iSPP) = rCR(:,iSPP) - ff_RD(:) * cDieDeltaC
!           end where
!
!           rCL(:,iSPP)   = MAX(0.0, rCL(:,iSPP))
!           rCR(:,iSPP)   = MAX(0.0, rCR(:,iSPP))
!
!           rLIVE(:,iSPP) = rLIVE(:,iSPP) + 1.0
!
! !         * decide to give up?
!
!           where(rCA(:,iSPP) .le. 0.0)
!             rCA(:,iSPP)    = 0.0
!             rCL(:,iSPP)    = 0.0
!             rCR(:,iSPP)    = 0.0
!             rCWL(:,iSPP)   = 0.0
!             rCWR(:,iSPP)   = 0.0
!             rGrowW(:,iSPP) = 0.0
!             rGrowT(:,iSPP) = 0.0
!             rDie(:,iSPP)   = 0.0
!             rLIVE(:,iSPP)  = -1.0
!             dGDeadCnt(:) = dGDeadCnt(:) + 1.0
!           end where
!
! !         * transfer water from sublayer to rooting zone layer because of
! !           root growth
!
! !         zDW(:)   = 0.0
!
! !         zDW(:)   = MAX(0.0, pC_WSMAX * rCWR(:,iSPP) - 0.5 * pjedi_paw(:)) - MAX(0.0, pC_WSMAX * rCWR(:,iSPP) - 0.5 * pjedi_paw(:))
!           zDW(:)   = MAX(0.0, pC_WSMAX * SQRT(rCWR(:,iSPP)) * pjedi_paw(:) - pPAW0) - MAX(0.0, zWMAX(:) - pPAW0)
!           zDW(:)   = rWSUB(:,iSPP) / cWSUB_MAX * zDW(:)
!           zDW(:)   = MIN(zDW(:), rWSUB(:,iSPP))
!           zDW(:)   = MAX(0.0, zDW(:))
!
!           rW(:,iSPP)    = rW(:,iSPP)    + zDW(:)
!           rWSUB(:,iSPP) = rWSUB(:,iSPP) - zDW(:)
!
!         end where
!
! !       ----------------------------------------------------------------
! !       update accumulated variables
! !       ----------------------------------------------------------------
!
!         if (nonepoint .eq. 1) then      !! extra output fields for single point runs
!
! !         * species-based accumulation for one point
!           dSARAbd(:,iSPP)  = dSARAbd(:,iSPP) + dRAbd(:,iSPP)
!           dSAGPP(:,iSPP)   = dSAGPP(:,iSPP) + zGPP(:)
!           dSACS(:,iSPP)    = dSACS(:,iSPP)  + rCS(:,iSPP)
!           dSACA(:,iSPP)    = dSACA(:,iSPP)  + rCA(:,iSPP)
!           dSACL(:,iSPP)    = dSACL(:,iSPP)  + rCL(:,iSPP)
!           dSACR(:,iSPP)    = dSACR(:,iSPP)  + rCR(:,iSPP)
!           dSACWL(:,iSPP)   = dSACWL(:,iSPP) + rCWL(:,iSPP)
!           dSACWR(:,iSPP)   = dSACWR(:,iSPP) + rCWR(:,iSPP)
!
!           dSALAI(:,iSPP)   = dSALAI(:,iSPP)  + zLAI(:)
!           dSAWMAX(:,iSPP)  = dSAWMAX(:,iSPP) + zWMAX(:)
!           dSANPP(:,iSPP)   = dSANPP(:,iSPP)  + zNPP(:)
!           dSARES(:,iSPP)   = dSARES(:,iSPP)  + zRES(:)
!           dSALIT(:,iSPP)   = dSALIT(:,iSPP)  + zLIT(:)
!           dSAFVEG(:,iSPP)  = dSAFVEG(:,iSPP) + zFVEG(:)
!           dSAFH2O(:,iSPP)  = dSAFH2O(:,iSPP) + zFH2O(:)
!           dSAFT(:,iSPP)    = dSAFT(:,iSPP)   + zFT(:)
!
!           dSAET(:,iSPP)    = dSAET(:,iSPP)    + (zSEVAP(:) + zLEVAP(:) + zBEVAP(:) + zTRANS(:))
!           dSASEVAP(:,iSPP) = dSASEVAP(:,iSPP) + zSEVAP(:)
!           dSALEVAP(:,iSPP) = dSALEVAP(:,iSPP) + zLEVAP(:)
!           dSABEVAP(:,iSPP) = dSABEVAP(:,iSPP) + zBEVAP(:)
!           dSATRANS(:,iSPP) = dSATRANS(:,iSPP) + zTRANS(:)
!
! !         * grid-based accumuation
! !           using relative abundances of previous time step
!
!           dGAGPP(:)  = dGAGPP(:)  +  dRAbd(:,iSPP) * zGPP(:)
!           dGANPP(:)  = dGANPP(:)  +  dRAbd(:,iSPP) * zNPP(:)
!           dGARES(:)  = dGARES(:)  +  dRAbd(:,iSPP) * zRES(:)
!           dGALIT(:)  = dGALIT(:)  +  dRAbd(:,iSPP) * zLIT(:)
!           dGAFH2O(:) = dGAFH2O(:) +  dRAbd(:,iSPP) * zFH2O(:)
!           dGAFT(:)   = dGAFT(:)   +  dRAbd(:,iSPP) * zFT(:)
!
!           dGAET(:)   = dGAET(:)   +  dRAbd(:,iSPP) * (zSEVAP(:) + zLEVAP(:) + zBEVAP(:) + zTRANS(:))
!
!           dGAALB(:)  = dGAALB(:)  +  dRAbd(:,iSPP) * zALB(:)
!           dGALAI(:)  = dGALAI(:)  +  dRAbd(:,iSPP) * zLAI(:)
!           dGAFVEG(:) = dGAFVEG(:) +  dRAbd(:,iSPP) * zFVEG(:)
!           dGAFFOR(:) = dGAFFOR(:) +  dRAbd(:,iSPP) * zFFOR(:)
!           dGAWMAX(:) = dGAWMAX(:) +  dRAbd(:,iSPP) * zWMAX(:)
!         else
!
! !         * grid-based accumuation
! !           using relative abundances of previous time step
!
!           dGAGPP(:)    = dGAGPP(:)  +  dRAbd(:,iSPP) * zGPP(:)
!           dGANPP(:)    = dGANPP(:)  +  dRAbd(:,iSPP) * zNPP(:)
!           dGARES(:)    = dGARES(:)  +  dRAbd(:,iSPP) * zRES(:)
!           dGALIT(:)    = dGALIT(:)  +  dRAbd(:,iSPP) * zLIT(:)
!           dGAFH2O(:)   = dGAFH2O(:) +  dRAbd(:,iSPP) * zFH2O(:)
!           dGAFT(:)     = dGAFT(:)   +  dRAbd(:,iSPP) * zFT(:)
!
!           dGAET(:)     = dGAET(:)   +  dRAbd(:,iSPP) * (zSEVAP(:) + zLEVAP(:) + zBEVAP(:) + zTRANS(:))
!
!           dGAALB(:)    = dGAALB(:)  +  dRAbd(:,iSPP) * zALB(:)
!           dGALAI(:)    = dGALAI(:)  +  dRAbd(:,iSPP) * zLAI(:)
!           dGAFVEG(:)   = dGAFVEG(:) +  dRAbd(:,iSPP) * zFVEG(:)
!           dGAFFOR(:)   = dGAFFOR(:) +  dRAbd(:,iSPP) * zFFOR(:)
!           dGAWMAX(:)   = dGAWMAX(:) +  dRAbd(:,iSPP) * zWMAX(:)
!
! !         * species-based accumuation
!
!           dSARAbd(:,iSPP) = dSARAbd(:,iSPP) + dRAbd(:,iSPP)
!           dSAGPP(:,iSPP)  = dSAGPP(:,iSPP)  + zGPP(:)
!           dSANPP(:,iSPP)  = dSANPP(:,iSPP)  + zNPP(:)
!           dSACS(:,iSPP)   = dSACS(:,iSPP)   + rCS(:,iSPP)
!           dSACA(:,iSPP)   = dSACA(:,iSPP)   + rCA(:,iSPP)
!           dSACL(:,iSPP)   = dSACL(:,iSPP)   + rCL(:,iSPP)
!           dSACR(:,iSPP)   = dSACR(:,iSPP)   + rCR(:,iSPP)
!           dSACWL(:,iSPP)  = dSACWL(:,iSPP)  + rCWL(:,iSPP)
!           dSACWR(:,iSPP)  = dSACWR(:,iSPP)  + rCWR(:,iSPP)
!         endif
!
!       enddo
!
!       nAccuCount  = nAccuCount + 1
!       nSACount    = nSACount   + 1
!
! !     ------------------------------------------------------------------
!
! !     * deallocate temporary fields
!
!       __deallocate(zWMAX)
!       __deallocate(zALB)
!       __deallocate(zFVEG)
!       __deallocate(zLAI)
!       __deallocate(zFFOR)
!       __deallocate(zFSNOW)
!       __deallocate(zWLMAX)
!       __deallocate(zFW)
!       __deallocate(zFT2)
!       __deallocate(zFT3)
!       __deallocate(zFT)
!       __deallocate(zA)
!       __deallocate(zB)
!       __deallocate(zThFall)
!       __deallocate(zTRANS)
!       __deallocate(zRunoff)
!       __deallocate(zInfil)
!       __deallocate(zMelt)
!       __deallocate(zDrain)
!       __deallocate(zSET)
!       __deallocate(zFH2O)
!       __deallocate(zPAR)
!       __deallocate(zGPP)
!       __deallocate(zNPP)
!       __deallocate(zRES)
!       __deallocate(zLIT)
!       __deallocate(ff_GrowT)
!       __deallocate(ff_GrowW)
!       __deallocate(ff_Grow)
!       __deallocate(ff_Die)
!       __deallocate(ff_Seed)
!       __deallocate(ff_A)
!       __deallocate(ff_AL)
!       __deallocate(ff_AS)
!       __deallocate(ff_AR)
!       __deallocate(ff_GrowL)
!       __deallocate(ff_GrowR)
!       __deallocate(ff_LW)
!       __deallocate(ff_RW)
!       __deallocate(ff_LD)
!       __deallocate(ff_RD)
!       __deallocate(ZALLOC)
!       __deallocate(fff_)
!       __deallocate(ffx_)
!       __deallocate(zDW)
!
!       __deallocate(zPET)
!       __deallocate(zLEVAP)
!       __deallocate(zRAIN)
!       __deallocate(zBEVAP)
!       __deallocate(zSSG)
!       __deallocate(zzSSG)
!       __deallocate(zSALB)
!       __deallocate(zSALBMAX)
!       __deallocate(zSALBMIN)
!
!       __deallocate(zSNOW)
!       __deallocate(zSNOWALB)
!       __deallocate(zRNET)
!       __deallocate(zSEVAP)
!       __deallocate(zWL)
!       __deallocate(zSubDrain)
!
!       __deallocate(zdt)
!       __deallocate(zdlen)
!
!       __deallocate(zLifePlants)
!       __deallocate(zSuccessPlants)
!
!       return
!       end subroutine jedi_step_model_ind
