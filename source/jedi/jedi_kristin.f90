      subroutine jedi_init_fields_fromGP (gp,iSPP)
      use jedi_mod
      implicit none

      integer :: iSPP, gp

      if (kmypid == kroot) then
        grCA(:,iSPP)    = 0.0
        grCL(:,iSPP)    = 0.0
        grCR(:,iSPP)    = 0.0
        grCWL(:,iSPP)   = 0.0
        grCWR(:,iSPP)   = 0.0
        grLIVE(:,iSPP)  = -1.0
        grCS(:,iSPP)    = 0.0
        if (kPop_dyn .eq. 2) grArea(:,iSPP) = 0.0
        if (kPop_dyn .eq. 2) grAreaBare(:)  = 1.0
        if (kPop_dyn .eq. 0) grCS(:,:)      = 0.0
        kSPPID(iSPP)    = iSPP
      endif

      rW(:,iSPP)     = 0.0
      rWS(:,iSPP)    = 0.0
      rWL(:,iSPP)    = 0.0
      rWSUB(:,iSPP)  = 0.0
      rCA(:,iSPP)    = 0.0
      rCL(:,iSPP)    = 0.0
      rCR(:,iSPP)    = 0.0
      rCWL(:,iSPP)   = 0.0
      rCWR(:,iSPP)   = 0.0
      rLIVE(:,iSPP)  = -1.0
      rGrowW(:,iSPP) = 0.0
      rGrowG(:,iSPP) = 0.0
      rGrowT(:,iSPP) = 0.0
      rDie(:,iSPP)   = 0.0
      rCS(:,iSPP)    = 0.0
      if (kPop_dyn .eq. 2) rArea(:,iSPP) = 0.0
      if (kPop_dyn .eq. 2) rAreaBare(:)  = 1.0
      if (kPop_dyn .eq. 0) rCS(:,:)      = 0.0


      rW(gp,iSPP)     = 0.0
      rWS(gp,iSPP)    = 0.0
      rWL(gp,iSPP)    = 0.0
      rWSUB(gp,iSPP)  = 0.0
      rCA(gp,iSPP)    = 0.0
      rCL(gp,iSPP)    = 0.0
      rCR(gp,iSPP)    = 0.0
      rCWL(gp,iSPP)   = 0.0
      rCWR(gp,iSPP)   = 0.0
      rLIVE(gp,iSPP)  = 0.0
      rGrowW(gp,iSPP) = 0.49
      rGrowG(gp,iSPP) = 0.49
      rGrowT(gp,iSPP) = 0.49
      rDie(gp,iSPP)   = 0.0
      rCS(gp,iSPP)    = pA0
      if (kPop_dyn .eq. 2) rArea(gp,iSPP) = 0.0
      if (kPop_dyn .eq. 2) rAreaBare(gp)  = 1.0

      return
      end subroutine jedi_init_fields_fromGP

!     ******************************************************************
!     JEDI_DYN_AREASTEP
!     ******************************************************************

      subroutine jedi_dyn_areastep ()
      use jedi_dyn_mod
      use jedi_mod
      implicit none

      integer :: iSPP, jSPP, i, j
      character(len=135) :: msg ! Only for output into jedi_diag

      rCtot(:,:) = rCWL(:,:) + rCL(:,:) + rCR(:,:) + rCWR(:,:) + rCA(:,:)

      do i = 1, NHOR

!       * update bare area
        rAreaBare(i) = max (0.0, 1.0 - sum(rArea(i,FirstSPP:kMaxSPP)) )

!       * biomass dominance
        sumCab(i) = MAX(cTiny, SUM(rCtot(i,FirstSPP:kMaxSPP)) )
        dominance(i,:) = rCtot(i,:)/sumCab(i)

!       * ---------------------------------------
!       * Establishment: bare ground colonization
!       * ---------------------------------------

        dfAreaGerm(:) = 0.0
        fG(:)         = 0.0
        sumfG(:)      = 0.0
        NPP0(:)       = 0.0
        tau(:)        = 0.0
        sumdfAreaG    = 0.0

        if (pfGerm < 10000 ) then
!         * Seed Competition
          fG(FirstSPP:kMaxSPP)  = 1.0 - EXP(-1.0 * pfGerm * fGerm(i,FirstSPP:kMaxSPP))
        else
!         * No seed competition
          fG(FirstSPP:kMaxSPP)  = 1.0
        endif

!       fG(FirstSPP:kMaxSPP)  = 1.0 - EXP(-1.0 * pfGerm * fGerm(i,FirstSPP:kMaxSPP))

        sumfG(i) = MAX(cTiny, SUM(fG(FirstSPP:kMaxSPP)))

        do iSPP = FirstSPP, kMaxSPP
          if (fG(iSPP) .gt. 0.0 ) then
            NPP0(iSPP)   = pC_GPP * p15(iSPP) * zFT(i) * 0.55 * pjedi_swdown(i)
!           tau(iSPP)    = MIN(1.0,NPP0(iSPP)/max(NPP0(iSPP)*10,rCtot(i,iSPP)))
            if (rCtot(i,iSPP) .ne. 0.0) then
              tau(iSPP)  = MIN(1.0,NPP0(iSPP)/rCtot(i,iSPP))
            else
              tau(iSPP)  = 1.0
            endif
            dfAreaGerm(iSPP) = fG(iSPP) / sumfG(i) * rAreaBare(i) * tau(iSPP)
          endif
        enddo

        sumdfAreaG = SUM(dfAreaGerm(FirstSPP:kMaxSPP))

!       * should never happen, just if tau
        if (sumdfAreaG .gt. rAreaBare(i)) then
          dfAreaGerm(FirstSPP:kMaxSPP) = dfAreaGerm(FirstSPP:kMaxSPP)  &
     &                                 / sumdfAreaG * rAreaBare(i)
          sumdfAreaG = SUM(dfAreaGerm(FirstSPP:kMaxSPP))
        endif

!       * ------------------------------------------------
!       * competition: invasion and exclusion by dominance
!       * ------------------------------------------------

        comp(:,:)     = 0.0
        dfAreaComp(:) = 0.0
        dfAreaExcl(:) = 0.0

        if (pfComp < 10000 ) then
          do iSPP = FirstSPP, kMaxSPP
            if (dominance(i,iSPP) .gt. 0.0 ) then
              do jSPP = FirstSPP, kMaxSPP
                if (rArea(i,jSPP) .gt. 0.0 ) then
                  comp(iSPP,jSPP) = max(0.0, ((dominance(i,iSPP) - dominance(i,jSPP))) &
     &                            ** pfComp * rArea(i,jSPP) * tau(iSPP) * rArea(i,iSPP))
                endif
              enddo
            endif
          enddo
        endif

!       do iSPP = FirstSPP, kMaxSPP
!         do jSPP = FirstSPP, kMaxSPP
!           comp(iSPP,jSPP) = MAX(0.0, (dominance(i,iSPP) - dominance(i,jSPP)) &
!    &                      * pfComp * rArea(i,jSPP) * tau(iSPP))
!         enddo
!       enddo

        do iSPP = FirstSPP, kMaxSPP
          dfAreaComp(iSPP) = SUM(comp(iSPP,:))
          dfAreaExcl(iSPP) = SUM(comp(:,iSPP))
        enddo

!       * --------------------------------------
!       * mortality and perturbation/disturbance
!       * --------------------------------------

        dfAreaM(:) = 0.0
        do iSPP = FirstSPP, kMaxSPP
          if (rCtot(i,iSPP) .le. cTiny) then
            mTau(i,iSPP) = max(1/365.0, fLIT(i,iSPP) + zRES(i,iSPP))
          else
            mTau(i,iSPP) = (fLIT(i,iSPP) + zRES(i,iSPP)) / rCtot(i,iSPP)
          endif
        enddo
        dfAreaM(:) = MIN(rArea(i,:), pMortTau * mTau(i,:) * rArea(i,:))

!       * ------------
!       * Accumulation
!       * ------------

        dGAtau(i)    = dGAtau(i)    + SUM(tau(FirstSPP:kMaxSPP))
        dGAtauM(i)   = dGAtauM(i)   + SUM(mTau(i,FirstSPP:kMaxSPP))
        dSAtau(i,:)  = dSAtau(i,:)  + tau(:)
        dSAtauM(i,:) = dSAtauM(i,:) + mTau(i,:)

        if (nonepoint .eq. 1) then
          dSAMort(i,:) = dSAMort(i,:) + dfAreaM(:)
          dSAGerm(i,:) = dSAGerm(i,:) + dfAreaGerm(:)
          dSACol(i,:)  = dSACol(i,:)  + dfAreaComp(:)
          dSAExcl(i,:) = dSAExcl(i,:) + dfAreaExcl(:)
        endif

        dGAMort(i)   = dGAMort(i)   + SUM(dfAreaM(:))
        dGAGerm(i)   = dGAGerm(i)   + SUM(dfAreaGerm(:))
        dGACol(i)    = dGACol(i)    + SUM(dfAreaComp(:))
        dGAExcl(i)   = dGAExcl(i)   + SUM(dfAreaExcl(:))

!       * -----------
!       * Update Area
        rArea(i,:) = rArea(i,:) + dfAreaGerm(:) + dfAreaComp(:) - dfAreaExcl(:) - dfAreaM(:)
        rAreaBare(i) = rAreaBare(i) - SUM(dfAreaGerm(:)) + SUM(dfAreaM(:))

!       * death population
        do iSPP = FirstSPP, kMaxSPP
          if (rArea(i,iSPP) .le. cTiny .or. rCtot(i,iSPP) .lt. cTiny ) then
            rAreaBare(i)   = rAreaBare(i) + rArea(i,iSPP)
            rArea(i,iSPP)  = 0.0
            fLIT(i,iSPP)   = fLIT(i,iSPP) + rCA(i,iSPP) + rCL(i,iSPP) + rCR(i,iSPP) + rCWR(i,iSPP) + rCWL(i,iSPP)
            rCA(i,iSPP)    = 0.0
            rCL(i,iSPP)    = 0.0
            rCR(i,iSPP)    = 0.0
            rCWL(i,iSPP)   = 0.0
            rCWR(i,iSPP)   = 0.0
            rDie(i,iSPP)   = 0.0
            rLive(i,iSPP)  = 0.0
            if (rCS(i,iSPP) .le. 0.0) rLive(i,iSPP) = -1.0  ! locally extinct
          endif
        enddo
      enddo

!     __diag(kFile_Diag,'areastep')
!     __diag_num(kFile_Diag,'   AREASTEP with  ')

      return
      end subroutine jedi_dyn_areastep
