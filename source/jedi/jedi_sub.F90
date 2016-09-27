#define __ACTIVATE
#include "../globe/globe_macros.f90"

!     ******************************************************************
!     JEDI_INIT_FIELDS
!     ******************************************************************

      subroutine jedi_init_fields (iSPP)
      use jedi_mod
      implicit none

      integer :: iSPP

      if (kmypid == kroot) then
        if (iSPP .le. kNumSPP + FirstSPP - 1) then
          grCA(:,iSPP)    = 0.0
          grCL(:,iSPP)    = 0.0
          grCR(:,iSPP)    = 0.0
          grCWL(:,iSPP)   = 0.0
          grCWR(:,iSPP)   = 0.0
          grLIVE(:,iSPP)  = 0.0
          grCS(:,iSPP)    = pA0
          if (kPop_dyn .ge. 1) grArea(:,iSPP) = 0.0
          if (kPop_dyn .ge. 1) grAreaBare(:)  = 1.0
          if (kPop_dyn .le. 1) grCS(:,iSPP)   = 0.0
          kSPPID(iSPP) = iSPP
        else
          grCA(:,iSPP)    = 0.0
          grCL(:,iSPP)    = 0.0
          grCR(:,iSPP)    = 0.0
          grCWL(:,iSPP)   = 0.0
          grCWR(:,iSPP)   = 0.0
          grLIVE(:,iSPP)  = -1.0
          grCS(:,iSPP)    = 0.0
          if (kPop_dyn .ge. 1) grArea(:,iSPP) = 0.0
          if (kPop_dyn .ge. 1) grAreaBare(:)  = 1.0
          if (kPop_dyn .le. 0) grCS(:,iSPP)   = 0.0
          kSPPID(iSPP) = 0
        endif
      endif

      if (iSPP .le. kNumSPP + FirstSPP - 1) then
        rW(:,iSPP)     = 0.0
        rWS(:,iSPP)    = 0.0
        rWL(:,iSPP)    = 0.0
        rWSUB(:,iSPP)  = 0.0
        rCA(:,iSPP)    = 0.0
        rCL(:,iSPP)    = 0.0
        rCR(:,iSPP)    = 0.0
        rCWL(:,iSPP)   = 0.0
        rCWR(:,iSPP)   = 0.0
        rLIVE(:,iSPP)  = 0.0
        rGrowW(:,iSPP) = 0.49
        rGrowG(:,iSPP) = 0.49
        rGrowT(:,iSPP) = 0.49
        rDie(:,iSPP)   = 0.0
        rCS(:,iSPP)    = pA0
        rCtot(:,iSPP)  = 0.0
        if (kPop_dyn .ge. 1) rArea(:,iSPP) = 0.0
        if (kPop_dyn .ge. 1) rAreaBare(:)  = 1.0
        if (kPop_dyn .le. 1) rCS(:,iSPP)   = 0.0
      else
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
        rCtot(:,iSPP)  = 0.0
        if (kPop_dyn .ge. 1) rArea(:,iSPP) = 0.0
        if (kPop_dyn .ge. 1) rAreaBare(:)  = 1.0
        if (kPop_dyn .le. 1) rCS(:,iSPP)   = 0.0
      endif

      return
      end subroutine jedi_init_fields

!     ******************************************************************
!     JEDI_OUTPUT_INIT
!     ******************************************************************

      subroutine jedi_output_init ()
      use jedi_mod
      implicit none

      integer :: iSPP, kFile_onepoint, iostat
      logical :: exist
      character (40) :: znoneout

      __diag(kFile_Diag,'jedi_output_init: start')

!     * allocate output variables

      __allocate(dRAbd,(NHOR,kMaxSPP))
      __allocate(dGGermCnt,(NHOR))
      __allocate(dGDeadCnt,(NHOR))

      __allocate(dGARAbd,(NHOR))
      __allocate(dGAGPP,(NHOR))
      __allocate(dGANPP,(NHOR))
      __allocate(dGARES,(NHOR))
      __allocate(dGARESH,(NHOR))
      __allocate(dGALIT,(NHOR))

      __allocate(dGACA,(NHOR))
      __allocate(dGACL,(NHOR))
      __allocate(dGACR,(NHOR))
      __allocate(dGACWL,(NHOR))
      __allocate(dGACWR,(NHOR))
      __allocate(dGACVEG,(NHOR))

      __allocate(dGALOSS_LIT_CS,(NHOR))
      __allocate(dGALOSS_LIT_CL,(NHOR))
      __allocate(dGALOSS_LIT_CR,(NHOR))
      __allocate(dGALOSS_LIT_CWL,(NHOR))
      __allocate(dGALOSS_LIT_CWR,(NHOR))
      __allocate(dGALOSS_LIT_CSLOW,(NHOR))

      __allocate(dGASOLRAD,(NHOR))
      __allocate(dGALH,(NHOR))
      __allocate(dGASH,(NHOR))

      __allocate(dGACALOSS,(NHOR))
      __allocate(dGACLLOSS,(NHOR))
      __allocate(dGACRLOSS,(NHOR))
      __allocate(dGACWLLOSS,(NHOR))
      __allocate(dGACWRLOSS,(NHOR))

      __allocate(dGACLALLOC,(NHOR))
      __allocate(dGACRALLOC,(NHOR))
      __allocate(dGACWLALLOC,(NHOR))
      __allocate(dGACWRALLOC,(NHOR))
      __allocate(dGACSALLOC,(NHOR))

      __allocate(dGAFH2O,(NHOR))
      __allocate(dGAFT,(NHOR))

      __allocate(dGAET,(NHOR))
      __allocate(dGAQSURF,(NHOR))
      __allocate(dGAQTOT,(NHOR))
      __allocate(dGATRANS,(NHOR))
      __allocate(dGAESOIL,(NHOR))
      __allocate(dGALEVAP,(NHOR))

      __allocate(dGAALB,(NHOR))
      __allocate(dGALAI,(NHOR))
      __allocate(dGAFVEG,(NHOR))
      __allocate(dGAFFOR,(NHOR))
      __allocate(dGAWMAX,(NHOR))

      __allocate(dGAtau,(NHOR))
      __allocate(dGAtauM,(NHOR))
      __allocate(dGAGerm,(NHOR))
      __allocate(dGACol,(NHOR))
      __allocate(dGAExcl,(NHOR))
      __allocate(dGAMort,(NHOR))
      __allocate(dGADist,(NHOR))
      __allocate(dGASpec,(NHOR))
      __allocate(dGAMig,(NHOR))
      __allocate(dGAExt,(NHOR))

      __allocate(gdGARAbd,(gNHOR))
      __allocate(gdGAGPP,(gNHOR))
      __allocate(gdGANPP,(gNHOR))
      __allocate(gdGANEE,(gNHOR))
      __allocate(gdGARES,(gNHOR))
      __allocate(gdGARESH,(gNHOR))
      __allocate(gdGALIT,(gNHOR))
      __allocate(gdGARESE,(gNHOR))

      __allocate(gdGACA,(gNHOR))
      __allocate(gdGACL,(gNHOR))
      __allocate(gdGACR,(gNHOR))
      __allocate(gdGACWL,(gNHOR))
      __allocate(gdGACWR,(gNHOR))

      __allocate(gdGACVEG,(gNHOR))

      __allocate(gdGASOLRAD,(gNHOR))
      __allocate(gdGALH,(gNHOR))
      __allocate(gdGASH,(gNHOR))

      __allocate(gdGACALOSS,(gNHOR))
      __allocate(gdGACLLOSS,(gNHOR))
      __allocate(gdGACRLOSS,(gNHOR))
      __allocate(gdGACWLLOSS,(gNHOR))
      __allocate(gdGACWRLOSS,(gNHOR))

      __allocate(gdGALOSS_LIT_CS,(gNHOR))
      __allocate(gdGALOSS_LIT_CL,(gNHOR))
      __allocate(gdGALOSS_LIT_CR,(gNHOR))
      __allocate(gdGALOSS_LIT_CWL,(gNHOR))
      __allocate(gdGALOSS_LIT_CWR,(gNHOR))
      __allocate(gdGALOSS_LIT_CSLOW,(gNHOR))

      __allocate(gdGACLALLOC,(gNHOR))
      __allocate(gdGACRALLOC,(gNHOR))
      __allocate(gdGACWLALLOC,(gNHOR))
      __allocate(gdGACWRALLOC,(gNHOR))
      __allocate(gdGACSALLOC,(gNHOR))

      __allocate(gdGAFT,(gNHOR))
      __allocate(gdGAFH2O,(gNHOR))

      __allocate(gdGAET,(gNHOR))
      __allocate(gdGAQTOT,(gNHOR))
      __allocate(gdGAQSURF,(gNHOR))
      __allocate(gdGAESOIL,(gNHOR))
      __allocate(gdGALEVAP,(gNHOR))
      __allocate(gdGATRANS,(gNHOR))

      __allocate(gdGAALB,(gNHOR))
      __allocate(gdGALAI,(gNHOR))
      __allocate(gdGAFVEG,(gNHOR))
      __allocate(gdGAFFOR,(gNHOR))
      __allocate(gdGAWMAX,(gNHOR))

      __allocate(gdGAtau,(gNHOR))
      __allocate(gdGAtauM,(gNHOR))
      __allocate(gdGAGerm,(gNHOR))
      __allocate(gdGACol,(gNHOR))
      __allocate(gdGAExcl,(gNHOR))
      __allocate(gdGAMort,(gNHOR))
      __allocate(gdGADist,(gNHOR))
      __allocate(gdGASpec,(gNHOR))
      __allocate(gdGAMig,(gNHOR))
      __allocate(gdGAExt,(gNHOR))

      __allocate(BM_total,(gNHOR))
      __allocate(Richness,(gNHOR))

      __allocate(gdSARAbd,(gNHOR,kMaxSPP))
      __allocate(gdSAGPP,(gNHOR,kMaxSPP))
      __allocate(gdSANPP,(gNHOR,kMaxSPP))
      __allocate(gdSACS,(gNHOR,kMaxSPP))
      __allocate(gdSACA,(gNHOR,kMaxSPP))
      __allocate(gdSACL,(gNHOR,kMaxSPP))
      __allocate(gdSACR,(gNHOR,kMaxSPP))
      __allocate(gdSACWL,(gNHOR,kMaxSPP))
      __allocate(gdSACWR,(gNHOR,kMaxSPP))

      __allocate(gdSACTOT,(gNHOR,kMaxSPP))
      __allocate(gdSACVEG,(gNHOR,kMaxSPP))

      __allocate(gdSAtau,(gNHOR,kMaxSPP))
      __allocate(gdSAtauM,(gNHOR,kMaxSPP))
      __allocate(gdSACol,(gNHOR,kMaxSPP))
      __allocate(gdSAGerm,(gNHOR,kMaxSPP))
      __allocate(gdSAMort,(gNHOR,kMaxSPP))
      __allocate(gdSAExcl,(gNHOR,kMaxSPP))

      __allocate(dSARAbd,(NHOR,kMaxSPP))
      __allocate(dSAGPP,(NHOR,kMaxSPP))
      __allocate(dSANPP,(NHOR,kMaxSPP))
      __allocate(dSACS,(NHOR,kMaxSPP))
      __allocate(dSACA,(NHOR,kMaxSPP))
      __allocate(dSACL,(NHOR,kMaxSPP))
      __allocate(dSACR,(NHOR,kMaxSPP))
      __allocate(dSACWL,(NHOR,kMaxSPP))
      __allocate(dSACWR,(NHOR,kMaxSPP))

      __allocate(dSACTOT,(NHOR,kMaxSPP))
      __allocate(dSACVEG,(NHOR,kMaxSPP))

      __allocate(dSAtau,(NHOR,kMaxSPP))
      __allocate(dSAtauM,(NHOR,kMaxSPP))

      if (nonepoint .eq. 1) then
        __allocate(dSALAI,(NHOR,kMaxSPP))
        __allocate(dSAWMAX,(NHOR,kMaxSPP))
        __allocate(dSARES,(NHOR,kMaxSPP))
        __allocate(dSALIT,(NHOR,kMaxSPP))
        __allocate(dSAFVEG,(NHOR,kMaxSPP))
        __allocate(dSAFH2O,(NHOR,kMaxSPP))
        __allocate(dSAFT,(NHOR,kMaxSPP))

        __allocate(dSAET,(NHOR,kMaxSPP))
        __allocate(dSASEVAP,(NHOR,kMaxSPP))
        __allocate(dSALEVAP,(NHOR,kMaxSPP))
        __allocate(dSABEVAP,(NHOR,kMaxSPP))
        __allocate(dSATRANS,(NHOR,kMaxSPP))

        __allocate(dSAMort,(NHOR,kMaxSPP))
        __allocate(dSAGerm,(NHOR,kMaxSPP))
        __allocate(dSACol,(NHOR,kMaxSPP))
        __allocate(dSAExcl,(NHOR,kMaxSPP))
      endif

!     * open files

      if (koutput == 1) call globe_open_output(sfile_Grid, kFile_Grid, kFile_Diag)

      call jedi_output_reset

      call jedi_output_species_reset


      if (kmypid == kroot .and. nonepoint .eq. 1) then
        __diag(kFile_Diag,'jedi_output_init: creating single species files')
        do iSPP = FirstSPP, kMaxSPP
          write(znoneout,'("species",i5.5)') iSPP + kFrPl - 2
          kFile_onepoint = 100 + iSPP
          open(kFile_onepoint,file=trim(znoneout), status='replace',   &
     &                        form='formatted', iostat=iostat)
          call error_open_new(iostat,znoneout)
        enddo
      endif

      if (kmypid == kroot .and. kSpec_dyn .gt. 0) then
        __diag(kFile_Diag,'jedi_output_init: creating extinct species file')
        open(kFile_SPPExtinct, FILE=sfile_SPPExtinct, STATUS='replace', &
     &                         FORM='formatted', IOSTAT=iostat)
        call error_open_new(iostat,sfile_SPPExtinct)
      endif

!     * end

      __diag(kFile_Diag,'jedi_output_init: end')

      return
      end subroutine jedi_output_init

!     ******************************************************************
!     JEDI_OUTPUT_STOP
!     ******************************************************************

      subroutine jedi_output_stop ()
      use jedi_mod
      implicit none

      integer :: iSPP, kFile_onepoint

      if (koutput == 1) call globe_close_output(kFile_Grid, kFile_Diag)

      if (kmypid == kroot .and. nonepoint .eq. 1) then
        do iSPP = FirstSPP, kMaxSPP
          kFile_onepoint = 100 + iSPP
          close(kFile_onepoint)
        enddo
      endif

      if (kmypid == kroot .and. kSpec_dyn .gt. 0) then
        close(kFile_SPPExtinct)
      endif

      return
      end subroutine jedi_output_stop

!     ******************************************************************
!     JEDI_OUTPUT_RESET (reset monthly accumulations)
!     ******************************************************************

      subroutine jedi_output_reset ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none
      __diag(kFile_Diag,'jedi_output_reset: start')

!     * reset accumulated variables
      dGARAbd    = 0.0
      dGAGPP     = 0.0
      dGANPP     = 0.0
      dGARES     = 0.0
      dGARESH    = 0.0
      dGALIT     = 0.0

      dGACA     = 0.0
      dGACL     = 0.0
      dGACR     = 0.0
      dGACWL     = 0.0
      dGACWR     = 0.0
      dGACVEG    = 0.0

      dGALOSS_LIT_CS(:)    = 0.0
      dGALOSS_LIT_CL(:)    = 0.0
      dGALOSS_LIT_CR(:)    = 0.0
      dGALOSS_LIT_CWL(:)   = 0.0
      dGALOSS_LIT_CWR(:)   = 0.0
      dGALOSS_LIT_CSLOW(:) = 0.0

      dGASOLRAD = 0.0
      dGALH     = 0.0
      dGASH     = 0.0

      dGACALOSS(:)  = 0.0
      dGACLLOSS(:)  = 0.0
      dGACRLOSS(:)  = 0.0
      dGACWLLOSS(:) = 0.0
      dGACWRLOSS(:) = 0.0

      dGACLALLOC(:)  = 0.0
      dGACRALLOC(:)  = 0.0
      dGACWLALLOC(:) = 0.0
      dGACWRALLOC(:) = 0.0
      dGACSALLOC(:)  = 0.0

      dGAFH2O    = 0.0
      dGAFT      = 0.0
      dGAET      = 0.0
      dGATRANS      = 0.0
      dGAESOIL      = 0.0
      dGALEVAP      = 0.0
      dGAQSURF      = 0.0
      dGAQTOT      = 0.0
      dGAALB     = 0.0
      dGALAI     = 0.0
      dGAFVEG    = 0.0
      dGAFFOR    = 0.0
      dGAWMAX    = 0.0
      dGGermCnt  = 0.0
      dGDeadCnt  = 0.0

      dGAtau(:)  = 0.0
      dGAGerm(:) = 0.0
      dGACol(:)  = 0.0
      dGAExcl(:) = 0.0
      dGAMort(:) = 0.0
      dGADist(:) = 0.0
      dGASpec(:) = 0.0
      dGAMig(:)  = 0.0
      dGAExt(:)  = 0.0

      nAccuCount = 0

      if (kSpec_dyn .gt. 0) then
        Ext(:,:)     = 0.0
        countExt(:)  = 0.0
        countSpec(:) = 0.0
      endif

      __diag(kFile_Diag,'jedi_output_reset: end')

      return
      end subroutine jedi_output_reset

!     ******************************************************************
!     JEDI_OUTPUT_SPECIES_RESET
!     ******************************************************************

      subroutine jedi_output_species_reset ()
      use jedi_mod
      implicit none

      __diag(kFile_Diag,'jedi_output_species_reset: start')

!     * reset accumulated variables

      if (kmypid == kroot) then
        gdSARAbd(:,:)    = 0.0
        gdSAGPP(:,:)     = 0.0
        gdSANPP(:,:)     = 0.0
        gdSACS(:,:)      = 0.0
        gdSACA(:,:)      = 0.0
        gdSACL(:,:)      = 0.0
        gdSACR(:,:)      = 0.0
        gdSACWL(:,:)     = 0.0
        gdSACWR(:,:)     = 0.0
        gdSAtau(:,:)     = 0.0
        
        gdSACVEG(:,:)    = 0.0
        gdSACTOT(:,:)    = 0.0
       endif

      dSARAbd(:,:)    = 0.0
      dSAGPP(:,:)     = 0.0
      dSANPP(:,:)     = 0.0
      dSACS(:,:)      = 0.0
      dSACA(:,:)      = 0.0
      dSACL(:,:)      = 0.0
      dSACR(:,:)      = 0.0
      dSACWL(:,:)     = 0.0
      dSACWR(:,:)     = 0.0

      dSACVEG(:,:)      = 0.0
      dSACTOT(:,:)      = 0.0

      dSAtau(:,:)     = 0.0

      nSACount = 0

      __diag(kFile_Diag,'jedi_output_species_reset: end')

      return
      end subroutine jedi_output_species_reset

!     ******************************************************************
!     JEDI_RELABD
!     ******************************************************************

      subroutine jedi_relabd ()
      use jedi_mod
      implicit none

      integer :: iSPP, j

      REAL,allocatable,dimension(:) :: zCpts
      REAL,allocatable,dimension(:) :: zCtot

!     * allocate temporary fields

      allocate(zCpts(NHOR))
      allocate(zCtot(NHOR))

!     * calculate relative abundances

      zCtot(:) = 0.0
      dRAbd(:,:) = 0.0

      if (kPop_dyn .eq. 0) then
        do iSPP = FirstSPP, kMaxSPP
          where (rCS(:,iSPP) .ge. p14(iSPP) * pA0)
            zCpts(:) = rCA(:,iSPP) + rCR(:,iSPP) + rCL(:,iSPP) + rCWR(:,iSPP) + rCWL(:,iSPP)
            zCtot(:) = zCtot(:) + zCpts(:)
          end where
        enddo

        do iSPP = FirstSPP, kMaxSPP
          where (rCS(:, iSPP) .ge. p14(iSPP) * pA0)
            zCpts(:) = rCA(:,iSPP) + rCR(:,iSPP) + rCL(:,iSPP) + rCWR(:,iSPP) + rCWL(:,iSPP)
            dRAbd(:,iSPP) = zCpts(:) / MAX(1.E-06, zCtot(:))
          end where
        enddo
        dRAbd(:,BareSPP) = 0.0
      endif

      if (kPop_dyn .eq. 2) then
        if (YEAR < 69 .and. kspinup .eq. 1) then
          do iSPP = FirstSPP, kMaxSPP
            zCpts(:) = rCA(:,iSPP) + rCR(:,iSPP) + rCL(:,iSPP) + rCWR(:,iSPP) + rCWL(:,iSPP)
            zCtot(:) = zCtot(:) + zCpts(:)
          enddo
          do iSPP = FirstSPP, kMaxSPP
            zCpts(:) = rCA(:,iSPP) + rCR(:,iSPP) + rCL(:,iSPP) + rCWR(:,iSPP) + rCWL(:,iSPP)
            where (zCpts(:) .gt. cTiny)
              dRAbd(:,iSPP) = zCpts(:) / MAX(1.E-06, zCtot(:))
            end where
          enddo
          dRAbd(:,BareSPP) = 0.0
        else
          dRAbd(:,:) = rArea(:,:)
          dRAbd(:,BareSPP) = rAreaBare(:)
        endif
      endif

!     * deallocate temporary fields

      deallocate(zCpts)
      deallocate(zCtot)

      return
      end subroutine jedi_relabd

!     ******************************************************************
!     JEDI_TRANSFORM_PARMS
!     ******************************************************************

      subroutine jedi_transform_parms (iSPP)
      use jedi_mod
      implicit none

      integer :: iSPP
      real    :: zA

      if (kPop_dyn .eq. 0) then
        p01(iSPP) = 10.0**(4.0 * p01(iSPP) - 2.0)
        p02(iSPP) = 10.0**(4.0 * p02(iSPP) - 2.0)
        p13(iSPP) = 10.0**(4.0 * p13(iSPP) - 2.0)
        p14(iSPP) = 10.0**(8.0 * p14(iSPP) - 7.0)
      else
        zA = (p05(iSPP) + p06(iSPP) + p07(iSPP) + p08(iSPP))
        p01(iSPP) = 10.0**(4.0 * p01(iSPP) - 2.0)
        p02(iSPP) = 10.0**(4.0 * p02(iSPP) - 2.0)
        p03(iSPP) = 15.00 * p03(iSPP)
        p04(iSPP) = 10.0**(3.0 * p04(iSPP) - 3.0)
        p05(iSPP) = p05(iSPP) / zA
        p06(iSPP) = p06(iSPP) / zA
        p07(iSPP) = p07(iSPP) / zA
        p08(iSPP) = p08(iSPP) / zA
        p11(iSPP) = 1.0 / (80.0 * 365.0 * p11(iSPP))
        p12(iSPP) = 1.0 / (10.0**(1.48 * p12(iSPP) + 0.30) * (365.0/12.0))
        p13(iSPP) = 10.0**(4.0 * p13(iSPP) - 1.0)
        p15(iSPP) = 0.01 + pCN_Leaf * p15(iSPP)
      endif

      return
      end subroutine jedi_transform_parms

!     ******************************************************************
!     JEDI_TRANSFORM_PARMS_BACK
!     ******************************************************************

      subroutine jedi_transform_parms_back (iSPP)
      use jedi_mod
      implicit none

      integer :: iSPP
      real :: zA

      if (kPop_dyn .eq. 0) then
        p01(iSPP) = (LOG10(p01(iSPP)) + 2.0) / 4.0
        p02(iSPP) = (LOG10(p02(iSPP)) + 2.0) / 4.0
        p13(iSPP) = (LOG10(p13(iSPP)) + 2.0) / 4.0
        p14(iSPP) = (LOG10(p14(iSPP)) + 7.0) / 8.0
      else
        zA = (p05(iSPP) + p06(iSPP) + p07(iSPP) + p08(iSPP))
        p01(iSPP) = (LOG10(p01(iSPP)) + 2.0) / 4.0
        p02(iSPP) = (LOG10(p02(iSPP)) + 2.0) / 4.0
        p03(iSPP) = p03(iSPP) / 15.00
        p04(iSPP) = (LOG10(p04(iSPP)) + 3.0) / 3.0
        p05(iSPP) = p05(iSPP) * zA
        p06(iSPP) = p06(iSPP) * zA
        p07(iSPP) = p07(iSPP) * zA
        p08(iSPP) = p08(iSPP) * zA
        p11(iSPP) = 1.0 / (80.0 * 365.0 * p11(iSPP))
        p12(iSPP) = (LOG10((p12(iSPP) * (365.0/12.0))) + 0.30) / (-1.48)
        p13(iSPP) = (LOG10(p13(iSPP)) + 1.0) / 4.0
        p15(iSPP) = (p15(iSPP) - 0.01) / pCN_Leaf
      endif

      return
      end subroutine jedi_transform_parms_back
