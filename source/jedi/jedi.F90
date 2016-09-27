#define __ACTIVATE
#include "../globe/globe_macros.f90"
!> \file jedi.F90
!> \brief Main routines of JEDI

!
!     JeDi - PopDyn
!     =============
!
!> @author code development and testing by (in alphabetical order):
!       Kristin Bohn
!       Axel Kleidon
!       Ryan Pavlick
!       Bjoern Reu
!       Steffen Richter
!       Kerstin Sickel
!
!     contact:
!       Axel Kleidon (akleidon@bgc-jena.mpg.de)
!       Biospheric Theory and Modelling Group
!       Max-Planck-Institute for Biogeochemistry
!       Hans-Knoell-Str. 10
!       07745 Jena, Germany
!
!     Jena, September 2008

!     ******************************************************************
!     JEDI_INIT
!     ******************************************************************

      subroutine jedi_init ()
      use jedi_mod
      implicit none

      integer :: iostat, i
      logical :: exist
      integer :: iSPP
      integer :: ihead(8)
      integer :: kmapinput = 1

      real, allocatable :: tmp_input_g(:)
      real, allocatable :: tmp_input_ll(:)

!    * display version information

      if (kmypid == kroot) call globe_show_version

!    * set values for tiny and epsilon

      r8_tiny = tiny(r8_tiny)
      r8_epsilon = epsilon(r8_epsilon)

!    * open diagnostic output file

      call globe_open_diag(sfile_Diag, kFile_Diag)

      __diag(kFile_Diag,'jedi_init: start')

!    * read namelist parameters

      call jedi_read_namelist

!    * generating a seed for each process
      rand_seed=((kmypid + 1) * (-1) * seed)

      if (krestart .le. 0) then

!    * allocate fields

        call jedi_alloc
        if (kPop_dyn .gt. 0) call jedi_dyn_alloc

!    * initialization of fields

        rCSLOW(:)     = 0.0
        rLIT_CL(:)    = 0.0
        rLIT_CR(:)    = 0.0
        rLIT_CWL(:)   = 0.0
        rLIT_CWR(:)   = 0.0

        grCSLOW(:)     = 0.0
        grLIT_CL(:)    = 0.0
        grLIT_CR(:)    = 0.0
        grLIT_CWL(:)   = 0.0
        grLIT_CWR(:)   = 0.0

        do iSPP = BareSPP, kMaxSPP
          if (initGP .gt. 0 .and. kMig_dyn .eq. 1) then
            call jedi_init_fields_fromGP(initGP, iSPP)
          else
            call jedi_init_fields(iSPP)
          endif
        enddo

        rCS(:,BareSPP)   = 0.0
        grCS(:,BareSPP)  = 0.0

!       * read species parameter

        call jedi_read_species
        __globe_mpbc(kSPPID)

!       * initialize output files

        call jedi_output_init

      else
        call jedi_read_restart
      endif

      if (kdynamic .eq. 0) dRAbd(:,:) = rArea(:,:)

!    * optional overwrite rLive with values from file

      if (ktopspecies == 1) then
        __diag(kFile_Diag,'jedi_init: read sfile_SPPTopGrid')
        __allocate(tmp_input_g,(gNHOR))
        __allocate(tmp_input_ll,(nlon*nlat))
!       * open species input file
        call globe_open_input(sfile_SPPTopGrid, kFile_SPPTopGrid, kFile_Diag)
        if (kmypid == kroot) then
!         * check the file format
          __globe_read_srvnc(kFile_SPPTopGrid,ihead,iostat=iostat)
          if (ihead(5)*ihead(6) == kNumGPts) kmapinput=0
!         * skip land-sea-mask if present
          if (ihead(1) == 172) then
            __globe_read_srvnc(kFile_SPPTopGrid,ihead,tmp_input_g(1:kNumGPts))
          endif
        endif
        do i = 1, kMaxSPP
          if (kmypid == kroot) then
!           * read the data for the level i
            if (kmapinput == 0) then
              __globe_read_srvnc(kFile_SPPTopGrid,ihead,tmp_input_g(1:kNumGPts))
            else
              __globe_read_srvnc(kFile_SPPTopGrid,ihead,tmp_input_ll)
              tmp_input_g(1:kNumGPts) = tmp_input_ll(gkGPID(:))
            endif
!           * fill the extra points
            tmp_input_g(kNumGPts:gNHOR) = tmp_input_g(kNumGPts)
          endif
!         * scatter the field to all CPUs
          __globe_mpsc_from_to(tmp_input_g,rLive(:,i))
        enddo
!       * close species input file
        call globe_close_input(kFile_SPPTopGrid, kFile_Diag)
        __deallocate(tmp_input_g)
        __deallocate(tmp_input_ll)
      endif

!     * optimizing version

      if (kOpti == 1) call jedi_opti_init

!     * end of initialization

      __diag(kFile_Diag,'jedi_init: end')

      return
      end subroutine jedi_init

!     ******************************************************************
!     JEDI_STEP
!     ******************************************************************

      subroutine jedi_step ()
      use jedi_mod
      implicit none

      select case(kPop_dyn)
      ! case(0)
      !   call jedi_relabd
      !   call jedi_step_model_ind
      case(1)
        call jedi_step_model_dyn
      case(2)
        call jedi_relabd
        call jedi_step_model_dyn
      end select

      return
      end subroutine jedi_step

!     ******************************************************************
!     JEDI_STOP
!     ******************************************************************

      subroutine jedi_stop ()
      use jedi_mod
      implicit none

      __diag(kFile_Diag,'jedi_stop: start')

!     * write successfull species

      call jedi_success

!     * write species parameter

      if (kspec_yrout == 0) then
        call jedi_output_sppgrid
      !  call jedi_output_species !fixme
        call jedi_output_gaspp
      endif

!     * close output

      call jedi_output_stop

!     * optimizing version

      if (kOpti == 1) call jedi_opti_stop

!     * write restart file

      if (krestart .ge. 0) call jedi_write_restart

!     deallocate the fields to clean free the memory

!     call jedi_dealloc

!     __diag(kFile_Diag,'jedi_stop: dealloc end')
!     if (kPop_dyn .gt. 0) call jedi_dyn_dealloc

      __diag(kFile_Diag,'jedi_stop: end')

      call globe_close_diag(kFile_Diag)

      return
      end subroutine jedi_stop

!     ******************************************************************
!     JEDI_OUTPUT
!     ******************************************************************

      subroutine jedi_output ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, j, k, iLON, jLAT

      REAL,allocatable,dimension(:) :: zCVEG, zW, zCS, zCA, zCL, zCR
      REAL,allocatable,dimension(:) :: zCWL, zCWR, zArea
      REAL,allocatable,dimension(:) :: zSuccess, zShannon
      REAL,allocatable,dimension(:) :: zgEvenness, zgShannon
      REAL,allocatable,dimension(:) :: zgCVEG, zgSuccess
      REAL,allocatable,dimension(:) :: zdMigSeed, zdMigOut, zdMigIn
      REAL,allocatable,dimension(:) :: zExt ,abdN

!     * allocate temporary fields

      __allocate(zCVEG,(NHOR))
      __allocate(zW,(NHOR))
      __allocate(zCS,(NHOR))
      __allocate(zCA,(NHOR))
      __allocate(zCL,(NHOR))
      __allocate(zCR,(NHOR))
      __allocate(zCWL,(NHOR))
      __allocate(zCWR,(NHOR))
      __allocate(zArea,(NHOR))
      __allocate(zSuccess,(NHOR))
      __allocate(zShannon,(NHOR))

      __allocate(abdN,(NHOR))

      __allocate(zgEvenness,(gNHOR))
      __allocate(zgShannon,(gNHOR))
      __allocate(zgCVEG,(gNHOR))
      __allocate(zgSuccess,(gNHOR))
      __allocate(zdMigSeed,(gNHOR))
      __allocate(zdMigOut,(gNHOR))
      __allocate(zdMigIn,(gNHOR))
      __allocate(zExt,(gNHOR))

!     * initialize

      zCVEG(:) = 0.0
      zW(:)    = 0.0
      zCS(:)   = 0.0
      zCA(:)   = 0.0
      zCL(:)   = 0.0
      zCR(:)   = 0.0
      zCWL(:)  = 0.0
      zCWR(:)  = 0.0
      zArea(:) = 0.0

      zSuccess(:)  = 0.0
      zShannon(:)  = 0.0

      zgEvenness(:) = 0.0

      zdMigSeed(:) = 0.0
      zExt(:) = 0

      __diag(kFile_Diag,'jedi_output: start')

!     ------------------------------------------------------------------
!     * gather fields
!     ------------------------------------------------------------------

      __globe_mpga_to_from(gdGARAbd,dGARAbd)
      __globe_mpga_to_from(gdGAGPP,dGAGPP)
      __globe_mpga_to_from(gdGANPP,dGANPP)
      __globe_mpga_to_from(gdGARES,dGARES)
      __globe_mpga_to_from(gdGARESH,dGARESH)

      __globe_mpga_to_from(gdGACR,dGACR)
      __globe_mpga_to_from(gdGACA,dGACA)
      __globe_mpga_to_from(gdGACL,dGACL)
      __globe_mpga_to_from(gdGACWL,dGACWL)
      __globe_mpga_to_from(gdGACWR,dGACWR)

      __globe_mpga_to_from(gdGACVEG,dGACVEG)

      __globe_mpga_to_from(gdGASOLRAD,dGASOLRAD)
      __globe_mpga_to_from(gdGAFH2O,dGAFH2O)
      __globe_mpga_to_from(gdGAFT,dGAFT)

      __globe_mpga_to_from(gdGALOSS_LIT_CS,dGALOSS_LIT_CS)
      __globe_mpga_to_from(gdGALOSS_LIT_CL,dGALOSS_LIT_CL)
      __globe_mpga_to_from(gdGALOSS_LIT_CR,dGALOSS_LIT_CR)
      __globe_mpga_to_from(gdGALOSS_LIT_CWL,dGALOSS_LIT_CWL)
      __globe_mpga_to_from(gdGALOSS_LIT_CWR,dGALOSS_LIT_CWR)
      __globe_mpga_to_from(gdGALOSS_LIT_CSLOW,dGALOSS_LIT_CSLOW)

      __globe_mpga_to_from(grCSLOW,rCSLOW)
      __globe_mpga_to_from(grLIT_CL,rLIT_CL)
      __globe_mpga_to_from(grLIT_CR,rLIT_CR)
      __globe_mpga_to_from(grLIT_CWL,rLIT_CWL)
      __globe_mpga_to_from(grLIT_CWR,rLIT_CWR)

      __globe_mpga_to_from(gdGAET,dGAET)
      __globe_mpga_to_from(gdGAESOIL,dGAESOIL)
      __globe_mpga_to_from(gdGALEVAP,dGALEVAP)
      __globe_mpga_to_from(gdGAQSURF,dGAQSURF)
      __globe_mpga_to_from(gdGATRANS,dGATRANS)
      __globe_mpga_to_from(gdGAQTOT,dGAQTOT)

      __globe_mpga_to_from(gdGAALB,dGAALB)
      __globe_mpga_to_from(gdGALAI,dGALAI)
      __globe_mpga_to_from(gdGAFVEG,dGAFVEG)
      __globe_mpga_to_from(gdGAWMAX,dGAWMAX)
      __globe_mpga_to_from(gdGAtau,dGAtau)
      __globe_mpga_to_from(gdGAGerm,dGAGerm)
      __globe_mpga_to_from(gdGACol,dGACol)
      __globe_mpga_to_from(gdGAExcl,dGAExcl)
      __globe_mpga_to_from(gdGAMort,dGAMort)
      __globe_mpga_to_from(gdGADist,dGADist)
      __globe_mpga_to_from(gdGASpec,dGASpec)
      __globe_mpga_to_from(gdGAMig,dGAMig)
      __globe_mpga_to_from(gdGAExt,dGAExt)

      __globe_mpga_to_from(gdGACALOSS,dGACALOSS)
      __globe_mpga_to_from(gdGACLLOSS,dGACLLOSS)
      __globe_mpga_to_from(gdGACRLOSS,dGACRLOSS)
      __globe_mpga_to_from(gdGACWLLOSS,dGACWLLOSS)
      __globe_mpga_to_from(gdGACWRLOSS,dGACWRLOSS)

      __globe_mpga_to_from(gdGACLALLOC,dGACLALLOC)
      __globe_mpga_to_from(gdGACRALLOC,dGACRALLOC)
      __globe_mpga_to_from(gdGACWLALLOC,dGACWLALLOC)
      __globe_mpga_to_from(gdGACWRALLOC,dGACWRALLOC)
      __globe_mpga_to_from(gdGACSALLOC,dGACSALLOC)


!     * calculate grid-summed seed carbon
      zCS(:)   = zCS(:) + SUM(rCS(:,:), DIM=2)

!     ------------------------------------------------------------------
!     * calculate diversity indices
!     ------------------------------------------------------------------

      select case(kPop_dyn)
      ! case(0)
      !   do iSPP = FirstSPP, kMaxSPP
      !     where (rCS(:,iSPP) .ge. p14(iSPP) * pA0)
      !       zSuccess(:) = zSuccess(:) + 1.0
      !     end where
      !     where (dRAbd(:,iSPP) .gt. 0.0)
      !       zShannon(:) = zShannon(:) - dRAbd(:,iSPP) * log(dRAbd(:,iSPP))
      !     end where
      !   enddo

      case(1)
        do iSPP = FirstSPP, kMaxSPP
          where (rCtot(:,iSPP) .gt. pA0)
            zSuccess(:) = zSuccess(:) + 1.0
            zShannon(:) = zShannon(:) - dRAbd(:,iSPP) * log(dRAbd(:,iSPP))
          end where
        enddo

      case(2)
        zArea(:) = SUM(dRAbd(:,FirstSPP:kMaxSPP), dim=2)
        if (nonepoint .eq. 1) then
          do iSPP = FirstSPP, kMaxSPP
            where (rArea(:,iSPP) .gt. cTiny)
              zSuccess(:) = zSuccess(:) + 1.0
            end where
            where (dRAbd(:,iSPP) .gt. 0.0 .and. zArea(:) .gt. 0.0 )
              abdN(:) = dRAbd(:,iSPP) / zArea(:)
              zShannon(:) = zShannon(:) - abdN(:) * log(abdN(:))
            end where
          enddo
        else
          do iSPP = FirstSPP, kMaxSPP
            where (rArea(:,iSPP) .gt. cTiny)
              zSuccess(:) = zSuccess(:) + 1.0
            end where
            where (dRAbd(:,iSPP) .gt. 0.0)
              zShannon(:) = zShannon(:) - dRAbd(:,iSPP) * log(dRAbd(:,iSPP))
            end where
          enddo
        endif
        if (kSpec_dyn .gt. 0) countExt(:) = SUM(Ext(:,:), DIM=2)
      end select

      __globe_mpga_to_from(zgShannon,zShannon)
      __globe_mpga_to_from(BM_total,zCVEG)
      __globe_mpga_to_from(Richness,zSuccess)
      if (kSpec_dyn .gt. 0) then
        __globe_mpga_to_from(zExt,countExt)
      endif

      if (kmypid == kroot) then

        if (MAXVAL(zgShannon(:)) .ne. 0.0) zgEvenness(:) = zgShannon(:) / MAXVAL(zgShannon(:))

        if (kMig_dyn .eq. 1) then
          zdMigOut(:)  = SUM(gfseedout(:,:), dim=2)
          zdMigIn(:)   = SUM(gfseedin(:,:), dim=2)
          zdMigSeed(:) = SUM(gfseedin(:,:), dim=2) - SUM(gfseedout(:,:), dim=2)
        endif

!     ------------------------------------------------------------------
!       * averaging
!     ------------------------------------------------------------------

        if (nAccuCount .gt. 0) then
          gdGARAbd(:)    = gdGARAbd(:) / REAL(nAccuCount)
          gdGAGPP(:)     = gdGAGPP(:)  / REAL(nAccuCount)
          gdGANPP(:)     = gdGANPP(:)  / REAL(nAccuCount)
          gdGARES(:)     = gdGARES(:)  / REAL(nAccuCount)
          gdGARESH(:)    = gdGARESH(:) / REAL(nAccuCount)
          gdGALIT(:)     = gdGALIT(:)  / REAL(nAccuCount)

          gdGACR(:)     = gdGACR(:)  / REAL(nAccuCount)
          gdGACA(:)     = gdGACA(:)  / REAL(nAccuCount)
          gdGACL(:)     = gdGACL(:)  / REAL(nAccuCount)
          gdGACWL(:)    = gdGACWL(:)  / REAL(nAccuCount)
          gdGACWR(:)    = gdGACWR(:)  / REAL(nAccuCount)
          gdGACVEG(:)    = gdGACVEG(:)  / REAL(nAccuCount)

          gdGALOSS_LIT_CS(:)    = gdGALOSS_LIT_CS(:)    / REAL(nAccuCount)
          gdGALOSS_LIT_CL(:)    = gdGALOSS_LIT_CL(:)    / REAL(nAccuCount)
          gdGALOSS_LIT_CR(:)    = gdGALOSS_LIT_CR(:)    / REAL(nAccuCount)
          gdGALOSS_LIT_CWL(:)   = gdGALOSS_LIT_CWL(:)   / REAL(nAccuCount)
          gdGALOSS_LIT_CWR(:)   = gdGALOSS_LIT_CWR(:)   / REAL(nAccuCount)
          gdGALOSS_LIT_CSLOW(:) = gdGALOSS_LIT_CSLOW(:) / REAL(nAccuCount)

          gdGACALOSS(:)  = gdGACALOSS(:)  / REAL(nAccuCount)
          gdGACLLOSS(:)  = gdGACLLOSS(:)  / REAL(nAccuCount)
          gdGACRLOSS(:)  = gdGACRLOSS(:)  / REAL(nAccuCount)
          gdGACWLLOSS(:) = gdGACWLLOSS(:) / REAL(nAccuCount)
          gdGACWRLOSS(:) = gdGACWRLOSS(:) / REAL(nAccuCount)

          gdGACLALLOC(:)  = gdGACLALLOC(:)  / REAL(nAccuCount)
          gdGACRALLOC(:)  = gdGACRALLOC(:)  / REAL(nAccuCount)
          gdGACWLALLOC(:) = gdGACWLALLOC(:) / REAL(nAccuCount)
          gdGACWRALLOC(:) = gdGACWRALLOC(:) / REAL(nAccuCount)
          gdGACSALLOC(:)  = gdGACSALLOC(:)  / REAL(nAccuCount)

          gdGAFH2O(:)    = gdGAFH2O(:)   / REAL(nAccuCount)
          gdGASOLRAD(:)  = gdGASOLRAD(:) / REAL(nAccuCount)
          gdGAFT(:)      = gdGAFT(:)     / REAL(nAccuCount)

          gdGAQTOT(:)      = gdGAQTOT(:)     / REAL(nAccuCount)
          gdGAET(:)      = gdGAET(:)     / REAL(nAccuCount)
          gdGAESOIL(:)      = gdGAESOIL(:)     / REAL(nAccuCount)
          gdGALEVAP(:)      = gdGALEVAP(:)     / REAL(nAccuCount)
          gdGAQSURF(:)      = gdGAQSURF(:)     / REAL(nAccuCount)
          gdGATRANS(:)      = gdGATRANS(:)     / REAL(nAccuCount)


          gdGAALB(:)     = gdGAALB(:)    / REAL(nAccuCount)
          gdGALAI(:)     = gdGALAI(:)    / REAL(nAccuCount)
          gdGAFVEG(:)    = gdGAFVEG(:)   / REAL(nAccuCount)
          gdGAWMAX(:)    = gdGAWMAX(:)   / REAL(nAccuCount)

          gdGAtau(:)     = gdGAtau(:)  / REAL(nAccuCount)
          gdGAGerm(:)    = gdGAGerm(:) / REAL(nAccuCount)
          gdGACol(:)     = gdGACol(:)  / REAL(nAccuCount)
          gdGAExcl(:)    = gdGAExcl(:) / REAL(nAccuCount)
          gdGAMort(:)    = gdGAMort(:) / REAL(nAccuCount)
          gdGADist(:)    = gdGADist(:) / REAL(nAccuCount)
          gdGASpec(:)    = gdGASpec(:) / REAL(nAccuCount)
          gdGAMig(:)     = gdGAMig(:)  / REAL(nAccuCount)
          gdGAExt(:)     = gdGAExt(:)  / REAL(nAccuCount)

        endif

!       ----------------------------------------------------------------
!       * change some units for CLAMP
!       ----------------------------------------------------------------

        gdGANPP(:)     = gdGAGPP(:) - gdGARES(:)  ! net primary productivity
        gdGANEE(:)     = gdGARESH(:) - gdGANPP(:) ! net ecosystem exchange
        gdGARESE(:)    = gdGARESH(:) + gdGARES(:) ! ecosystem respiration
        gdGALH(:)      = gdGAET(:) * cLambda      ! latent heat
        gdGASH(:)      = gdGAFT(:) - gdGALH(:)    ! sensible heat

!       ----------------------------------------------------------------
!       * output fields
!       ----------------------------------------------------------------

        if (nonepoint .eq. 1) then
          __globe_writeoutput(kFile_Grid,zgShannon,5371)
        else
          __globe_writeoutput(kFile_Grid,zgEvenness,5371)
        endif

        __globe_writeoutput(kFile_Grid,gdGAGPP,   kcode_5300)
        __globe_writeoutput(kFile_Grid,gdGANPP,   kcode_5301)
        __globe_writeoutput(kFile_Grid,gdGANEE,   kcode_5302)
        __globe_writeoutput(kFile_Grid,gdGARES,   kcode_5310)
        __globe_writeoutput(kFile_Grid,gdGARESH,  kcode_5311)
        __globe_writeoutput(kFile_Grid,gdGARESE,  5312)

        __globe_writeoutput(kFile_Grid,gdGAET,    kcode_5182)

        __globe_writeoutput(kFile_Grid,gdGAESOIL,    5183)
        __globe_writeoutput(kFile_Grid,gdGALEVAP,    5184)
        __globe_writeoutput(kFile_Grid,gdGAQSURF,    5185)
        __globe_writeoutput(kFile_Grid,gdGATRANS,    5186)
        __globe_writeoutput(kFile_Grid,gdGAQTOT,    5187)

        __globe_writeoutput(kFile_Grid,gdGAALB,   kcode_5174)
        __globe_writeoutput(kFile_Grid,gdGALAI,   kcode_5200)
        __globe_writeoutput(kFile_Grid,gdGAFVEG,  kcode_5199)
        __globe_writeoutput(kFile_Grid,gdGAWMAX,  kcode_5229)
        __globe_writeoutput(kFile_Grid,gdGAFH2O,  kcode_5379)
        __globe_writeoutput(kFile_Grid,gdGALH,    5380)
        __globe_writeoutput(kFile_Grid,gdGASH,    5381)
        __globe_writeoutput(kFile_Grid,gdGASOLRAD,5382)
        __globe_writeoutput(kFile_Grid,gdGAFT,    5383)

        __globe_writeoutput(kFile_Grid,grCSLOW,   5600)
        __globe_writeoutput(kFile_Grid,grLIT_CL,  5601)
        __globe_writeoutput(kFile_Grid,grLIT_CR,  5602)
        __globe_writeoutput(kFile_Grid,grLIT_CWL, 5603)
        __globe_writeoutput(kFile_Grid,grLIT_CWR, 5604)

        __globe_writeoutput(kFile_Grid,gdGACLLOSS,   5611)
        __globe_writeoutput(kFile_Grid,gdGACRLOSS,   5612)
        __globe_writeoutput(kFile_Grid,gdGACWLLOSS,  5613)
        __globe_writeoutput(kFile_Grid,gdGACWRLOSS,  5614)
        __globe_writeoutput(kFile_Grid,gdGACALOSS,   5616)

        __globe_writeoutput(kFile_Grid,gdGACLALLOC,  5621)
        __globe_writeoutput(kFile_Grid,gdGACRALLOC,  5622)
        __globe_writeoutput(kFile_Grid,gdGACWLALLOC, 5623)
        __globe_writeoutput(kFile_Grid,gdGACWRALLOC, 5624)
        __globe_writeoutput(kFile_Grid,gdGACSALLOC,  5625)

        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CSLOW,  5630)
        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CL,     5631)
        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CR,     5632)
        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CWL,    5633)
        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CWR,    5634)
        __globe_writeoutput(kFile_Grid,gdGALOSS_LIT_CS,     5635)

        if (kPop_dyn .eq. 2 ) then
          __globe_writeoutput(kFile_Grid,gdGACol,    kcode_5411)
          __globe_writeoutput(kFile_Grid,gdGAtau,    5409)
          __globe_writeoutput(kFile_Grid,gdGAGerm,   kcode_5410)
          __globe_writeoutput(kFile_Grid,gdGAExcl,   kcode_5412)
          __globe_writeoutput(kFile_Grid,gdGAMort,   kcode_5413)
          __globe_writeoutput(kFile_Grid,gdGADist,   kcode_5414)
          __globe_writeoutput(kFile_Grid,gdGAExt,    kcode_5415)
          if (kMig_dyn .eq. 1) then
            __globe_writeoutput(kFile_Grid,zdMigSeed,kcode_5402)
            __globe_writeoutput(kFile_Grid,zdMigOut, kcode_5403)
            __globe_writeoutput(kFile_Grid,zdMigIn,  kcode_5404)
            __globe_writeoutput(kFile_Grid,gdGAMig,  kcode_5416)
          endif
          if (kSpec_dyn .eq. 1) __globe_writeoutput(kFile_Grid,gdGASpec,5417)
        endif
      endif

      __globe_writeoutput(kFile_Grid,zW,      5140)
      __globe_writeoutput(kFile_Grid,zCS,     5372)
      __globe_writeoutput(kFile_Grid,gdGACA,     5373)
      __globe_writeoutput(kFile_Grid,gdGACL,     5374)
      __globe_writeoutput(kFile_Grid,gdGACR,     5375)
      __globe_writeoutput(kFile_Grid,gdGACWL,    5376)
      __globe_writeoutput(kFile_Grid,gdGACWR,    5377)
      __globe_writeoutput(kFile_Grid,gdGACVEG,kcode_5304)
      __globe_writeoutput(kFile_Grid,Richness,kcode_5370)

      if (kPop_dyn .ge. 1) then
        __globe_writeoutput(kFile_Grid,gdGARAbd,5400)
      endif

      if (kPop_dyn .eq. 2) then
        __globe_writeoutput(kFile_Grid,zExt,5401)
        if (kSpec_dyn .eq. 1) __globe_writeoutput(kFile_Grid,countSpec,5405)
      endif

!     * free the temporary space

      __deallocate(zCVEG)
      __deallocate(zW)
      __deallocate(zCS)
      __deallocate(zCA)
      __deallocate(zCL)
      __deallocate(zCR)
      __deallocate(zCWL)
      __deallocate(zCWR)
      __deallocate(zArea)
      __deallocate(zSuccess)
      __deallocate(zShannon)

      __deallocate(abdN)

      __deallocate(zgEvenness)
      __deallocate(zgShannon)
      __deallocate(zgCVEG)
      __deallocate(zgSuccess)
      __deallocate(zdMigSeed)
      __deallocate(zdMigIn)
      __deallocate(zdMigOut)
      __deallocate(zExt)

!     ------------------------------------------------------------------
!     * reset fields
!     ------------------------------------------------------------------

      call jedi_output_reset

!     ------------------------------------------------------------------
!     * end of output
!     ------------------------------------------------------------------

      __diag(kFile_Diag,'jedi_output: end')

      return
      end subroutine jedi_output

!     ******************************************************************
!     JEDI_DIAG
!     ******************************************************************

      subroutine jedi_diag ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP
      character(len=135) :: msg

      REAL,allocatable,dimension(:) :: zCVEG, zSuccess, zArea, zShannon
      REAL,allocatable,dimension(:) :: zgCVEG, zgSuccess, zgdGAGPP, zgArea, zgShannon

      REAL,allocatable,dimension(:) :: abdN

      __diag(kFile_Diag,'jedi_diag: start')

!     * allocate temporary fields

      __allocate(zCVEG,(NHOR))
      __allocate(zSuccess,(NHOR))
      __allocate(zArea,(NHOR))
      __allocate(zShannon,(NHOR))
      __allocate(zgCVEG,(gNHOR))
      __allocate(zgSuccess,(gNHOR))
      __allocate(zgdGAGPP,(gNHOR))
      __allocate(zgArea,(gNHOR))
      __allocate(zgShannon,(gNHOR))

      __allocate(abdN,(gNHOR))

      zCVEG(:)    = 0.0
      zSuccess(:) = 0.0
      zArea(:)    = 0.0
      zShannon(:) = 0.0

      if (nonepoint .eq. 1) then   ! single point version

        if (kPop_dyn .eq. 2) then
          zArea(:) = SUM(dRAbd(:,FirstSPP:kMaxSPP),dim=2)
          do iSPP = FirstSPP, kMaxSPP
            where (dRAbd(:,iSPP) .gt. cTiny)
              zCVEG(:) = zCVEG(:) + dRAbd(:,iSPP) * (rCA(:,iSPP) + rCWL(:,iSPP) + rCWR(:,iSPP) + rCL(:,iSPP) + rCR(:,iSPP))
            end where
            where (rArea(:,iSPP) .gt. cTiny)
              zSuccess(:) = zSuccess(:) + 1.0
            end where
            where (dRAbd(:,iSPP) .gt. 0.0 .and. zArea(:) .gt. 0.0)
              abdN(:) = dRAbd(:,iSPP) / zArea(:)
              zShannon(:) = zShannon(:) - abdN(:) * log(abdN(:))
            end where
          enddo
        endif

        if (nAccuCount .gt. 0) dGAGPP(:) = dGAGPP(:) / REAL(nAccuCount)

        write(msg,*) 'START: YEAR = ', OUTYEAR, 'MONTH = ', MONTH
        __diag(kFile_Diag,trim(msg))
        __diag_num(kFile_Diag,'   gdGAGPP  = ',dGAGPP(1))
        __diag_num(kFile_Diag,'   zCVEG    = ',zCVEG(1))
        __diag_num(kFile_Diag,'   gd rcs   = ',SUM(dSACS)/REAL(nAccuCount))
        __diag_num(kFile_Diag,'   zdRABD   = ',zArea(1))
        __diag_num(kFile_Diag,'   zbare    = ',rAreaBare(1))
        __diag_num(kFile_Diag,'   zSuccess = ',zSuccess(1))
        __diag_num(kFile_Diag,'   zShannon = ',zShannon(1))
        write(msg,*) 'END: YEAR = ', OUTYEAR, 'MONTH = ', MONTH
        __diag(kFile_Diag,trim(msg))
        if (kdiag == 1) flush(kFile_Diag)

      else     ! global version

        select case(kPop_dyn)
        ! case(0)
        !   do iSPP = FirstSPP, kMaxSPP
        !     where (dRAbd(:,iSPP) .gt. 0.0)
        !       zCVEG(:) = zCVEG(:) + dRAbd(:,iSPP) * (rCA(:,iSPP) + rCWL(:,iSPP) + rCWR(:,iSPP) + rCL(:,iSPP) + rCR(:,iSPP))
        !     end where
        !     where (rCS(:,iSPP) .ge. p14(iSPP) * pA0)
        !       zSuccess(:) = zSuccess(:) + 1.0
        !     end where
        !   enddo

        case(1)
          zCVEG(:) = SUM(dRAbd(:,:) * (rCA(:,:) + rCWL(:,:) + rCWR(:,:) + rCL(:,:) + rCR(:,:)), DIM=2)
          do iSPP = FirstSPP, kMaxSPP
            where (rCTot(:,iSPP) .gt. pA0)
              zSuccess(:) = zSuccess(:) + 1.0
            end where
          enddo
          zArea(:) = SUM(dRAbd(:,FirstSPP:kMaxSPP), dim=2)

        case(2)
          do iSPP = FirstSPP, kMaxSPP
            where (dRAbd(:,iSPP) .gt. 0.0)
              zCVEG(:) = zCVEG(:) + dRAbd(:,iSPP) * (rCA(:,iSPP) + rCWL(:,iSPP) + rCWR(:,iSPP) + rCL(:,iSPP) + rCR(:,iSPP))
            end where
            where (rArea(:,iSPP) .gt. 0.0)
              zSuccess(:) = zSuccess(:) + 1.0
            end where
          enddo
          zArea(:) = SUM(dRAbd(:,:), dim=2)
        end select

        __globe_mpga_to_from(zgCVEG,zCVEG)
        __globe_mpga_to_from(zgSuccess,zSuccess)
        __globe_mpga_to_from(zgArea,dGALAI)
        __globe_mpga_to_from(zgdGAGPP,dGAGPP)

        if (nAccuCount .gt. 0) zgdGAGPP(:) = zgdGAGPP(:) / REAL(nAccuCount)
        if (nAccuCount .gt. 0) zgArea(:)   = zgArea(:)   / REAL(nAccuCount)

        if (kmypid == kroot) then
          write(msg,*) 'START: YEAR = ', OUTYEAR, 'MONTH = ', MONTH
          __diag(kFile_Diag,trim(msg))
          __diag(kFile_Diag,'gdGAGPP =')
          call globe_plot_map(kFile_Diag, zgdGAGPP)
          __diag(kFile_Diag,'zCVEG =')
          call globe_plot_map(kFile_Diag, zgCVEG)
          __diag(kFile_Diag,'zRABD =')
          call globe_plot_map(kFile_Diag, zgArea)
          __diag(kFile_Diag,'zSuccess =')
          call globe_plot_map(kFile_Diag, zgSuccess)
          write(msg,*) 'END: YEAR = ', OUTYEAR, 'MONTH = ', MONTH, 'pPCO2 =', pPCO2
          __diag(kFile_Diag,trim(msg))
          flush(kFile_Diag)
        endif
      endif

!     * free temporary space

      __deallocate(zCVEG)
      __deallocate(zSuccess)
      __deallocate(zArea)
      __deallocate(zShannon)
      __deallocate(zgCVEG)
      __deallocate(zgSuccess)
      __deallocate(zgdGAGPP)
      __deallocate(zgArea)
      __deallocate(zgShannon)

      __deallocate(abdN)

!     * onepoint output

      if (kmypid == kroot .and. nonepoint .eq. 1) call jedi_output_spt

!     * end of output

      __diag(kFile_Diag,'jedi_diag: end')

      return
      end subroutine jedi_diag

