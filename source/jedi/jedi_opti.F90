#define __ACTIVATE
#include "../globe/globe_macros.f90"

!     extension to JEDI to optimize species parameters to maximize
!     productivity
!
!     Axel Kleidon, 28.06.2008
!
!     VARIABLES:
!     - kOpti:         switch for optimization run
!     - pOptiDropFrac: lower fraction of ranked spp to drop
!     - kOptiYears:    number of years to average over for sorting of productivity
!     - dOptiGoalAcc:  goal function accumulated for each species
!

!     ******************************************************************
!     JEDI_OPTI_INIT
!     ******************************************************************

      subroutine jedi_opti_init ()
      use jedi_mod
      implicit none

      integer :: iostat

      __diag(kFile_Diag,'jedi_opti_init')

      dOptiGoalAcc(:,:) = 0.

      if (kmypid == kroot) then

        open(kFile_Opti, FILE=sfile_Opti, STATUS='replace',            &
     &                   FORM='formatted', IOSTAT=iostat)
        call error_open_new(iostat,sfile_Opti)
        write(kFile_Opti,'(5A15)') 'nOptiCnt', 'zGoalSum', 'zGoalMax', 'nBestCnt', 'nReplaceCnt'
        close(kFile_Opti)

        call globe_open_output(sfile_OptiOut, kFile_OptiOut, kFile_Diag)

        open(kFile_OptiSPP, FILE=sfile_OptiSPP, STATUS='replace',      &
     &                      FORM='formatted', IOSTAT=iostat)
        call error_open_new(iostat,sfile_OptiSPP)
        close(kFile_OptiSPP)
      endif

      return
      end subroutine jedi_opti_init

!     ******************************************************************
!     JEDI_OPTI_FILTER
!     ******************************************************************

      subroutine jedi_opti_filter ()
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      real       :: zGoalSum, zGoalMax, zGoal

      integer    :: i, j, k, iSPP, iSPP2, iSPPexist
      integer    :: kReplace
      integer    :: nBestCnt, nReplaceCnt

      real       :: jedi_opti_ran1
      integer    :: iostat
      character(len=135) :: msg
      
      integer, allocatable :: inUsedPosSPP(:)

      allocate(inUsedPosSPP(kMaxSPP))

      __diag(kFile_Diag,'jedi_opti_filter')

      nOptiCnt = nOptiCnt + 1

!     * gather accumulated goal function value
      __globe_mpga_to_from(gdOptiGoalAcc,dOptiGoalAcc)

!     * run evaluation on root only
      if (kmypid == kroot) then

!       * sort species by their productivities
        write(*,*) 'jedi_opti_filter: sorting...'

        i = FirstSPP
        do iSPP = FirstSPP, kMaxSPP
          inUsedPosSPP(i) = iSPP
          i = i + 1
          if (kSpec_dyn .gt. 0) then
            if (lastFree .ge. firstFree) then
              if (lastFree .ne. 0) then
                do j = firstFree, lastFree
                  if (inFreePosSPP(j) == iSPP) then
                    i = i - 1
                  endif
                enddo
              endif
            else
              do j = 1, lastFree
                if (inFreePosSPP(j) == iSPP) then
                  i = i - 1
                endif
              enddo
              do j = firstFree, kMaxSPP
                if (inFreePosSPP(j) == iSPP) then
                  i = i - 1
                endif
              enddo
            endif
          endif
        enddo

        if (kSpec_dyn .gt. 0) then
          if (inUsedPosSPP(kNumSPP) .eq. kMaxSPP .and. i .ne. kNumSPP) then
            write(*,*) "ERROR: Wrong number of living plants"
            stop
          endif
          if (inUsedPosSPP(kNumSPP) .ne. kMaxSPP .and. i-1 .ne. kNumSPP) then
            write(*,*) "ERROR: Wrong number of living plants"
            stop
          endif
        endif

        do i = FirstSPP, kNumSPP
          kOptiSorted(:,i) = inUsedPosSPP(i)
        enddo

        do i = FirstSPP, kNumSPP
          do j = i+1, kNumSPP
            do k = 1, gNHOR
              iSPP  = kOptiSorted(k, i)
              iSPP2 = kOptiSorted(k, j)
              if (gdOptiGoalAcc(k,iSPP) .lt. gdOptiGoalAcc(k,iSPP2)) then
                kOptiSorted(k,i) = iSPP2
                kOptiSorted(k,j) = iSPP
              endif
            enddo
          enddo
        enddo

!       * find maximum

        gdOptiGoalMax(:) = 0.0

        do k = 1, gNHOR
          i = FirstSPP
          do iSPP = FirstSPP, kMaxSPP
            if (inUsedPosSPP(i) == iSPP) then
              gdOptiGoalMax(k) = MAX(gdOptiGoalMax(k), gdOptiGoalAcc(k,iSPP))
              i = i + 1
            endif
          enddo
        enddo

!       * diagnostics of goal function
        zGoalSum = SUM(gdOptiGoalAcc(:,inUsedPosSPP(FirstSPP:kNumSPP))) / real(gNHOR)
        __diag_num(kFile_Diag,'zGoalSum = ',zGoalSum/real(kOptiYears))

        zGoalMax = SUM(gdOptiGoalMax(:)) / real(gNHOR) / real(kOptiYears) / 365.0
        __diag_num(kFile_Diag,'zGoalMax = ',zGoalMax)

!       * generate list which marks those species that are best at at least
!         one grid cell to save these from zeroing

        __diag(kFile_Diag,'jedi_opti_filter: create list with best species...')

        kOptiBest(:) = 0
        do k = 1, gNHOR
          iSPP            = kOptiSorted(k,FirstSPP)
          kOptiBest(iSPP) = kOptiBest(iSPP) + 1
          gOptiBestMap(k) = real(iSPP)
        enddo

        nBestCnt     = 0
        do iSPP = FirstSPP, kMaxSPP
          if (kOptiBest(iSPP) .gt. 0) then
            nBestCnt = nBestCnt + 1
          endif
        enddo

        __diag_num(kFile_Diag,'nBestCnt = ',nBestCnt)

        do iSPP = FirstSPP, kMaxSPP
          if (kOptiBest(iSPP) .gt. 0) then
            write(msg,*) 'spp ', iSPP, ' best at ', kOptiBest(iSPP), &
     &                   ' gpts', ' with goal=', SUM(gdOptiGoalAcc(:,iSPP))
            __diag(kFile_Diag,trim(msg))
          endif
        enddo

!       * go through sorted list and randomly set productivity to zero
!         with probability = goal/max(goal) of the lowest pOptiFrac fraction.

        write(*,*) 'jedi_opti_filter: zeroing...'

        if (pOptiFrac .gt. 0.) then
          do i = FirstSPP, kNumSPP
            do k = 1, gNHOR
              iSPP  = kOptiSorted(k, i)
!              if ((real(i)/real(kMaxSPP) .gt. jedi_opti_ran1(kOptiSeed)) .and. (kOptiBest(iSPP) == 0)) then
              if (((gdOptiGoalAcc(k, iSPP)/gdOptiGoalMax(k)/pOptiFrac) .lt. jedi_opti_ran1(kOptiSeed))) then
                gdOptiGoalAcc(k, iSPP) = 0.0
              endif
            enddo
          enddo
        endif

!       * sort species by their global productivities.  this list is used to
!         take attributes for new species

        write(*,*) 'jedi_opti_filter: global sorting...'

        do i = FirstSPP, kNumSPP
          kOptiSortedGlobal(i) = inUsedPosSPP(i)
        enddo

        do i = FirstSPP, kNumSPP
          do j = i+1, kNumSPP
            iSPP  = kOptiSortedGlobal(i)
            iSPP2 = kOptiSortedGlobal(j)
            if (SUM(gdOptiGoalAcc(:, iSPP)) .lt. SUM(gdOptiGoalAcc(:, iSPP2))) then
              kOptiSortedGlobal(i) = iSPP2
              kOptiSortedGlobal(j) = iSPP
            endif
          enddo
        enddo

!       * if species has zero productivity at all grid points, replace by
!         new species

        nReplaceCnt  = 0
        do i = FirstSPP, kNumSPP
          iSPP = inUsedPosSPP(i)
          kReplace = 1
          do k = 1, gNHOR
            if ((gdOptiGoalAcc(k, iSPP) .gt. 0.0) .or. (pOptiFrac .eq. 0.0)) kReplace = 0
          enddo

          if (kReplace == 1) then

!           * use namelist parameter pOptiFracRan as probability to create
!             a new species by (a) entirely random, or (b) from one of the best species

            call jedi_transform_parms_back(iSPP)

            if (pOptiFracRan .lt. jedi_opti_ran1(kOptiSeed)) then

!             * replace by random new species

              p01(iSPP) = jedi_opti_ran1(kOptiSeed)
              p02(iSPP) = jedi_opti_ran1(kOptiSeed)
              p05(iSPP) = jedi_opti_ran1(kOptiSeed)
              p06(iSPP) = jedi_opti_ran1(kOptiSeed)
              p07(iSPP) = jedi_opti_ran1(kOptiSeed)
              p08(iSPP) = jedi_opti_ran1(kOptiSeed)
              p09(iSPP) = jedi_opti_ran1(kOptiSeed)
              p10(iSPP) = jedi_opti_ran1(kOptiSeed)
              p13(iSPP) = jedi_opti_ran1(kOptiSeed)
              p16(iSPP) = jedi_opti_ran1(kOptiSeed)
              p14(iSPP) = jedi_opti_ran1(kOptiSeed)
              p15(iSPP) = jedi_opti_ran1(kOptiSeed)
              p04(iSPP) = jedi_opti_ran1(kOptiSeed)
              p17(iSPP) = jedi_opti_ran1(kOptiSeed)
              p12(iSPP) = jedi_opti_ran1(kOptiSeed)
              p11(iSPP) = jedi_opti_ran1(kOptiSeed)
              p03(iSPP) = jedi_opti_ran1(kOptiSeed)
              p18(iSPP) = jedi_opti_ran1(kOptiSeed)
              p19(iSPP) = jedi_opti_ran1(kOptiSeed)
              p20(iSPP) = jedi_opti_ran1(kOptiSeed)

            else

!             * vary existing species

              iSPPexist    = (nBestCnt-1) * jedi_opti_ran1(kOptiSeed) + 1 + FirstSPP
              iSPPexist    = kOptiSortedGlobal(iSPPexist)

              p01(iSPP)  = p01(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p02(iSPP)  = p02(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p05(iSPP)  = p05(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p06(iSPP)  = p06(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p07(iSPP)  = p07(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p08(iSPP)  = p08(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p09(iSPP)  = p09(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p10(iSPP)  = p10(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p13(iSPP)  = p13(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p16(iSPP)  = p16(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p14(iSPP)  = p14(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p15(iSPP)  = p15(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p04(iSPP)  = p04(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p17(iSPP)  = p17(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p12(iSPP)  = p12(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p11(iSPP)  = p11(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p03(iSPP)  = p03(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p18(iSPP)  = p18(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p19(iSPP)  = p19(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))
              p20(iSPP)  = p20(iSPPexist) * (1. + 2. * pOptiVar * (jedi_opti_ran1(kOptiSeed) - 0.5))

!             * make sure parameters do not turn negative!

              p01(iSPP)  = MAX(0., p01(iSPP))
              p02(iSPP)  = MAX(0., p02(iSPP))
              p05(iSPP)  = MAX(0., p05(iSPP))
              p06(iSPP)  = MAX(0., p06(iSPP))
              p07(iSPP)  = MAX(0., p07(iSPP))
              p08(iSPP)  = MAX(0., p08(iSPP))
              p09(iSPP)  = MAX(0., p09(iSPP))
              p10(iSPP)  = MAX(0., p10(iSPP))
              p13(iSPP)  = MAX(0., p13(iSPP))
              p16(iSPP)  = MAX(0., p16(iSPP))
              p14(iSPP)  = MAX(0., p14(iSPP))
              p15(iSPP)  = MAX(0., p15(iSPP))
              p04(iSPP)  = MAX(0., p04(iSPP))
              p17(iSPP)  = MAX(0., p17(iSPP))
              p12(iSPP)  = MAX(0., p12(iSPP))
              p11(iSPP)  = MAX(0., p11(iSPP))
              p03(iSPP)  = MAX(0., p03(iSPP))
              p18(iSPP)  = MAX(0., p18(iSPP))
              p19(iSPP)  = MAX(0., p19(iSPP))
              p20(iSPP)  = MAX(0., p20(iSPP))

            endif

            call jedi_transform_parms(iSPP)

!           * reset pools

            call jedi_init_fields(iSPP)

            nReplaceCnt = nReplaceCnt + 1

          endif
        enddo

        write(*,*) 'nReplaceCnt = ', nReplaceCnt

!       * write optimization log

        open(kFile_Opti, FILE=sfile_Opti, STATUS='old', IOSTAT=iostat, &
     &                   FORM='formatted', POSITION='append')
        write(kFile_Opti,'(I15,2E15.6,2I15)') nOptiCnt, zGoalSum, zGoalMax, nBestCnt, nReplaceCnt
        close(kFile_Opti)

!       * write best species and value of goal function

        __globe_writeoutput(kFile_OptiOut,gdOptiGoalMax,kcode_goalfct)
        __globe_writeoutput(kFile_OptiOut,gOptiBestMap,kcode_optispp)
        flush(kFile_OptiOut)

!       * write species properties file

        open(kFile_OptiSPP, FILE=sfile_OptiSPP, STATUS='old', IOSTAT=iostat, &
     &                      FORM='formatted', POSITION='append')
        write(kFile_OptiSPP,*) '=== species of optimization step ', nOptiCnt, ' ==='
        write(kFile_OptiSPP,*) nOptiCnt, nBestCnt

        do iSPP = FirstSPP, kMaxSPP
          if (kOptiBest(iSPP) .gt. 0) then
            call jedi_transform_parms_back(iSPP)
            write (kFile_OptiSPP,'(i8,20f6.3)') iSPP,            &
     &            p01(iSPP), p02(iSPP), p05(iSPP), p06(iSPP), p07(iSPP), &
     &            p08(iSPP), p09(iSPP), p10(iSPP), p13(iSPP), p16(iSPP), &
     &            p14(iSPP), p15(iSPP), p04(iSPP), p17(iSPP), p12(iSPP), &
     &            p11(iSPP), p03(iSPP), p18(iSPP), p19(iSPP), p20(iSPP)
            call jedi_transform_parms(iSPP)
          endif
        enddo
        close(kFile_OptiSPP)

      endif

!     * reset accumulation field

      if (kmypid == kroot) write(*,*) 'jedi_opti_filter: reset accumulation field'

      dOptiGoalAcc(:,:) = 0.

      if (kmypid == kroot) write(*,*) 'jedi_opti_filter: done'

      deallocate(inUsedPosSPP)

      return
      end subroutine jedi_opti_filter

!     ******************************************************************
!     JEDI_OPTI_STOP
!     ******************************************************************

      subroutine jedi_opti_stop ()
      use jedi_mod
      implicit none

      __diag(kFile_Diag,'jedi_opti_stop')

      if (kmypid == kroot) then
        close(kFile_Opti)
        call globe_close_output(kFile_OptiOut, kFile_Diag)
        close(kFile_OptiSPP)
      endif

      return
      end subroutine jedi_opti_stop

!     ==================================================================
!     JEDI_OPTI_RAN1
!     ==================================================================

!     * random number generator, from Numerical Recipes

      FUNCTION jedi_opti_ran1(idum)
      INTEGER      idum, IA, IM, IQ, IR, NTAB, NDIV
      REAL         jedi_opti_ran1, AM, EPS, RNMX

      PARAMETER   (IA=16807, IM=2147483647, AM=1./IM, &
     &             IQ=127773, IR=2836, NTAB=32, &
     &             NDIV=1+(IM-1)/NTAB, EPS=1.2e-7, &
     &             RNMX=1.-EPS)

      INTEGER      j, k, iv(NTAB), iy
      SAVE         iv, iy
      DATA         iv /NTAB*0/, iy /0/

      IF (idum .LE.0 .OR. iy .EQ. 0) THEN

        idum = MAX(-idum, 1)

        DO j = NTAB+8, 1, -1
          k    = idum/IQ
          idum = IA * (idum - k * IQ) - IR * k
          IF (idum .lt. 0) idum = idum + IM
          IF (j .LE. NTAB) iv(j)= idum
        ENDDO

        iy   = iv(1)

      ENDIF

      k    = idum/IQ
      idum = IA * (idum - k * IQ) - IR * k
      IF (idum .lt. 0) idum = idum + IM
      j    = 1 + iy/NDIV
      iy   = iv(j)
      iv(j)= idum
      jedi_opti_ran1 = MIN(AM*iy, RNMX)

      RETURN
      END
