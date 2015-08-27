#define __ACTIVATE
#include "../globe/globe_macros.f90"

!
!     ------------------------------------------------------------------
!     SUBROUTINES FOR MIGRATION, SPECIATION & EXTINCTION IN JEDI V3
!     ------------------------------------------------------------------

!     ******************************************************************
!     JEDI_COLLECTSEEDS --- sum up seedflux for migration
!     ******************************************************************

      subroutine jedi_collectSeeds
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP

!     * seedbank decay for seeds to migrate
!     fseedout = MAX(0.0, fseedout * (1.0 - pSeedTau))

!     * calculate flux from landpoint to outside

      where (rLive .gt. 0.0)                      ! outflow just if species are living
        fseedout = fseedout + fseed * pmig !*p19  ! fseedout collects daily produced seeds that disperse
                                                  ! pmig: part of seedflux that migrate
        rCS      = rCS - (fseed * pmig)           ! seed pools rCS are actualiyed by leaving seeds
      end where

      return
      end subroutine jedi_collectSeeds

!     ******************************************************************
!     JEDI_MIGRATION --- disperse collected seeds to neigbors
!     ******************************************************************
!     subroutine for calculating migration from a cell to the 8 neigbored cells
!     collect information, do migration just once a month

      subroutine jedi_migration
      use jedi_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP, i, j, nnLpt

!     * once a month disperse collected seeds

      if (kmypid == kroot) gfseedout(:,:) = 0.0

      __globe_mpga_to_from(gfseedout,fseedout)  ! root collects outflow of seeds
      __globe_mpga_to_from(grCS,rCS)  ! root collects all seedpools rCS

      fseedout(:,:) = 0.0
      gfseedin(:,:) = 0.0

!     * root actualizes seed pool by adding the outflow of the neigbors

      if (kmypid == kroot) then
        do i = 1, gNHOR                ! landpoint i
          do j = 1, 8
            nnLpt = klpt_nn(i,j)       ! NN landpoint of landpoint i
            if (nnLpt .gt. 0) then     ! if neighbor is a landpoint
              do iSPP = FirstSPP, kMaxSPP
                grCS(i,iSPP) = grCS(i,iSPP) + 1./8. * gfseedout(nnLpt,iSPP) ! actualize seed pools, same seeds for 8 neigbors
                gfseedin(i,iSPP) = gfseedin(i,iSPP) + 1./8. * gfseedout(nnLpt,iSPP) ! store seedinflow
              enddo
            endif
          enddo
        enddo
      endif
      __globe_mpsc_from_to(grCS,rCS)  ! scatter actualized rC

      return
      end subroutine jedi_migration

!     ******************************************************************
!     JEDI_EXTINCTION
!     ******************************************************************

      subroutine jedi_extinction
      use jedi_mod
      use globe_mod
      use jedi_dyn_mod
      implicit none

      integer :: iSPP

      __globe_mpga_to_from(grLive,rLive)  ! rLive is criteria for extinction

      if (kmypid == kroot) then
        do iSPP = FirstSPP, kMaxSPP

          if ( (.not. (MAXVAL(grLive(:,iSPP)) .gt. -1.0 )) .and. (kSPPID(iSPP) .gt. 0)) then     !Extinct Species
            call jedi_write_extinct (iSPP)             ! write out extinct species to a file
            kSPPID(iSPP) = 0                           ! tag that extinct, not again write it out!

!           * add stored position of iSPP to FreePosition array
            if ((lastFree .eq. kMaxSPP)) then
              lastFree = 1
            else
              lastFree = lastFree + 1
            endif
            inFreePosSPP(lastFree) = iSPP
            count = count + 1                           ! number of free positions
            kNumSPP = kNumSPP - 1                       ! number of species alive
          endif
        enddo
      endif

      __globe_mpbc(kSPPID)

      return
      end subroutine jedi_extinction

!     ******************************************************************
!     JEDI_SPECIATION
!     ******************************************************************

      subroutine jedi_speciation
      use jedi_mod
      use jedi_dyn_mod
      use globe_mod
      implicit none

      integer :: iSPP, i
      integer :: winner, countSPP,fromS, toS
      integer :: newSPP
      real    :: gasdev, ran1, rand, var
      real    :: probsum, sumfseed

      var = pSpec
      countSPP = 0

      Afseed = Afseed + fseed

      if (mod(kdiy,30) .eq. 0) then ! Speciation
        __globe_mpga_to_from(gAfseed,Afseed)
        __globe_mpga_to_from(gfseed,fseed)
        __globe_mpga_to_from(grCS,rCS)
        __globe_mpga_to_from(grLive,rLive)

!       * set accumulated seedflux to zero
        Afseed = 0.0

        if (kmypid == kroot) then     ! root makes the speciation
          if (count .ge. gNHOR) then  ! if enought space for a new species in every grispoint
            do i=1, gNHOR
              winner = -1
              sumfseed = SUM(gfseed(i,FirstSPP:kMaxSPP)) ! to norm seedflux to 1.0
              if (sumfseed .gt. 0.0) then         ! avoid division by zero
                probsum = 0.0
                rand = ran1(rand_seed)
                do iSPP = FirstSPP, kMaxSPP       ! winner chosen by roulett
                  if (gfseed(i,iSPP) .gt. 0.0) then
                    probsum = probsum + gfseed(i,iSPP) / sumfseed
                    if (probsum .gt. rand) then
                      winner=iSPP
                      exit
                    endif
                  endif
                enddo
              endif

              if ( winner .gt. 0) then
!             if ( gfseed(i,winner) .gt. 0.0 ) then
                countSPP = countSPP + 1

!               * look up where to store parameter information of new Species
                  newSPP = inFreePosSPP(firstFree) ! index of where to information about new species
                  if (newSPP .gt. kMaxSPP) exit    ! species array is full, should never happen

!               * remove this index from list of free positions
                if ((firstFree .eq. kMaxSPP)) then
                  firstFree = 1
                else
                  firstFree = firstFree + 1
                endif
                kNumSPPID = kNumSPPID + 1                       ! ID of this Species, and number of total number of species until now
                kSPPID(newSPP) = kNumSPPID                      ! set ID to newSPP
                count = count - 1                               ! one less free position for new species
                kNumSPP = kNumSPP + 1                           ! one more existing species at the moment
                call jedi_transform_parms_back(winner)

                p01(newSPP) = MIN( 1.0, (p01(winner) + gasdev(rand_seed,0.0,var)) )
                if (p01(newSPP) .le. 0.0 ) p01(newSPP) = cTiny
                p02(newSPP) = MIN( 1.0, (p02(winner) + gasdev(rand_seed,0.0,var)) )
                if (p02(newSPP) .le. 0.0 ) p02(newSPP) = cTiny
                p05(newSPP) = MIN( 1.0, (p05(winner) + gasdev(rand_seed,0.0,var)) )
                if (p05(newSPP) .le. 0.0 ) p05(newSPP) = cTiny
                p06(newSPP) = MIN( 1.0, (p06(winner) + gasdev(rand_seed,0.0,var)) )
                if (p06(newSPP) .le. 0.0 ) p06(newSPP) = cTiny
                p07(newSPP) = MIN( 1.0, (p07(winner) + gasdev(rand_seed,0.0,var)) )
                if (p07(newSPP) .le. 0.0 ) p07(newSPP) = cTiny
                p08(newSPP) = MIN( 1.0, (p08(winner) + gasdev(rand_seed,0.0,var)) )
                if (p08(newSPP) .le. 0.0 ) p08(newSPP) = cTiny
                p09(newSPP) = MIN( 1.0, (p09(winner) + gasdev(rand_seed,0.0,var)) )
                if (p09(newSPP) .le. 0.0 ) p09(newSPP) = cTiny
                p10(newSPP) = MIN( 1.0, (p10(winner) + gasdev(rand_seed,0.0,var)) )
                if (p10(newSPP) .le. 0.0 ) p10(newSPP) = cTiny
                p13(newSPP) = MIN( 1.0, (p13(winner) + gasdev(rand_seed,0.0,var)) )
                if (p13(newSPP) .le. 0.0 ) p13(newSPP) = cTiny
                p16(newSPP) = MIN( 1.0, (p16(winner) + gasdev(rand_seed,0.0,var)) )
                if (p16(newSPP) .le. 0.0 ) p16(newSPP) = cTiny
                p14(newSPP) = MIN( 1.0, (p14(winner) + gasdev(rand_seed,0.0,var)) )
                if (p14(newSPP) .le. 0.0 ) p14(newSPP) = cTiny
                p15(newSPP) = MIN( 1.0, (p15(winner) + gasdev(rand_seed,0.0,var)) )
                if (p15(newSPP) .le. 0.0 ) p15(newSPP) = cTiny
                p04(newSPP) = MIN( 1.0, (p04(winner) + gasdev(rand_seed,0.0,var)) )
                if (p04(newSPP) .le. 0.0 ) p04(newSPP) = cTiny
                p17(newSPP) = MIN( 1.0, (p17(winner) + gasdev(rand_seed,0.0,var)) )
                if (p17(newSPP) .le. 0.0 ) p17(newSPP) = cTiny
                p12(newSPP) = MIN( 1.0, (p12(winner) + gasdev(rand_seed,0.0,var)) )
                if (p12(newSPP) .le. 0.0 ) p12(newSPP) = cTiny
                p11(newSPP) = MIN( 1.0, (p11(winner) + gasdev(rand_seed,0.0,var)) )
                if (p11(newSPP) .le. 0.0 ) p11(newSPP) = cTiny
                p03(newSPP) = MIN( 1.0, (p03(winner) + gasdev(rand_seed,0.0,var)) )
                if (p03(newSPP) .le. 0.0 ) p03(newSPP) = cTiny
                p18(newSPP) = MIN( 1.0, (p18(winner) + gasdev(rand_seed,0.0,var)) )
                if (p18(newSPP) .le. 0.0 ) p18(newSPP) = cTiny
                p19(newSPP) = MIN( 1.0, (p19(winner) + gasdev(rand_seed,0.0,var)) )
                if (p19(newSPP) .le. 0.0 ) p19(newSPP) = cTiny
                p20(newSPP) = MIN( 1.0, (p20(winner) + gasdev(rand_seed,0.0,var)) )
                if (p20(newSPP) .le. 0.0 ) p20(newSPP) = cTiny

                call jedi_transform_parms(newSPP)
                call jedi_transform_parms(winner)

!               * give the new mutated SPP some initial seedmass ; need to give inital conditions?
                grCS(i,winner) = grCS(i,winner) - 0.5 * gfseed(i,winner)
                grCS(i,newSPP) = 0.5 * gfseed(i,winner)
                grLive(i,newSPP) = 0.0
                countSpec(i) = countSpec(i) + 1.0
              endif
            enddo
          endif
        endif

!       * broadcast spe array to all procs
!       __globe_mpbc(kNumSPP)
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

        __globe_mpsc_from_to(grCS,rCS)   ! scatter actualized rCS
        __globe_mpsc_from_to(grLive,rLive)
      endif

      return
      end subroutine jedi_speciation

!     ******************************************************************
!     JEDI_DISTURBANCE --- just root does, roots knows NN
!     ******************************************************************

      subroutine jedi_disturbance
      use jedi_mod
      use jedi_dyn_mod
      use globe_mod
      implicit none

      real    :: ran1

      integer :: i, iSPP
      real    :: fdist, dfAD, sumdfAD, wkdist

      wkdist = 1.0 / (pDistTau * 12.0)
      dfAD = 0.0

      do i = 1, NHOR
        dfAreaM(:) = 0.0
        sumdfAD = 0.0
        fdist   = ran1(rand_seed)
        if (fdist .lt. wkdist) then
!          fLIT(i,:)    = fLIT(i,:) + rCA(i,:) + rCL(i,:) + rCR(i,:) + rCWR(i,:) + rCWL(i,:)
!          rCA(i,:)     = 0.0
!          rCL(i,:)     = 0.0
!          rCR(i,:)     = 0.0
!          rCWL(i,:)    = 0.0
!          rCWR(i,:)    = 0.0
!          rDie(i,:)    = 0.0
!          rLive(i,:)   = 0.0
          dfAreaM(:)   = rArea(i,:) * 1.0

          rArea(i,:)   = rArea(i,:) - dfAreaM(:)
          rAreaBare(i) = rAreaBare(i) + SUM(dfAreaM(:))

          dGADist(i)   = dGADist(i) + 1.0
          __diag(kFile_Diag,'Dist')
        endif
      enddo

      return
      end subroutine jedi_disturbance

!     ==================================================================
!     FUNCTION FOR RANDOM from Numerical Recipies
!     ==================================================================

      function ran1(idum)
      implicit none

      integer, parameter :: K4B = selected_int_kind(9)
      integer(K4B), intent(INOUT) :: idum
      real :: ran1

      ! minimal random number generator of park and miller ...
      ! returns a uniform random deviate between 0.0 und 1.0
      ! call with idum a negetive integr to initialize; thereafter, do not alter idum
      ! except to reinitialize. Period: 3.1 x 18^18

      integer(K4B),parameter ::ia=16807, im=2147483647, iq=127773, ir=2836
      real, save ::am
      integer(K4B),save :: ix=-1, iy=-1,k

      if (idum <= 0 .or. iy < 0) then             ! Init
        am = nearest(1.0, -1.0) / im
        iy = ior(ieor(888889999, ABS(idum)), 1)
        ix = ieor(777755555, ABS(idum))
        idum = ABS(idum) + 1                      ! set idum positive
      endif

      ix = ieor(ix, ishft(ix, 13))                ! Marsaglia shift sequence with period 2^32 -1
      ix = ieor(ix, ishft(ix, -17))
      ix = ieor(ix, ishft(ix, 5))
      k = iy / iq                                 ! Park-Miller sequence by method of Schrage, period 2^31 -2
      iy = ia * (iy - k * iq) - ir * k
      if (iy < 0) iy = iy + im
      ran1 = am * ior(iand(im, ieor(ix, iy)), 1)  ! combine two generators wth masking to ensure nonzero value

      end function ran1

!     ==================================================================
!     FUNCTION FOR GAUSSDISTRIBUTED RANDOM from Numerical Recipies
!     ==================================================================

      function gasdev(idum, erw, var)

      integer :: idum
      real    :: gasdev, erw, var
      integer :: iset
      real    :: fac, gset, rsq, v1, v2, ran
      save    :: iset, gset
      data    iset/0/

      if (idum .eq. 0) iset = 0
      if (iset .eq. 0) then
1       v1 = 2.0 * ran1(idum) - 1
        v2 = 2.0 * ran1(idum) - 1
        rsq = v1**2 + v2**2
        if (rsq .ge. 1 .or. rsq .eq. 0 ) goto 1

        fac = sqrt(-2.0 * log(rsq) / rsq)
        gset = v1 * fac
        gasdev = v2 * fac
        iset = 1
      else
        gasdev = gset
        iset = 0
      endif
      gasdev = (gasdev * sqrt(var)) + erw

      return
      end function gasdev
