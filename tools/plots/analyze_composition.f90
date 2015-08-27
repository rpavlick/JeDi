!
!     A N A L Y Z E _ C O M P O S I T I O N
!
!     purpose:
!     ========
!
!     this program takes JeDi model output of successful species and uses it
!     to derive:
!
!     - richness-evenness relationship
!     - relative abundance distribution
!     - ranked abundance distribution
!     - relative productivity distribution
!     - ranked productivity distribution
!
!     for a given number of richness classes.  
!
!     files:
!     ======
!
!     jedi_species.txt           input file from jedi.
!                                text file, fixed-format, unit 10
!
!     analyze-evenness.txt       output file containing richness and evenness 
!                                for each grid point.  
!                                text file, tab-delimited, unit 20
!
!     analyze-relabd_spp_dX.txt  output file containing binned relative 
!                                abundances for a given richness class X.  The
!                                first column contains the upper bin boundary
!                                value (relative abundance), the second the 
!                                normalized number of species with the given
!                                abundance. 
!                                text file, tab-delimited, unit 21
!
!     analyze-rank_spp_dX.txt    output file containing ranked relative 
!                                abundances for a given richness class X.  The
!                                first column contains the absolute rank value,
!                                the second its abundance.  
!                                text file, tab-delimited, unit 23
!
!     analyze-rank_mostabd.txt   output file containing the most abundant 
!                                species and its relative abundances for all
!                                grid points.  The first column contains
!                                the grid point index, second contains the
!                                species id, the third contains the abundance
!                                text file, tab-delimited, unit 24
!                                
!     analyze-relabd_gpp_dX.txt  output file containing binned frequencies of
!                                productivity values for a given richness class
!                                X.  Not sure if this is really useful, the low
!                                diversity values screw up the histogram...                    
!                                text file, tab-delimited, unit 25
!
!     analyze-rank_gpp_dX.txt    output file containing species ranked according
!                                to their normalized productivity for a given
!                                richness class X.  The first column contains 
!                                the rank, the second the normalized 
!                                productivity relative to the maximum found at
!                                each grid cell.
!                                text file, tab-delimited, unit 26
!
!     history:
!     ========
!     version 1         28.03.2008, Axel Kleidon
!

      program analyze_composition
      
      implicit     none

!     * potential array sizes

      parameter    nMaxSPP     = 4000
      parameter    nMaxGPts    = 11000

!     * number of bins for binning and bin limits

      parameter    nNumBins    = 30
      real         pBinHigh(0:nNumBins)

      parameter    nNumGPPBins = 120
      real         pGPPBinHigh(0:nNumGPPBins)

!     * number of classes for richness classification
!       and class limits

      parameter    nNumClasses = 4
      real         pClassHigh(0:nNumClasses)

!     * actual array sizes

      integer      nNumSPP      
      integer      nNumGPts
      integer      nSamples

!     * data read in from jedi output

      integer      kLon, kLat
      integer      kSPPIDtemp
      real         zRAbd, zGPP, zCS, zCA, zCL, zCR, zCWL, zCWR
      integer      kLonOld, kLatOld
      
      real         pRelAbd(nMaxGPts, nMaxSPP)
      real         pGPP(nMaxGPts, nMaxSPP)
      real         pCS(nMaxGPts, nMaxSPP)
      real         pCTOT(nMaxGPts, nMaxSPP)
      real         pCA(nMaxGPts, nMaxSPP)
      integer      kSPPID(nMaxGPts, nMaxSPP)
      integer      kSPPCnt(nMaxGPts)
      
!     * loop variables and indices

      integer      i, j, k, i2, j2
      integer      kGridPt

!     * fields of computed richness and evenness
!       used to create file unit 20, and for evaluations below

      real         pRichness(nMaxGPts)
      real         pEvenness(nMaxGPts)

!     * fields to record abundance histograms
!       (normalized by richness of class)
!       used to create file unit 21

      real         pFreqDist(nNumBins)
      real         pClassRichness

!     * fields for ranking of species by abundance
!       used to create file unit 23

      real         pRankedAbd(nMaxSPP)
      real         pRankTemp
      integer      kRankSPP(nMaxSPP)
      real         kRankSPPTemp
      real         pSPPCount2(nMaxGPts)
      real         pAccRankAbd(nMaxSPP)

!     * fields to find most abundant species and its abundance
!       used to create file unit 24

      integer      kMostAbdSPP(nMaxGPts)
      real         pMostAbdSPPAbd(nMaxGPts)

!     * fields to record productivity histograms
!       (normalized by max. productivity of each grid point)
!       used to create file unit 25

      real         pGPPMax

!     * fields for ranking according to productivity
!       used to create file unit 26

      real         pRankedGPP(nMaxSPP)
      real         pAccRankGPP(nMaxSPP)

!     * fields for carbon allocation analysis

      real         pResTime(nMaxGPts, nMaxSPP)
      real         pCSTime(nMaxGPts, nMaxSPP)

      real         pSTMax, pRTMax

      real         pGPGPPMean
      real         pGPSTMean
      real         pGPRTMean
      
      integer      kCnt      
      
!     * variables for file name handling

      character*1  sText(0:9)
      data         sText/'0','1','2','3','4','5','6','7','8','9'/
      
!     ================================================================      

      write(*,*) 'analyze_composition: start'

!     ----------------------------------------------------------------
!     * initialization
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: create classes and bins'

!     * create classes

      pClassHigh(0)    = 0.0
      pClassHigh(1)    = 0.25
      pClassHigh(2)    = 0.5
      pClassHigh(3)    = 0.75
      pClassHigh(4)    = 1.0

!     * create bins

      pBinHigh(0)      = 10**(-4. + 4. * real(-1)/(real(nNumBins)-1.))
      do i2 = 1, nNumBins
        pBinHigh(i2)   = 10**(-4. + 4. * real(i2-1)/(real(nNumBins)-1.))
      end do

      pGPPBinHigh(0)    = 0.
!      pGPPBinHigh(0)   = 10**(-3. + 3. * real(-1)/(real(nNumGPPBins)-1.))
      do i2 = 1, nNumGPPBins
        pGPPBinHigh(i2) = real(i2)/real(nNumGPPBins)
!        pGPPBinHigh(i2)   = 10**(-3. + 3. * real(i2-1)/(real(nNumBins)-1.))
      end do
      
!     ----------------------------------------------------------------
!     * read-in KM2000/JeDi 1 output from Bjoern's ISOMAP matrix
!       modified with an integer header value with the number of 
!       entries in one row
!     ----------------------------------------------------------------

!      open(10, file='diva_SpecAbd.txt')
!      read(10,*) nNumSPP
!      write(*,*) 'analyze_composition: nNumSPP in diva_SpecAbd.txt = ',nNumSPP
!      read(10,*) (kSPPID(i), i = 1, nNumSPP)
!      kGridPt        = 1
!  100 continue
!      read(10,*,end=149) (pRelAbd(kGridPt,j), j=1, nNumSPP)
!      kGridPt        = kGridPt + 1
!      goto 100
!  149 continue
!      nNumGPts       = kGridPt
!      close(10)
!
!      open(11, file='diva_GPP.txt')
!      read(11,*) nNumSPP
!      write(*,*) 'analyze_composition: nNumSPP in diva_GPP.txt = ',nNumSPP
!      read(11,*)
!      kGridPt        = 1
!  150 continue
!      read(11,*,end=199) (pGPP(kGridPt,j), j=1, nNumSPP)
!      kGridPt        = kGridPt + 1
!      goto 150
!  199 continue
!      nNumGPts       = kGridPt
!      close(11)

!     ----------------------------------------------------------------
!     * read-in jedi 1 output 
!     ----------------------------------------------------------------

      kLonOld        = -999
      kLatOld        = -999
      kGridPt        = 0
      kSPPCnt(:)     = 0
      nNumSPP        = 0
      pRelAbd(:,:)   = 0.
      pGPP(:,:)      = 0.
      pCS(:,:)       = 0.
      pCA(:,:)       = 0.
      pCTot(:,:)     = 0.
      
      open(10, file='jedi_species.txt', form='formatted')

  200 continue

      read(10,'(2I5,I10,8E12.4)', end=249) kLon, kLat, kSPPIDtemp,      &
     &                   zRAbd, zGPP, zCS, zCA, zCL, zCR, zCWL, zCWR
!      write(*,'(2I5,I10,8E12.4)') kLon, kLat, kSPPIDtemp,      &
!     &                   zRAbd, zGPP, zCS, zCA, zCL, zCR, zCWL, zCWR

!     * check for new grid point
      if ((kLon .ne. kLonOld) .or. (kLat .ne. kLatOld)) then
        write(*,*) 'analyze_composition: ', kGridPt, kSPPCnt(kGridPt), kLonOld, kLatOld
        kGridPt          = kGridPt + 1
        kSPPCnt(kGridPt) = 0
        kLonOld          = kLon
        kLatOld          = kLat
      end if

      kSPPCnt(kGridPt)                   = kSPPCnt(kGridPt) + 1
      pRelAbd(kGridPt, kSPPCnt(kGridPt)) = zRAbd
      pGPP(kGridPt, kSPPCnt(kGridPt))    = zGPP
      pCS(kGridPt, kSPPCnt(kGridPt))     = zCS
      pCA(kGridPt, kSPPCnt(kGridPt))     = zCA
      pCTot(kGridPt, kSPPCnt(kGridPt))   = zCA + zCL + zCR + zCWL + zCWR      
      kSPPID(kGridPt, kSPPCnt(kGridPt))  = kSPPIDtemp
      nNumSPP                            = max(nNumSPP, kSPPCnt(kGridPt))
      goto 200

  249 continue
      close(10)
 
      nNumGPts       = kGridPt
      
      write(*,*) 'analyze_composition: nNumGPts = ', nNumGPts
      write(*,*) 'analyze_composition: nNumSPP  = ', nNumSPP
      
!     ----------------------------------------------------------------
!     * analyze: richness versus evenness
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: richness vs. evennness'

      pRichness(:)   = 0.
      pEvenness(:)   = 0.
      
      do i = 1, nNumGPts
      do j = 1, nNumSPP
        if (pRelAbd(i, j) .gt. 0.) then
          pRichness(i) = pRichness(i) + 1.
          pEvenness(i) = pEvenness(i) - pRelAbd(i, j) * log10(pRelAbd(i, j))
        end if      
      end do
      end do
      
      pEvenness(:)  = pEvenness(:)/MAXVAL(pEvenness(:))
      pRichness(:)  = pRichness(:)/MAXVAL(pRichness(:))
      
      open(20, file='analyze-evenness.txt')
      do i = 1, nNumGPts
        write(20, *) pRichness(i), char(9), pEvenness(i)  
      end do
      close(20)

!     ----------------------------------------------------------------
!     * analyze: frequency distribution of species abundance
!       for given diversity values using log abundance bins
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: relative abundance distribution'
      
!     * do binning

      do k = 1, nNumClasses
        pFreqDist(:)   = 0.
        nSamples       = 0
        pClassRichness  = 0.

        open(20, file='analyze-rawabd_d'//sText(k)//'.txt')
        
        do i = 1, nNumGPts
          if (pRichness(i) .gt. pClassHigh(k-1) .and. pRichness(i) .le. pClassHigh(k)) then
            nSamples   = nSamples + 1
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pClassRichness = pClassRichness + 1
                do i2 = 1, nNumBins
                  if (pRelAbd(i, j) .gt. pBinHigh(i2-1) .and. pRelAbd(i, j) .le. pBinHigh(i2)) then
                    pFreqDist(i2) = pFreqDist(i2) + 1.
                  end if
                end do
                write(20, *) i, pRelAbd(i, j)            
              end if
            end do
          end if
        end do
        close(20)
        
        pClassRichness = pClassRichness/real(nSamples)
        
        open(21, file='analyze-relabd_spp_d'//sText(k)//'.txt')
        write(21, *) 'relabd_spp'//char(9)//'d'//sText(k)
        do i2 = 1, nNumBins
          write(21, *) pBinHigh(i2), char(9), pFreqDist(i2)/real(nSamples)/pClassRichness
        end do
        close(21)
        
      end do

!     ----------------------------------------------------------------
!     * analyze: rank abundance plot of species
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: ranked abundance'

      kMostAbdSPP(:)   = 0
      pMostAbdSPPAbd(:)= 0.
      
!     * do ranking

      do k = 1, nNumClasses
        pAccRankAbd(:)            = 0.
        nSamples                  = 0    
        do i = 1, nNumGPts
          if (pRichness(i) .gt. pClassHigh(k-1) .and. pRichness(i) .le. pClassHigh(k)) then
            nSamples              = nSamples + 1
            do i2 = 1, nMaxSPP
              pRankedAbd(i2)      = pRelAbd(i, i2)
              kRankSPP(i2)        = i2
            end do
!           * count total number of species for normalization            
            pSPPCount2(i)         = 0.
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pSPPCount2(i)     = pSPPCount2(i) + 1.
              end if
            end do
!           * sort species by abundance
            do i2 = 1, nNumSPP
              do j2 = i2+1, nNumSPP
                if (pRankedAbd(i2) .le. pRankedAbd(j2)) then
                  pRankTemp       = pRankedAbd(i2)
                  pRankedAbd(i2)  = pRankedAbd(j2)
                  pRankedAbd(j2)  = pRankTemp
                  kRankSPPTemp    = kRankSPP(i2)
                  kRankSPP(i2)    = kRankSPP(j2)
                  kRankSPP(j2)    = kRankSPPTemp
                end if
              end do
            end do
!           * accumulate
            do i2 = 1, pSPPCount2(i)
              pAccRankAbd(i2)     = pAccRankAbd(i2) + pRankedAbd(i2)
            end do
!           * save most abundant species
            kMostAbdSPP(i)        = kRankSPP(1)
            pMostAbdSPPAbd(i)     = pRankedAbd(1)
          end if
        end do
        
        open(23, file='analyze-rank_spp_d'//sText(k)//'.txt')
        write(23, *) 'rank_spp'//char(9)//'d'//sText(k)
        do i2 = 1, nNumSPP
           write(23, *) i2, char(9), pAccRankAbd(i2)/real(nSamples)
        end do
        close(23)
        
      end do

      open(24, file='analyze-rank_mostabd.txt')
      write(24, *) 'grdpt'//char(9)//'kMostAbdSPP'//char(9)//'pMostAbdSPPAbd'
      do i = 1, nNumGPts
        write(24, *) i, char(9), kSPPID(i, kMostAbdSPP(i)), char(9), pMostAbdSPPAbd(i)
      end do
      close(24)

!     ----------------------------------------------------------------
!     * analyze: frequency distribution of species productivity
!       for given diversity values using log abundance bins
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: relative productivity distribution'
      
!     * do binning

      do k = 1, nNumClasses
        pFreqDist(:)   = 0.
        nSamples       = 0
        pClassRichness  = 0.        
        do i = 1, nNumGPts
          if (pRichness(i) .gt. pClassHigh(k-1) .and. pRichness(i) .le. pClassHigh(k)) then
            nSamples   = nSamples + 1
!           * find maximum productivity
            pGPPMax               = 0.
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pGPPMax           = max(pGPPMax, pGPP(i, j))
              end if
            end do
!           * bin relative productivity
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pClassRichness = pClassRichness + 1
                do i2 = 1, nNumGPPBins
                  if (pGPP(i, j)/pGPPMax .gt. pGPPBinHigh(i2-1) .and. pGPP(i, j)/pGPPMax .le. pGPPBinHigh(i2)) then
                    pFreqDist(i2) = pFreqDist(i2) + 1.
                  end if
                end do
              end if
            end do
          end if
        end do

        pClassRichness = pClassRichness/real(nSamples)
                
        open(25, file='analyze-relabd_gpp_d'//sText(k)//'.txt')
        write(25, *) 'relabd_gpp'//char(9)//'gpp_d'//sText(k)
        do i2 = 1, nNumGPPBins
          write(25, *) pGPPBinHigh(i2), char(9), pFreqDist(i2)/real(nSamples)/pClassRichness
        end do
        close(25)
        
      end do
      
!     ----------------------------------------------------------------
!     * analyze: rank abundance plot for productivity
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: ranked productivity'

!     * do ranking

      do k = 1, nNumClasses
        pAccRankGPP(:)            = 0.
        nSamples                  = 0
        do i = 1, nNumGPts
          if (pRichness(i) .gt. pClassHigh(k-1) .and. pRichness(i) .le. pClassHigh(k)) then
            nSamples              = nSamples + 1
            do i2 = 1, nMaxSPP
              pRankedGPP(i2)      = pGPP(i, i2)
              kRankSPP(i2)        = i2
            end do
!           * find maximum GPP of grid cell
            pGPPMax               = 0.
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pGPPMax           = max(pGPPMax, pGPP(i, j))
              end if
            end do
!           * count total number of species for normalization            
            pSPPCount2(i)         = 0.
            do j = 1, nNumSPP
              if (pRelAbd(i, j) .gt. 0.) then
                pSPPCount2(i)     = pSPPCount2(i) + 1.
              end if
            end do
!           * sort species by gpp
            do i2 = 1, nNumSPP
              do j2 = i2+1, nNumSPP
                if (pRankedGPP(i2) .le. pRankedGPP(j2)) then
                  pRankTemp       = pRankedGPP(i2)
                  pRankedGPP(i2)  = pRankedGPP(j2)
                  pRankedGPP(j2)  = pRankTemp
                  kRankSPPTemp    = kRankSPP(i2)
                  kRankSPP(i2)    = kRankSPP(j2)
                  kRankSPP(j2)    = kRankSPPTemp
                end if
              end do
            end do
!           * accumulate relative GPP
            do i2 = 1, pSPPCount2(i)
              pAccRankGPP(i2)     = pAccRankGPP(i2) + pGPP(i, kRankSPP(i2))/pGPPMax
            end do
          end if
        end do
        
        open(26, file='analyze-rank_gpp_d'//sText(k)//'.txt')
        write(26, *) 'rank_spp'//char(9)//'gpp_d'//sText(k)
        do i2 = 1, nNumSPP
          write(26, *) real(i2), char(9), pAccRankGPP(i2)/real(nSamples)
        end do
        close(26)
        
      end do

!     ----------------------------------------------------------------
!     * analyze: carbon partitioning in richness classes
!     ----------------------------------------------------------------

!      write(*,*) 'analyze_composition: carbon partitioning'
!
!      where (pGPP(:,:) .gt. 0.)
!        pResTime(:,:)  = pCTot(:,:)/pGPP(:,:)
!      end where
!
!      where (pGPP(:,:) .gt. 0.)
!        pCSTime(:,:)   = pCS(:,:)/pGPP(:,:)
!      end where
!
!!      pGPGPPMean(:)    = 0.
!!      pGPSTMean(:)     = 0.
!!      pGPRTMean(:)     = 0.
!      
!      pGPPMax          = 0.
!      pSTMax           = 0.
!      pRTMax           = 0.
!                 
!      do k = 1, nNumClasses
!        open(27, file='analyze-carbon_d'//sText(k)//'.txt')
!        write(27, *) 'rich_d'//sText(k)//char(9)//'even_d'//sText(k)//char(9)  &
!     &             //'gpp_d'//sText(k)//char(9)//'st_d'//sText(k)//char(9)     &
!     &             //'rt_d'//sText(k)//char(9)
!        do i = 1, nNumGPts
!          if (pRichness(i) .gt. pClassHigh(k-1) .and. pRichness(i) .le. pClassHigh(k)) then
!            pGPGPPMean  = 0.
!            pGPSTMean   = 0.
!            pGPRTMean   = 0.
!            kCnt        = 0
!            do j = 1, nNumSPP
!              if (pRelAbd(i, j) .gt. 0.) then
!                pGPGPPMean = pGPGPPMean + pGPP(i, j)
!                pGPSTMean  = pGPSTMean  + pCSTime(i, j)
!                pGPRTMean  = pGPRTMean  + pResTime(i, j)
!                kCnt       = kCnt       + 1
!              end if
!            end do
!            pGPGPPMean  = pGPGPPMean/real(kCnt)
!            pGPSTMean   = pGPSTMean/real(kCnt)
!            pGPRTMean   = pGPRTMean/real(kCnt)
!            write(27, '(E12.4,A1,E12.4,A1,E12.4,A1,E12.4,A1,E12.4)')           &
!     &                   pRichness(i), char(9), pEvenness(i), char(9),         &
!     &                   pGPGPPMean, char(9), pGPSTMean, char(9), pGPRTMean
!          end if
!        end do
!        close(27)
!      end do

!     ----------------------------------------------------------------
!     * finish
!     ----------------------------------------------------------------

      write(*,*) 'analyze_composition: end'
      
      end
      