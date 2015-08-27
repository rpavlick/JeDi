; ***********************************************
; traits.ncl
;
; - plot mean trait values weighted by rel. abd.
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************

; read grid description
fm2 = addfile("/Net/Groups/C-Side/BTM/scratch/rpavlick/clamp_data/surface_model/lnd_T42.nc","r")
lat = fm2->lat
lon = fm2->lon
area = fm2->area
landfrac = fm2->landfrac
landmask = fm2->landmask

; cdo -f nc copy jedi_species.srv jedi_species.nc
; read in relative abundances from jedi_species.nc
fm = addfile("jedi_species.nc","r")
abd0 = fm->var5400
grabd = (/abd0(0,:,::-1,:)/)

abdmask = dim_sum_n(grabd,0)
abdmask@_FillValue = -9e20
abdmask = where(abdmask.gt.0.0,abdmask,abdmask@_FillValue)

abd = abd / conform_dims((/kMaxSPP,kNumLat,kNumLon/),abdmask,(/1,2/))
delete(abdmask)

; define dimension sizes
grabdsizes = dimsizes(grabd)
kMaxSPP = grabdsizes(0)
kNumLat = grabdsizes(1)
kNumLon = grabdsizes(2)
kNumTraits = 15
kTotalTraits = 20
spp = fspan(1,kMaxSPP,kMaxSPP)

; write trait matrix to file
grabd!0 = "lev"  
grabd!1 = "lat"  
grabd!2 = "lon"
grabd&lev = spp
grabd&lat = lat
grabd&lon = lon

; read in jedi_success file
success = readAsciiTable("jedi_success.txt",kTotalTraits+1,"float",0)

; get index of successful species from first column of success
sppidx0  = success(:,0)
sppidx = floattointeger(sppidx0) - 1

; create trait matrix
traits = new((/kMaxSPP,kNumTraits/),"float")
traits(sppidx,:) = success(:,1:kNumTraits)

; define bin sizes
kNumBins = 20
pBinMin = 0.
pBinMax = 1.
pBinDelta = (pBinMax - pBinMin) / (1. * kNumBins)

; bin minimums and maximums
z0 = fspan(pBinMin,pBinMax - pBinDelta,kNumBins)
z1 = fspan(pBinDelta,pBinMax,kNumBins)
kTraits = ispan(1,15,1)

gbins = new((/kNumBins,kNumTraits,kNumLat,kNumLon/),"float")
gbins!0 = "bin"
gbins!1 = "trait"  
gbins!2 = "lat"  
gbins!3 = "lon"
gbins&bin = z0
gbins&trait = kTraits
gbins&lat = lat
gbins&lon = lon

do j = 0,kNumTraits-1
do i = 0,kNumBins-1
zidx = ind(traits(:,j) .ge. z0(i) .and. traits(:,j) .lt. z1(i))
if(.not. any(ismissing(zidx))) then
if (dimsizes(zidx) .ge. 2) then
gbins(i,j,:,:) = (/dim_sum_n(grabd(zidx,:,:),0)/)
else
gbins(i,j,:,:) = (/grabd(zidx,:,:)/)
end if
end if
delete(zidx)
end do
end do

geven = new((/kNumTraits,kNumLat,kNumLon/),"float")
geven!0 = "lev"  
geven!1 = "lat"  
geven!2 = "lon"
geven&lev = kTraits
geven&lat = lat
geven&lon = lon

even = new((/kNumLat,kNumLon/),"float")
even!0 = "lat"  
even!1 = "lon"
even&lat = lat
even&lon = lon


ww = where(grabd .ne. 0., grabd * -log(grabd), 0.0)
ww0 = grabd
ww0@_FillValue = -999.
ww0 = where(grabd .ne. 0., 1.0, ww0@_FillValue)
ww1 = dim_sum_n_Wrap(ww0,0)
ww2 = where(ww1 .gt. 1, ww1, ww0@_FillValue)
even = dim_sum_n(ww,0) / log(ww2)

yy = where(gbins .ne. 0., gbins * -log(gbins), 0.0)
geven = (/dim_sum_n(yy,0)/) / log(kNumBins)
beven = geven / conform_dims((/kNumTraits,kNumLat,kNumLon/),even,(/1,2/))

beven!0 = "lev"  
beven!1 = "lat"  
beven!2 = "lon"
beven&lev = kTraits
beven&lat = lat
beven&lon = lon


h = addfile("gbins.nc","c")
h->gbins = gbins
h->geven = geven
h->even = even
h->beven = beven

