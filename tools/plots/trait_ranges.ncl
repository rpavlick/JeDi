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
fm = addfile("jedi_species.nc","r")
abd0 = fm->var5400
abd = (/abd0(0,:,:,:)/)

abdsizes = dimsizes(abd)

kNumLat = abdsizes(1)
kNumLon = abdsizes(2)
kMaxSPP = abdsizes(0)
kNumTraits = 20

; read in jedi_success file
success = readAsciiTable("jedi_success.txt",kNumTraits+1,"float",0)

; get index of successful species from first column of success
sppidx0  = success(:,0)
sppidx = floattointeger(sppidx0) - 1

; create trait matrix
traits = new((/kMaxSPP,kNumTraits/),"float")
traits(sppidx,:) = success(:,1:kNumTraits)

gtrait = conform_dims( (/kMaxSPP,kNumTraits,kNumLat,kNumLon/), traits , (/0,1/) ) 

alloctraitsum = (gtrait(:,2,:,:) + gtrait(:,3,:,:) + gtrait(:,4,:,:) + gtrait(:,5,:,:))
alloctraitsum@_FillValue = -9e20
alloctraitsum = where(alloctraitsum.gt.0,alloctraitsum,alloctraitsum@_FillValue)

i = 2
do while(i.le.5)
 gtrait(:,i,:,:) = gtrait(:,i,:,:) / alloctraitsum 
 i = i + 1
end do


gabd = conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),abd,(/0,2,3/))
