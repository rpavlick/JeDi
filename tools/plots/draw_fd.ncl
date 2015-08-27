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


mtrait = dim_sum_n(gtrait *  conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),abd,(/0,2,3/)),0) 

gtrait = (gtrait - conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),mtrait,(/1,2,3/)))^2.0
qtrait = dim_sum_n(gtrait *  conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),abd,(/0,2,3/)),0) 


wks = gsn_open_wks("ps","panel")               ; open a ps file
plot = new(17,graphic)                          ; create a plot array

gsn_define_colormap(wks,"testcmap")        ; choose colormap
res          = True
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance fram
res@cnFillOn              = True               ; turn on color fill
res@cnLinesOn             = False              ; turn off contour lines
res@gsnSpreadColors       = True               ; use full range of color map
res@mpMinLatF             = -60

abdmask = dim_sum_n(abd(1:kMaxSPP-1,:,:),0)
abdmask@_FillValue = -9e20
abdmask = where(abdmask.gt.0.0,abdmask,abdmask@_FillValue)



i = 0
do while(i.le.16)
  ctrait = (/qtrait(i,:,:)/) 
  ctrait@_FillValue = -9e20
  ctrait = mask(ctrait,landmask(::-1,:),1) / abdmask
  ctrait = ctrait(::-1,:)
  ctrait!0 = "lat"  
  ctrait!1 = "lon"
  ctrait&lat = lat
  ctrait&lon = lon
  ctrait@units  = "Unitless"  
  j = i + 1
  ctrait@long_name = "Trait " + j
  plot(i) = gsn_csm_contour_map_ce(wks,ctrait,res)     
  i=i+1
end do 


pres = True
pres@gsnMaximize           = True

gsn_panel(wks,plot,(/4,5/),pres)             ; now draw as one plot

res@cnLevelSelectionMode  = "ManualLevels"
res@cnMinLevelValF        = 0.2              ; Min level
res@cnMaxLevelValF        = 0.8           ; Max level
res@cnLevelSpacingF       = 0.05            ; interval

i = 0
do while(i.le.16)
  ctrait = (/mtrait(i,:,:)/) 
  ctrait@_FillValue = -9e20
  ctrait = mask(ctrait,landmask(::-1,:),1) / abdmask
  ctrait = ctrait(::-1,:)
  ctrait!0 = "lat"  
  ctrait!1 = "lon"
  ctrait&lat = lat
  ctrait&lon = lon
  ctrait@units  = "Unitless"  
  j = i + 1
  ctrait@long_name = "Trait " + j
  plot(i) = gsn_csm_contour_map_ce(wks,ctrait,res)     
  i=i+1
end do 

pres = True
pres@gsnMaximize           = True

gsn_panel(wks,plot,(/4,5/),pres)             ; now draw as one plot

