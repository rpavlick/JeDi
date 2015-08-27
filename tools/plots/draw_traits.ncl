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
abd(0,:,:) = 0.0

abdsizes = dimsizes(abd)

kNumLat = abdsizes(1)
kNumLon = abdsizes(2)
kMaxSPP = abdsizes(0)
kNumTraits = 15
kTotalTraits = 20

; read in jedi_success file
success = readAsciiTable("jedi_success.txt",kTotalTraits+1,"float",0)

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

abdmask = dim_sum_n(abd,0)
abdmask@_FillValue = -9e20
abdmask = where(abdmask.gt.0.0,abdmask,abdmask@_FillValue)

abd = abd / conform_dims((/kMaxSPP,kNumLat,kNumLon/),abdmask,(/1,2/))
delete(abdmask)

mtrait = dim_sum_n(gtrait *  conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),abd,(/0,2,3/)),0) 

gtrait = (gtrait - conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),mtrait,(/1,2,3/)))^2.0 
qtrait = dim_sum_n(sqrt(gtrait) *  conform_dims((/kMaxSPP,kNumTraits,kNumLat,kNumLon/),abd,(/0,2,3/)),0) 

fdis = sqrt(dim_sum_n(dim_sum_n(gtrait,1) * abd,0))

  mtrait = mask(mtrait,conform_dims((/kNumTraits,kNumLat,kNumLon/),landmask,(/1,2/),1)
  mtrait = mtrait(:,::-1,:)

  mtrait!0 = "trait"
  mtrait!1 = "lat"
  mtrait!2 = "lon"
  mtrait&trait = ispan(1,20,1)
  mtrait&lat = lat
  mtrait&lon = lon

  qtrait = mask(qtrait,conform_dims((/kNumTraits,kNumLat,kNumLon/),landmask,(/1,2/),1)
  qtrait = qtrait(:,::-1,:)

  qtrait!0 = "trait"
  qtrait!1 = "lat"
  qtrait!2 = "lon"
  qtrait&trait = ispan(1,20,1)
  qtrait&lat = lat
  qtrait&lon = lon

;  fdis = mask(fdis,landmask,1)
;  fdis = fdis(::-1,:)

  fdis!0 = "lat"
  fdis!1 = "lon"
  fdis&lat = lat
  fdis&lon = lon

fm = addfile("traits.nc","c")
fm->CWT = mtrait
fm->FD  = qtrait
fm->Fdis = fdis

; define biodiversity colormap
  cmap = (/(/255., 255., 255./), \
 (/0., 0., 0./), \
 (/255., 240., 150./), \
 (/220., 240., 84./),   \
 (/191., 255., 100./),  \
 (/153., 255., 71./),  \
 (/26., 255., 26./),  \
 (/0., 240., 179./),  \
 (/0., 127., 240./),  \
 (/170., 0., 240./),  \
 (/215., 0., 0./)/)
 cmap = cmap/256.0

wks = gsn_open_wks("ps","panel")               ; open a ps file
plot = new(15,graphic)                          ; create a plot array

gsn_define_colormap(wks,"testcmap")        ; choose colormap
res          = True
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance fram
res@cnFillMode            = "RasterFill"
res@cnFillOn              = True               ; turn on color fill
res@cnLinesOn             = False              ; turn off contour lines
res@gsnSpreadColors       = True               ; use full range of color map
res@mpMinLatF             = -60

res@cnLevelSelectionMode  = "ManualLevels"
res@cnMinLevelValF        = 0.2              ; Min level
res@cnMaxLevelValF        = 0.8           ; Max level
res@cnLevelSpacingF       = 0.05            ; interval

i = 0
do while(i.le.14)
  ctrait = (/mtrait(i,:,:)/) 
  ctrait@_FillValue = -9e20
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

gsn_panel(wks,plot,(/5,3/),pres)             ; now draw as one plot


  gsn_define_colormap(wks,cmap)        ; choose colormap
  delete(res@cnMinLevelValF)
  delete(res@cnMaxLevelValF)
  delete(res@cnLevelSpacingF)
  res@cnLevelSelectionMode = "AutomaticLevels" ; set explicit contour levels
;  maxval = 0.28
;  levels = maxval *(/ 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8 /)
;  res@cnLevels = levels

i = 0
do while(i.le.14)
  ctrait = (/qtrait(i,:,:)/) 
  ctrait@_FillValue = -9e20
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

gsn_panel(wks,plot,(/5,3/),pres)             ; now draw as one plot


delete(plot)


  gsn_define_colormap(wks,cmap)        ; choose colormap
;  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels
;  maxval = 0.50
;  levels = maxval *(/ 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8 /)
;  res@cnLevels = levels

  ctrait = (/fdis(:,:)/)
  ctrait@_FillValue = -9e20
  ctrait!0 = "lat"
  ctrait!1 = "lon"
  ctrait&lat = lat
  ctrait&lon = lon
  ctrait@units  = "Unitless"
  ctrait@long_name = "Functional Dispersion"
  plot = gsn_csm_contour_map_ce(wks,ctrait,res)

  draw(plot)
  frame(wks)
