load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


kNumTraits=15
kNumBins=20
kNumLat=64
kNumLon=128
kNumClass=3

f=addfile("gbins.nc","r")
gbins = f->gbins

g=addfile("class_koeppen.nc","r")
class0 = g->var1(0,:,:)
;g=addfile("land_class_T42.nc","r")
;class0 = tofloat(g->LAND_CLASS)

class = conform_dims((/kNumBins,kNumLat,kNumLon/),class0,(/1,2/))

poo = new((/kNumClass,kNumBins/),"float")
wks = gsn_open_wks("pdf","class")

; trait number
do t=0,kNumTraits-1
mask0 = class
mask0@_FillValue=-999
mask0(:,:,:) = 0.


i=(/4/)
 do j=0,0

 mask0 = where(class .eq. i(j),1.0,mask0@_FillValue)
 masksum = sum(mask0)

 abc = mask0 * gbins(:,t,:,:)
 poo(j,:) = dim_avg_n(dim_avg_n(abc,2),1) / dim_sum_n(dim_sum_n(abc,2),1)


end do

;*********************************
; create plot
;*********************************
 res   = True
 ;res@gsnScale = True

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

res@vpHeightF = 0.35        ; Changes the aspect ratio
res@vpWidthF  = 0.68

 res@tiXAxisString = "Normalized trait value"
 res@tiYAxisString = "Frequency"
res@trXMinF = 0.0
res@trXMaxF = 1.0

res@xyLineThicknesses = (/4., 4., 4./)
res@xyDashPatterns = (/0, 0, 0/)
res@xyLineColors      = (/"blue","green","red"/)          ; change line color



; res@gsnMaximize       = True

 res@gsnXYBarChart = True
 res@gsnXYBarChartOutlineOnly = True
 p = t + 1
 res@tiMainString = "Trait " + p

 plot = gsn_csm_xy (wks,gbins&bin,poo,res)

end do

