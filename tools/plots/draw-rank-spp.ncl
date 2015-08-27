
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

file1="analyze-rank_spp_d1.txt"
file2="analyze-rank_spp_d2.txt"
file3="analyze-rank_spp_d3.txt"
file4="analyze-rank_spp_d4.txt"

data01 = readAsciiTable(file1, 2, "float", 1)
data02 = readAsciiTable(file2, 2, "float", 1)
data03 = readAsciiTable(file3, 2, "float", 1)
data04 = readAsciiTable(file4, 2, "float", 1)

nspp = numAsciiRow(file1)
nspp = nspp - 1

data      = new((/4,nspp/),float)
sppidx = data01(:,0)

data(0,:) = data01(:,1)
data(1,:) = data02(:,1)
data(2,:) = data03(:,1)
data(3,:) = data04(:,1)

 wks   = gsn_open_wks ("ps","rankspp")                  ; open workstation
 res = True
 res@gsnScale = True
; res@gsnMaximize = True
 res@vpHeightF = 0.3        ; Changes the aspect ratio
 res@vpWidthF  = 0.8 
 res@tiMainString = "(b)"
 res@xyLineThicknesses = (/4., 4., 4., 4./)
 res@xyDashPatterns = (/0, 0, 0, 0/)
 res@xyLineColors      = (/"green","blue","purple","red"/)          ; change line color

 res@tiXAxisString = "Ranked growth strategy"
 res@trXMinF = 1
 res@trXMaxF = nspp
  
 res@tiYAxisString = "Relative abundance"
 res@trYLog        = True
 res@trYMinF = 0.0000001
 res@trYMaxF = 0.1
 data = where(data .le. 0.0 ,data@_FillValue,data)

; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend

 res@pmLegendSide           = "Top"                 ; Change location of
 res@pmLegendParallelPosF   = 0.80                   ; move units right
 res@pmLegendOrthogonalPosF = -0.8                  ; more neg = down
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.25                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@xyExplicitLegendLabels = (/" Q1"," Q2"," Q3"," Q4"/)
res@pmLegendHeightF = 0.12
 plot  = gsn_csm_xy(wks,sppidx,data,res) ; create plot


end
