
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

file1="analyze-relabd_spp_d1.txt"
file2="analyze-relabd_spp_d2.txt"
file3="analyze-relabd_spp_d3.txt"
file4="analyze-relabd_spp_d4.txt"

data01 = readAsciiTable(file1, 2, "float", 1)
data02 = readAsciiTable(file2, 2, "float", 1)
data03 = readAsciiTable(file3, 2, "float", 1)
data04 = readAsciiTable(file4, 2, "float", 1)

nspp = numAsciiRow(file1)
nspp = nspp - 1

data      = new((/4,nspp,2/),float)
sppidx = data01(:,0)
printVarSummary(data01)
data(0,:,:) = data01(:,:)
data(1,:,:) = data02(:,:)
data(2,:,:) = data03(:,:)
data(3,:,:) = data04(:,:)

;*********************************
; create plot
;K********************************
 wks   = gsn_open_wks ("ps", "relabd" )
 res   = True
 res@gsnScale = True               

; these four resources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

; res@gsnMaximize = True
 res@gsnScale = True
 res@vpHeightF = 0.4        ; Changes the aspect ratio
 res@vpWidthF  = 0.6

 res@tiXAxisString = "Relative abundance of growth strategy"
 res@tiYAxisString = "Frequency"
 res@trXLog = True
res@trXMinF = 0.0001

; res@gsnMaximize       = True

 res@xyLineThicknesses = (/4., 4., 4., 4./)
 res@xyDashPatterns = (/0, 0, 0, 0/)
 res@xyLineColors      = (/"green","blue","purple","red"/)          ; change line color

 res@gsnXYBarChart = True
 res@gsnXYBarChartOutlineOnly = True

; add a legend
 res@pmLegendDisplayMode    = "Always"              ; turn on legend
 
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .85                   ; move units right
 res@pmLegendOrthogonalPosF = -0.8                  ; more neg = down
 
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.25                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
; labels for the legend
 res@xyExplicitLegendLabels = (/" Q1"," Q2"," Q3"," Q4"/)
 res@tiMainString = "(a)"
res@pmLegendHeightF = 0.12
 plot = gsn_csm_xy (wks,data(:,:,0),data(:,:,1),res)


end
