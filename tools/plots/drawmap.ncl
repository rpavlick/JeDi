load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
;************************************************

;************************************************
; open file and read in data
;************************************************
  in = addfile("jedi_mean.nc","r")
  rich  = in->RICH
  even  = in->EVEN
;************************************************
; create plot
;************************************************
  wks = gsn_open_wks("pdf" ,"ce")                ; open a ps file
biodiv = (/ (/255,255,255/), (/0,0,0/), (/255, 240, 150/), \
(/220, 240, 84/), \
(/191, 255, 100/),\
(/153, 255, 71/), \
(/26, 255, 26/), \
(/0, 240, 179/), \
(/0, 127, 240/), \
(/170, 0, 240/), \
(/215, 0, 0/)/)
biodiv_cmap = biodiv / 255.0
  gsn_define_colormap(wks,biodiv_cmap)        ; choose colormap
  plot = new(2,graphic)                          ; create a plot array

  res                       = True     ; plot mods desired
  res@cnFillMode            = "RasterFill"
  res@cnFillOn              = True     ; turn on color fill
  res@cnLinesOn             = False    ; turn of contour lines
  res@gsnSpreadColors  = True
  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels 
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame
  res@mpMinLatF            = -60.      ; range to zoom in on
  res@mpMaxLatF            =  90.

  minval = 3
  maxval = max(rich)
  res@cnLevelSelectionMode = "ExplicitLevels" ; set explicit contour levels

  levels = (/ 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8 /)
  res@cnLevels = levels
  rich = where(rich .ge. minval, rich/maxval, rich@_FillValue) 
  res@tiMainString = "(a)"
  res@gsnLeftString = "Functional richness"
  res@gsnRightString = "unitless"

  plot(0) = gsn_csm_contour_map_ce(wks,rich(0,:,:), res)

  levels = (/ 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8 /)
  res@cnLevels = levels
  res@tiMainString = "(b)"

  res@gsnLeftString = "Functional evenness"
  res@gsnRightString = "unitless"
  even = where((maxval * rich) .ge. minval, even, even@_FillValue) 
  plot(1) = gsn_csm_contour_map_ce(wks,even(0,:,:), res)

  panres = True ; panel resource list
  panres@gsnMaximize                = True

;************************************************
; create panel
;************************************************
  gsn_panel(wks,plot,(/2,1/),panres)             ; now draw as one plot
