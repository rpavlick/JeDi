; ###############################################################################
; # create_paw.ncl
; # Copyright (c) 2015 Ryan Pavlick
; #
; # Permission is hereby granted, free of charge, to any person obtaining a copy
; # of this software and associated documentation files (the "Software"), to deal
; # in the Software without restriction, including without limitation the rights
; # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
; # copies of the Software, and to permit persons to whom the Software is
; # furnished to do so, subject to the following conditions:
; #
; # The above copyright notice and this permission notice shall be included in all
; # copies or substantial portions of the Software.
; #
; # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
; # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
; # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
; # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
; # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
; # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
; # SOFTWARE.
; ###############################################################################

; this script creates the unit available water input file (paw.srv) for JeDi-DGVM
; it uses a global map of soil texture classes from the HWSD
; and applies average AWC values from the USDA
; this file can be regridded as necessary


infile = addfile("hwsd.final.hlf.nc4","r")
outfile = addfile("paw.nc","c")

lat = infile->lat
lon = infile->lon
original = toint(infile->soiltexture)
paw = tofloat(original)
paw!0 = "lat"
paw!1 = "lon"
paw&lat = lat
paw&lon = lon
paw@_FillValue = default_fillvalue("float")
paw = 0.167 ; assign default to loam


paw = where(original .eq. 1,0.042,paw) ; Sand
paw = where(original .eq. 2,0.083,paw) ; Loamy Sand
paw = where(original .eq. 3,0.125,paw) ; Sandy Loam
paw = where(original .eq. 4,0.167,paw) ; Loam
paw = where(original .eq. 5,0.167,paw) ; Silt Loam
paw = where(original .eq. 6,0.167,paw) ; Silt
paw = where(original .eq. 7,0.183,paw) ; Sandy Clay Loam
paw = where(original .eq. 8,0.183,paw) ; Clay Loam
paw = where(original .eq. 9,0.183,paw) ; Silty Clay Loam
paw = where(original .eq. 10,0.167,paw) ; Sandy Clay
paw = where(original .eq. 11,0.167,paw) ; Silty Clay
paw = where(original .eq. 12,0.167,paw) ; Clay


paw@title = "Unit available water capacity"
paw@units = "mm water per mm soil"
outfile->paw = paw
system("cdo -b F64 -f srv -sellonlatbox,-180,180,90,-90 paw.nc paw.srv")
system("cdo -sellonlatbox,-180,180,90,-90 paw.nc p.nc; mv p.nc paw.nc")

