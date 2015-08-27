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
; the data from Batjes 1996 is only provided in ranges
; these are converted to specific values following Kleidon and Heimann 2000
; this file can be regridded as necessary


infile = addfile("original_batjes1996.nc","r")
outfile = addfile("paw.nc","c")

lat = infile->lat
lon = infile->lon
original = tofloat(infile->Band1)
original = original

paw = original
paw!0 = "lat"
paw!1 = "lon"
paw&lat = lat
paw&lon = lon
paw@_FillValue = default_fillvalue("float")
paw = 0.105 ; default value for grid cells misaligned with other forcing datasets

paw = where(original .eq. 1,0.040,paw)
paw = where(original .eq. 2,0.075,paw)
paw = where(original .eq. 3,0.105,paw)
paw = where(original .eq. 4,0.135,paw)
paw = where(original .eq. 5,0.175,paw)
paw = where(original .eq. 6,0.350,paw)
paw = where(original .eq. 7,0.060,paw)
paw = where(original .eq. 8,0.120,paw)
paw = where(original .eq. 9,0.200,paw)

paw@title = "Unit available water capacity"
paw@units = "mm water per mm soil"
paw@source = "ISRIC-WISE Soil Moisture Retention"
outfile->paw = paw
system("cdo -invertlat -sellonlatbox,-180,180,90,-90 paw.nc p.nc; mv p.nc paw.nc")
system("cdo -b F64 -f srv -sellonlatbox,-180,180,90,-90 paw.nc paw.srv")
