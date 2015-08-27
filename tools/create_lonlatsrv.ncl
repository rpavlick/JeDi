
system("rm latitude.nc longitude.nc latitude.srv longitude.srv")

infile = addfile("landsea.nc","r")
lat = infile->lat
lon = infile->lon
nlat=dimsizes(lat)
nlon=dimsizes(lon)
lat2d = conform_dims((/nlat,nlon/),lat,0)
lon2d = conform_dims((/nlat,nlon/),lon,1)

lat2d!0	 					  =	"lat"
lat2d!1							=	"lon"
lat2d&lat 					=	lat
lat2d&lon 					=	lon

lat2dfile = addfile("latitude.nc","c")
lat2dfile->lat2d= lat2d

lon2d!0	 					  =	"lat"
lon2d!1						  =	"lon"
lon2d&lat 					=	lat
lon2d&lon 					=	lon

lon2dfile = addfile("longitude.nc","c")
lon2dfile->lon2d = lon2d

system("cdo -f srv -b F64 copy latitude.nc latitude.srv")
system("cdo -f srv -b F64 copy longitude.nc longitude.srv")
