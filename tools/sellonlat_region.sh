MIN_LON="-80"
MAX_LON="-50"
MAX_LAT="5.5"
MIN_LAT="-20"

IN_DIR="/Users/rpavlick/github/JeDi/data/0.5degree"

for NCFILE in "elevation.nc" "landsea.nc" "paw.nc"; do
cdo -sellonlatbox,${MIN_LON},${MAX_LON},${MAX_LAT},${MIN_LAT} ${IN_DIR}/${NCFILE} ${NCFILE}

cdo -f srv -b F64 copy ${NCFILE} ${NCFILE%.nc}.srv

done
