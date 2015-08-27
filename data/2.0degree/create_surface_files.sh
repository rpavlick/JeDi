mv fractional_land.2-deg.nc landfraction.nc
cdo griddes -sellonlatbox,-180,180,90,-90 landfraction.nc > griddes
cdo remapbil,griddes ../0.5degree/paw.nc paw.nc
cdo remapnn,griddes  ../0.5degree/glacier.nc glacier.nc
cdo -mul glacier.nc -sellonlatbox,-180,180,90,-90 -setclonlatbox,0,0,360,-57,-90 -gtc,0.2 landfraction.nc landsea.nc
cdo -ifthenc,1.0 -gtc,0 landsea.nc elevation.nc

for NCFILE in "elevation.nc" "landsea.nc" "paw.nc"; do
 cdo -f srv -b F64 copy ${NCFILE} ${NCFILE%.nc}.srv
done
