# script for converting ISI-MIP forcing data to SRV format for JeDi input

oldres="original"
newres="0.5degree"
basedir=/home/rpavlick/nobackup0/FORCING/ISIMIP
oldmaskfile=${basedir}/${oldres}/landsea.nc
newmaskfile=${basedir}/${newres}/landsea.nc
gridfile=${basedir}/${newres}/griddes
#weightsfile=${basedir}/${newres}/remapweights_${oldres}_to_${newres}.nc
#cdo gencon,${gridfile} ${oldmaskfile} ${weightsfile}
stefanboltzmann="5.670373E-8"
emissivity="0.98"

for model in "HadGEM2-ES"; do
for scenario in "spinup"; do # "historical" "rcp8p5"
for var in "pr" "rlds" "rsds" "tas"; do
olddir=${basedir}/${oldres}/${model}/${scenario}

newdir=${basedir}/${newres}/${model}/${scenario}
mkdir -p ${newdir}

for ifile in ${olddir}/${var}*.nc4; do
ofile=${newdir}/${var}.srv
cdo -f srv -b F64 cat -ifthen ${newmaskfile} ${ifile} ${ofile}
done
done
cdo -add ${newdir}/rlds.srv -mulc,-1 -mulc,${emissivity} -mulc,${stefanboltzmann} -pow,4 ${newdir}/tas.srv ${newdir}/rlns.srv
rm ${newdir}/rlds.srv
done
done
