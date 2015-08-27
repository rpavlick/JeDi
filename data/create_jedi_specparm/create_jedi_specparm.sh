ENSEMBLE_SIZE=20
TRAITS=15
for SPECIES in "10" "20" "50" "100" "200" "500" "1000" "2000"
do

for ((a=1; a <= ENSEMBLE_SIZE ; a++))
do
SEED=`expr $a + $SPECIES`
./latin_random.x $TRAITS $SPECIES $SEED
done
done
