file title  : ISRIC-WISE Soil Moisture Retention (0-1 m)
file type   : ascii GRID format
columns     : 720
rows        : 360
ref. system : lat/long
ref. units  : deg
unit dist.  : 1.0000000
min. X      : -180.0000000
max. X      : 180.0000000
min. Y      : -90.0000000
max. Y      : 90.0000000
pos'n error : unknown
resolution  : 0.5000000
min. value  : 0
max. value  : 11
value units : classes
value error : unknown
flag value  : none
flag def'n  : none
legend cats : 12
category  0 : Background
category  1 : T<=60 mm/m    (> 66% of grid; mm m-1)
category  2 : 60<T<=90      (> 66% of grid; mm m-1)
category  3 : 90<T<=120     (> 66% of grid; mm m-1)
category  4 : 120<T<=150    (> 66% of grid; mm m-1)
category  5 : 150<T<=200    (> 66% of grid; mm m-1)
category  6 : 200<T<500     (> 66% of grid; mm m-1)
category  7 : T<=90     Complex (> 50% of grid; mm m-1)
category  8 : 90<T<=150 Complex (> 50% of grid; mm m-1)
category  9 : 150<=T    Complex  (> 50% of grid; mm m-1)
category 10 : Glaciers      (> 66% of grid; mm m-1)
category 11 : Oceans        (> 66% of grid; mm m-1)
comment     : T is abbreviation for Total Available Water Capacity (TAWC) in mm
comment     : to a depth of 100 cm (except for shallow Lithosols, Rankers and
comment     : Rendzinas).
lineage     : Source: 
lineage     : Batjes, N.H., 1996. Development of a world data set of soil water
lineage     : retention properties using pedotransfer rules. Geoderma (71), 31-5
lineage     : Derived from World Inventory of Soil Emission Potentials Database 
lineage     : Copyright, ISRIC, Wageningen (1996) 


To read data into ARC INFO or GRID, use the following commands:

ARC:   ASCIIGRID <in_ascii_file> <out_grid> {INT | FLOAT}
or
GRID:  newgrid =  ASCIIGRID (<in_ascii_file>, {INT | FLOAT})
