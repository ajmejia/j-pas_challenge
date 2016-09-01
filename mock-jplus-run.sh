# $1: script
# $2: SED flavor (e.g., spec, jpas)
# $3: list name (e.g., danzon, chapi)

setenv  DYNBAS     ./DynBaS
setenv  DB_MODELS  ./models
setenv  DB_FILTERS ./jplus-filters/JPLUS_FILTERBIN.RES

nohup ./$1 $2/$3_00 > outs/$2/$3_output_00.log &
nohup ./$1 $2/$3_01 > outs/$2/$3_output_01.log &
nohup ./$1 $2/$3_02 > outs/$2/$3_output_02.log &
nohup ./$1 $2/$3_03 > outs/$2/$3_output_03.log &
nohup ./$1 $2/$3_04 > outs/$2/$3_output_04.log &
nohup ./$1 $2/$3_05 > outs/$2/$3_output_05.log &
nohup ./$1 $2/$3_06 > outs/$2/$3_output_06.log &
nohup ./$1 $2/$3_07 > outs/$2/$3_output_07.log &
nohup ./$1 $2/$3_08 > outs/$2/$3_output_08.log &
nohup ./$1 $2/$3_09 > outs/$2/$3_output_09.log &

nohup ./$1 $2/$3_10 > outs/$2/$3_output_10.log &
nohup ./$1 $2/$3_11 > outs/$2/$3_output_11.log &
nohup ./$1 $2/$3_12 > outs/$2/$3_output_12.log &
nohup ./$1 $2/$3_13 > outs/$2/$3_output_13.log &
nohup ./$1 $2/$3_14 > outs/$2/$3_output_14.log &
nohup ./$1 $2/$3_15 > outs/$2/$3_output_15.log &
nohup ./$1 $2/$3_16 > outs/$2/$3_output_16.log &
nohup ./$1 $2/$3_17 > outs/$2/$3_output_17.log &
nohup ./$1 $2/$3_18 > outs/$2/$3_output_18.log &
nohup ./$1 $2/$3_19 > outs/$2/$3_output_19.log &

nohup ./$1 $2/$3_20 > outs/$2/$3_output_20.log &
nohup ./$1 $2/$3_21 > outs/$2/$3_output_21.log &
nohup ./$1 $2/$3_22 > outs/$2/$3_output_22.log &
nohup ./$1 $2/$3_23 > outs/$2/$3_output_23.log &
nohup ./$1 $2/$3_24 > outs/$2/$3_output_24.log &
nohup ./$1 $2/$3_25 > outs/$2/$3_output_25.log &
nohup ./$1 $2/$3_26 > outs/$2/$3_output_26.log &
nohup ./$1 $2/$3_27 > outs/$2/$3_output_27.log &
nohup ./$1 $2/$3_28 > outs/$2/$3_output_28.log &
nohup ./$1 $2/$3_29 > outs/$2/$3_output_29.log &
