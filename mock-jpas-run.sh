# $1: script
# $2: SED flavor (e.g., spec, jpas)
# $3: list name (e.g., danzon, chapi)

nohup ./$1 $2/$3_00 $2 > outs/$2/$3_output_00.log &
nohup ./$1 $2/$3_01 $2 > outs/$2/$3_output_01.log &
