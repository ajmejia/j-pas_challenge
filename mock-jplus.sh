#!/usr/bin/bash

export DYNBAS="$PWD/DynBaS"
export DB_MODELS="$PWD/models"
export DB_FILTERS="$PWD/filters/FILTERBIN.RES"

while read -r file
do
  if [ ! -f `echo "outs/$1/dynbasfit_$file" | sed "s/.txt/.log/"` ]
  then
    echo "fitting SED in $file"
    time $DYNBAS/bin/dynbas --Zsun=0.02 --models!=Z0.0000 --output-dir="outs/$1" -v --L-passband=66 --passbands=60,...,71 "$1/$file"
  fi
done < $1/$2
