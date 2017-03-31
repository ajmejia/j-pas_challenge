#!/usr/bin/csh

setenv  DYNBAS     $PWD/DynBaS
setenv  DB_MODELS  $PWD/models
setenv  DB_FILTERS $PWD/filters/FILTERBIN.RES

foreach file ( `cat $1/$2` )
  if ( ! -e `echo "outs/$1/dynbasfit_$file" | sed "s/.txt/.log/"` ) then
    echo "fitting SED in $file"
    $DYNBAS/bin/dynbas --Zsun=0.02 --models!=Z0.0000 --output-dir="outs/$1" -v --L-passband=66 --passbands=61,...,64,66,68,70 "$1/$file"
  endif
end
