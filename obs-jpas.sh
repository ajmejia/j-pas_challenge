#!/usr/bin/csh

setenv  DYNBAS     $PWD/DynBaS
setenv  DB_MODELS  $PWD/models
setenv  DB_FILTERS $PWD/filters/FILTERBIN.RES

foreach file ( `cat $1` )
  if ( ! -e `echo "outs/$2/dynbasfit_$file" | sed "s/.txt/.log/"` ) then
    echo "fitting SED in $file"
    $DYNBAS/bin/dynbas --Zsun=0.02 --models!=Z0.0000 --output-dir="outs/$2" -v --L-passband=17 --passbands=1,...,58 "$2/$file"
  endif
end
