#!/usr/bin/csh

setenv  DYNBAS     $PWD/DynBaS
setenv  DB_MODELS  $PWD/models
setenv  DB_FILTERS $PWD/jpas-filters/JPAS_FILTERBIN.RES

foreach file ( `cat $1` )
  if ( ! -e `echo "outs/$2/dynbasfit_$file" | sed "s/.txt/.log/"` ) then
    echo "fitting SED in $file"
    $DYNBAS/bin/dynbas --Zsun=0.0199 --models!=Z0.0000 --output-dir="outs/$2" -v --L-passband=15 --passbands=1,...,55 "$2/$file"
  endif
end
