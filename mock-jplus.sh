#!/usr/bin/csh

foreach file ( `cat $1` )
  if ( ! -e `echo "outs/$2/dynbasfit_$file" | sed "s/.txt/.log/"` ) then
    echo "fitting SED in $file"
    $DYNBAS/bin/dynbas --Zsun=0.0199 --models!=Z0.0000 --output-dir=outs/$2 -v --L-passband=9 --passbands=2,...,6,8,9,10,12,14 $2/$file
  endif
end
