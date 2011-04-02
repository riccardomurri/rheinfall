#! /bin/sh

wget -nv --recursive -l 2 --continue \
  --no-host-directories --cut-dirs=2 \
  -D ljk.imag.fr -I 'membres/Jean-Guillaume.Dumas/Matrices' \
  -A '.sms.gz',BIBD,Homology,Grobner,Forest,G5,Kocay,Relat,CAG,Margulies,Taha,Franz,Trefethen,GL7d,Cellules/diffSL6,Cellules/diffGL6,SPG,Mgn,Smooshed \
  http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/Matrices/ 
