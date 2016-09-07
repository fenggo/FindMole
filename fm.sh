#!/bin/bash


for i in `ls hug*.lammpstrj`

do 

echo $i >  in.fm
echo 1 >> in.fm
echo 5.0 >> in.fm
echo 1 >> in.fm


findmole<in.fm

j=`echo $i | sed 's/hug//g'`
jj=`echo $j | sed 's/.lammpstrj//g'`

mv molecular_structure.txt fm/mole$jj.txt
mv product.txt fm/pd$jj.txt

done

exit 0 

