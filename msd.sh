#!/bin/bash

for i in `ls hug*.lammpstrj`

do
   j=`echo $i | sed 's/hug//g'`
   jj=`echo $j | sed 's/.lammpstrj//g'`
   echo $i > in.t~
   echo 1 >> in.t~
   echo 5.0 >> in.t~
   echo 3 >> in.t~
   echo 5 >> in.t~
   findmole<in.t~
   mv msd.txt msd$jj
done

exit 0
