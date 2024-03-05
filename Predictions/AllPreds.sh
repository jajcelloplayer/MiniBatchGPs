#!/bin/bash
# set -x

for i in `seq 1 52`;
	do echo This server has run $i datasets;

next="$(head -1 $1)" # Find the next data set to run

echo $next

if [ $next -gt 50 ]
then
  echo Done
  exit
fi

sed -i '1d' $1 

Rscript Final_Prediction_Script.R $next

done
