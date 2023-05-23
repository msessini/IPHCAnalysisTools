#!/usr/bin/bash

for i in {1..503}
do
  sed -i 's/que=.*/que="'$1'"/g' Set_$i/Qsub_Set_$i
done

