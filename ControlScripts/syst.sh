#!/bin/bash

for i in {1..503}
do
  sed -i '/UncertType: */d' Set_$i/Input.txt
  sed -i '/Analysis: HCPMuTau/a UncertType: default' Set_$i/Input.txt
  sed -i -e '/UncertType/{r syst.txt' -e 'd;}' Set_$i/Input.txt
done

