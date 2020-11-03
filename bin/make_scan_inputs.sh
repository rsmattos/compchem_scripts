#!/bin/bash

rm geom.trj.xyz

num_folders=$(ls -l DISPLACEMENT | grep -c ^d -1)

for X in $(eval echo {0..$((num_folders-1))})
do
    geom2xyz.pl DISPLACEMENT/CALC.c1.d$X/geom
    mkdir -p scan/step_0$X
    cp orca42.pbs unbv.inp scan/step_0$X
    tail -n +3 DISPLACEMENT/CALC.c1.d$X/geom.xyz >> scan/step_0$X/unbv.inp
    cat DISPLACEMENT/CALC.c1.d$X/geom.xyz >> geom.trj.xyz 
    echo "*" >> scan/step_0$X/unbv.inp
    echo "" >> scan/step_0$X/unbv.inp
done
