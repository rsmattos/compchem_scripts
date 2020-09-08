#!/bin/bash

rm geom.trj.xyz

for X in $(eval echo {1..$1})
do
    geom2xyz.pl DISPLACEMENT/CALC.c1.d$X/geom
    mv geom.xyz DISPLACEMENT/CALC.c1.d$X/geom.xyz

    mkdir -p scan/step_0$X
    cp orca42.pbs scan/step_0$X
    cp unbv.inp scan/step_0$X
    tail -n +3 DISPLACEMENT/CALC.c1.d$X/geom.xyz >> scan/step_0$X/unbv.inp
    cat DISPLACEMENT/CALC.c1.d$X/geom.xyz >> geom.trj.xyz 
    echo "*" >> scan/step_0$X/unbv.inp
done
