#!/bin/bash

for X in 0 1 2 3 4 5 6 7
do
    cp dens_ana.in step_0$X/
    cd step_0$X/
    orca_2mkl unbv -molden
    cp unbv.cis orca.cis
    analyze_tden.py
    cd ..
done

~/Scripts/compchem_scripts/bin/scan.py . --descriptors -g -s svg --noshow
