#!/bin/bash

basis="sv(p)"
scfiter=300
functional="wb97x-d"
grid="m3"
scfconv=7
exstates=10
exopt=""
memory=5300
use_symm=false
use_dft=true
use_ri=true

# do ground state optimization
# excited state optimization

while test $# -gt 0
do
  if [[ "$1" == "basis" ]]; then
    shift; basis="$1"; shift
  elif [[ "$1" == "functional" ]]; then
    shift; functional="$1"; shift
  elif [[ "$1" == "grid" ]]; then
    shift; grid="$1"; shift
  elif [[ "$1" == "scfconv" ]]; then
    shift; scfconv=$1; shift
  elif [[ "$1" == "memory" ]]; then
    shift; memory=$1; shift
  elif [[ "$1" == "scfiter" ]]; then
    shift; scfiter=$1; shift
  elif [[ "$1" == "states" ]]; then
    shift; exstates=$1; shift
  elif [[ "$1" == "exopt" ]]; then
    shift; exopt=$1; exstates=$(($1+3)); shift
  else shift
  fi
done

declare -a arr=("" "" 
                "a coord" 
                "ired"
                "desy"
                "*"

                "b" "all $basis" "" ""
                "*" 

                "eht" "" "" ""

                "scf" "iter" "$scfiter" "conv" "$scfconv" "" 
                "dft" "on" "func $functional" "grid $grid" "q" 
                "ri" "on" "q"
                "ex" "rpas" "q" "a $exstates" "q" "q" ""
                "cc" "memory $(echo "$memory*0.7 / 1" | bc)" "*"
                "q")

rm control

input=("" "" "a coord")
if [ "$use_symm" = true ]; then
    iniput+=("desy")
fi
input+=("ired" "*")
input+=("b" "all $basis" "" "" "*"
input+=("eht" "" "" "")
input+=("scf" "iter" "$scfiter" "conv" "$scfconv" "")
if [ "$use_dft" = true ]; then
    input+=("dft" "on" "func $functional" "grid $grid" "q")
fi
if [ "$use_ri" = true ]; then
    input+=("ri" "on" "q")
fi
if [ "$vexc" = true ]; then
    input+=("ex" "rpas" "q" "a $exstates" "q" "q" "")
fi
input+=("cc" "memory $(echo "$memory*0.7 / 1" | bc)" "*")
input+=("q")

for x in "${input[@]}"
do
    echo "$x"
    sleep 0.1
done
#| define

if [ ! -z "$exopt" ]
  then
    head -n -1 control > temp; mv temp control
    echo "\$exopt $exopt" >> control
    echo "\$end" >> control
fi
