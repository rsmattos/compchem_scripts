#!/bin/sh

basis="sv(p)"
scfiter=300
functional="wb97x-d"
grid="m3"
scfconv=7
exstates=10
exopt=""
memory=5300

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

for x in "${arr[@]}"
do
    echo "$x"
    sleep 0.1
done | define

if [ ! -z "$exopt" ]
  then
    head -n -1 control > temp; mv temp control
    echo "\$exopt $exopt" >> control
    echo "\$end" >> control
fi
