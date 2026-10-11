# run all 8 problems on all 4 boundary types
# assumes Gmsh and GetDP are in the path
# script with coarse mesh (s=1.5), can be changed below

#!/bin/bash

set -euo pipefail

OPT="-setnumber s 1.5"
gmsh main.geo -3 $OPT

getdp_run() {
    PROB="-setnumber prob $1"
    for ((b=1; b<=4; b++)); do
        getdp main.pro -solve ResMain -pos PostMain $OPT $PROB -setnumber bound $b
    done
}

for ((p=1; p<=8; p++)); do
    getdp_run $p
done

