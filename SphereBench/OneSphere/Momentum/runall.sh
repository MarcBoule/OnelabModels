# run all 7 problems on all 2 boundary types
# assumes Gmsh and GetDP are in the path
# script with coarse mesh (s=1.5), can be changed below

#!/bin/bash

set -euo pipefail

OPT="-setnumber s 1.5"
gmsh main.geo -3 $OPT

getdp_run() {
    PROB="-setnumber prob $1"
    getdp main.pro -solve ResMain -pos PostMain $OPT $PROB -setnumber iabc 1
    getdp main.pro -solve ResMain -pos PostMain $OPT $PROB -setnumber iabc 0
}

for ((p=1; p<=7; p++)); do
    getdp_run $p
done

