# run all 3 problems on all 4 boundary types
# problem 3 has 3 separate cases, two of which require quarters=2
# assumes Gmsh and GetDP are in the path
# script with coarse mesh (s=1.5), can be changed below

#!/bin/bash

set -euo pipefail

getdp_run() {
    PROB="-setnumber prob $1"
    for ((b=1; b<=4; b++)); do
        getdp main.pro -solve ResMain -pos PostMain $OPT $PROB -setnumber bound $b
    done
}

OPT="-setnumber s 1.5 -setnumber quarters 1 -setnumber mags 1"
gmsh main.geo -3 $OPT
getdp_run 1
getdp_run 2
getdp_run 3
OPT="-setnumber s 1.5 -setnumber quarters 2 -setnumber mags 2"
gmsh main.geo -3 $OPT
getdp_run 3
OPT="-setnumber s 1.5 -setnumber quarters 2 -setnumber mags 3"
getdp_run 3

