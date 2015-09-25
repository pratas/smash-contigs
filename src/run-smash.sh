#!/bin/bash
#
# Run: run-smash $1 $2
#
# $1 = contigs   filename
# $2 = reference filename
#
./smash-map    -v -n 4 -r 100 $1 $2
./smash-reduce -v -e 10 match.map
./smash-visual -v -l 3 -c 1000 match.map.red 
