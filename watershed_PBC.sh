#!/bin/bash

gcc -O3 -o watershed_PBC watershed_PBC.c

filename='GEBCO_2014_2D_90N_90S_f8'
hcutoff='0000'

./watershed_PBC $filename $hcutoff

./plot.py $filename $hcutoff
