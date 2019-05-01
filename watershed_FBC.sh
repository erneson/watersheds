#!/bin/bash

gcc -O3 -o watershed_FBC watershed_FBC.c

filename='fbm2d_l0032_h0.7_s1231'
hcutoff='0000'

./watershed_FBC $filename $hcutoff

./plot.py $filename $hcutoff
