#!/bin/bash

gfortran -O3 -o fbm2d fbm2d.f90

l=0032
h=0.7
s=1231

./fbm2d "$l" "$h" "$s" l"$l"_h"$h"_s"$s"
