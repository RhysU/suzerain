#!/bin/bash
# Given textual output from ./apps/perfect/perfect_mean,
# generate a PNG of boundary layer behavior on standard output.
title=$(awk '2 == NR { print $1; }' $1)
gnuplot << EOF
set title "$1, t = $title"
set term png noenhanced
set key autotitle columnheader
set logscale x
set xlabel 'distance from surface'
plot '$1' using 'y':'bar_rho' axis x1y1 with linespoints \
   , ''   using 'y':'bar_u'   axis x1y1 with linespoints \
   , ''   using 'y':'bar_v'   axis x1y1 with linespoints \
   , ''   using 'y':'bar_w'   axis x1y1 with linespoints \
   , ''   using 'y':'bar_T'   axis x1y1 with linespoints \
   , ''   using 'y':'bar_p'   axis x1y1 with linespoints \
   , ''   using 'y':'bar_H0'  axis x1y1 with linespoints
EOF
