#!/usr/bin/gnuplot
set autoscale
set format y "%s"
set title "Linear fit"
set key reverse Left outside
set grid
#set style data linespoints
plot "data.dat" using 1:2:3 with yerrorbars title "Data", \
     "fit.dat"  using 1:2 title "Fit" with dots

exit
