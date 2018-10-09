#!/usr/bin/gnuplot
set autoscale
set format y "%s"
set title "Linear fit"
set key reverse Left outside
set grid
#set style data linespoints
f(x) = a*x + b
fit f(x) "data.dat" using 1:2 via a,b
plot "data.dat" using 1:2:3 with yerrorbars, f(x)

exit
