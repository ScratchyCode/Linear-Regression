#!/bin/bash
DIR=$(dirname "$0")
cd $DIR

gcc -c stat.c
gcc Least_squares.c stat.o -lm -o fit.out
./fit.out
gnuplot script.gnuplot -p 2>/dev/null

sleep 2

rm fit.out
rm stat.o

exit
