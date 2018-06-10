#!/bin/bash
DIR=$(dirname "$0")
cd $DIR
gcc Least_squares.c -lm -o fit.out
./fit.out
gnuplot script.gnuplot -p

sleep 1
rm fit.dat
rm data.dat
rm fit.out
exit
