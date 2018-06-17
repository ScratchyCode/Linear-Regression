#!/bin/bash
DIR=$(dirname "$0")
cd $DIR
gcc -c stat.c
gcc Least_squares.c stat.o -lm -o fit.out
./fit.out
gnuplot script.gnuplot -p

sleep 1
rm fit.dat
rm data.dat
rm fit.out
rm stat.o
exit
