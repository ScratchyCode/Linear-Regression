# Linear-regression
Performs the linear fit of a dataset with the least squares method.

Create a file with following format: x y sigma(y)

CAUTION: the file's lines must be equal to the number of points!

# Usage
To run program under GNU/Linux simply execute bash script with: bash start.sh

Under Windows first compile stat.c library with: gcc -c stat.c -lm

After compile Least_Squares.c with: gcc Least_Squares.c stat.o -lm
