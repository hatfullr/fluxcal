# FluxCal

dimen*.dat files contain details about the 2D driving grid in cgs units.
Each file contains one line with the minimum x position, width of the
cells, number of cells in the x-direction, minimum y position, height of
the cells, and number of cells in the y-direction.

# Installing

To install FluxCal, run the following commands from the terminal
```
git clone https://github.com/hatfullr/fluxcal.git
cd fluxcal
./install
```

# Testing

When running the "standard_check" test, you may get files that differ
from the original. This could be due to a difference in precision, which
is controlled by the compiler. Check to make sure you have the same
version compiler installed when running the test. The version of GNU
Fortran used to run the "original" directory can be found in the file
"GNUFortran.txt" there.

