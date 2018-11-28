# FluxCal

dimen*.dat files contain details about the 2D driving grid in cgs units.
Each file contains one line with the minimum x position, width of the
cells, number of cells in the x-direction, minimum y position, height of
the cells, and number of cells in the y-direction.


# Testing

When running the "standard_check" test, you may get files that differ
from the original. This could be due to a difference in precision, which
is controlled by the compiler. Check to make sure you have the same
version compiler installed when running the test. The "original"
directory was run with

GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.10) 5.4.0 20160609

