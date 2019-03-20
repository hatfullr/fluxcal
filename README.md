# Installing

To install FluxCal, run the following commands from the terminal
```
git clone https://github.com/hatfullr/fluxcal.git
cd fluxcal
./install
```


# Using FluxCal

To run FluxCal, try the following:
```
mkdir fluxcaltest
cd fluxcaltest
ln -s <fluxcal-directory>/flux_cal .
cp <fluxcal-directory>/defaults/flux_cal.input .
cp <datafile> .
```
Now edit `flux_cal.input` appropriately, save the file, and run `flux_cal`
with
```
./flux_cal
```
You can find a premade input file to test with in `<fluxcal-directory>/test/standard_check/test/fluxcal_0000.dat`.




# Troubleshooting

When running the "standard_check" test, you may get files that differ
from the original. This could be due to a difference in precision, which
is controlled by the compiler. Check to make sure you have the same
version compiler installed when running the test. The version of GNU
Fortran used to run the "original" directory can be found in the file
"GNUFortran.txt" there.

