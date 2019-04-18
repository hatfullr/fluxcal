# Using FluxCal

## Installation

### Dependencies
FluxCal uses `gfortran` to compile its code. Install using,
```
sudo apt install gfortran
```

### Linux
To install FluxCal on Linux, run the following commands from the terminal
```
git clone https://github.com/hatfullr/fluxcal.git
cd fluxcal
./install
```
DO NOT RUN AS ROOT. This script will unzip the required included data files,
build `flux_cal`, and link the plotting suite to your home directory. We strongly
recommend you perform a file integrity check by running
```
./integrity_check
```
This may take up to 10 minutes or so depending on your system. If you run into
issues here, a solution may be available in the "Troubleshooting" section below.

### Windows
FluxCal has not been tested on Windows machines, but it might still work. Try
[installing the Ubuntu terminal on Windows 10](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0), then continue installation as though you were using
Linux.

## Uninstallation
To uninstall FluxCal, run the script `uninstall`, then do
```
echo "cd .."
echo "rm -rf fluxcal"
```
The `uninstall` script only removes the links to the plotting utilities in your
`$HOME/bin` directory. All other files pertaining to FluxCal are located in the
directory `fluxcal`. Thus, to finish installation,
```
cd ..
rm -rf fluxcal
```


## Running FluxCal

To run FluxCal, try the following
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
You can find a premade input file to test with in `<fluxcal-directory>/test/standard_check/original/fluxcal_0000.dat`. Do not move this file.


## Plotting utilities

FluxCal comes prepackaged with a suite of plotting utilities which can be found
at `<fluxcal-directory>/tools/`,
```
plot_Teff
plot_closest
plot_lc
```
Upon installation, these are linked to `$HOME/bin/.` and thus, as long as this
directory is in your `PATH` variable, you can run any of these simply as
`plot_Teff`, `plot_closest`, and `plot_lc`. Please note that `plot_closest` and
`plot_lc` are not as sophisticated as `plot_Teff`.


# Troubleshooting

## Failing compilation on install

There is a known bug with `gfortran-7.3.0` where compilation fails upon trying
to build `mathlib/odeint.o`,
```
/usr/bin/gfortran -O4 -ffixed-line-length-132 -mcmodel=large  -c -o mathlib/odeint.o mathlib/odeint.f
mathlib/odeint.f:186:0:

       END
 
internal compiler error: Segmentation fault
Please submit a full bug report,
with preprocessed source if appropriate.
See <file:///usr/share/doc/gcc-7/README.Bugs> for instructions.
<builtin>: recipe for target 'mathlib/odeint.o' failed
make: *** [mathlib/odeint.o] Error 1

Failed
```
If you encounter this error, try compiling again,
```
./install
```
This may fail up to 10 times before finally succeeding. Hopefully this error
disappears with future `gfortran` releases. We do not recommend downgrading
your GNU compiler suite, as this is generally harmful.


## Failing standard_check test

When running the "standard_check" test, you may get files that differ
from the original. This could be due to a difference in precision, which
is controlled by the `gfortran` compiler. Check to make sure you have the same
version compiler installed when running the test. The version of `gfortran`
used to run the original test is located at
`<fluxcal-directory>/test/standard_check/original/GNUFortran.txt`.


