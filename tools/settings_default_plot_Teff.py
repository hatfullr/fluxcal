# default settings_plot_Teff.py
# DO NOT MOVE THIS FILE.
# Setting usedefault = True will only work when running plot_Teff.py from this directory.

usedefault = True # Set to True to use this file. All other settings files will be omitted.

runit=6.9598e10
#runit = 1.e0
tunit = 60.*60.*24. # days default
tempunit = 1e0
# Axis limits activate when NOT using adaptive plotting
# Setting limits to (0.,0.) will tell plot_Teff.py to find its own
# limits and stick to those. The limits it uses are the absolute
# min and max of all data fed to the program.
xlim = (None,None)
ylim = (None,None)
Tefflim = (None,None)
logteff = False
cmapname = 'nipy_spectral' # Ignored when usebb is False
usebb = False              # Use blackbody spectrum
spatialaxesadapt = False    # Use adaptive spatial limits
colorbaradapt = False      # Use adaptive colorbar limits
filepattern = "" # Search pattern for data files. Set to "" to prompt every time
paperfriendly = False      # Paper friendly plot [Experimental]

fontsize = 12
figsize = (8.0,8.0) # In inches
xlabel = "$x$"
ylabel = "$y$"
cbarlabel = "$T_{\\mathrm{eff}}\\ \\left[K\\right]$"
timelabel = "t = {:15.3f}"
rotationslabel = "Rotations (x,y,z) = ({:06.2f}$^o$,{:06.2f}$^o$,{:06.2f}$^o$)"
