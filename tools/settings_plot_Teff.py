# settings_plot_Teff.py
# You can copy this file to another directory to use it by default there.
# If this doesn't exist in the directory from which plot_Teff was called,
# plot_Teff will use this file instead.

runit=6.9598e10
#runit = 1.e0
tunit = 60.*60.*24. # seconds -> days
tempunit = 1e0

# Axis limits activate when NOT using adaptive plotting
# Setting limits to (None,None) will tell plot_Teff.py to find its own
# limits and stick to those. The limits it uses are the absolute
# min and max of all data fed to the program.
xlim = (None,None)
ylim = (None,None)
Tefflim = (0.,10000.)

# Set the following to None to prompt the user
usebb = False              # Use blackbody spectrum.
spatialaxesadapt = None   # Use adaptive spatial limits.
colorbaradapt = None      # Use adaptive colorbar limits.
filepattern = None        # Search pattern for data files.
paperfriendly = False      # Paper friendly plot [Experimental].

# Miscellaneous plotting options
logteff = False
cmapname = 'nipy_spectral' # Ignored when usebb is False
fontsize = 16
figsize = (8.0,8.0) # In inches
xlabel = "$x\\ \\left[R_\\odot\\right]$"
ylabel = "$y\\ \\left[R_\\odot\\right]$"
cbarlabel = "$T_{\\mathrm{eff}}\\ \\left[K\\right]$"
timelabel = "t = {:15.3f} days"
rotationslabel = "Rotations (x,y,z) = ({:06.2f}$^o$,{:06.2f}$^o$,{:06.2f}$^o$)"
