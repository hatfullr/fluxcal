# settings_plot_Teff.py

#runit=6.9598e10
runit = 1.e0
tunit = 1.e0/60./60./24.
# Axis limits activate when NOT using adaptive plotting
# Setting limits to (0.,0.) will tell plot_Teff.py to find its own
# limits and stick to those. The limits it uses are the absolute
# min and max of all data fed to the program.
xlim = (0.,0.)
ylim = (0.,0.)
Tefflim = (4300.,10000.)
cmapname = 'nipy_spectral' # Ignored when usebb is False
usebb = False              # Use blackbody spectrum
spatialaxesadapt = True    # Use adaptive spatial limits
colorbaradapt = False      # Use adaptive colorbar limits
filepattern = "teffs*.dat" # Search pattern for data files
paperfriendly = False      # Paper friendly plot
