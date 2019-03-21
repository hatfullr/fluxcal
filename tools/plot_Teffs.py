import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import matplotlib.colors as cols
import glob
import sys
import matplotlib as mpl
import os.path

# This allows the program to be used by both python2 and python3
try:
    input = raw_input
except NameError:
    pass

#runit=6.9598e10
runit=1.e0
tunit=1.e0
xlim = (None,None)
ylim = (None,None)
Tefflim = (None,None)
cmapdefault = 'nipy_spectral'
usebb = None
spatialaxesadapt = None
colorbaradapt = None
filepattern = ""
paperfriendly = False

if os.path.isfile('settings_plot_Teff.py'):
    print("Importing settings_plot_Teff.py")
    from settings_plot_Teff import *

files = []
while len(files) == 0:
    if filepattern == "":
        filepattern = input("Enter file name(s) or patterns: ")
    filepattern = filepattern.split(" ")
    #if not isinstance(filepattern,list): filepattern = [filepattern]
    for pattern in filepattern:
        fs = sorted(glob.glob(pattern))
        if len(fs) > 0:
            for i in fs:
                if ((i[:5] != "teffs") or (i[-4:] != ".dat") or (len(i)<=0)):
                    print("ERROR: You must provide teffs files from FluxCal output.")
                    print("Error on file",i)
                else:
                    files.append(i)
        else:
            print("ERROR: File pattern '"+pattern+"' not found")
    

if usebb == None:
    loop = True
    while loop:
        user_input = input("Use blackbody color spectrum? (y/n): ")
        if user_input == "y":
            usebb = True
            loop = False
        elif user_input == "n":
            loop = False
            
if((spatialaxesadapt == None) and (colorbaradapt == None)):
    loop = True
    while loop:
        user_input = input("Use adaptive plot limits? (y/n): ")
        if user_input == "y":
            #adaptivelimits = True
            loop2 = True
            while loop2:
                user_input = input("On spatial axes? (y/n): ")
                if user_input == "y":
                    spatialaxesadapt = True
                    loop2 = False
                elif user_input == "n":
                    loop2 = False
                else:
                    print("ERROR: Please respond with 'y' or 'n'.")
            loop2 = True
            while loop2:
                user_input = input("On colorbar? (y/n): ")
                if user_input == "y":
                    colorbaradapt = True
                    loop2 = False
                elif user_input == "n":
                    loop2 = False
                else:
                    print("ERROR: Please respond with 'y' or 'n'.")
            loop = False
        elif user_input != "n":
            print("ERROR: Please respond with 'y' or 'n'.")

if paperfriendly == None:
    loop = True
    paperfriendly = False
    while loop:
        user_input = input("Paper friendly? (y/n): ")
        if user_input == "y":
            paperfriendly = True
            loop = False
        elif user_input == "n":
            loop = False
        else:
            print("ERROR: Please respond with 'y' or 'n'.")

mpl.rcParams['font.size'] = 12
mpl.rcParams['figure.figsize'] = (8.0,8.0)
if paperfriendly:
    # This format is appropriate for ApJ style papers that have
    # two columns of text per page. Each column has a width of
    # 3.5225 inches.
    #mpl.rcParams['font.size'] = 10
    #mpl.rcParams['figure.figsize'] = (3.35225,2.5)
    mpl.rcParams['text.color'] = 'black'
    mpl.rcParams['axes.facecolor'] = '0.95'
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['xtick.color'] = 'black'
    mpl.rcParams['ytick.color'] = 'black'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['figure.edgecolor'] = 'white'
    mpl.rcParams['figure.subplot.left'] = 0.15
    mpl.rcParams['figure.subplot.right'] = 0.80
    mpl.rcParams['figure.subplot.bottom'] = 0.12
    mpl.rcParams['figure.subplot.top'] = 1.0
    mpl.rcParams['axes.titlepad'] = 0.0
else:
    #mpl.rcParams['font.size'] = 12
    #mpl.rcParams['figure.figsize'] = (8.0,8.0)
    mpl.rcParams['text.color'] = 'white'
    mpl.rcParams['axes.facecolor'] = 'black'
    mpl.rcParams['axes.edgecolor'] = 'white'
    mpl.rcParams['axes.labelcolor'] = 'white'
    mpl.rcParams['xtick.color'] = 'white'
    mpl.rcParams['ytick.color'] = 'white'
    mpl.rcParams['figure.facecolor'] = 'black'
    mpl.rcParams['figure.edgecolor'] = 'black'
    mpl.rcParams['figure.subplot.right'] = 0.85
    mpl.rcParams['figure.subplot.bottom'] = 0.10
    mpl.rcParams['figure.subplot.top'] = 0.95
    mpl.rcParams['figure.titleweight'] = 300


# Setup the matplotlib rcParams
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'monospace'
mpl.rcParams['font.weight'] = 300
mpl.rcParams['font.monospace'] = 'DejaVu Sans'
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.major.size'] = 9
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.width'] = 0.5
mpl.rcParams['xtick.labelsize'] = 'medium'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.minor.visible'] = True

mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.major.size'] = 9
mpl.rcParams['ytick.minor.size'] = 4
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.width'] = 0.5
mpl.rcParams['ytick.labelsize'] = 'medium'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.minor.visible'] = True

#mpl.rcParams['font.size'] = 12
#mpl.rcParams['figure.titleweight'] = 300
#mpl.rcParams['figure.figsize'] = (8.0, 8.0)
#mpl.rcParams['figure.subplot.right'] = 0.85
#mpl.rcParams['figure.subplot.bottom'] = 0.10
#mpl.rcParams['figure.subplot.top'] = 0.95

x0win = 0.125
y0win = 0.18273809523809526
x1win = 0.8095238095238095
y1win = 0.8672619047619047

        
# Create the plot

fig = plt.figure()
ax = plt.gca()
ax.set_position([x0win,y0win,x1win-x0win,y1win-y0win])

        
# Create the colorbar

vmin = 0.
vmax = 10000.
norm = cols.Normalize(vmin,vmax)

if Tefflim[0] != None:
    vmin = Tefflim[0]
if Tefflim[1] != None:
    vmax = Tefflim[1]

mlibcolorbar=False
if usebb:
    try:
        print("Getting live website data for blackbody spectrum colors.")
        import urllib2
        website = 'http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html'
    
        page = urllib2.urlopen(website)
        for i in range(0,20): page.readline() # Skip the header
    
        content = page.readlines()[1:-4][::2]
        content = np.asarray(content)
    
        #content = np.delete(content,0,axis=0)
        bbcols = np.zeros((len(content),4))
        for i in range(0,len(content)):
            content[i] = content[i][34:-17]
            bbcols[i][0] = float(content[i][:5])
            bbcols[i][1] = float(content[i][44:50])
            bbcols[i][2] = float(content[i][51:57])
            bbcols[i][3] = float(content[i][58:64])
    
        Nbbcols = len(bbcols)
        colors = np.reshape(bbcols[:,1:],(Nbbcols,3))
        colors = colors[np.where((bbcols[:,0] > vmin) & (bbcols[:,0] < vmax))[0]]
        my_cm = cols.LinearSegmentedColormap.from_list(
            'blackbody', colors, N=Nbbcols)
        print("SUCCESS")

    except Exception as e:
        print("FAILED:", e)
        try:
            print("Getting blackbody spectrum colors from blackbody_temps.dat in the working directory.")
            bbcols = np.genfromtxt("blackbody_temps.dat",usecols=(0,6,7,8))
            Nbbcols = len(bbcols)
            colors = np.reshape(bbcols[:,1:],(Nbbcols,3))
            colors = colors[np.where((bbcols[:,0] > vmin) & (bbcols[:,0] < vmax))[0]]
            my_cm = cols.LinearSegmentedColormap.from_list(
                'blackbody', colors, N=Nbbcols)
            print("SUCCESS")
        except Exception as e:
            print("FAILED:",e)
            print("Cannot use the blackbody color spectrum. Please fix one of the above errors or choose 'n' at the blackbody color spectrum prompt.")
            sys.exit()
else:
    print("Using the default colorbar '"+cmapdefault+"'.")
    try:
        mlibcolorbar = True
        my_cm = plt.cm.get_cmap(cmapdefault)
        print("SUCCESS")
    except Exception as e:
        print("FAILED:",e)
        print("ERROR: Please fix the above error and try again.")
        sys.exit()

if paperfriendly:
    my_cm.set_bad('white',1.)
else:
    my_cm.set_bad('black',1.)
    
im = plt.imshow([[np.nan],[np.nan]],extent=(ax.get_xlim()[0],ax.get_xlim()[1],ax.get_ylim()[0],ax.get_ylim()[1]),cmap=my_cm)
pos = ax.get_position()
cbwidth = 0.03
cax = fig.add_axes([pos.xmax+0.01,pos.ymin,cbwidth,(pos.ymax-pos.ymin)])
caxposorig = cax.get_position()
cb = plt.colorbar(im,cax=cax)
plt.clim(vmin,vmax)
cb.ax.minorticks_off()

cbar_label = "$T_{\\mathrm{eff}}\\ \\left[K\\right]$"
cbstate=0

if not spatialaxesadapt or not colorbaradapt:
    # Read all the datafiles and find the maximum viewport
    xmaxx = -1e30
    ymaxx = -1e30
    xminn = 1e30
    yminn = 1e30
    aminn = 1.e30
    amaxx = 0.
    for datafile in files:
        grid = np.genfromtxt(datafile,max_rows=1)

        xmin = grid[0]*runit
        hxmap = grid[1]*runit
        nx = int(grid[2])
        ymin = grid[3]*runit
        hymap = grid[4]*runit
        ny = int(grid[5])
        
        data = 10.**(np.genfromtxt(datafile,skip_header=1,max_rows=ny))
        data[np.isinf(data)] = np.nan # Set all inf to NaN
        data[np.isneginf(data)] = np.nan # Set all -inf to NaN
        data[np.where(data == 0)] = np.nan # Set all T=0 to NaN

        aminn = min(np.amin(data[np.isfinite(data)]),aminn)
        amaxx = max(np.amax(data[np.isfinite(data)]),amaxx)
        
        xmax = xmin + hxmap * nx
        ymax = ymin + hymap * ny

        if xmin < xminn: xminn = xmin
        if ymin < yminn: yminn = ymin
        if xmax > xmaxx: xmaxx = xmax
        if ymax > ymaxx: ymaxx = ymax


    dx = xmaxx - xminn
    dy = ymaxx - yminn
    if dx > dy:
        ymaxx = yminn + (dy+dx)/2.
        yminn = yminn + (dy-dx)/2.
    elif dx < dy:
        xmaxx = xminn + (dx+dy)/2.
        xminn = xminn + (dx-dy)/2.

triggered = [ False, False, False ]
pos=[]
for datafile in files:
    
    # Get the grid data

    hdr = np.genfromtxt(datafile,max_rows=1)
    
    xmin = hdr[0]*runit
    hxmap = hdr[1]*runit
    nx = int(hdr[2])
    ymin = hdr[3]*runit
    hymap = hdr[4]*runit
    ny = int(hdr[5])
    t = hdr[6] #in seconds
    t = t*tunit # in days
    
    xmax = xmin + hxmap * nx
    ymax = ymin + hymap * ny
    
    # Get the Teff data
    
    data = 10.**(np.genfromtxt(datafile,skip_header=1,max_rows=ny))

    data[np.isinf(data)] = np.nan # Set all inf to NaN
    data[np.isneginf(data)] = np.nan # Set all -inf to NaN
    data[np.where(data == 0.)] = np.nan # Set all T=0 to NaN

    # Maintain plot axis dimensions
    if type(pos) == 'matplotlib.transforms.Bbox':
        ax.set_position([pos.x0,pos.y0,pos.width,pos.height])
        
    amax = np.amax(data[np.isfinite(data)])
    amin = np.amin(data[np.isfinite(data)])
    if spatialaxesadapt:
    
        # Center and square the plotted region
    
        dx = xmax - xmin
        dy = ymax - ymin
        
        if dx > dy:
            ax.set_ylim(ymin + (dy-dx)/2., ymin + (dy+dx)/2.)
        elif dx < dy:
            ax.set_xlim(xmin + (dx-dy)/2., xmin + (dx+dy)/2.)

        # Set temperature bounds
    else:
        myxmin = xminn
        myxmax = xmaxx
        if xlim[0] != None:
            myxmin = xlim[0]
        if xlim[1] != None:
            myxmax = xlim[1]

        myymin = yminn
        myymax = ymaxx
        if ylim[0] != None:
            myymin = ylim[0]
        if ylim[1] != None:
            myymax = ylim[1]
            
        ax.set_xlim(myxmin,myxmax)
        ax.set_ylim(myymin,myymax)

        
    if colorbaradapt:
        Tmax = amax
        Tmin = amin
    else:
        if Tefflim[0] != None:
            Tmin = Tefflim[0]
        else:
            Tmin = aminn

        if Tefflim[1] != None:
            Tmax = Tefflim[1]
        else:
            Tmax = amaxx

    # Plot the data
    im = ax.imshow(data,cmap=my_cm,vmin=Tmin,vmax=Tmax,
                    aspect='equal',extent=(xmin,xmax,ymin,ymax),norm=norm,origin='lower')
    
    if type(pos) != 'matplotlib.transforms.Bbox':
        pos = ax.get_position()

    if colorbaradapt:
        cb.set_clim(Tmin,Tmax)
        cb.update_normal(cax)

    # Resize colorbar and put arrow on it for stuff out of bounds
    if not colorbaradapt:
        cbstateold = cbstate
        # cbstate = 0 no arrows
        # cbstate = 1 top arrow
        # cbstate = 2 bottom arrow
        # cbstate = 3 both arrows
        cbstate = 0 # No arrows
        if amax > Tmax: # top arrow
            triggered[0] = True
            cbstate = 1
        if amin < Tmin: # bottom arrow
            triggered[1] = True
            cbstate = 2
        if((amax > Tmax) and (amin < Tmin)): # both arrows
            triggered[2] = True
            cbstate = 3
        height = caxposorig.height
        x0 = caxposorig.x0 + 0.0007
        y0 = caxposorig.y0
        width = caxposorig.width
        if triggered[0]:
            height = height + 0.034
        if triggered[1]:
            y0 = y0 - 0.034
            height = height + 0.034

        if cbstate == 0:
            cb.extend = 'neither'
        elif cbstate == 1:
            cb.extend = 'max'
        elif cbstate == 2:
            cb.extend = 'min'
        elif cbstate == 3:
            cb.extend = 'both'

        if cbstateold != cbstate:
            cax.set_position([x0,y0,width,height])
            cb._inside = cb._slice_dict[cb.extend]
            cb.update_normal(cax)


                

    # Re-plot the axis labels

    ax.annotate("$x$",((pos.x0+pos.x1)/2.,0.05),xycoords='figure fraction', ha='center')
    ax.annotate("$y$",(0.0,(pos.y1+pos.y0)/2.),xycoords='figure fraction', ha='left',va='center',rotation='vertical')
    ax.annotate(cbar_label,(1.0,(pos.y1+pos.y0)/2.),xycoords='figure fraction', va='center',ha='right',rotation='vertical')

    
    # Update the plot title
    
    rotations = datafile.split("/")[-1][12:-4].split("_")
    zdeg = int(rotations[0])
    ydeg = int(rotations[1])
    xdeg = int(rotations[2])
    if not paperfriendly:
        pos = ax.get_position()
        ax.annotate("t = {:15.3f}".format(t),(pos.x0,0.95),xycoords='figure fraction',ha='left')
        ax.annotate("Rotations (x,y,z) = ({:06.2f}$^o$,{:06.2f}$^o$,{:06.2f}$^o$)".format(xdeg,ydeg,zdeg),(pos.x0,0.92),xycoords='figure fraction', ha='left')

    # Maintain plot axis dimensions
    if type(pos) == 'matplotlib.transforms.Bbox':
        ax.set_position([pos.x0,pos.y0,pos.width,pos.height])

    # Save the image

    savename = "teffs"+datafile.split("/")[-1][5:-4]
    if paperfriendly:
        savename = savename+".eps"
    else:
        savename = savename+".png"
    print("Saving", savename)
    plt.savefig(savename,facecolor=fig.get_facecolor())

    # Maintain the same axis size no matter what
    pos = ax.get_position()
    
    # Clear the axis
    ax.clear()

print("Finished.")
if len(files) > 1:
    print("")
    print("Use this command to make a movie:")
    print("convert -delay 10 -loop 0 teffs*.png teffs.gif")

