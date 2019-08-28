# Makefile for flux_cal
# Use gfortran because it's more standard
FC = $(shell which gfortran)
ifeq (, $(shell which gfortran))
$(error ERROR: gfortran is not installed. Try 'sudo apt install gfortran')
endif

FFLAGS = -O4 -ffixed-line-length-132 -mcmodel=large

LIB = -L/usr/lib
BUILD_DIR = ./build

executable = flux_cal

mathlib = mathlib
llib = lib
opdeplib = optical_depth
findTefflib = find_teff

# These need to be built first
math_modules = $(mathlib)/bilinear_interpolation.o

# Library objects
lib_obj = $(llib)/getTemperature.o      \
          $(llib)/derivs2.o             \
          $(llib)/kernels.o             \
          $(llib)/init.o                \
          $(llib)/getOpacity.o          \
          $(llib)/read_fluxcal.o        \
          $(llib)/getLocalQuantities.o  \
          $(llib)/getLocalvz.o          \
          $(llib)/useeostable.o         \
          $(llib)/readineostable.o      \
          $(llib)/trackParticles.o      \
          $(llib)/getOpacitySub.o       \
          $(llib)/quicksort.o           \
          $(llib)/makeOutputFile.o      \
          $(llib)/output.o              \
          $(llib)/getLocalAngle.o       \
	  $(llib)/opacity.o             \

# Math objects
math_obj = $(mathlib)/odeint.o               \
           $(mathlib)/rkck.o                 \
           $(mathlib)/rkqs.o                 \
           $(mathlib)/fourPointArea.o        \

# Optical depth objects
opdep_obj = $(opdeplib)/optical_depth.o            \
            $(opdeplib)/createGrid.o               \
            $(opdeplib)/getFlux.o                  \
            $(opdeplib)/prepareIntegration.o       \
            $(opdeplib)/integrateTau.o             \
            $(opdeplib)/peakWavelengths.o          \
            $(opdeplib)/setViewingAngle.o          \
            $(opdeplib)/useDimenFile.o             \
            $(opdeplib)/writeDimenFile.o           \
            $(opdeplib)/particlesAtPos.o           \
            $(opdeplib)/integrationAtPos.o         \
            $(opdeplib)/getTpractical.o            \
            $(opdeplib)/getClosestParticles.o      \
            $(opdeplib)/integrationAtAllPos.o      \
            $(opdeplib)/writeTempsFile.o           \
            $(opdeplib)/init_grid.o                \
	    $(opdeplib)/integrateSpherical.o	   \

# find_teff objects
findTeff_obj = $(findTefflib)/get_slop.o                \
               $(findTefflib)/get_teff.o                \
               $(findTefflib)/get_kappa.o               \
               $(findTefflib)/ini_opacity_photosphere.o \
               $(findTefflib)/rho_out_pt.o              \
               $(findTefflib)/ini_slops.o               \

# All modules
modules = $(math_modules)

# All objects
process_obj = $(modules) $(lib_obj) $(math_obj) $(opdep_obj) $(findTeff_obj)

$(executable): $(executable).o $(process_obj)
	$(FC) -I$(BUILD_DIR) -o $(executable) $(executable).o $(process_obj) $(LIB)
	@\mv $(process_obj) $(BUILD_DIR)/.
	@\mv *.mod $(BUILD_DIR)/.

# Delete all made stuff
clean:
	@\rm -f $(BUILD_DIR)/* $(executable)
#	@\rm -rf $(executable).o $(process_obj) $(modules) $(executable)
