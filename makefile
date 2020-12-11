# Makefile for flux_cal
# Use gfortran because it's more standard
FC = $(shell which gfortran)
ifeq (, $(shell which gfortran))
$(error ERROR: gfortran is not installed. Try 'sudo apt install gfortran')
endif

FFLAGS = -O4 -ffixed-line-length-132 -mcmodel=large

LIB = /usr/lib
BUILD_DIR = build

executable = flux_cal

mathlib = mathlib
llib = lib
opdeplib = optical_depth
findTefflib = find_teff

# These need to be built first
math_modules = $(mathlib)/bilinear_interpolation.o \
	       $(mathlib)/toms760.o                \

lib_modules = $(llib)/iso_fortran_env.o

# Library objects
lib_obj = $(llib)/getTemperature.o      \
          $(llib)/derivs.o              \
          $(llib)/kernels.o             \
          $(llib)/init.o                \
          $(llib)/getOpacity.o          \
          $(llib)/read_fluxcal.o        \
          $(llib)/getLocalQuantities.o  \
          $(llib)/getLocalvz.o          \
          $(llib)/useeostable.o         \
          $(llib)/readineostable.o      \
          $(llib)/trackParticles.o      \
          $(llib)/quicksort.o           \
          $(llib)/makeOutputFile.o      \
          $(llib)/output.o              \
          $(llib)/getLocalAngle.o       \
	  $(llib)/opacity.o             \
	  $(llib)/opacityTables.o       \

# Math objects
math_obj = $(mathlib)/odeint.o           \
           $(mathlib)/rkck.o             \
           $(mathlib)/rkqs.o             \
	   $(mathlib)/interp.o           \
	   $(mathlib)/toms526.o          \
	   $(mathlib)/tseval_f90.o       \
	   $(mathlib)/simps.o            \

# Optical depth objects
opdep_obj = $(opdeplib)/createGrid.o               \
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
modules = $(lib_modules) $(math_modules)

# All objects
process_obj = $(modules) $(lib_obj) $(math_obj) $(opdep_obj) $(findTeff_obj)

$(executable): $(executable).o $(process_obj)
	$(FC) $(FFLAGS) -o $(executable) $(executable).o $(process_obj) -L$(LIB) -I$(BUILD_DIR) 
	@\mv $(executable).o $(process_obj) $(BUILD_DIR)/.
	@\mv $(shell find . -name "*.mod") $(BUILD_DIR)/.

# Delete all made stuff
clean:
	@\rm -f $(BUILD_DIR)/* $(executable) $(process_obj)

