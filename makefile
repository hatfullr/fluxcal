# Makefile for flux_cal
# Use gfortran because it's more standard
FC = $(shell which gfortran)
FFLAGS = -O4 -ffixed-line-length-132 -mcmodel=large

LIB = -L/usr/lib

executable = flux_cal

mathlib = mathlib
llib = lib
opdeplib = optical_depth
findTefflib = find_teff

# Library objects
lib_obj = $(llib)/usekappatable.o      \
          $(llib)/getTemperature.o     \
          $(llib)/derivs2.o            \
          $(llib)/kernels.o            \
          $(llib)/init.o               \
          $(llib)/getOpacity.o         \
          $(llib)/read_fluxcal.o       \
          $(llib)/getLocalQuantities.o \
          $(llib)/getLocalvz.o         \
          $(llib)/useeostable.o        \
          $(llib)/readineostable.o     \
          $(llib)/trackParticles.o     \
          $(llib)/getOpacitySub.o      \
	  $(llib)/opacities_cold.o     \
          $(llib)/quicksort.o          \
          $(llib)/makeOutputFile.o     \
          $(llib)/output.o             \

# Math objects
math_obj = $(mathlib)/odeint.o \
           $(mathlib)/rkck.o   \
           $(mathlib)/rkqs.o   \

# Optical depth objects
opdep_obj = $(opdeplib)/optical_depth.o        \
            $(opdeplib)/createGrid.o           \
            $(opdeplib)/getFlux.o              \
            $(opdeplib)/prepareIntegration.o   \
            $(opdeplib)/integrateTau.o         \
            $(opdeplib)/peakWavelengths.o      \
            $(opdeplib)/setViewingAngle.o      \
            $(opdeplib)/useDimenFile.o         \
            $(opdeplib)/writeDimenFile.o       \
            $(opdeplib)/particlesAtPos.o       \
            $(opdeplib)/integrationAtPos.o     \
            $(opdeplib)/getParticleInfo.o      \
            $(opdeplib)/getTpractical.o        \
            $(opdeplib)/getClosestParticles.o  \

# find_teff objects
findTeff_obj = $(findTefflib)/get_slop.o                \
               $(findTefflib)/get_teff.o                \
               $(findTefflib)/get_kappa.o               \
               $(findTefflib)/ini_opacity_photosphere.o \
               $(findTefflib)/rho_out_pt.o              \
               $(findTefflib)/ini_slops.o               \

# All objects
process_obj = $(executable).o $(lib_obj) $(math_obj) $(opdep_obj) $(findTeff_obj)

$(executable): $(process_obj)
	$(FC) -o $(executable) $(process_obj) $(LIB)

# Delete all made stuff
clean:
	@\rm -rf $(process_obj) $(executable)
