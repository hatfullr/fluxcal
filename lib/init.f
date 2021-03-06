      subroutine init
!     This is meant to initialize flux_cal with user's input
      include 'flux_cal.h'
      integer i
      logical inputexists,baseunitsexists,fileexists
      character*255 eosfilepath,inputfile, opacfile1, opacfile2
      character*255 opacityfile_temp
      character*255 info_particle_string
      logical use_rosseland, use_planck
      logical throw_error

      write(o,*)" ___________________________________________________  "
      write(o,*)"/         ___  _              ___        _          \ "
      write(o,*)"|        | __|| | _  _ __ __ / __| __ _ | |         | "
      write(o,*)"|        | _| | || || |\ \ /| (__ / _` || |         | "
      write(o,*)"|        |_|  |_| \_,_|/_\_\ \___|\__,_||_|         | "
      write(o,*)"\___________________________________________________/ "
      write(o,*)"|                                                   | "
      write(o,*)"|  FluxCal Copyright (C) 2020. Written by,          | "
      write(o,*)"|                                                   | "
      write(o,*)"|   Roger Hatfull      [University of Alberta]      | "
      write(o,*)"|   Natasha Ivanova    [University of Alberta]      | "
      write(o,*)"|   James C. Lombardi  [Allegheny College]          | "
      write(o,*)"|                                                   | "
      write(o,*)"|  This program comes with ABSOLUTELY NO WARRANTY.  | "
      write(o,*)"|  This is free software, and you are welcome to    | "
      write(o,*)"|  redistribute it under certain conditions. See    | "
      write(o,*)"|  LICENSE file for details, or visit               | "
      write(o,*)"|  https://www.gnu.org/licenses                     | "
      write(o,*)"\___________________________________________________/ "
      write(o,*)"|                   |                                 "
      write(o,*)"|   Version 1.0.0   |                                 "
      write(o,*)"\___________________/                                 "
      write(o,*)"                                                      "

      call flush(o)
      
      call sleep(1)
      
      write(o,*) "** Initializing **"



!-----------------------------------------------------------------------



      
      
      !### Files to read
      ! flux_cal starts at file number 'start' and ends at file number
      ! 'finish', taking steps of 'step'. All values must be integers.
      start=100
      finish=100
      step=0

      ! The directory flux_cal was installed to. This is used to find the
      ! included tabulated equation of state 'sph.eos' and other data
      ! tables required for flux_cal to run.
      flux_cal_dir='/home/hatfull/fluxcal_hatfull/'

      

      
      
      !### Runtime options
      ! Calculate the flux at each grid cell and record the data.
      get_fluxes=.true.

      ! Set this to true to create a data file containing all the
      ! particles whose surfaces are closest to the observer  along each
      ! driving grid cell. In other words, the particles that flux_cal
      ! first interacts with at each grid cell. The file will be called
      ! 'closest_particles_XXXX.dat' where 'XXXX' is a file number.
      get_closest_particles=.false.

      ! Finds all the particles the ray-tracer comes in contact with at
      ! posx, posy. Units must be in cgs.
      get_particles_at_pos=.false.
      posx=0.d0
      posy=0.d0

      ! Finds the data being used at each integration step at posx, posy
      get_integration_at_pos=.false.

      ! Finds the data being used at each integration step for all grid
      ! cells. Creates one output file for every grid cell, so use with
      ! extreme caution! You must also set get_fluxes to .true. to use this
      ! routine.
      get_integration_at_all_pos=.false.

      ! Track a list of particles through each input file and create
      ! output files containing their data. Specify the file name that
      ! contains the list of particle IDs in 'trackfile', and whether
      ! the resulting output file should be in binary format or in
      ! human readable format by setting 'binary_tracking_file' to
      ! '.true.' or '.false.', respectively.
      ! 'track_all' set to '.true.' tracks all particles without reading
      ! the trackfile.
      track_particles=.false.
      track_all=.false.
      trackfile='fluxcal.particles'
      binary_tracking_file=.false.

      ! Get the information of one particle printed to the standard
      ! output
      get_info_of_particle=.false.
      info_particle=0




      
      !### SPH simulation parameters
      ! nkernel determines what kind of kernel to use. You should set
      ! this to the kernel you used to calculate your data.
      ! 0 is cubic spline
      ! 1 is Wendland 3,3
      ! 2 is Wendland C4
      nkernel=2




      
      !### Integration
      ! Set rotation angles about the x, y, and z axes (degrees)
      anglexdeg=0.d0
      angleydeg=0.d0
      anglezdeg=0.d0
      
      ! yscale proportionality constant
      yscalconst=0.06d0

      ! When the detected flux from the simulation changes by less than
      ! fracaccuracy, grid refinement stops.
      fracaccuracy=0.01d0

      ! These control the resolution of the driving and integration grids.
      ! We vary Nx and Ny to achieve convergence in the flux. You cannot
      ! set max_Nx and max_Ny greater than 512 unless you change the
      ! value of "NXMAPMAX" and "NYMAPMAX" in
      ! optical_depth/optical_depth.h (not recommended).
      min_Nx=3                  ! (default=3)
      min_Ny=3                  ! (default=3)
      max_Nx=512                ! (default=512)
      max_Ny=512                ! (default=512)
      
      ! Choose which integrator to use.
      ! 0 = fourth order adaptive stepsize Runge-Kutta
      ! 1 = Simpson's rule (default)
      ! 2 = Pseudo Simpson's rule (deprecated)
      ! 0 is probably fastest but may encounter errors. 1 is probably
      ! slower, but you can directly control its accuracy using the
      ! variables below.
      integrator=1

      ! These are for the Runge--Kutta integrator (integrator=0).
      ! step1 is the step size the integrator should take when it is
      ! near tau of 1
      step1=1d0
      step2=1d30
      step3=1d4
      step4=1d15

      ! These control the precision of the Simpson's Rule integrator
      ! (integrator=1). simps_max_frac_dF is the maximum by which the
      ! integrated flux is allowed to change with each step, and the
      ! step size is continually refined until dF is less than or equal
      ! to simps_max_frac_dF. The same is true for optical depth, but
      ! with simps_max_dtau as the controller instead. If you set any of
      ! these to 0.d0, no maximum is imposed (not recommended for
      ! simps_max_frac_dF and simps_max_dtau). All the above is equally
      ! true for the minimum values, but it is okay to set the minimum
      ! to zero, as there are safeguards in place against stepsize
      ! underflows.
      !
      ! Each limit must be given in cgs units.
      simps_min_dz=0            ! (default=0)
      simps_max_dz=0            ! (default=0)

      simps_min_dtau=0          ! (default=0)
      simps_max_dtau=1.d0       ! (default=1.d0)
      
      simps_min_frac_dF=0       ! (default=0)
      simps_max_frac_dF=0       ! (default=0)
      
      ! This controls the threshold for the value that is checked
      ! against simps_min_frac_dF and simps_max_frac_dF. We calculate
      ! that value as |dF/(F+simps_F_cutoff)|. This effectively sets a
      ! "cutoff" value for the flux, such that contributions to the
      ! integrated flux are ignored when they are below simps_F_cutoff.
      ! We advise against setting this to some very small number like
      ! 1.d-30 because the stepsize will become arbitrarily small and
      ! the integrator=1 routine is not good at escaping extremely small
      ! stepsizes.
      simps_F_cutoff=1.d-4      ! (default=1.d-4)
            
      ! This is the fractional amount by which the stepsize is changed
      ! for the Simpson's Rule integrator (integrator=1). That is, when
      ! the stepsize must be decreased, we take dz=simps_alpha*dz, and
      ! when it must be increased, we take dz=(1.d0+simps_alpha)*dz. We
      ! choose 0.5d0 by default such that the stepsize is halved when
      ! it is too large and increased by 1.5 when too small.
      simps_alpha=0.5d0
      
      ! The maximum number of steps the integrator will take before
      ! failing. We recommend MAXSTP=10000 for integrator=0. For
      ! integrator=1, MAXSTP is used to make an initial guess on the
      ! first stepsize. The first stepsize depends on if simps_min_dz
      ! and/or simps_max_dz are non-zero (turned on). We first take the
      ! initial guess of dz=(z2-z1)/abs(MAXSTP), and then make sure
      ! that guess fits within simps_min_dz and simps_max_dz. If it is
      ! outside this domain, we set dz to either simps_min_dz or
      ! simps_max_dz depending on which is closer.
      !
      ! When MAXSTP<0, there is no upper limit on the number of
      ! integration steps for each line of sight. If you do want to
      ! impose a limit, set MAXSTP to a positive value and set
      ! simps_max_step_error=.true. (default). We suggest MAXSTP=-2
      ! (default) for integrator=1. Neither an error nor a warning will
      ! be thrown when MAXSTP<0, so simps_max_step_error and
      ! simps_max_step_warning can be set to anything.
      MAXSTP=-2
      simps_max_step_error=.true.
      simps_max_step_warning=.false.
            
      ! Optical depth at which to stop the integration. Only applies if
      ! get_fluxes or get_integration_at_pos are .true.. Set to 1.d30 to
      ! integrate through the entirety of the fluid on each line of
      ! sight. The value of taulimit_threshold determines how close to
      ! taulimit the integrator must be to stop integration (1.d-6
      ! default, and only works with integrator=1).
      taulimit=1.d1
      taulimit_threshold=1.d-6

      ! Optical depth by which the integrator will detect that a region
      ! is optically thick.
      tau_thick_integrator=1.d1

      ! The value of optical depth that a particle must exceed to be
      ! handled by the envelope fitting routine. 1.d1 by default.
      tau_thick_envfit=1.d1
      
      ! Decide to use the envelope fitting routine or not. Set this to
      ! .false. to use only the SPH temperatures.
      envfit=.true.
      




      
      !### Opacities
      ! When opacityfiles are not present or when a given (rho,T) lies
      ! outside the domain of the opacityfiles, FluxCal uses analytic
      ! approximations to calculate the boundfree, freefree, electron
      ! scatter, and negative H ion opacities, as done in Onno Pols notes
      ! on Stellar Evolution.
      ! astro.ru.nl/~onnop/education/stev_utrecht_notes/chapter5-6.pdf
      
      ! Decide to throw either an error, warning, both, or neither when
      ! (rho,T) is out of bounds for both the opacityfiles and any
      ! available analytic approximations in the code. When using the
      ! analytic approximations, a warning will always be thrown. By
      ! default, the opacity is set to 0.d0 when out of bounds of both the
      ! opacityfiles and analytic approximations.
      !
      ! Setting opacity_analytic_warning to .false. suppresses the warnings
      ! given when FluxCal uses analytic approximations to calculate the 
      ! opacity. This should only be set to .false. for debugging purposes.
      opacity_oob_error=.false.
      opacity_oob_warning=.true.
      opacity_analytic_warning=.true.

      ! Define the rho for which any rho <= this value forces opacity to
      ! equal 0.d0. Default is 1.d-19, which is the minimum density allowed
      ! by the cop subroutine.
      opacity_rho_cutoff=1.d-19

      ! The smoothing window in the temperature direction for opacities. A
      ! window of 0.d0 completely disables smoothing. Units in Kelvin. Do
      ! not make this larger than 1.d4, as this is the minimum distance in
      ! T required to calculate e- scattering and kramer opacities in the
      ! range 1.d4 < T < 1.d8.
      smoothing_window_T = 1000.d0
      smoothing_window_rho = 10.d0
      
      ! Define the opacity files you would like to be read in, as well as
      ! how you would like the code to transition between them. You must
      ! give the maximum and minimum logT and logR values (logTmins,
      ! logTmaxs, logRmins, and logRmaxs). Set these to -1.d30 to use up to
      ! the maximum and minimum in the file. You must also give the point
      ! at which you want opacity blending to begin (logT_blend1,
      ! logT_blend2, logR_blend1, logR_blend2). Blending will always occur
      ! in the region between mins and blend1 and between blend2 and maxs.
      ! No blending is done between blend1 and blend2. Set blend values to
      ! -1.d30 to turn off the blending region.
      !
      ! logR = logRho - 3*logT + 18
      !
      ! You can only use up to 9 opacity files.
      !
      ! The first two values in opacityfiles array are reserved for the
      ! included lowT rosseland and planck opacity files. These files
      ! require the included code "opacity.f" to be read.
      ! For the first two files, 1.0D-19 < rho < 1.0D-7 and
      ! 500 < T < 10000.
      ! logR = logrho - 3logT + 18
      ! for rhomin,Tmin, logR = -9.0969
      ! for rhomin,Tmax, logR = -13
      ! for rhomax,Tmin, logR = 2.90309
      ! for rhomax,Tmax, logR = -1
      ! Thus, the low temperature opacity routine can handle logR in the
      ! following ranges:
      ! -9.0969 < logR < 2.90309  for rhomin < rho < rhomax and T=Tmin
      !     -13 < logR < -1       for rhomin < rho < rhomax and T=Tmax
      ! -9.0969 < logR < -13      for Tmin < T < Tmax and rho=rhomin
      !      -1 < logR < 2.90309  for Tmin < T < Tmax and rho=rhomax
      ! Thus, for any constant temperature within the temperature range,
      ! varying density over the density range gives -13 < logR < 2.90309.
      ! However, for any constant density, varying the temperature over
      ! the temperature range reveals a gap in -9.0969 < logR < -1.
      !
      ! FluxCal always favors results from opacityfiles(1) and
      ! opacityfiles(2) over all other opacityfiles.

      opacityfiles(1)='kP_h2001.dat' ! Planck opacities (lowT)
      logTmins(1)=-1.d30        ! This does nothing. Limits are set in code
      logTmaxs(1)=3.041393      ! Not using Planck opacities above this T
      logRmins(1)=-1.d30        ! This does nothing. Limits are set in code
      logRmaxs(1)=-1.d30        ! This does nothing. Limits are set in code
      logT_blend1(1)=-1.d30     ! No blending at lower T bound
      logT_blend2(1)=2.995635d0 ! Start blending from Planck to Rosseland 
      logR_blend1(1)=-1.d30     ! No blending at lower R bound
      logR_blend2(1)=-1.d30     ! No blending at upper R bound
   
      opacityfiles(2)='kR_h2001.dat' ! Rosseland opacities (lowT)
      logTmins(2)=2.995635d0    ! Not using Rosseland opacities below this
      logTmaxs(2)=-1.d30       ! This does nothing. Limits are set in code
      logRmins(2)=-1.d30        ! This does nothing. Limits are set in code
      logRmaxs(2)=-1.d30        ! This does nothing. Limits are set in code
      logT_blend1(2)=3.041393d0 ! Stop blending from Planck to Rosseland 
      logT_blend2(2)=-1.d30     ! Start blending from Ross. to user files
      logR_blend1(2)=-1.d30     ! No blending at lower R bound
      logR_blend2(2)=-1.d30     ! No blending at upper R bound
      
      opacityfiles(3)='gs98_z0.02_x0.7.data'
      logTmins(3)=3.75d0
      logTmaxs(3)=8.7d0
      logRmins(3)=-8.d0
      logRmaxs(3)=1.d0
      logT_blend1(3)=-1.d30
      logT_blend2(3)=-1.d30
      logR_blend1(3)=-1.d30
      logR_blend2(3)=-1.d30

      opacityfiles(4)=''
      logTmins(4)=-1.d30
      logTmaxs(4)=-1.d30
      logRmins(4)=-1.d30
      logRmaxs(4)=-1.d30
      logT_blend1(4)=-1.d30
      logT_blend2(4)=-1.d30
      logR_blend1(4)=-1.d30
      logR_blend2(4)=-1.d30

      opacityfiles(5)=''
      logTmins(5)=-1.d30
      logTmaxs(5)=-1.d30
      logRmins(5)=-1.d30
      logRmaxs(5)=-1.d30
      logT_blend1(5)=-1.d30
      logT_blend2(5)=-1.d30
      logR_blend1(5)=-1.d30
      logR_blend2(5)=-1.d30

      opacityfiles(6)=''
      logTmins(6)=-1.d30
      logTmaxs(6)=-1.d30
      logRmins(6)=-1.d30
      logRmaxs(6)=-1.d30
      logT_blend1(6)=-1.d30
      logT_blend2(6)=-1.d30
      logR_blend1(6)=-1.d30
      logR_blend2(6)=-1.d30

      opacityfiles(7)=''
      logTmins(7)=-1.d30
      logTmaxs(7)=-1.d30
      logRmins(7)=-1.d30
      logRmaxs(7)=-1.d30
      logT_blend1(7)=-1.d30
      logT_blend2(7)=-1.d30
      logR_blend1(7)=-1.d30
      logR_blend2(7)=-1.d30

      opacityfiles(8)=''
      logTmins(8)=-1.d30
      logTmaxs(8)=-1.d30
      logRmins(8)=-1.d30
      logRmaxs(8)=-1.d30
      logT_blend1(8)=-1.d30
      logT_blend2(8)=-1.d30
      logR_blend1(8)=-1.d30
      logR_blend2(8)=-1.d30
      
      opacityfiles(9)=''
      logTmins(9)=-1.d30
      logTmaxs(9)=-1.d30
      logRmins(9)=-1.d30
      logRmaxs(9)=-1.d30
      logT_blend1(9)=-1.d30
      logT_blend2(9)=-1.d30
      logR_blend1(9)=-1.d30
      logR_blend2(9)=-1.d30
      

      

      ! This is for controlling the Planck and Rosseland opacities of
      ! the lowT opacity files 'kR_h2001.dat' and 'kP_h2001.dat'.
      !    dust_model controls the silicate type
      !       'nrm' - "normal" silicate dust model,    Fe/(Fe+Mg) is 0.3,
      !       'ips' - "iron-poor" silicate dust model, Fe/(Fe+Mg) is 0.0,
      !       'irs' - "iron-rich" silicate dust model, Fe/(Fe+Mg) is 0.4,
      !    dust_topology controls the topology of the grains
      !       'h' - homogeneous particles,
      !       'c' - composite particles,
      !       'p' - porous composite particles,
      !    dust_shape controls the shape of the particles
      !       's' - spherical dust,
      !       'a' - aggregate dust,
      !       '5' - 5-layered spherical dust
      ! For any dust_model, you can not use:
      !    dust_topology is 'h' AND dust_shape is '5'
      !    dust_topology is 'p' AND dust_shape is 'a'
      dust_model='nrm'
      dust_topology='h'
      dust_shape='s'

      
      
      !### Miscellaneous
      ! Solar metallicity
      metallicity=0.02d0
      


      !### Data units
      ! Give the conversion factors from your data's units to cgs.
      ! For example, set runit to 6.955d10 if your distances are measured
      ! in solar radii. This is intended for codes that use intricate
      ! unit systems, such as StarSmasher. If your data is in units that
      ! aren't that intricate, use the flux_cal.baseunits file instead.
      ! For example, StarSmasher uses vunit of sqrt(G*Msun/Rsun) cm/s and
      ! tunit of sqrt(Rsun^3/GMsun).
      !
      ! Be sure to give these in cgs. If you want to use other units,
      ! edit the flux_cal.baseunits file.
      runit=6.9598d10           ! distance [cm]
      munit=1.9892d33           ! mass [grams]
      tunit=1.59353d3           ! time [s]
      vunit=4.36754d7           ! velocity [cm/s]
      Eunit=3.79447d48          ! energy [ergs]
      rhounit=5.90049d0         ! density [grams/cm^3]
      muunit=1.d0               ! mean molecular weight [grams]
      gunit=2.74079d4           ! local g [cm/s^2]

      ! Give the conversion factor from your input data units to what
      ! you desire the output units to be in. Internal (cgs) values will
      ! be divided by these to produce your output.
      ! These values are affected by the flux_cal.baseunits file.
      runit_out=1.d0            ! distance
      munit_out=1.d0            ! mass
      tunit_out=1.d0            ! time
      vunit_out=1.d0            ! velocity
      Eunit_out=1.d0            ! energy
      rhounit_out=1.d0          ! density
      muunit_out=1.d0           ! mean molecular weight
      gunit_out=1.d0            ! gravitational acceleration
      tempunit_out=1.d0         ! temperature
      punit_out=1.d0            ! pressure
      Lunit_out=1.d0            ! luminosity
      kunit_out=1.d0            ! opacity
      sunit_out=1.d0		! specific entropy



      
      !### File names
      ! File name for the EOS, opacity files, filters file,
      ! and output file. For opacity files, set to '' to not
      ! use any file (turns off that part of the code).
      eosfile='sph.eos'
      filtersfile='filters.dat'
      outfile='flux_cal.output'


      
      !### Internal uses
      ! The debug flag is used internally to make debugging easier.
      debug=.false.

      
!-----------------------------------------------------------------------


      inquire(file='flux_cal.input',exist=inputexists)
      if(inputexists) then
         open(44,file='flux_cal.input')
         read(44,input)
         close(44)
      else
         write(e,*)"No 'flux_cal.input' file found."
         error stop "init.f"
      end if

c      if(tau_thick.ne.-1.d30) then
c         tau_thick_integrator = taulimit
c         tau_thick_envfit = tau_thick
c      end if
      
      
      ! Base units (cgs). These are 
      gram=1.d0
      sec=1.d0
      cm=1.d0
      kelvin=1.d0

      ! Get the base units
      inquire(file='flux_cal.baseunits',exist=baseunitsexists)
      if(baseunitsexists) then
         write(o,*)"'flux_cal.baseunits' file found. Using"//
     $        " user-specified physical base units."
         open(42,file='flux_cal.baseunits')
         read(42,baseunits)
         close(42)
      else
         inquire(file=trim(adjustl(flux_cal_dir))//
     $        "/defaults/flux_cal.baseunits",exist=baseunitsexists)
         if(baseunitsexists) then
            write(o,*)"Using default baseunits file '",
     $           trim(adjustl(flux_cal_dir)),"/defaults/",
     $           "flux_cal.baseunits'."
            open(42,file=trim(adjustl(flux_cal_dir))//
     $           "/defaults/flux_cal.baseunits")
            read(42,baseunits)
            close(42)
         else
            write(e,*)"WARNING: No 'flux_cal.baseunits' default file "//
     $           "found in either the current working directory or "//
     $           "'"//trim(adjustl(flux_cal_dir))//"/defaults'. "//
     $           "Assuming cgs units for all variables."
         end if
      end if
      
      ! Constants
      erg=gram*cm**2.d0/sec**2.d0
      boltz=1.380658d-16*erg/kelvin
      crad = 2.997924580d+10*cm/sec
      planck = 6.6260755d-27*gram*cm**2.d0/sec
      crad2 = crad**2.d0
      sigma = pi**2.d0*boltz*(boltz*2.0d0*pi/planck)**3.d0/60.0d0/crad2
      arad  = 4.d0*sigma/crad
      qconst = 1.5d0*boltz/arad
      gravconst = 6.67408d-8*cm**3.d0/gram/sec**2.d0
      coeff=planck*crad/(boltz*1d6)
      pc=3.08567758d18*cm       ! number of cm in a pc
      distance=10.d0*pc         ! distance to get absolute magnitude


      runit=runit*cm
      munit=munit*gram
      tunit=tunit*sec
      vunit=vunit*cm/sec
      Eunit=Eunit*gram*cm**2.d0/sec**2.d0
      rhounit=rhounit*gram/cm**3.d0
      muunit=muunit*gram
      gunit=gunit*cm/sec**2.d0
      
      runit_out=runit_out*cm
      munit_out=munit_out*gram
      tunit_out=tunit_out*sec
      vunit_out=vunit_out*cm/sec
      Eunit_out=Eunit_out*gram*cm**2.d0/sec**2.d0
      rhounit_out=rhounit_out*gram/cm**3.d0
      muunit_out=muunit_out*gram
      gunit_out=gunit_out*cm/sec**2.d0
      tempunit_out=tempunit_out*kelvin
      Lunit_out=Lunit_out*gram*cm**2.d0/sec**2.d0
      kunit_out=kunit_out*cm**2.d0/gram
      sunit_out=sunit_out*cm**2.d0/sec**2.d0/kelvin

      
 97   format(A7," = ",I6,"   ",A10," = ",g10.4)
 98   format(A10," = ",I10,"   ",A10," = ",I10)
 99   format(A13," = ",E10.4,"   ",A11," = ",E10.4)
 999  format(A13," = ",E10.4)
 100  format(A23," = ",E10.4)
 101  format(A23," = ",I10)
 102  format(A27," = ",L1)
 103  format(A22," = '",A,"'")
 104  format("   ",A," = ",L1)
 105  format("   ",A," = ",E10.4)
 106  format("   ",A," = '",A,"'")
 107  format("   ",A," = ",I8)
 108  format(A23," = ",L10)
 109  format(A13," = '",A,"'")
 110  format(A13,"(",I1,") = '",A,"'")
 111  format(A15,"(",I1,") = ",E10.4)
 112  format(A27," = ",L1)
 113  format(A27," = ",E10.4)

c     Write everything to the terminal
      write(o,*) ""
      write(o,99) "gram  ",gram  ,"arad       ",arad
      write(o,99) "sec   ",sec   ,"gravconst  ",gravconst
      write(o,99) "cm    ",cm    ,"pc         ",pc
      write(o,99) "kelvin",kelvin,"distance   ",distance
      write(o,99) "erg   ",erg   ,"boltz      ",boltz
      write(o,99) "crad  ",crad  ,"planck     ",planck
      write(o,99) "sigma ",sigma
      write(o,*) ""
      write(o,99) "runit       ",runit       ,"munit      ",munit
      write(o,99) "tunit       ",tunit       ,"vunit      ",vunit
      write(o,99) "Eunit       ",Eunit       ,"rhounit    ",
     $     rhounit
      write(o,99) "muunit      ",muunit      ,"gunit      ",gunit
      write(o,*) ""
      write(o,99) "runit_out   ",runit_out   ,"munit_out  ",
     $     munit_out
      write(o,99) "tunit_out   ",tunit_out   ,"vunit_out  ",
     $     vunit_out
      write(o,99) "Eunit_out   ",Eunit_out   ,"rhounit_out",
     $     rhounit_out
      write(o,99) "muunit_out  ",muunit_out  ,"gunit_out  ",
     $     gunit_out
      write(o,99) "tempunit_out",tempunit_out,"punit_out  ",
     $     punit_out
      write(o,99) "Lunit_out   ",Lunit_out   ,"kunit_out  ",
     $     kunit_out
      write(o,999)"sunit_out   ",sunit_out
      write(o,*) ""
      
      write(o,97) "start ",start ,"anglexdeg",anglexdeg
      write(o,97) "finish",finish,"angleydeg",angleydeg
      write(o,97) "step  ",step  ,"anglezdeg",anglezdeg

      write(o,109) "flux_cal_dir",trim(adjustl(flux_cal_dir))
      
      write(o,*) ""
      write(o,101) "nkernel               ",nkernel
      write(o,100) "yscalconst            ",yscalconst
      write(o,100) "fracaccuracy          ",fracaccuracy
      write(o,101) "MAXSTP                ",MAXSTP
      write(o,101) "min_Nx                ",min_Nx
      write(o,101) "min_Ny                ",min_Ny
      write(o,101) "max_Nx                ",max_Nx
      write(o,101) "max_Ny                ",max_Ny
      write(o,*) ""
      if ( integrator .eq. 0 ) then
         write(o,100) "step1                 ",step1
         write(o,100) "step2                 ",step2
         write(o,100) "step3                 ",step3
         write(o,100) "step4                 ",step4
      else if(integrator.eq.1) then
         write(o,100) "simps_min_dz          ",simps_min_dz
         write(o,100) "simps_max_dz          ",simps_max_dz
         write(o,100) "simps_min_dtau        ",simps_min_dtau
         write(o,100) "simps_min_dtau        ",simps_max_dtau
         write(o,100) "simps_min_frac_dF     ",simps_min_frac_dF
         write(o,100) "simps_max_frac_dF     ",simps_max_frac_dF
         write(o,100) "simps_F_cutoff        ",simps_f_cutoff
         write(o,108) "simps_max_step_error  ",simps_max_step_error
         write(o,108) "simps_max_step_warning",
     $        simps_max_step_warning
         write(o,100) "taulimit_threshold    ",taulimit_threshold
      end if
      write(o,*) ""
      write(o,100) "taulimit              ",taulimit
      write(o,100) "tau_thick_integrator  ",tau_thick_integrator
      write(o,100) "tau_thick_envfit      ",tau_thick_envfit
      write(o,108) "envfit                ",envfit
      write(o,*) ""
      write(o,100) "metallicity           ",metallicity
      
      write(o,*) ""
      write(o,102) "get_fluxes                ",get_fluxes
      write(o,102) "get_true_luminosity       ",get_true_luminosity
      write(o,102) "get_particles_at_pos      ",get_particles_at_pos
      write(o,102) "get_integration_at_pos    ",
     $     get_integration_at_pos
      write(o,102) "get_integration_at_all_pos",
     $     get_integration_at_all_pos
      
      if(get_particles_at_pos.or.get_integration_at_pos) then
         write(o,105) "posx",posx
         write(o,105) "posy",posy
      end if
      
      write(o,102) "get_closest_particles     ",
     $     get_closest_particles
      write(o,102) "track_particles           ",track_particles
      
      if(track_particles) then
         write(o,104) "track_all           ",track_all
         write(o,104) "binary_tracking_file",binary_tracking_file
         write(o,106) "trackfile           ",
     $        trim(adjustl(trackfile))
      end if
      
      write(o,102) "get_info_of_particle      ",get_info_of_particle
      if(get_info_of_particle) then
         write(o,107) "info_particle", info_particle
      end if

      write(o,*)
      write(o,112) "opacity_oob_error         ",opacity_oob_error
      write(o,112) "opacity_oob_warning       ",opacity_oob_warning
      write(o,112) "opacity_analytic_warning  ",
     $     opacity_analytic_warning
      write(o,113) "opacity_rho_cutoff        ",opacity_rho_cutoff
      write(o,113) "smoothing_window_T        ",smoothing_window_T
      

      do i=1,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).gt.0) then
            write(o,*)
            write(o,110) "opacityfiles",i,
     $           trim(adjustl(opacityfiles(i)))
            write(o,111) "      logTmins",i,logTmins(i)
            write(o,111) "      logTmaxs",i,logTmaxs(i)
            write(o,111) "      logRmins",i,logRmins(i)
            write(o,111) "      logRmaxs",i,logRmaxs(i)
            write(o,111) "   logT_blend1",i,logT_blend1(i)
            write(o,111) "   logT_blend2",i,logT_blend2(i)
            write(o,111) "   logR_blend1",i,logR_blend1(i)
            write(o,111) "   logR_blend2",i,logR_blend2(i)
         end if
      end do
      
      write(o,*) ""
      write(o,103) "eosfile              ",trim(adjustl(eosfile))
      write(o,103) "filtersfile          ",
     $     trim(adjustl(filtersfile))
      write(o,103) "outfile              ",trim(adjustl(outfile))
      write(o,*) ""

      call flush(o)

c     Catch some runtime errors
      if(tau_thick_integrator .gt. taulimit) then
         write(o,*) "tau_thick_integrator > taulimit. Integration "//
     $        "will always stop before reaching optically thick "//
     $        "material."
         error stop "init.f"
      end if

      if((tau_thick_integrator .gt. tau_thick_envfit).and.
     $     (envfit)) then
         write(o,*) "tau_thick_integrator > tau_thick_envfit. "//
     $        "Integrator will not return physical values."
         error stop "init.f"
      end if
      
      if((.not.envfit).and.(taulimit.lt.tau_thick).and.
     $     (tau_thick.ne.-1.d30)) then ! LEGACY
         write(o,*) "envfit is off and taulimit < tau_thick. "//
     $        "Integration will always stop before reaching "//
     $        "optically thick material."
         error stop "init.f"
      end if
      
      if ((start.gt.finish).and.(step.gt.0)) then
         write(o,*) "Must have start < finish when step > 0."
         write(o,*) "start=",start,"finish=",finish
         error stop "init.f"
      end if
      if ((start.lt.finish).and.(step.lt.0)) then
         write(o,*) "Must have start > finish when step < 0."
         write(o,*) "start=",start,"finish=",finish
         error stop "init.f"
      end if

      if(start.lt.0) then
         write(o,*) "Starting file index must be > 0."
         write(o,*) "start=",start
         error stop "init.f"
      end if

      if((step.eq.0).or.(start.eq.finish)) then
         finish=start
         step=1
      end if

      if (MAXSTP.eq.0) then
         write(o,*) "Cannot have MAXSTP = 0"
         error stop "init.f"
      end if
      if (abs(MAXSTP).eq.1) then
         write(o,*) "Cannot have |MAXSTP| = 1"
         error stop "init.f"
      end if
      if ( integrator.eq.1 ) then
         if ( (simps_min_dz.ne.0 .and. simps_max_dz.ne.0) .and.
     $        (dabs(simps_min_dz).gt.dabs(simps_max_dz)) ) then
            write(o,*) "|simps_min_dz| > |simps_max_dz|"
            write(o,*) "simps_min_dz = ",simps_min_dz
            write(o,*) "simps_max_dz = ",simps_max_dz
            error stop "init.f"
         end if
         if ( (simps_min_dtau.ne.0 .and. simps_max_dtau.ne.0) .and.
     $        (dabs(simps_min_dtau).gt.dabs(simps_max_dtau)) ) then
            write(o,*) "|simps_min_dtau| > |simps_max_dtau|"
            write(o,*) "simps_min_dtau = ",simps_min_dtau
            write(o,*) "simps_max_dtau = ",simps_max_dtau
            error stop "init.f"
         end if
         if ( (simps_min_frac_dF.ne.0.and.simps_max_frac_dF.ne.0).and.
     $        (dabs(simps_min_frac_dF).gt.dabs(simps_max_frac_dF)))then
            write(o,*) "|simps_min_frac_dF| > |simps_max_frac_dF|"
            write(o,*) "simps_min_frac_dF = ",simps_min_frac_dF
            write(o,*) "simps_max_frac_dF = ",simps_max_frac_dF
            error stop "init.f"
         end if
         if ( simps_min_frac_dF.ne.0 .and. simps_max_frac_dF.ne.0 .and.
     $        simps_F_cutoff.eq.0 ) then
            write(o,*) "Must give a non-zero value for simps_F_cutoff"//
     $           " when simps_min_frac_dF != 0 and simps_max_frac_dF "//
     $           "!= 0"
            error stop "init.f"
         end if
         if ( taulimit.ne.1.d30 .and. taulimit_threshold.le.0 ) then
            write(o,*) "When taulimit=1.d30, you must give a positive"//
     $           ", non-zero value for taulimit_threshold"
            error stop "init.f"
         end if
      end if

      if (min_Nx.lt.3) then
         write(o,*) "min_Nx = ",min_Nx
         write(o,*) "Must have min_Nx >= 3"
         error stop "init.f"
      end if
      if (min_Ny.lt.3) then
         write(o,*) "min_Ny = ",min_Ny
         write(o,*) "Must have min_Ny >= 3"
         error stop "init.f"
      end if
      if (max_Nx.gt.NXMAPMAX) then
         write(o,*) "max_Nx = ",max_Nx
         write(o,*) "Must have max_Nx <= ",NXMAPMAX
         error stop "init.f"
      end if
      if (max_Ny.gt.NYMAPMAX) then
         write(o,*) "max_Ny = ",max_Ny
         write(o,*) "Must have max_Ny <= ",NYMAPMAX
         error stop "init.f"
      end if

      call flush(o)
      
 200  format(A,"/defaults/",A)
      
      write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $     trim(adjustl(eosfile))
      inquire(file=trim(adjustl(inputfile)),exist=fileexists)
      if(.not. fileexists) then
         write(o,*) "Could not find EOS file '"//
     $        trim(adjustl(eosfile))//"' in flux_cal_dir/defaults."
         error stop "init.f"
      end if
      call readineostable(inputfile)

      call tabulinit


      if ( len(trim(adjustl(opacityfiles(1)))).eq.0 .and.
     $     len(trim(adjustl(opacityfiles(2)))).eq.0 ) then
         write(o,*) "opacityfiles(1) = '"//
     $        trim(adjustl(opacityfiles(1)))//"'"
         write(o,*) "opacityfiles(2) = '"//
     $        trim(adjustl(opacityfiles(2)))//"'"
         write(o,*) "There must be valid arguments for one"//
     $        "of or both opacityfiles(1) and opacityfiles(2)"//
     $        " (Planck and Rosseland opacities, respectively)."
         error stop "init.f"
      end if
      
      use_planck = .false.
      use_rosseland = .false.
      if (len(trim(adjustl(opacityfiles(1)))).ne.0) use_planck=.true.
      if (len(trim(adjustl(opacityfiles(2)))).ne.0) use_rosseland=.true.

      if((use_planck.eqv..false.).and.(use_rosseland.eqv..false.)) then
         write(o,*) "No Planck or Rosseland file found. There should"//
     $        " at least be 'kP_h2001.dat' and 'kR_h2001.dat' in the"//
     $        " defaults directory of your flux_cal installation." //
     $        " Please check that flux_cal installed properly."
         error stop "init.f"
      end if
      
      throw_error = .false.
      do i=1,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).gt.0) then
            ! Check current working directory for the file
            inquire(file=trim(adjustl(opacityfiles(i))),
     $           exist=fileexists)
            if (fileexists.eqv..false.) then
               write(o,*) "WARNING: Could not find '"//
     $              trim(adjustl(opacityfiles(i)))//"'. "//
     $              "Checking the defaults directory in '"//
     $              trim(adjustl(flux_cal_dir))//"/defaults'."
               write(opacityfile_temp,200) trim(adjustl(flux_cal_dir)),
     $              trim(adjustl(opacityfiles(i)))
               opacityfiles(i) = trim(adjustl(opacityfile_temp))
               inquire(file=trim(adjustl(opacityfiles(i))),
     $              exist=fileexists)
               if(fileexists) then
                  write(o,*) "Found '"//trim(adjustl(opacityfiles(i)))//
     $                 "'."
               else
                  write(o,*)"Could not find '"//
     $                 trim(adjustl(opacityfiles(i)))//"' in "//
     $                 "flux_cal_dir/defaults. make sure flux_cal_dir"//
     $                 " is correct, and if flux_cal was properly "//
     $                 "installed."
                  throw_error = .true.
               end if
            end if

            ! Check to make sure the blending regions are valid
            if ( logT_blend1(i).ne.-1.d30
     $           .and.
     $           ((logTmins(i).ne.-1.d30 .and.
     $           logT_blend1(i).lt.logTmins(i))
     $           .or.
     $           (logTmaxs(i).ne.-1.d30 .and.
     $           logT_blend1(i).gt.logTmaxs(i)))) then
               write(o,*) "Your logT_blend1(i) is invalid for"//
     $              " opacityfiles(i) = '"//
     $              trim(adjustl(opacityfiles(i)))//"'."
               throw_error = .true.
            end if
            if ( logT_blend2(i).ne.-1.d30
     $           .and.
     $           ((logTmaxs(i).ne.-1.d30 .and.
     $           logT_blend2(i).gt.logTmaxs(i))
     $           .or.
     $           (logTmins(i).ne.-1.d30 .and.
     $           logT_blend2(i).lt.logTmins(i)))) then
               write(o,*) "Your logT_blend2(i) is invalid for"//
     $              " opacityfiles(i) = '"//
     $              trim(adjustl(opacityfiles(i)))//"'."
               throw_error = .true.
            end if
            if ( logR_blend1(i).ne.-1.d30
     $           .and.
     $           ((logRmins(i).ne.-1.d30 .and.
     $           logR_blend1(i).lt.logRmins(i))
     $           .or.
     $           (logRmaxs(i).ne.-1.d30 .and.
     $           logR_blend1(i).gt.logRmaxs(i)))) then
               write(o,*) "Your logR_blend1(i) is invalid for"//
     $              " opacityfiles(i) = '"//
     $              trim(adjustl(opacityfiles(i)))//"'."
               throw_error = .true.
            end if
            if ( logR_blend2(i).ne.-1.d30
     $           .and.
     $           ((logRmaxs(i).ne.-1.d30 .and.
     $           logR_blend2(i).gt.logRmaxs(i))
     $           .or.
     $           (logRmins(i).ne.-1.d30 .and.
     $           logR_blend2(i).lt.logRmins(i)))) then
               write(o,*) "Your logR_blend2(i) is invalid for"//
     $              " opacityfiles(i) = '"//
     $              trim(adjustl(opacityfiles(i)))//"'."
               throw_error = .true.
            end if
         end if
      end do

      if(throw_error) error stop "init.f"
      

      if (use_planck .or. use_rosseland) then
         if(use_planck) then
            write(o,*)"Using Planck opacities from '"//
     $           trim(adjustl(opacityfiles(1)))//"'"
         else
            write(o,*) "WARNING (init.f): Not using any Planck"//
     $           " opacities."
         end if
         if(use_rosseland) then
            write(o,*)"Using Rosseland opacities from '"//
     $           trim(adjustl(opacityfiles(2)))//"'"
         else
            write(o,*) "WARNING (init.f): Not using any "//
     $           "Rosseland opacities."
         end if
         call init_lowT_opac(use_planck,use_rosseland,
     $        opacityfiles(1),opacityfiles(2),0) ! in opacity.f
      end if

      call init_highT_opac(0) ! in opacity_highT.f




      
      if(envfit) then
         ! Initialize find_teff dependencies
         write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $        "lowT_fa05_gs98_z0.02_x0.7.data"
         inquire(file=trim(adjustl(inputfile)), exist=fileexists)
         if(.not. fileexists) then
            write(o,*) "Could not find '"//trim(adjustl(inputfile))//
     $           "' in"
            write(o,*)"flux_cal_dir/defaults. make sure flux_cal_dir "//
     $           "is correct, and if flux_cal was properly installed."
            error stop "init.f"
         end if
         
         call ini_opacity_photospehre_lowT(inputfile)
         call create_density_array
      
         write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $        "table_nabla.dat"
         inquire(file=trim(adjustl(inputfile)),exist=fileexists)
         if(.not. fileexists) then
            write(o,*) "Could not find '"//trim(adjustl(inputfile))//
     $           "' in"
            write(o,*)"flux_cal_dir/defaults. make sure flux_cal_dir "//
     $           "is correct, and if flux_cal was properly installed."
            error stop "init.f"
         end if
      
         call ini_slops(inputfile)
      end if
         
c     Check to see if the filtersfile exists and read it if it does
      if(len(trim(adjustl(filtersfile))).gt.0) then
         write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $        trim(adjustl(filtersfile))
         inquire(file=trim(adjustl(inputfile)),exist=fileexists)
         if(.not.fileexists) then
            write(o,*) "Could not find '"//trim(adjustl(inputfile))//
     $           "'"//
     $           " in flux_cal_dir/defaults. make sure flux_cal_dir"//
     $           " is correct, and if flux_cal was properly installed."
            error stop "init.f"
         end if
      
         write(o,*) "Reading filtersfile"
         open(16,file=trim(adjustl(inputfile)),status='old')
         read (16,*) numfilters
      
         do i=1,numfilters
            read(16,*) filtername(i), wavelength(i),absoluteflux(i)
            yscalfactor(i)=yscalconst*wavelength(i)
         enddo
         close(16)
      end if

      
      if(track_particles) then
         if(.not.track_all) then ! Use trackfile
            inquire(file=trim(adjustl(trackfile)),exist=fileexists)
            if(.not.fileexists) then
               write(o,*) "Could not find trackfile '"//
     $              trim(adjustl(trackfile))//"'."
               error stop "init.f"
            end if
         end if                 ! If track_all, don't use trackfile
      end if

c      if(taulimit.lt.tau_thick) then
c         write(o,*) "WARNING: taulimit < tau_thick. This will result in"
c         write(o,*) "integrations that stop part way through optically"
c         write(o,*) "thin material, and as such, results should not be"
c         write(o,*) "expected to be physical."
c      end if

c     Initialize some variables
      !posx = posx*runit
      !posy = posy*runit
      printIntegrationSteps = .false.

      if(get_info_of_particle) then
         write(info_particle_string,*) info_particle
         pinfo_file = trim("particle_" //
     $        trim(adjustl(info_particle_string))) // ".dat"
         call makeOutputFile(pinfo_file)
      end if
      
      write(o,*) "** Complete **"
      write(o,*) ""

      call flush(o)
      
      return
      end subroutine
