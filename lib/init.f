      subroutine init
!     This is meant to initialize flux_cal with user's input
      include 'flux_cal.h'
      integer i
      logical inputexists,baseunitsexists,fileexists
      character*255 eosfilepath,inputfile, opacfile1, opacfile2
c      logical baseunitsdefaultexists
c      character*255 baseunits_default


      write(*,*)" ___________________________________________________  "
      write(*,*)"/         ___  _              ___        _          \ "
      write(*,*)"|        | __|| | _  _ __ __ / __| __ _ | |         | "
      write(*,*)"|        | _| | || || |\ \ /| (__ / _` || |         | "
      write(*,*)"|        |_|  |_| \_,_|/_\_\ \___|\__,_||_|         | "
      write(*,*)"\___________________________________________________/ "
      write(*,*)"|                                                   | "
      write(*,*)"|  FluxCal Copyright (C) 2018. Written by,          | "
      write(*,*)"|                                                   | "
      write(*,*)"|   Roger Hatfull      [University of Alberta]      | "
      write(*,*)"|   Natasha Ivanova    [University of Alberta]      | "
      write(*,*)"|   James C. Lombardi  [Allegheny College]          | "
c      write(*,*)"|   Jose Nandez        [Uni. of West. Ontario]      | "
      write(*,*)"|                                                   | "
      write(*,*)"|  This program comes with ABSOLUTELY NO WARRANTY.  | "
      write(*,*)"|  This is free software, and you are welcome to    | "
      write(*,*)"|  redistribute it under certain conditions. See    | "
      write(*,*)"|  LICENSE file for details, or visit               | "
      write(*,*)"|  https://www.gnu.org/licenses                     | "
      write(*,*)"\___________________________________________________/ "
      write(*,*)"|                   |                                 "
      write(*,*)"|   Version 0.0.1   |                                 "
      write(*,*)"\___________________/                                 "
      write(*,*)"                                                      "

      call sleep(1)
      
      write(*,*) "** Initializing **"
      
!-----------------------------------------------------------------------



      
      
      !### Files to read
      ! flux_cal starts at file number 'start' and ends at file number
      ! 'finish', taking steps of 'step'. All values must be integers.
      start=0
      finish=0
      step=0

      ! The directory flux_cal was installed to. This is used to find the
      ! included tabulated equation of state 'sph.eos' and other data
      ! tables required for flux_cal to run.
      flux_cal_dir = ''

      

      
      
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
      ! posx, posy
      get_particles_at_pos=.false.
      posx=0.d0
      posy=0.d0

      ! Finds the data being used at each integration step at posx, posy
      get_integration_at_pos=.false.

      ! Track a list of particles through each input file and create
      ! output files containing their data. Specify the file name that
      ! contains the list of particle IDs in 'trackfile', and whether
      ! the resulting output file should be in binary format or in
      ! human readable format by setting 'binary_tracking_file' to
      ! '.true.' or '.false.', respectively.
      track_particles=.false.
      trackfile='fluxcal.particles'
      binary_tracking_file=.false.

      ! Get the information of one particle printed to the standard
      ! output
      get_info_of_particle=.false.
      info_particle=0




      
      !### Simulation parameters
      ! nkernel determines what kind of kernel to use. You should set
      ! this to the kernel you used to calculate your data.
      ! 0 = cubic spline
      ! 1 = Wendland 3,3
      ! 2 = Wendland C4
      nkernel=0




      
      !### Integration
      ! yscale proportionality constant
      yscalconst=0.06d0
      
      fracaccuracy=0.01d0

      ! step1 is the step size the integrator should take when it is
      ! near tau=1
      step1=1d0
      step2=1d30
      step3=1d4
      step4=1d15
      
      ! Optical depth at which to stop the integration. Only applies to
      ! get_fluxes and get_integration_at_pos.
      taulimit=1.d0

      ! Optical depth value that signifies optically thick material. If
      ! the integrator finds a region where tau is equal to this value,
      ! it will tell flux_cal to use the envelope fitting routine.
      ! Otherwise, flux_cal will use the SPH temperature. Set to 10.d0 by
      ! default, because the envelope fitting routine works mostly for
      ! tau > 10.
      tau_thick=1.d1
      
      ! Decide to use the envelope fitting routine or not. Set this to
      ! .false. to use only the SPH temperatures.
      envfit=.true.




      
      !### Data units
      ! Give the conversion factors from your data's units to cgs.
      ! For example, set runit=6.955d10 if your distances are measured
      ! in solar radii. This is intended for codes that use intricate
      ! unit systems, such as StarSmasher. If your data is in units that
      ! aren't that intricate, use the flux_cal.baseunits file instead.
      ! For example, StarSmasher uses vunit=sqrt(G*Msun/Rsun) cm/s and
      ! tunit=sqrt(Rsun^3/GMsun).
      !
      ! Be sure to give these in cgs. If you want to use other units,
      ! edit the flux_cal.baseunits file.
      runit=1.d0                ! distance [cm]
      munit=1.d0                ! mass [grams]
      tunit=1.d0                ! time [s]
      vunit=1.d0                ! velocity [cm/s]
      Eunit=1.d0                ! energy [ergs]
      rhounit=1.d0              ! density [grams/cm^3]
      muunit=1.d0               ! mean molecular weight [grams]
      gunit=1.d0                ! local g [cm/s^2]

      ! Give the conversion factor from your input data units to what
      ! you desire the output units to be in. Internal (cgs) values will
      ! be divided by these to produce your output.
      ! These values are affected by the flux_cal.baseunits file.
      runit_out=1.d0
      munit_out=1.d0
      tunit_out=1.d0
      vunit_out=1.d0
      Eunit_out=1.d0
      rhounit_out=1.d0
      muunit_out=1.d0
      gunit_out=1.d0
      tempunit_out=1.d0
      punit_out=1.d0
      Lunit_out=1.d0
      
      
      
      !### Physics
      ! If true, use Rosseland opacities only. If false, will use
      ! Rosseland and Planck opacities with a smoothing function.
      rossonly=.true.
      
      ! Solar metallicity
      metallicity=0.02d0
      
      ! Rotation angle about x axis (degrees)
      anglexdeg=0.d0

      ! Rotation angle about y axis (degrees)
      angleydeg=0.d0

      ! Rotation angle about z axis (degrees)
      anglezdeg=0.d0

      ! Dust formation radius (set to 1.d30 for no dust)
      Rform=1.d30




      
      !### File names
      ! File name for the opacity files, filters file, and output file.
      ! Set opacitydustfile='' for no dust.
      eosfile='sph.eos'
      opacityfile='sph.opacities_ferguson_yes_grains_and_molecules'
      opacitydustfile='sph.opacities_yes_grains_and_molecules_ncooling2'
      filtersfile='filters.dat'
      outfile='flux_cal.output'




      
!-----------------------------------------------------------------------


      inquire(file="flux_cal.baseunits",exist=baseunitsexists)
      
      ! Base units (cgs). These are 
      gram=1.d0
      sec=1.d0
      cm=1.d0
      kelvin=1.d0

      ! Get the base units
      inquire(file='flux_cal.baseunits',exist=baseunitsexists)
      if(baseunitsexists) then
         write(*,*) "'flux_cal.baseunits' file found. Using"//
     $        " user-specified"
         write(*,*) "physical base units."
         open(42,file='flux_cal.baseunits')
         read(42,baseunits)
         close(42)
      else
         inquire(file=trim(adjustl(flux_cal_dir))//
     $        "/defaults/flux_cal.baseunits",exist=baseunitsexists)
         if(baseunitsexists) then
            write(*,*) "Using default baseunits file from "//
     $           "flux_cal_dir/defaults."
            open(42,file=trim(adjustl(flux_cal_dir))//
     $           "/defaults/flux_cal.baseunits")
            read(42,baseunits)
            close(42)
         else
            write(*,*)"WARNING: No 'flux_cal.baseunits' default file "//
     $           "found"
            write(*,*)"in either the current working directory or the"//
     $           " flux_cal_dir/defaults directory. Defaulting to cgs"//
     $           " units."
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


      inquire(file='flux_cal.input',exist=inputexists)
      if(inputexists) then
         open(44,file='flux_cal.input')
         read(44,input)
         close(44)
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
      else
         write(*,*) "No 'flux_cal.input' file found."
         error stop "init.f line 241"
      end if
 97   format(A7," = ",I6,"   ",A10," = ",g10.4)
 98   format(A10," = ",I10,"   ",A10," = ",I10)
 99   format(A10," = ",E10.4,"   ",A10," = ",E10.4)
 999  format(A10," = ",E10.4)
 100  format(A13," = ",E10.4)
 101  format(A13," = ",I10)
 102  format(A23," = ",L1)
 103  format(A16," = '",A,"'")
 104  format("   ",A," = ",L1)
 105  format("   ",A," = ",E10.4)
 106  format("   ",A," = '",A,"'")
 107  format("   ",A," = ",I8)

c     Write everything to the terminal
      write(*,*) ""
      write(*,99) "gram  ",gram  ,"arad     ",arad
      write(*,99) "sec   ",sec   ,"gravconst",gravconst
      write(*,99) "cm    ",cm    ,"pc       ",pc
      write(*,99) "kelvin",kelvin,"distance ",distance
      write(*,99) "erg   ",erg   ,"boltz    ",boltz
      write(*,99) "crad  ",crad  ,"planck   ",planck
      write(*,99) "sigma ",sigma
      write(*,*) ""
      write(*,99) "runit       ",runit       , "munit      ",munit
      write(*,99) "tunit       ",tunit       , "vunit      ",vunit
      write(*,99) "Eunit       ",Eunit       , "rhounit    ",rhounit
      write(*,99) "muunit      ",muunit      , "gunit      ",gunit
      write(*,*) ""
      write(*,99) "runit_out   ",runit_out   , "munit_out  ",munit_out
      write(*,99) "tunit_out   ",tunit_out   , "vunit_out  ",vunit_out
      write(*,99) "Eunit_out   ",Eunit_out   , "rhounit_out",rhounit_out
      write(*,99) "muunit_out  ",muunit_out  , "gunit_out  ",gunit_out
      write(*,99) "tempunit_out",tempunit_out, "punit_out  ",punit_out
      write(*,99) "Lunit_out   ",Lunit_out
      write(*,*) ""
      
      write(*,97) "start ",start ,"anglexdeg",anglexdeg
      write(*,97) "finish",finish,"angleydeg",angleydeg
      write(*,97) "step  ",step  ,"anglezdeg",anglezdeg

      write(*,103) "flux_cal_dir",trim(adjustl(flux_cal_dir))
      
      write(*,*) ""
      write(*,101) "n           ",n
      write(*,101) "nkernel     ",nkernel
      write(*,100) "yscalconst  ",yscalconst
      write(*,100) "fracaccuracy",fracaccuracy
      write(*,100) "step1       ",step1
      write(*,100) "step2       ",step2
      write(*,100) "step3       ",step3
      write(*,100) "step4       ",step4
      write(*,100) "taulimit    ",taulimit
      write(*,100) "tau_thick   ",tau_thick
      write(*,102) "envfit      ",envfit
      write(*,*) ""
c      write(*,102) "rossonly    ",rossonly
      write(*,100) "metallicity ",metallicity
      write(*,100) "Rform       ",Rform
      
      write(*,*) ""
      write(*,102) "get_fluxes            ",get_fluxes

      write(*,102) "get_particles_at_pos  ",get_particles_at_pos
      write(*,102) "get_integration_at_pos",get_integration_at_pos
      
      if(get_particles_at_pos.or.get_integration_at_pos) then
         write(*,105) "posx",posx
         write(*,105) "posy",posy
      end if
      
      write(*,102) "get_closest_particles ",get_closest_particles
      write(*,102) "track_particles       ",track_particles
      
      if(track_particles) then
         write(*,104) "binary_tracking_file",binary_tracking_file
         write(*,106) "trackfile           ",trim(adjustl(trackfile))
      end if
      
      write(*,102) "get_info_of_particle  ",get_info_of_particle
      if(get_info_of_particle) then
         write(*,107) "info_particle", info_particle
      end if
      
      write(*,102) "rossonly              ",rossonly

      write(*,*) ""
      write(*,103) "eosfile        ",trim(adjustl(eosfile))
      write(*,103) "opacityfile    ",trim(adjustl(opacityfile))
      write(*,103) "opacitydustfile",trim(adjustl(opacitydustfile))
      write(*,103) "filtersfile    ",trim(adjustl(filtersfile))
      write(*,103) "outfile        ",trim(adjustl(outfile))
      write(*,*) ""

 200  format(A,"/defaults/",A)
      
      write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $     trim(adjustl(eosfile))
      inquire(file=trim(adjustl(inputfile)),exist=fileexists)
      if(.not. fileexists) then
         write(*,*) "Could not find EOS file '",
     $        trim(adjustl(eosfile)),"' in flux_cal_dir/defaults."
         error stop "init.f line 358"
      end if
      call readineostable(inputfile)

      call tabulinit

      write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $     trim(adjustl(opacityfile))
      inquire(file=trim(adjustl(inputfile)),exist=fileexists)
      if (.not. fileexists) then
         write(*,*) "Could not find opacity file '",
     $        trim(adjustl(opacityfile)),"'"
         write(*,*) "in flux_cal_dir/defaults. Make sure flux_cal_dir"
         write(*,*)"is correct, and if flux_cal was properly installed."
         error stop "init.f line 392"
      end if
      
      call readinkappatable(inputfile)

      if(len(trim(adjustl(opacitydustfile))).gt.0) then
         write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $        trim(adjustl(opacitydustfile))
         inquire(file=trim(adjustl(inputfile)),exist=fileexists)
         if(fileexists) call readinkappatabledust(inputfile)
      else
         write(*,*) "No dust opacity table found. Setting the radius at"
         write(*,*) "which dust begins to form to Rform=1.d30 to avoid"
         write(*,*) "any dust calculations."
         Rform=1.d30
      end if
      
c      inquire(file=opacitydustfile,exist=fileexists)
c      if(fileexists .and. (Rform .gt. 0.d0)) then
c         call readinkappatabledust
c      end if

      call ini_opac_dust_and_gas   ! cold T opacities initialization

      
      write(inputfile,200) trim(adjustl(flux_cal_dir)),"lowT_fa05_gs98_z0.02_x0.7.data"
      inquire(file=trim(adjustl(inputfile)), exist=fileexists)
      if(.not. fileexists) then
         write(*,*) "Could not find 'opac.dat' in flux_cal_dir/defaults"
         write(*,*) "Make sure flux_cal_dir is correct, and if flux_cal"
         write(*,*) "was properly installed."
         error stop "init.f line 402"
      end if

      call ini_opacity_photospehre_lowT(inputfile)
      call create_density_array
      
      write(inputfile,200) trim(adjustl(flux_cal_dir)),"table_nabla.dat"
      inquire(file=trim(adjustl(inputfile)),exist=fileexists)
      if(.not. fileexists) then
         write(*,*) "Could not find 'table_nable.dat' in "
         write(*,*) "flux_cal_dir/defaults. Make sure flux_cal_dir is"
         write(*,*) "correct, and if flux_cal was properly installed."
         error stop "init.f line 413"
      end if
      
      call ini_slops(inputfile)

c     Check to see if the filtersfile exists and read it if it does
      write(inputfile,200) trim(adjustl(flux_cal_dir)),
     $     trim(adjustl(filtersfile))
      inquire(file=trim(adjustl(inputfile)),exist=fileexists)
      if(.not.fileexists) then
         write(*,*) "Could not find '",trim(adjustl(filtersfile)),"' in"
         write(*,*) "flux_cal_dir/defaults. make sure flux_cal_dir is"
         write(*,*) "correct, and if flux_cal was properly installed."
         error stop "init.f line 426"
      end if
      
      write(*,*) "Reading filtersfile"
      open(16,file=trim(adjustl(inputfile)),status='old')
      read (16,*) numfilters
      
      do i=1,numfilters
         read(16,*) filtername(i), wavelength(i),absoluteflux(i)
         yscalfactor(i)=yscalconst*wavelength(i)
      enddo
      close(16)

      
      if(track_particles) then
         inquire(file=trim(adjustl(trackfile)),exist=fileexists)
         if(.not.fileexists) then
            write(*,*) "Could not find trackfile '",
     $           trim(adjustl(trackfile)),"'."
            error stop "init.f line 445"
         end if
         
         write(*,*) "Reading trackfile"
         open(17,file=trim(adjustl(trackfile)),status='old')
         do i=1,nmax
            read(17,*,end=120) track(i)
         enddo
 120     close(17)
      end if

      
      
c     Catch some runtime errors
      if ((start.gt.finish).and.(step.gt.0)) then
         write(*,*) "Must have start < finish when step > 0."
         write(*,*) "start=",start,"finish=",finish
         error stop "init.f line 387"
      end if
      if ((start.lt.finish).and.(step.lt.0)) then
         write(*,*) "Must have start > finish when step < 0."
         write(*,*) "start=",start,"finish=",finish
         error stop "init.f line 392"
      end if

      if(start.lt.0) then
         write(*,*) "Starting file index must be > 0."
         write(*,*) "start=",start
         error stop "init.f line 398"
      end if

      if((step.eq.0).or.(start.eq.finish)) then
         finish=start
         step=1
      end if

c      if(taulimit.lt.tau_thick) then
c         write(*,*) "WARNING: taulimit < tau_thick. This will result in"
c         write(*,*) "integrations that stop part way through optically"
c         write(*,*) "thin material, and as such, results should not be"
c         write(*,*) "expected to be physical."
c      end if
      
      write(*,*) "** Complete **"
      write(*,*) ""



c      write(*,*)" ___________________________________________________  "
c      write(*,*)"/         ___  _              ___        _          \ "
c      write(*,*)"|        | __|| | _  _ __ __ / __| __ _ | |         | "
c      write(*,*)"|        | _| | || || |\ \ /| (__ / _` || |         | "
c      write(*,*)"|        |_|  |_| \_,_|/_\_\ \___|\__,_||_|         | "
c      write(*,*)"\___________________________________________________/ "
c      write(*,*)"|                                                   | "
c      write(*,*)"|  FluxCal Copyright (C) 2018. Written by,          | "
c      write(*,*)"|                                                   | "
c      write(*,*)"|   Roger Hatfull      [Uni. of Alberta]            | "
c      write(*,*)"|   Natasha Ivanova    [Uni. of Alberta]            | "
c      write(*,*)"|   James C. Lombardi  [Allegheny College]          | "
cc      write(*,*)"|   Jose Nandez        [Uni. of West. Ontario]      | "
c      write(*,*)"|                                                   | "
c      write(*,*)"|  This program comes with ABSOLUTELY NO WARRANTY.  | "
c      write(*,*)"|  This is free software, and you are welcome to    | "
c      write(*,*)"|  redistribute it under certain conditions. See    | "
c      write(*,*)"|  LICENSE file for details, or visit               | "
c      write(*,*)"|  https://www.gnu.org/licenses                     | "
c      write(*,*)"\___________________________________________________/ "
c      write(*,*)"|                   |                                 "
c      write(*,*)"|   Version 0.0.1   |                                 "
c      write(*,*)"\___________________/                                 "
c      write(*,*)"                                                      "





      
      return
      end subroutine
