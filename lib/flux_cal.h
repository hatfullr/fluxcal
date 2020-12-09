      implicit none
      integer nmax,nnmax,ntab
      PARAMETER (NMAX=420000,NNMAX=128) ! Also in odeint.f
c      PARAMETER (NMAX=301000,NNMAX=128)
c      PARAMETER (NMAX=221000,NNMAX=1000)
c      PARAMETER (NMAX=180000,NNMAX=100)
      PARAMETER (NTAB=100000)
      real*8 hmin,hmax
      integer n,nout,nit
      common/intpar/ hmin,hmax,n
      common/out/ nout,nit
      real*8 x(nmax),y(nmax),z(nmax),am(nmax),hp(nmax),rho(nmax),
     $      vx(nmax),vy(nmax),vz(nmax),a(nmax),wmeanmolecular(nmax),
     $      localg(nmax),metal(nmax),tempp(nmax),pp(nmax),entropy(nmax)
      common/part/ x,y,z,am,hp,rho,vx,vy,vz,a,wmeanmolecular,localg,
     $      metal,tempp,pp,entropy
	   
      real*8 pi,kelvin,gram,sec,cm,erg,boltz,crad,planck,crad2,
     $      sigma,arad,qconst,coeff,pc,distance,gravconst
      parameter (pi=3.1415926535897932384626d0)
      common/consts/ kelvin,gram,sec,cm,erg,boltz,crad,planck,crad2,
     $      sigma,arad,qconst,coeff,pc,distance,!munit,runit,
     $      gravconst
      namelist/baseunits/ gram,sec,cm,kelvin
      real*8 wtab(ntab),ctab
      common/wtabul/ wtab,ctab

      integer neighbors(nmax) ! Remove later
      common/neigh/ neighbors ! Remove later
	
      real*8 mu
      common/surfaceangle/ mu
	
      CHARACTER*12 outfilename
      integer outfilenum,iform,iout,ni1,ni2,ni3,innit,nnit
      common/counter/ innit
      real*8 anglexdeg,angleydeg,anglezdeg
      common/azimuth/anglexdeg,angleydeg,anglezdeg
      integer numfilters,numfiltersmax
      parameter (numfiltersmax=15)
      character*2 filtername(numfiltersmax)
      real*8 wavelength(numfiltersmax),absoluteflux(numfiltersmax)
      common /splotfilters/ wavelength,absoluteflux,numfilters
      common /filternames/ filtername

      logical opacity_oob_error,opacity_oob_warning,opacity_analytic_warning

      real*8 boundfree,freefree,escattering,negativeHion

      common/OOB_opacity_behavior/ opacity_oob_error,opacity_oob_warning,
     $      opacity_analytic_warning
	
      integer numopacityfiles
      parameter (numopacityfiles=9)
      character*255 opacityfiles(numopacityfiles)
      common/common_opacityfiles/ opacityfiles
      real*8 logTmins(numopacityfiles), logTmaxs(numopacityfiles)
      common/common_Tminsmaxs/ logTmins,logTmaxs
      real*8 logRmins(numopacityfiles), logRmaxs(numopacityfiles)
      common/common_Rminsmaxs/ logRmins, logRmaxs
      real*8 logT_blend1(numopacityfiles),
     $      logT_blend2(numopacityfiles),
     $      logR_blend1(numopacityfiles),
     $      logR_blend2(numopacityfiles)
      common/blending/ logT_blend1, logT_blend2,
     $      logR_blend1, logR_blend2
      integer nrows_opac_table,ncols_opac_table
      parameter (nrows_opac_table=300,ncols_opac_table=300)
      integer nrows_tables(numopacityfiles), ncols_tables(numopacityfiles)
      real*8 logopacitytables(numopacityfiles,
     $      nrows_opac_table,ncols_opac_table)
      real*8 logTs(numopacityfiles,nrows_opac_table),
     $     logRs(numopacityfiles,ncols_opac_table)

      real*8 opacity_rho_cutoff

      common/rhocutoff/ opacity_rho_cutoff
	
      common/TandR/ logTs, logRs
      common/tabledimens/ nrows_tables, ncols_tables
      common/opacity_hiT_tables/ logopacitytables
      
      character*255 opacityfile_rosseland,opacityfile_planck,
     $      filtersfile,trackfile,eosfile
      real*8 yscalconst,fracaccuracy,metallicity
      integer nkernel
      common/kernels/ nkernel
      common/simulation/ yscalconst,fracaccuracy,metallicity
      real*8 munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $      runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $      rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $      Lunit_out,kunit_out,sunit_out
      common/units/ munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $      runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $      rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $      Lunit_out,kunit_out,sunit_out

      integer start,finish,step
      common/filevars/ start,finish,step
      logical get_particles_at_pos,
     $      track_particles,binary_tracking_file,get_fluxes,
     $      get_integration_at_pos,get_info_of_particle,
     $      get_closest_particles,get_integration_at_all_pos,
     $      get_true_luminosity
      real*8 step1,step2,step3,step4,simps_alpha,
     $      simps_max_frac_dF,simps_max_dtau,simps_max_dz,
     $      simps_min_frac_dF,simps_min_dtau,simps_min_dz,simps_F_cutoff
      logical simps_max_step_error,simps_max_step_warning
      integer MAXSTP
      common/steps/ step1,step2,step3,step4,MAXSTP
      common/simpsonrule/simps_alpha,simps_min_dz,simps_max_dz,
     $      simps_min_frac_dF,simps_max_frac_dF,simps_min_dtau,
     $      simps_max_dtau,simps_max_step_error,simps_max_step_warning
      real*8 Tplanck
      common/opacity_controls/ Tplanck
      real*8 taulimit,posx,posy,taulimit_threshold
      real*8 tau_thick_envfit,tau_thick_integrator,tau_thick
      common/tauthick/ tau_thick_integrator,tau_thick_envfit,tau_thick
      integer info_particle
      logical envfit
      common/whatever/ info_particle
      character*255 flux_cal_dir
      common/inputstuff/ taulimit,posx,posy,get_particles_at_pos,
     $      track_particles,binary_tracking_file,get_fluxes,
     $      get_info_of_particle,envfit,taulimit_threshold,
     $      get_closest_particles,flux_cal_dir
      common/getint/ get_integration_at_pos,get_integration_at_all_pos,
     $      get_true_luminosity
      character*255 outfile
      common/outputfile/ outfile
      real*8 t
      common/time/ t
      character*3 dust_model
      character*1 dust_topology, dust_shape
      common/dust/ dust_model,dust_topology,dust_shape
      logical debug
      common/dbg/ debug
      logical track_all
      common/trackall/ track_all
      real*8 smoothing_window_T,smoothing_window_rho
      common/smoothingwindow/ smoothing_window_T,smoothing_window_rho
      namelist/input/ yscalconst,munit,runit,tunit,vunit,
     $      fracaccuracy,Eunit,rhounit,muunit,gunit,
     $      runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $      rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $      Lunit_out,dust_model,dust_topology,dust_shape,
     $      nkernel,opacityfile_rosseland,opacityfile_planck,filtersfile,
     $      metallicity,anglexdeg,angleydeg,anglezdeg,start,finish,step,
     $      step1,step2,step3,step4,taulimit,
     $      get_particles_at_pos,posx,posy,track_particles,trackfile,
     $      binary_tracking_file,get_fluxes,envfit,
     $      get_integration_at_pos,get_integration_at_all_pos,
     $      get_info_of_particle,info_particle,
     $      outfile,tau_thick_envfit,eosfile,get_true_luminosity,
     $      get_closest_particles,flux_cal_dir,tau_thick_integrator,
     $      tau_thick,track_all,kunit_out,sunit_out,opacityfiles,logTmins,
     $      logTmaxs,logRmins,logRmaxs,logT_blend1,logT_blend2,
     $      logR_blend1,logR_blend2,MAXSTP,opacity_oob_error,
     $      opacity_oob_warning,opacity_analytic_warning,
     $      opacity_rho_cutoff,debug,smoothing_window_T,
     $      smoothing_window_rho,simps_alpha,simps_max_frac_dF,integrator,
     $      simps_max_step_error,simps_max_step_warning,simps_max_dtau,
     $      taulimit_threshold,simps_max_dz,simps_min_dz,
     $      simps_min_frac_dF,simps_min_dtau,simps_F_cutoff
	
      common/inputfilenames/ filtersfile,trackfile,eosfile

      real*8 yscalfactor(numfiltersmax)
      common /yscales/ yscalfactor

      real*8 wmeanmu
      common/meanmu/ wmeanmu

      integer track(nmax)

      common/tracking/ track

      character*255 infname
      common/inputfilename/ infname

      real*8 natalum,sphlum,natalumhi,natalumlo
      common/splitlums/ natalum,sphlum,natalumhi,natalumlo


      real*8 getOpacity
      real*8 useeostable

      real*8 get_teff,getLocalAngle
      real*8 fourPointArea, fourPointArea2

      real*8 taucoef
      common/taucoefficient/ taucoef

      real*8 tauA(nmax), opac_sph(nmax)
      common/tauAopticaldepth/ tauA, opac_sph
	   
      real*8 Avis(nmax)
      common/surfaceareas/ Avis

      real*8 avgt
      integer numcell
      common/averageT/ avgt
      common/ncell/ numcell

      external quicksort,quicksort2

      real*8 Teff(nmax)
      common/teff/ Teff

c Opacity tables for enevelope fitting
      real*8 gridR(100)
      real*8 gridT(200)
      real*8 kap(100,200)
      integer numo, numr
      common /opacity/ kap, gridT, gridR, numo, numr
	
c     slops for envelope fitting
      real*8 grid_gn(400)
      real*8 grid_tn(400)
      real*8 nabla(400,400)
      integer num_g, num_t, slop_type
      common /slops/ nabla, grid_gn, grid_tn, num_g, num_t, slop_type

      logical printIntegrationSteps
      common/intatpos/ printIntegrationSteps

      character*255 pinfo_file
      common/infofile/ pinfo_file

      logical isInitGrid
      common/isgridinit/ isInitGrid

      logical dimenFileAlreadyExists
      common/dimenfilecheck/ dimenFileAlreadyExists

      real*8 Atot
      common/totalarea/ Atot

      real*8 table_nabla_Tmin,table_nabla_Tmax,
     $      table_nabla_gmin,table_nabla_gmax

      common/tablerange/ table_nabla_Tmin,table_nabla_Tmax,
     $      table_nabla_gmin,table_nabla_gmax

      real*8 rhooutsideTEOS,aoutsideTEOS
      common/outsideTEOS/ rhooutsideTEOS,aoutsideTEOS

      real*8 max_step_size,min_step_size
      integer min_steps_taken,max_steps_taken,ndtau_increase,
     $      ndF_increase,ndz_increase,
     $      ndtau_decrease,ndF_decrease,ndz_decrease,ndtaulimit_decrease,
     $      ndtaulimit_increase
      integer nstp
      common/boundsteptaken/ max_step_size,min_step_size,min_steps_taken,
     $      max_steps_taken,nstp,ndtau_increase,ndF_increase,ndz_increase,
     $      ndtau_decrease,ndF_decrease,ndz_decrease,ndtaulimit_decrease,
     $      ndtaulimit_increase
      
      real*8 taylor_expansion
	
      real*8 depth
      common/depthtracker/ depth

      real*8 rhoPlanckRossmin, rhoPlanckRossmax
      parameter(rhoPlanckRossmin=2.d-18,rhoPlanckRossmax=1.d-7)

      integer integrator
      common/integrator_type/ integrator
