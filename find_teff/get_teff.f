      function get_teff(p_sph, t_sph, g_sph,do_debug)
c     Check out pg. 73 of "Stellar Structure and Evolution" by
c     R. Kippenhahn & A. Weigert
c     implicit none
      include '../lib/flux_cal.h'
      real*8 get_kappa
      real*8 p_sph, t_sph,  g_sph
      real*8 teff_high, teff_low, teff_mid
      real*8 p_high, p_low, p_mid
      real*8 rho_high, rho_low, rho_mid
      real*8 pg_high, pg_low, pg_mid
      real*8 kap_high, kap_low, kap_mid
      real*8 slop
      real*8 get_slop,tmin
      real*8 rho_tab
      
      integer i, count, do_debug
      real*8 get_density, rho_sph      

      if(log10(t_sph).ge.5.999) then
         if(do_debug.eq.1) write(*,*) "THIS PARTICLE IS TOO HOT TO BE OUTSIDE", t_sph
         get_teff=-10000.
         return
      end if
      
!     slop is 4.25 for radiative envelopes and 2.5 for adiabatic      
      slop=get_slop(t_sph,g_sph,tmin,slop_type)
      slop=1./slop

      if(log10(t_sph).le.tmin) then
         if(do_debug.eq.1) write(*,*) "THIS PARTICLE IS COLD HOT TO BE TRETAD", t_sph
         get_teff=-20000.
         return
      end if

      
      teff_high=20000.
      IF(teff_high.gt.t_sph) teff_high=t_sph
      teff_low=10.**(tmin-0.01) ! T can not be lower for this g
      teff_mid=(teff_high+4.*teff_low)/5.      
      
      count=0
      do i=1,40 
         p_high=p_sph*(teff_high/t_sph)**slop
         p_mid=p_sph*(teff_mid/t_sph)**slop
         p_low=p_sph*(teff_low/t_sph)**slop

         rho_high=get_density(p_high,teff_high)
         rho_mid=get_density(p_mid,teff_mid)
         rho_low=get_density(p_low,teff_low)

         
         kap_high=get_kappa(teff_high, rho_high)
         kap_mid=get_kappa(teff_mid, rho_mid)
         kap_low=get_kappa(teff_low, rho_low)
         
         pg_high=2./3.*g_sph/kap_high
         pg_mid=2./3.*g_sph/kap_mid
         pg_low=2./3.*g_sph/kap_low
         
         if(pg_high/p_high.le.1.and.pg_mid/p_mid.ge.1) then
            teff_low=teff_mid
            teff_mid=(teff_high+teff_low)/2.
         else
            teff_high=teff_mid
            teff_mid=(teff_high+teff_low)/2.
         end if

         count=count+1
         if( (teff_high-teff_low).le.10) exit
         
      end do      

      get_teff=teff_mid

      if(do_debug.eq.1) write(*,100)  1./slop,
     &     10.**tmin, tmin,g_sph, log10(g_sph),t_sph, p_sph, get_teff
 100  format(" DEBUG GETTEFF",F7.3," T_MIN=",
     &     F7.2,F7.3, " g=",F8.1, F7.3, " T_SPH=", E11.5,
     &     " P_SPH=", E11.5," T_EFF=",F8.1)

      return
      end
