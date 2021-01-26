      function get_teff(p_sph, t_sph, g_sph,do_debug)
c     Check out pg. 73 of "Stellar Structure and Evolution" by
c     R. Kippenhahn & A. Weigert (1994)
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
      
      integer i, do_debug
      real*8 get_density, rho_sph      

c      if(log10(t_sph).gt.6) then
cc      if(log10(t_sph).ge.5.999) then
c         if(do_debug.eq.1) write(o,*) "THIS PARTICLE IS TOO HOT TO BE"//
c     $        " OUTSIDE", t_sph
c         get_teff=-10000.
c         return
c      end if
      
!     slop is 4.25 for radiative envelopes and 2.5 for adiabatic
!     Use the table_nabla.dat to find nabla ("slop") for values of
!     t_sph and g_sph. tmin and slop_type are outputs.
      slop=get_slop(t_sph,g_sph,tmin,slop_type)

      if(slop .le. -10.d0) then
         if(do_debug.eq.1) write(o,*) "No proper slope found."
         get_teff = 0.d0
         return
      end if
      
      slop=1./slop

c      if(log10(t_sph).lt.tmin) then
c         if(do_debug.eq.1) write(o,*) "THIS PARTICLE IS TOO COLD"//
c     $        " TO BE TREATED", t_sph
c         get_teff=-20000.
c         return
c      end if

      
      teff_high=20000.
      IF(teff_high.gt.t_sph) teff_high=t_sph
      teff_low=10.**(tmin-0.01) ! T can not be lower for this g
      teff_mid=(teff_high+4.*teff_low)/5. ! Average is heavily weighted toward the lower temperature

      ! Only uncomment this and associated lines for extreme debugging
c     open(50,file="teff_trace.dat",status='unknown')

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

         ! pg comes from the pressure at the photosphere. The full physical term
         ! has the "Eddington factor" in it, which is related to the Eddington limit.
         ! To get Pg, we say the "Eddington factor" is small. This is not true for
         ! high opacity stars and luminous stars, but is true for the Sun.
         ! Derivation stars with HE. Assumption is made that g is constant and equal
         ! to GM/R^2 and that opacity is not a function of tau. The kappa written below
         ! is meant to be the kappa at the photosphere. Therefore, g_sph is g at the
         ! surface.
         ! This assumes that g_sph = g_sph throughout the kernel, which is okay.
         ! g shouldn't be allowed to change through the kernel because the
         ! kernel function is not physical. These pressures represent the
         ! minimum allowed pressure values, as they are at the photosphere.
         pg_high=2./3.*g_sph/kap_high
         pg_mid=2./3.*g_sph/kap_mid
         pg_low=2./3.*g_sph/kap_low

         ! If the pressure at the photosphere is <= envelope solution pressure
         ! (at the upper bound) AND the middle photospheric pressure is >=
         ! the envelope solution pressure (at the midpoint).
         ! If the highest pressure is higher than the highest photospheric pressure
         ! and the "actual" pressure is lower than the "actual" photospheric pressure
         if(pg_high/p_high.le.1.and.pg_mid/p_mid.ge.1) then
            teff_low=teff_mid   ! Increase
            teff_mid=(teff_high+teff_low)/2.
c            write(o,*) "I increased Teff_low"
         ! Otherwise, the highest pressure must be lower than the highest photospheric
         ! pressure AND the "actual" pressure must be higher than the "actual"
         ! photospheric pressure, so we need to take decreasing steps.
         else
c            write(o,*) "I decreased teff_high"
            teff_high=teff_mid  ! Decrease
            teff_mid=(teff_high+teff_low)/2.
         end if

c         write(50,*) p_mid,rho_mid,kap_mid,pg_mid,teff_mid,p_sph,g_sph,
c     $        t_sph
         if( (teff_high-teff_low).le.10) exit
         
      end do
c      close(50)
      get_teff=teff_mid

c      if(do_debug.eq.1) write(o,100)  1./slop,
c     &     10.**tmin, tmin,g_sph, log10(g_sph),t_sph, p_sph, get_teff
 100  format(" DEBUG GETTEFF",F7.3," T_MIN=",
     &     F7.2,F7.3, " g=",F8.1, F7.3, " T_SPH=", E11.5,
     &     " P_SPH=", E11.5," T_EFF=",F8.1)

      return
      end
