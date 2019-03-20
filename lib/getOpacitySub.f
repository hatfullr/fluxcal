      subroutine getOpacitySub(x,y,z,tem,rho,tau,Rform,opacit)
c      include 'flux_cal.h'

      real*8 tem,x,y,z,Rform,rho,tau

      real*8 getOpacity
      external usetable,usetable2
      external usetabledust,usetabledust2

      real*8 opacit,opacit1,opacit2,fcondense
      integer nn

      if(tem.gt.0.d0) then
         if((x**2.d0 + y**2.d0 + z**2.d0)**0.5d0 .lt. Rform) then
            opacit=getOpacity(tem,rho,tau,usetable,usetable2)
         else
            opacit1=getOpacity(tem,rho,tau,usetable,usetable2)
            opacit2=getOpacity(tem,rho,tau,usetabledust,
     $           usetabledust2)
c     here we are setting f=f_max*(1-(r_form/r)^2) which comes from combining
c     eqn 5 of Nanni et al. (2013, http://adsabs.harvard.edu/abs/2013MNRAS.434.2390N)
c     with eqn 4 of Kochanek (2014, http://adsabs.harvard.edu/abs/2014MNRAS.444.2043K).
c     Kochanek uses a velocity scaling exponent n=0, while we use n=1:
            nn=1
c     the 1 at the beginning is what we have set f_max to. 0<f_max<1 so we assume 1 for now         
            fcondense=1.d0*(1.d0-(Rform/(x**2.d0+y**2.d0
     $           +z**2.d0)**0.5d0)**(2+nn))**3.d0
            
c     here we put our reduced condensation function f into eqn4 of Nanni et al 2013
c     to find the opacity with partial dust and molecule formation
            opacit=(1.d0-fcondense)*opacit1+fcondense*opacit2
         end if
      else
         opacit=0.d0
      end if


      return      
      end subroutine
