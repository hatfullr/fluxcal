      subroutine getOpacitySub(tem,rho,tau,opacit)
c      include 'flux_cal.h'

      real*8 tem,x,y,z,rho,tau

      real*8 getOpacity

      real*8 opacit,opacit1,opacit2,fcondense
      integer nn

      if(tem.gt.0.d0) then
         opacit = getOpacity(tem,rho,tau)

         
c         if((x**2.d0 + y**2.d0 + z**2.d0)**0.5d0 .lt. Rform) then
c            opacit=getOpacity(tem,rho,tau,usetable,usetable2)
c         else
c            opacit1=getOpacity(tem,rho,tau,usetable,usetable2)
c            opacit2=getOpacity(tem,rho,tau,usetabledust,
c     $           usetabledust2)
cc     here we are setting f=f_max*(1-(r_form/r)^2) which comes from combining
cc     eqn 5 of Nanni et al. (2013, http://adsabs.harvard.edu/abs/2013MNRAS.434.2390N)
cc     with eqn 4 of Kochanek (2014, http://adsabs.harvard.edu/abs/2014MNRAS.444.2043K).
cc     Kochanek uses a velocity scaling exponent n=0, while we use n=1:
c            nn=1
cc     the 1 at the beginning is what we have set f_max to. 0<f_max<1 so we assume 1 for now         
c            fcondense=1.d0*(1.d0-(Rform/(x**2.d0+y**2.d0
c     $           +z**2.d0)**0.5d0)**(2+nn))**3.d0
c            
cc     here we put our reduced condensation function f into eqn4 of Nanni et al 2013
cc     to find the opacity with partial dust and molecule formation
c            opacit=(1.d0-fcondense)*opacit1+fcondense*opacit2
c         end if
c      else
c         opacit=0.d0
      end if


      return      
      end subroutine
