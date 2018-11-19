      subroutine getOpacitySub(x,y,z,tem,rho,tau,ncooling,Rform,opacit)
c      include 'flux_cal.h'

      real*8 tem,x,y,z,Rform,rho,tau
      integer ncooling

      real*8 getOpacity
      external usetable,usetable2
      external usetabledust,usetabledust2

      real*8 opacit,opacit1,opacit2,fcondense
      integer nn

      if(tem.gt.0.d0) then
         if((x**2.d0 + y**2.d0 + z**2.d0)**0.5d0 .lt. Rform) then
            opacit=getOpacity(ncooling,tem,rho,tau,usetable,usetable2)
         else
            opacit1=getOpacity(ncooling,tem,rho,tau,usetable,usetable2)
            opacit2=getOpacity(ncooling,tem,rho,tau,usetabledust,
     $           usetabledust2)
            nn=1
            fcondense=1.d0*(1.d0-(Rform/(x**2.d0+y**2.d0
     $           +z**2.d0)**0.5d0)**(2+nn))**3.d0
            
            opacit=(1.d0-fcondense)*opacit1+fcondense*opacit2
         end if
      else
         opacit=0.d0
      end if


      return      
      end subroutine
