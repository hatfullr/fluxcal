      function useeostable(ucgs,rhocgs,which)
c     which=1 gives temperature
c     which=2 gives mean molecular mass mu
c     which=3 gives pressure
c     which=4 gives entropy
c     which=5 gives Gamma1
c     which=6 gives Grad_ad
      include 'flux_cal.h'
      real*8 pgas,prad
      real*8 rhocgs,ucgs,beta1,temperature,gam1
      real*8 sgas
      integer iu,irho
      integer numrho,numu
      real*8 utable1,rhotable1
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 eostable(maxtablesize,maxtablesize,6)
      real*8 meanmu

      real*8 steprho,stepu
      common/eoscom/ numrho,numu,rhotable1,utable1,
     $     steprho,stepu,eostable

      real*8 f00,f01,f10,f11,log10rho,log10u,
     $     rholow,rhohigh,ulow,uhigh
      integer which

      rhooutsideTEOS = 0.d0
      aoutsideTEOS = 0.d0
      
c     utable(iu)=utable(1)+(iu-1)*stepu
c     so, iu= (utable(iu)-utable(1))/stepu + 1

      log10rho=log10(rhocgs)
      log10u=log10(ucgs)

c      irho = min(max(1,int((log10rho-rhotable1)/steprho +1)),numrho-1)
      irho = int((log10rho-rhotable1)/steprho +1)
      iu =   int((log10u-utable1)/stepu +1)

      if(irho.ge.1 .and. irho.le.numrho-1 .and. iu.ge.1 .and.
     $     iu.le.numu-1) then
c      if(irho.ge.1 .and. irho.le.numrho-1) then
         rholow=log10rho-(rhotable1+(irho-1)*steprho)
         rhohigh=rhotable1+irho*steprho-log10rho

         if(iu.ge.1 .and. iu.le.numu-1) then
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            
c     Use bi-linear interpolation among the four cartesian
c     grid points (irho,iu), (irho+1,iu), (irho,iu+1), and (irho+1,iu+1)
            f00=rholow*ulow
            f10=rhohigh*ulow
            f01=rholow*uhigh
            f11=rhohigh*uhigh
            
            useeostable=(f00*eostable(iu+1,irho+1,which)
     $           + f10*eostable(iu+1,      irho,  which)
     $           + f01*eostable(iu,        irho+1,which)
     $           + f11*eostable(iu,        irho,which))/(steprho*stepu)
         else if(iu.le.0) then
            
c     print *,'l iu=',log10u,log10rho,iu,irho
            
            if(which.ne.2) then
c     We are at very low specific internal energy u, where the pressure
c     or temperature should be nearly proportional to rho*u (at fixed
c     composition)
               
c     Use linear interpolation between the two cartesian
c     grid points (1,irho), (1,irho+1) and in addition extrapolate
c     to smaller u
               useeostable=ucgs/10d0**utable1*
     $              (rholow*eostable(1,irho+1,which)
     $              + rhohigh*eostable(1,irho,which))/steprho
            else
               useeostable=
     $              (rholow*eostable(1,irho+1,which)
     $              + rhohigh*eostable(1,irho,which))/steprho
            endif
            
         else if(iu.ge.numu) then
            
            if(which.eq.3) then
c     We are at very high specific internal energy u, where the pressure
c     nearly proportional to u (at fixed rho and composition)
               
c     Use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               useeostable=ucgs/10d0**(utable1+(numu-1)*stepu)*
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            else if(which.eq.1) then
c     We are at very high specific internal energy u, where the temperature
c     is nearly proportional to u^(1/4) (at fixed rho and composition)
               
c     Use linear interpolation between the two cartesian
c     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
c     to larger u
               useeostable=
     $              (ucgs/10d0**(utable1+(numu-1)*stepu))**0.25d0*
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            else
               useeostable=
     $              (rholow*eostable(numu,irho+1,which)
     $              + rhohigh*eostable(numu,irho,which))/steprho
            endif
            
         endif
      else
c     At extreme densities we will use ideal gas + radiation pressure
         if(irho.lt.1) irho=1
         if(irho.gt.numrho) irho=numrho

         if(iu.lt.1) then
            meanmu=eostable(1,irho,2)
         else if (iu.ge.numu) then
            meanmu=eostable(numu,irho,2)
         else
            ulow=log10u-(utable1+(iu-1)*stepu)
            uhigh=utable1+iu*stepu-log10u
            meanmu=(ulow*eostable(iu+1,irho,2)
     $           + uhigh*eostable(iu,  irho,2))/stepu
         endif

         
         if(which.eq.2) then
            useeostable=meanmu
            return
         endif
         call getTemperature(3.d0*boltz*rhocgs/meanmu/arad/2.0d0,
     $        -ucgs*rhocgs/arad,temperature)
         rhooutsideTEOS = rhocgs
         aoutsideTEOS = ucgs
         
         if(which.eq.1) then
            useeostable=temperature
            return
         endif
         pgas=rhocgs*boltz*temperature/meanmu
         prad=arad*temperature**4/3.d0
         beta1=pgas/(pgas+prad)
         gam1=(32.d0-24.d0*beta1-3.d0*beta1**2) /
     $        (24.d0-21.d0*beta1)
         if(which.eq.5) then 
            useeostable=gam1
            return
         endif
         sgas=1.5d0*boltz/meanmu*
     $        log(ucgs*rhocgs**(-2.d0/3.d0))
     $          +1.5d0*boltz/meanmu
     $          *(5.d0/3.d0+log(4*pi/(3*planck**2)))        

         if(which.eq.4) then
            useeostable=sgas
            return
         endif 
         if(gam1.lt.0.999*4.d0/3.d0 .or. gam1.gt.1.001*5.d0/3.d0) then
            write(o,*) "WARNING: useeostable found a bad gamma value."
            write(o,*) "4/3 < gamma < 5/3"
c            write(o,*) beta1,pgas,prad,temperature,
c     $           -ucgs*rhocgs/arad
c            stop 'GAM1 value does not make sense'
         endif
         useeostable=pgas+prad
         
      endif
         
      end
