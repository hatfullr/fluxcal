      subroutine derivs2(s,tau,dtauds)
c     L1055 to L1242
      include 'flux_cal.h'
      INTEGER nrhs
      real*8 s,tau(*),dtauds(*)
      real*8 opacit, opacit1, opacit2
      COMMON /countrhs/nrhs
      real*8 XPOS, YPOS, rhocgs,xhp,ZPOS
      COMMON/TAUGRID/XPOS,YPOS
      real xh,t6,r
      common/localQuantities/ rhocgs,xh,t6,xhp
      real*8 fcondense
      integer ifilter
      real*8 Bnu,exponent
      real*8 eexponent,denominator
      real*8 nn

c      real*8 getOpacity
c      external usetable,usetabledust
c      external usetable2,usetabledust2

      real*8 tauestimate

      integer lastpart
      common/lastparticle/ lastpart
      


c     1990: Table 3 of http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1990PASP..102.1181B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
c      parameter(wavelength(1)=436.3, ! central wavelength of B filter in nm
c     $          wavelength(2)=544.8) ! central wavelength of V filter in nm
c     http://www.annualreviews.org/doi/pdf/10.1146/annurev.astro.41.082801.100251
c      parameter(wavelength(1)=436.1, ! central wavelength of B filter in nm
c     $          wavelength(2)=544.8) ! central wavelength of V filter in nm

c     1998: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1998A%26A...333..231B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
c      data wavelength /438d-7, 545d-7/ ! central wavelength of B and V filters in cm

c     These data come from Table A2 of Bessell et al. 1998 (http://adsabs.harvard.edu/abs/1998A%26A...333..231B)
c      data wavelength /366d-7, 438d-7, 545d-7, 641d-7, 798d-7,
c     $     1220d-7, 1630d-7, 2190d-7, 3450d-7/ ! central wavelength of UBVRIJHKL filters in cm
!      data wavelength /2124d-7,3776d-7,4769d-7/ ! central wavelength of NIRC-2 Kp and Lp filters as well as M filter

c     The NIRC-2 Kp filter has central wavelength of 2.124d-7cm, while
c     the NIRC-2 Lp filter has central wavelength of 3.776d-7cm (http://www2.keck.hawaii.edu/inst/nirc2/filters.html, http://www2.keck.hawaii.edu/inst/nirc2old/Manual/ObserversManual.html)

c     Kp central wavelength is 2120d-7cm = 2.12 microns, while L* is 3800d-7cm = 3.8 microns (Table A2 of Bessell et al.)
c     Ks central wavelength is 2157d-7cm = 2.157 microns, while L' is 3761d-7cm = 3.761 microns, and M is 4769d-7cm = 4.769 microns (Table 7.5 of Tokunga 2000: see http://www.ifa.hawaii.edu/~tokunaga/)

c     Most observations of V838 Mon used Cousins filter system for UBVRI observations, and 
c     the SAAO .75m telescope for JHKL observations, so these data use Bessell et al. 1998
c     as a reference for the UBVRI central wavelengths,
c     and http://old.saao.ac.za/facilities/instrumentation/infrared-photometer-mk-ii/ as a reference
c     for JHKL central wavelengths:
c      data wavelength /366d-7, 438d-7, 545d-7, 641d-7, 798d-7,
c     $     1250d-7, 1600d-7, 2200d-7, 3400d-7/ ! central wavelength of UBVRIJHKL filters in cm


c     Note, for future reference:
c     magλ = -2.5 log (fλ) - 21.100 - zp(fλ)
c     magν = -2.5 log (fν) - 48.598 - zp(fν)
c     where zp() is a zero-point magnitude.
      nrhs=nrhs+1

      call getLocalQuantities(XPOS,YPOS,s)

      ZPOS=s
      
c     t6 is the temperature in millions of degrees Kelvin

      if(t6.gt.0) then

         if ((xpos**2+ypos**2+zpos**2)**0.5d0.lt.Rform) then
            
c     We're too close to the origin for there to be dust
            opacit=getOpacity(ncooling,dble(t6*1d6),rhocgs,tau(1),
     $           usetable,usetable2)
         else

c     There will be some dust at this location
c     opacit1 is local Rosseland opacity without dust or molecules
c     opacit2 is local Rosseland opacity with dust and molecules
            opacit1=getOpacity(ncooling,dble(t6*1d6),rhocgs,tau(1),
     $           usetable,usetable2)
            opacit2=getOpacity(ncooling,dble(t6*1d6),rhocgs,tau(1),
     $           usetabledust,usetabledust2)

c     here we are setting f=f_max*(1-(r_form/r)^2) which comes from combining
c     eqn 5 of Nanni et al. (2013, http://adsabs.harvard.edu/abs/2013MNRAS.434.2390N)
c     with eqn 4 of Kochanek (2014, http://adsabs.harvard.edu/abs/2014MNRAS.444.2043K).
c     Kochanek uses a velocity scaling exponent n=0, while we use n=1:
            nn=1
c     the 1 at the beginning is what we have set f_max to. 0<f_max<1 so we assume 1 for now
            fcondense=1*(1-(Rform/(xpos**2+ypos**2+
     $           zpos**2)**0.5d0)**(2+nn))**3
               
c     here we put our reduced condensation function f into eqn4 of Nanni et al 2013
c     to find the opacity with partial dust and molecule formation
            opacit=(1-fcondense)*opacit1+fcondense*opacit2
            
            
         endif

c         write(*,*) "opacit, t6 = ", opacit, t6
         
         if(opacit.le.1d30) then
            dtauds(1)=-opacit*rhocgs ! NEGATIVE SIGN IS BECAUSE STEP IS IN NEGATIVE DIRECTION
            dtauds(2)=dtauds(1)
            dtauds(4)=(t6*1d6)**4*exp(-tau(1))*dtauds(1) ! dtauds(1)and hence dtauds(4) are NEGATIVE
                                                         ! BECAUSE STEP IS IN NEGATIVE DIRECTION
            do ifilter=1,numfilters
               exponent=coeff/(wavelength(ifilter)*t6)
               eexponent=exp(exponent)
               denominator=eexponent-1
               Bnu=2*planck*crad/wavelength(ifilter)**3d0/
     $              denominator
               dtauds(4+ifilter)=Bnu*exp(-tau(1))*dtauds(1) ! dtauds(1)and hence dtauds(4+ifilter) are
                                                            ! NEGATIVE BECAUSE STEP IS IN NEGATIVE DIRECTION
            enddo
         else
            dtauds(1)=0.d0
            dtauds(2)=0.d0
            dtauds(4)=0.d0
            do ifilter=1,numfilters
               dtauds(4+ifilter)=0.d0
            enddo
         endif

c     call usetable(rhocgs,dble(t6*1d6),opacit)
c     opacit=opacity(rhocgs,t6*1d6)
         dtauds(3)=-1/xhp       ! NEGATIVE SIGN IS BECAUSE STEP IS IN NEGATIVE DIRECTION
      else
         dtauds(1)=0.d0
         dtauds(2)=0.d0
         dtauds(3)=0.d0
         dtauds(4)=0.d0
         do ifilter=1,numfilters
            dtauds(4+ifilter)=0.d0
         enddo
      endif
      
      return
      end subroutine
