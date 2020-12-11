      subroutine derivs(s,tau,dtauds)
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

      integer i
      real*8 r2


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
      depth=s
      
c     t6 is the temperature in millions of degrees Kelvin

      if(t6.gt.0) then
         opacit = getOpacity(t6*1d6,rhocgs,xh)
c         if(opacit.eq.-1.d30) then
c            write(o,*) "derivs.f"
c            write(o,*) "xpos,ypos,zpos = ",xpos/runit_out,
c     $           ypos/runit_out,
c     $           zpos/runit_out
c            
c            do i=1,n
c               r2=(x(i)-xpos)**2.d0+(y(i)-ypos)**2.d0+(z(i)-zpos)**2.d0
c               if ( r2.lt.4.d0*hp(i)**2.d0 ) then
c                  write(o,*) "i = ",i
c               end if
c            end do
c            
c            stop
c         end if

         if(opacit.le.1d30) then
c            write(o,*) "opacit, rhocgs, dtauds(1) = ",opacit,rhocgs,
c     $           -opacit*rhocgs
            dtauds(1)=-opacit*rhocgs ! NEGATIVE SIGN IS BECAUSE STEP IS IN NEGATIVE DIRECTION
            dtauds(2)=dtauds(1)
            dtauds(4)=(t6*1d6)**4*exp(-tau(1))*dtauds(1) ! dtauds(1)and hence dtauds(4) are NEGATIVE
                                                         ! BECAUSE STEP IS IN NEGATIVE DIRECTION
c     We might want this later
c           do ifilter=1,numfilters
c               exponent=coeff/(wavelength(ifilter)*t6)
c               eexponent=exp(exponent)
c               denominator=eexponent-1
c               Bnu=2*planck*crad/wavelength(ifilter)**3d0/
c     $              denominator
c               dtauds(4+ifilter)=Bnu*exp(-tau(1))*dtauds(1) ! dtauds(1)and hence dtauds(4+ifilter) are
c                                                            ! NEGATIVE BECAUSE STEP IS IN NEGATIVE DIRECTION
c            enddo
         else
            dtauds(1)=0.d0
            dtauds(2)=0.d0
            dtauds(4)=0.d0
c     We might want this later
c            do ifilter=1,numfilters
c               dtauds(4+ifilter)=0.d0
c            enddo
         endif
         dtauds(3)=-1/xhp       ! NEGATIVE SIGN IS BECAUSE STEP IS IN NEGATIVE DIRECTION
      else
         dtauds(1)=0.d0
         dtauds(2)=0.d0
         dtauds(3)=0.d0
         dtauds(4)=0.d0
c     We might want this later
c         do ifilter=1,numfilters
c            dtauds(4+ifilter)=0.d0
c         enddo
      endif


      return
      end subroutine
