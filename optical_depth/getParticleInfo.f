      subroutine getParticleInfo(m)
      include 'optical_depth.h'

      character*255 filename
      
c      integer file
      integer m

      real*8 central_tau
      real*8 opacit,scgs

      
c      call getOpacitySub(x(m),y(m),z(m),tempp(m),rho(m),0.d0,ncooling,
c     $     Rform,opacit)
c      tauA = taucoef*am(m)*opacit/hp(m)**2.d0

c      if(tauA(m).gt.tau_thick) then
c         Teff = get_teff(pp(m),tempp(m),localg(m))
c      else
c         Teff = 0.d0
c      end if

      if(a(m).gt.0.d0) then
         scgs = useeostable(a(m),rho(m),4)
      else
         scgs = 0.d0
      end if
      
c      open(33,file=trim(adjustl(filename)),action='write',
c     $           position='append')
c     100  format(18E22.14)
c      call makeOutputFile(filename)
      call output(pinfo_file,m)
c      write(file,100) t,x(m),y(m),z(m),am(m),hp(m),rho(m),a(m),
c     $     wmeanmolecular(m),localg(m),tempp(m),pp(m),scgs,opacit,
c     $     Teff(m),tauA(m),Avis(m),sigma*4d0*Avis(m)*Teff(m)**4.d0

      
      end subroutine
      
