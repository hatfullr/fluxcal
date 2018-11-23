      subroutine getTpractical(z0,z1,zthick,thick_p,myh1,
     $     Tthin, Tthick, tau_thin)
c     Input = z0,z1,zthick,thick_p,h1
c     Output = Tthin, Tthick, tau_thin
      
      include 'optical_depth.h'
      integer i

      real*8 zstop
      real*8 z0,z1,zthick
      real*8 tau_thin
      integer thick_p
      real*8 myh1
      real*8 Tthin, Tthick

      if(envfit) then           ! If user wants to use get_teff AND ray-trace
         zstop = max(z0,zthick)
         if(z1.ne.zstop) then   ! If there is optically thin material
            nrhs=0
            do ifilter=1,numfilters+4
               taustart(ifilter)=0.d0
            enddo
            call odeint(taustart,4+numfilters,z1,zstop,
     $           eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
     $           rkqs)
            tau_thin = taustart(1)
            ! We need to find the temperature contribution from the 
            ! optically thin material in front of the optically thick
            ! taustart(4) is already attenuated in derivs2.f:
            ! dtauds(4)=(t6*1d6)**4*exp(-tau(1))*dtauds(1)
            Tthin=taustart(4)**0.25d0

            nphoto = taustart(3)
            sphoto = zstop
         else                   ! There is no optically thin material
            tau_thin = 0.d0
            Tthin=0.d0
            nphoto=0.
            sphoto=z0
         end if
         
         if(thick_p.gt.0) then
c     T of an optically thick particle
            
            if(Teff(thick_p).lt.0.and.Teff(thick_p).ge.-15000) then
               write(*,*) "a really hot particle is outside"
               stop
            end if

c            if(Teff(thick_p).ge.7000) then
c               write(*,*) "a moderately hot particle is outside:",
c     &              tempp(thick_p),Teff(thick_p), thick_p
c            end if

            
            if(Teff(thick_p).lt.-15000.and.Teff(thick_p).ge.-25000) then
c               write(*,*) "a cold particle is outside",
c     &              tempp(thick_p), rho(thick_p),
c     &              tauA(thick_p), opac_sph(thick_p),
c     &              hp(thick_p)*2./6.06e10,thick_p
               Tthick=tempp(thick_p)
            end if
            
            if(Teff(thick_p).gt.0) then
               Tthick = Teff(thick_p) ! Teff is calculated in lib/read_fluxcal.f
            end if            

            ! Use the smallest temperature possible, because the SPH methods will be more accurate
            ! at low T than the envelope fitting routine
            if(tempp(thick_p) .lt. Tthick) then
               write(*,*) " strange place. Why I am here?",
     &              tempp(thick_p),Teff(thick_p), Tthick, thick_p
               Tthick = tempp(thick_p)
               stop              
            end if
            
         else
            Tthick = 0.d0
         end if
      else                      ! User wants to use ONLY ray-trace
         nrhs=0
         do ifilter=1,numfilters+4
            taustart(ifilter)=0.d0
         enddo
         call odeint(taustart,4+numfilters,z1,z0,
     $        eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
     $        rkqs)
         sa=s(kount1-1)         ! z-position w/ optical depth tau_tot<tau_thick
         sb=s(kount1)           ! z-position w/ optical depth tau_tot>=tau_thick
c     Interpolate to find position where optical tau_tot=tau_thick:
         sphoto=(sa*(tau(1,kount1)-1.d0)
     $        +sb*(1.d0-tau(1,kount1-1)))/
     $        (tau(1,kount1)-tau(1,kount1-1))
         na=tau(3,kount1-1)
         nb=tau(3,kount1)
c     Interpolate to find number of smoothing lengths in at tau_tot=tau_thick:
         nphoto=(na*(tau(1,kount1)-1.d0)
     $        +nb*(1.d0-tau(1,kount1-1)))/
     $        (tau(1,kount1)-tau(1,kount1-1))
         if(kount1.gt.1) then
            Tthick=taustart(4)**0.25d0
         else
            kount1=0
            Tthick=0.d0
         end if
         
         Tthin=0.d0

         tau_thin=0.d0
         
      end if

      if(tau_thin .lt. 0.d0) then
         write(*,*) "tauthin < 0! Something went wrong in"//
     $        " getTpractical.f"
         write(*,*) "tauthin = ",tau_thin
         write(*,*) "Tthin = ",Tthin
         write(*,*) "Tthick = ",Tthick
         error stop
      end if

      end subroutine
