      subroutine getTpractical(z0,z1,zthick,thick_p,myh1,
     $     Tpractical, tau_thin, sphoto)
c     Start at z1 and integrate to max(z0,zthick) to find Teff
c     Input: z0 = z location of last visible gas on LOS
c            z1 = z location of first visible gas on LOS
c            zthick = z location of optically thick particle's surface
c            thick_p = ID of optically thick particle
c            h1 = Not sure
c     Output: Tthin = Temperature component of the flux Tthick, tau_thin, sphoto
c
c     Tthick = Surface temperature of an optically thick particle.
c     Tthin  = Temperature calculated from the outgoing flux from an
c              optically thick region made by many optically thin
c              particles.
c     Tthick > 0 iff an optically thick particle is along the LOS
c     Tthin > 0 iff optically thin particles are along the LOS
c
c     tau_thin = Attenuation factor, as in e^(-tau_thin)
c     sphoto   = z-location of the 
      include 'optical_depth.h'
      real*8 zstop,z0,z1,zthick,tau_thin,myh1,Tthin4,Tthick4
      integer i,thick_p
      integer outerid
      real*8 rn,mymaxdz,maxz,rv,rc,deltar2
      real*8 Tpractical

      Tthick4 = 0.d0
      Tthin4 = 0.d0
      sphoto = 0.d0
      tau_thin = 0.d0
      Tpractical = 0.d0
      
      zstop = max(zthick,z0)
      
c     Integrate down to the optically thick particle or until
c     accumulated tau is equal to tau_thick_integrator.
      ! Initialize the integrator
      nrhs=0
c     We might want this later
c      do ifilter=1,numfilters+4
c         taustart(ifilter)=0.d0
c      enddo
      do i=1,size(taustart)
         taustart(i)=0.d0
      end do

cc     Find the particle kernel closest to xpos,ypos,z1
c      maxz = -1.d30
c      outerid = 0
c      do i=1,n
c         mymaxdz = (4.d0*hp(i)**2.d0 - (x(i)-xpos)**2.d0
c     $        - (y(i)-ypos)**2.d0)**0.5d0
c         if(mymaxdz+z(i) .gt. maxz) then
c            maxz = mymaxdz+z(i)
c            outerid = i
c         end if
c      end do
c
c      deltar2 = (x(outerid)-xpos)**2.d0 + (y(outerid)-ypos)**2.d0
c     $     + (z(outerid)+2.d0*hp(outerid) - maxz)**2.d0
c
cc     Done by law of cosines
c      mu = 1.d0 - 0.125d0*deltar2/(hp(outerid)**2.d0)
c
c
c      if(mu .eq. 0.d0) then
c         mu = 1.d-30
c      end if

      if(z1.ne.zstop)then
         call odeint(taustart,4+numfilters,z1,zstop,
     $        eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
     $        rkqs)
         
         Tthin4=taustart(4)     ! int_0^tau T^4 e^(-tau') dtau'
         tau_thin = taustart(1) ! Attenuation factor

         if(kount1.gt.1) then   ! Integrator didn't reach an optically thick particle
            sa=s(kount1-1)      ! z-position w/ optical depth tau_tot<tau_thick_integrator
            sb=s(kount1)        ! z-position w/ optical depth tau_tot>=tau_thick_integrator
         
            ! Interpolate to find position where optical tau_tot=tau_thick_integrator:
            sphoto=(sa*(tau(1,kount1)-1.d0)
     $           +sb*(1.d0-tau(1,kount1-1)))/
     $           (tau(1,kount1)-tau(1,kount1-1))
         else if (thick_p.gt.0) then ! Integrator reached an optically thick particle
            sphoto=zthick
            Tthick4 = Teff(thick_p)**4.d0 * exp(-tau_thin)
         end if
      else                      ! Integrator started and stopped at an optically thick particle
         Tthick4 = Teff(thick_p)**4.d0 ! No attenuation
         Tthin4 = 0.d0
      end if

         
      Tpractical = (Tthick4 + Tthin4)**0.25d0
      
cc     zstop = max(z0,zthick)
cc     write(*,*) "zthick,z0,z1 = ",zthick,z0,z1
c      if(thick_p.gt.0) then     ! There is a thick particle along the LOS
c         ! We need to try and integrate down to the surface of the thick
c         ! particle to check if thin particles create a thick region, and
c         ! also to attain tau_thin, the attenuation factor.
c         
cc     if(z1.ne.zstop) then ! If there is optically thin material
c         ! Initialize the integrator
c         nrhs=0
c         do ifilter=1,numfilters+4
c            taustart(ifilter)=0.d0
c         enddo
cc     call odeint(taustart,4+numfilters,z1,zstop,
cc     $           eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
cc     $           rkqs)
c         
c         ! Integrate with bounds set to the beginning to the surface of
c         ! the thick particle.
c         call odeint(taustart,4+numfilters,z1,zthick,
c     $        eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
c     $        rkqs)
c         Tthin=taustart(4)**0.25d0 ! int_0^tau T^4 e^(-tau') dtau'
c         tau_thin = taustart(1) ! Attenuation factor
c         nphoto = taustart(3)
c         if(kount1.gt.1) then   ! Integrator found a photosphere
c            sa=s(kount1-1)      ! z-position w/ optical depth tau_tot<tau_thick_integrator
c            sb=s(kount1)        ! z-position w/ optical depth tau_tot>=tau_thick_integrator
c            
c            ! Interpolate to find position where optical tau_tot=tau_thick_integrator:
c            sphoto=(sa*(tau(1,kount1)-1.d0)
c     $           +sb*(1.d0-tau(1,kount1-1)))/
c     $           (tau(1,kount1)-tau(1,kount1-1))
c            
c            call getLocalQuantities(xpos,ypos,sphoto) ! Get the temperature
c            Tthin = t6*1.d6
c         else                   ! Integrator did not find a photosphere
c            sphoto = zthick     ! Photosphere must be at surface of thick particle
c            if(Teff(thick_p).lt.0.and.Teff(thick_p).ge.-15000) then
c               write(*,*) "a really hot particle is outside"
c               stop
c            end if
c            
cc            if(Teff(thick_p).ge.7000) then
cc               write(*,*) "a moderately hot particle is outside:",
cc           &              tempp(thick_p),Teff(thick_p), thick_p
cc            end if
c            
c            
c            if(Teff(thick_p).lt.-15000.and.
c     $           Teff(thick_p).ge.-25000) then
cc               write(*,*) "a cold particle is outside",
cc           &              tempp(thick_p), rho(thick_p),
cc           &              tauA(thick_p), opac_sph(thick_p),
cc           &              hp(thick_p)*2./6.06e10,thick_p
c               Tthick=tempp(thick_p)
c            end if
c            
c            if(Teff(thick_p).gt.0) then
c               Tthick = Teff(thick_p) ! Teff is calculated in lib/read_fluxcal.f
c            end if
c            
c            ! Use the smallest temperature possible, because the SPH methods will be more accurate
c            ! at low T than the envelope fitting routine
c            if(tempp(thick_p) .lt. Tthick) then
c               write(*,*) " strange place. Why I am here?",
c     &              tempp(thick_p),Teff(thick_p), Tthick, thick_p
c               Tthick = tempp(thick_p)
c               stop
c            end if
c         end if
c         
c      else if((thick_p.eq.0).or.  ! There is no optically thick particle along the LOS OR
c     $        (.not.envfit)) then ! The user wants to use ONLY ray-trace
c         ! Initialize the integrator
c         nrhs=0
c         do ifilter=1,numfilters+4
c            taustart(ifilter)=0.d0
c         enddo
c
c         ! Integrate with bounds set across all fluid
c         call odeint(taustart,4+numfilters,z1,z0,
c     $        eps,0.25d0*myh1,myh1,nok,nbad,derivs2,
c     $        rkqs)
c         tau_thin = taustart(1) ! Attenuation factor
c         nphoto = taustart(3)
c         if(kount1.gt.1) then   ! Photosphere was found
c            sa=s(kount1-1)      ! z-position w/ optical depth tau_tot<tau_thick_integrator
c            sb=s(kount1)        ! z-position w/ optical depth tau_tot>=tau_thick_integrator
c            
c            ! Interpolate to find position where optical tau_tot=tau_thick_integrator:
c            sphoto=(sa*(tau(1,kount1)-1.d0)
c     $           +sb*(1.d0-tau(1,kount1-1)))/
c     $           (tau(1,kount1)-tau(1,kount1-1))
c            
c            na=tau(3,kount1-1)
c            nb=tau(3,kount1)
c            
c            ! Interpolate to find number of smoothing lengths in at tau_tot=tau_thick_integrator:
c            nphoto=(na*(tau(1,kount1)-1.d0)
c     $           +nb*(1.d0-tau(1,kount1-1)))/
c     $           (tau(1,kount1)-tau(1,kount1-1))
c            
c            call getLocalQuantities(xpos,ypos,sphoto) ! Get the local temperature
c            Tthin=t6*1.e6
c         else
c            sphoto = -1.d30     ! No photosphere here
c            nphoto = 0          ! No particles traveled through
c         end if
c         
c      end if
c      
      if(tau_thin .lt. 0.d0) then
         write(*,*) "tauthin < 0! Something went wrong in"//
     $        " getTpractical.f"
         write(*,*) "tau_thin = ",tau_thin
         write(*,*) "Tthin4 = ",Tthin4
         write(*,*) "Tthick4 = ",Tthick4
         error stop "getTpractical.f"
      end if
      
      end subroutine
      
