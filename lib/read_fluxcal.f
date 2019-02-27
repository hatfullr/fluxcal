      subroutine read_fluxcal(nnit)
c     This program will read the standard input format expected from the
c     user. The format is as follows:
c     
c     Header (all in one line):
c     t                 !time (real*8)
c     Body:
c     x(i)              !x position (real*8)
c     y(i)              !y position (real*8)
c     z(i)              !z position (real*8)
c     am(i)             !mass (real*8)
c     hp(i)             !smoothing length (real*8)
c     rho(i)            !density (real*8)
c     vx(i)             !x velocity (real*8)
c     vy(i)             !y velocity (real*8)
c     vz(i)             !z velocity (real*8)
c     a(i)              !specific internal energy (real*8)
c     wmeanmolecular(i) !mean molecular weight (real*8)
c     g(i)              !local g magnitude (real*8)
c
c     Note that this should be transposed, so that each row above is
c     really a column. For example, "x(i),y(i),z(i),am(i),..."
      include 'flux_cal.h'
      
c      character*255 fname
      logical fileexists
      real*8 opacit
      integer i, do_debug
      real*8 moutsideTEOS,eoutsideTEOS,mintempoutsideTEOS
      real*8 maxtempoutsideTEOS,minrhooutsideTEOS,maxrhooutsideTEOS
      
      if(nnit.le.9999) then
         write(infname,"('fluxcal_',i4.4,'.dat')") nnit
      else if(nnit.le.99999) then
         write(infname,"('fluxcal_',i5.5,'.dat')") nnit
      else
         write(infname,"('fluxcal_',i6.6,'.dat')") nnit
      end if

      inquire(file=trim(infname),exist=fileexists)
      if(.not.fileexists) then
         write(*,*) "Could not find the next output file"
         write(*,*) "Next file to read is '",trim(adjustl(infname)),"'"
         error stop "read_fluxcal.f"
      end if
      write(*,*) "Reading '",trim(adjustl(infname)),"'"
      
      open(12,file=trim(adjustl(infname)),form='unformatted')
      read(12) t
      t=t*tunit
      moutsideTEOS = 0.d0
      eoutsideTEOS = 0.d0
      mintempoutsideTEOS = 1.d30
      maxtempoutsideTEOS = -1.d30
      minrhooutsideTEOS = 1.d30
      maxrhooutsideTEOS = -1.d30
      do i=1,nmax
         read(12,end=200) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $        vx(i),vy(i),vz(i),a(i),wmeanmolecular(i),
     $        localg(i)         !,tempp(i),pp(i)
         x(i) = x(i)*runit
         y(i) = y(i)*runit
         z(i) = z(i)*runit
         am(i) = am(i)*munit
         hp(i) = hp(i)*runit
         rho(i) = rho(i)*rhounit
         vx(i) = vx(i)*vunit
         vy(i) = vy(i)*vunit
         vz(i) = vz(i)*vunit
         a(i) = a(i)*Eunit/munit
         wmeanmolecular(i)=wmeanmolecular(i)*muunit
         localg(i)=localg(i)*gunit
         if(a(i).gt.0.d0) then ! Core particle
            tempp(i) = useeostable(a(i),rho(i),1) * kelvin
            if(rhooutsideTEOS.gt.0.d0) then
               mintempoutsideTEOS = min(mintempoutsideTEOS,tempp(i))
               maxtempoutsideTEOS = max(maxtempoutsideTEOS,tempp(i))
               minrhooutsideTEOS = min(minrhooutsideTEOS,rho(i))
               maxrhooutsideTEOS = max(maxrhooutsideTEOS,rho(i))
               moutsideTEOS = moutsideTEOS + am(i)
               eoutsideTEOS = aoutsideTEOS + a(i)*am(i)
            end if
            pp(i) = useeostable(a(i),rho(i),3) * gram/cm/sec**2.d0
         else
            tempp(i) = 0.d0
            pp(i) = 0.d0
         end if

         call getOpacitySub(x(i),y(i),z(i),tempp(i),
     $        rho(i),0.d0,Rform,opacit)
         opac_sph(i)=opacit
c         if(tempp(i).le.8000) write(*,*) "OPACITY", tempp(i), opacit
         tauA(i) = taucoef*am(i)*opacit/hp(i)**2.d0

c     If the user wants to use get_teff, and particle i is optically
c     thick, store get_teff for use in getTpractical.f
         Teff(i) = 0.d0
         
         if(envfit .and. (tauA(i).gt.tau_thick_envfit)) then
c            if((table_nabla_Tmin.lt.tempp(i)).and.
c     $           (tempp(i).lt.table_nabla_Tmax)) then
c               if((table_nabla_gmin.lt.localg(i)).and.
c     $              (localg(i).lt.table_nabla_gmax)) then
c     Calculate Teff iff tempp(i) and localg(i) fall within the tabulated values
            if((3.51d0.le.log10(tempp(i))).and.
     $           (log10(tempp(i)).le.6.d0)) then
               if((-0.1d0.le.log10(localg(i))).and.
     $              (log10(localg(i)).le.4.2d0)) then
                  do_debug = 0
c                 if(i.eq.172258) do_debug=1
                  Teff(i) = get_teff(pp(i),tempp(i),localg(i),do_debug)
               end if
            end if
         end if
      end do
 200  close(12)

      if(i.eq.nmax) then
         write(*,*) "Too many particles in simulation"
         write(*,*) "Maximum number of particles allowed is ",nmax
         error stop "read_fluxcal.f"
      end if

      if((moutsideTEOS .gt. 0.d0).or.(eoutsideTEOS .gt. 0.d0)) then
         write(*,*) "WARNING: Some particles are outside the TEOS."//
     $        " Their cumulative properties are:"
         write(*,*) "         Mass = ",moutsideTEOS/munit_out
         write(*,*) "         Internal energy = ",eoutsideTEOS/eunit_out
         write(*,*) "         min(rho,T) = (",
     $        minrhooutsideTEOS/rhounit_out,",",
     $        mintempoutsideTEOS/tempunit_out,")"
         write(*,*) "         max(rho,T) = (",
     $        maxrhooutsideTEOS/rhounit_out,",",
     $        maxtempoutsideTEOS/tempunit_out,")"

      end if

      n=i-1
      
      end subroutine
