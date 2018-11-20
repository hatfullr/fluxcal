      include 'lib/flux_cal.h'
      real*8 TOTALpracticalLUM,avgrpos,avgxhp
      common/analysis2/ avgrpos,avgxhp
      common/foroutput/ TOTALpracticalLUM
      integer particle_info_file,holder,int_size
      character*255 pinfo_file,info_particle_string,papfname
      integer pid(nmax)
      integer nid,ip
      real*8 opacit

      call init                 ! Initialize

 230  format(3A15)
 231  format(3I15)
 232  format(A20,3E15.7)
 233  format(A20,3A15)
      
      
      write(*,*) "** Main loop **"
      write(*,*) "Clearing '",trim(adjustl(outfile)),"'"
      open(89,file=trim(adjustl(outfile)),status='replace')
      write(89,230) "start","finish","step"
      write(89,231) start,finish,step
      write(89,*) ""
      write(89,233) "filename","Time [days]","L [Lsun]","<T^4> [K^4]"
      write(89,233) "-------------------","--------------",
     $     "--------------","--------------"
      close(89)


      if(get_info_of_particle) then
c        Determine the character size of the integer
         holder = info_particle
         int_size = 0
         do while ( holder > 0 )
            holder = holder/10
            int_size = int_size + 1
         end do
         
         write(info_particle_string,*) info_particle
         pinfo_file = trim("particle_" //
     $        trim(adjustl(info_particle_string))) // ".dat"

         write(*,*) "Clearing '", trim(pinfo_file),"'"
         open(33,file=pinfo_file,status='replace')
 101     format(18A15)
         write(33,101) "t [s]","x [cm]","y [cm]",
     $        "z [cm]","am [g]","hp [cm]","rho [g/cm^3]","a [erg/g]",
     $        "mu*mH [g]","g [cm/s^2]","T [K]","P [g/cm/s^2]",
     $        "s [erg/K]","kappa [cm^2/g]","Teff [K]","tau",
     $        "A_vis [cm^2]","L_vis [erg/s]"
         close(33)
      end if
      
      DO innit=start,finish,step
         call read_fluxcal(innit)
         call setViewingAngle
         write(*,*) ""


c        ****************************************************************
         if(get_fluxes) then
            write(*,*) "Finding the fluxes at each grid point"
            call optical_depth
            write(*,*) "Writing to output file"
            open(89,file=trim(adjustl(outfile)),action='write',
     $           position='append')
            write(89,232) trim(adjustl(infname)), t,
     $           TOTALpracticalLUM,
     $           real(avgt)/real(numcell)
            close(89)
         end if
c        ****************************************************************



c        ****************************************************************
         if(get_closest_particles) then
            write(*,*) "Finding information on the closest particles"
            call getClosestParticles
         end if
c        ****************************************************************


         
c        ****************************************************************
         if(get_particles_at_pos) then
            write(*,*) "Finding all particles at posx, posy"
            call particlesAtPos(posx,posy,pid,nid)
         end if
c        ****************************************************************



         

c        ****************************************************************
         if(get_integration_at_pos) then
            write(*,*) "Finding data at integration steps for posx,posy"
            call integrationAtPos
         end if
c        ****************************************************************





c        ****************************************************************
         if(track_particles) then
            write(*,*) "Tracking particles"
            call trackParticles
         end if
c        ****************************************************************



         
         
c        ****************************************************************
         if(get_info_of_particle) then
            write(*,*) "Finding data of specific particle"
            open(33,file=trim(adjustl(pinfo_file)),action='write',
     $           position='append')
            call getParticleInfo(info_particle,33)
            close(33)
         end if
c        ****************************************************************


         
c         if(get_fluxes) then
c            write(*,*) "Finding the fluxes at each grid point"
c            call optical_depth
c            write(*,*) "Writing to output file"
c            open(89,file=trim(adjustl(outfile)),action='write',
c     $           position='append')
c            write(89,232) trim(adjustl(infname)), t*tunit,
c     $           TOTALpracticalLUM,
c     $           (avgt4/numcell)**0.25d0
c            close(89)
c         end if
         write(*,*) ""
         write(*,*) "**************************************************"
         write(*,*) ""
      enddo
c      close(89)
      write(*,*) "** Complete **"
      write(*,*) ""
 500  format(A20," ",f15.7," ",A4)

      end program
