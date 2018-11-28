      include 'lib/flux_cal.h'
      real*8 TOTALpracticalLUM,avgrpos,avgxhp
      common/analysis2/ avgrpos,avgxhp
      common/foroutput/ TOTALpracticalLUM
      integer particle_info_file,holder,int_size
      character*255 papfname
      integer pid(nmax)
      integer nid,ip
      real*8 opacit

      call init                 ! Initialize

 230  format(3A22)
 231  format(3I22)
 232  format(A20,3ES22.14)
 233  format(A20,3A22)
      
      
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


      DO innit=start,finish,step
         call read_fluxcal(innit)
         call setViewingAngle
         write(*,*) ""

         call init_grid


c        ****************************************************************
         if(get_fluxes) then
            ! write(*,*) "Finding the fluxes at each grid point"
!     call optical_depth
            call integrateTau
            call getFlux
            call writeTempsFile
            write(*,'(A,ES22.14)')' Total luminosity in enclosed area=',
     $           TOTALpracticalLUM/Lunit_out
            !write(*,'(A,ES22.14)')' Maximum Teff=',TMAX/tempunit_out
            write(*,'(A,ES22.14)')' Average Teff=',avgt/numcell/
     $           tempunit_out
            write(*,*) "Writing to '",trim(adjustl(outfile)),"'"
            open(89,file=trim(adjustl(outfile)),action='write',
     $           position='append')
            write(89,232) trim(adjustl(infname)), t,
     $           TOTALpracticalLUM,avgt/numcell
            close(89)
            write(*,*) ""
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
            write(*,*) "Finding data at integration steps at posx,posy"
            call integrationAtPos
         end if
c        ****************************************************************



         

c        ****************************************************************
         if(get_integration_at_all_pos) then
            write(*,*) "Finding data at integration steps at posx,posy"
            call integrationAtAllPos
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
            write(*,*) "Finding data of particle ",info_particle
            call output(pinfo_file,info_particle)
         end if
c        ****************************************************************


         

         write(*,*) ""
         write(*,*) "**************************************************"
         write(*,*) ""
      enddo
      write(*,*) "** Complete **"
      write(*,*) ""
 500  format(A20," ",f15.7," ",A4)

      end program
