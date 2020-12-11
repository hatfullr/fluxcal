
c      use iso_fortran_env
      include 'lib/flux_cal.h'
      real*8 TOTALpracticalLUM,avgrpos,avgxhp
      common/analysis2/ avgrpos,avgxhp
      common/foroutput/ TOTALpracticalLUM
      integer particle_info_file,holder,int_size
      character*255 papfname
      integer pid(nmax)
      integer nid,ip
      logical fileexists
      integer i, imax
      character*255 newoutfile
      real*8 dummy

      integer writefile         ! Remove later
      common/writefilee/writefile ! Remove later

      call init                 ! Initialize

c      writefile = 11            ! Remove later
c      open(writefile,file="opacity_calls.dat",status="unknown") ! Remove later
      
      imax = 999
      
 230  format(3A22)
 231  format(3I22)
 232  format(A20,3ES22.14)
 233  format(A20,3A22)

 234  format(A,"_",I3.3)
      
      write(o,*) "** Main loop **"
      inquire(file=trim(adjustl(outfile)),exist=fileexists)
      if(fileexists) then
         write(o,*) "File '",trim(adjustl(outfile)),"' "//
     $        "already exists"
         do i=1,imax
            write(newoutfile,234) trim(adjustl(outfile)),i
            inquire(file=trim(adjustl(newoutfile)),exist=fileexists)
            if(.not.fileexists) then
               exit
            end if
         end do
         outfile=newoutfile
      end if

      write(o,*) 
      
      write(o,*) "Creating '",trim(adjustl(outfile)),"'"
      open(89,file=trim(adjustl(outfile)),status='replace')
      write(89,230) "start","finish","step"
      write(89,231) start,finish,step
      write(89,*) ""
      write(89,233) "filename","Time","Liso","<T>"
      write(89,233) "-------------------","--------------",
     $     "--------------","--------------"
      close(89)

      
      DO innit=start,finish,step
         call read_fluxcal(innit)
         call setViewingAngle
         isInitGrid = .false.
         write(o,*) ""

         ! Remove if statement when finished testing
c         if(.not.get_true_luminosity) then
c            call init_grid
c         end if

c        ****************************************************************
         if(get_fluxes) then
            ! write(o,*) "Finding the fluxes at each grid point"
            call init_grid

            if(dimenFileAlreadyExists) then
               call integrateTau
               call getFlux
            end if
            
            call writeTempsFile
            write(o,*) ""
            write(o,'(A,ES22.14)')' Luminosity =',
     $           TOTALpracticalLUM/Lunit_out
!write(o,'(A,ES22.14)')' Maximum Teff=',TMAX/tempunit_out
            if(avgt.eq.0) then
               write(o,'(A,ES22.14)')' Average Teff =',0.d0
            else
               write(o,'(A,ES22.14)')' Average Teff =',avgt/numcell/
     $              tempunit_out
            end if
            write(o,'(A,ES22.14)') ' Effective visible area =',
     $           Atot/runit_out**2.d0
            write(o,*) ""
            write(o,*) "Writing to '",trim(adjustl(outfile)),"'"
            open(89,file=trim(adjustl(outfile)),action='write',
     $           position='append')
            write(89,232) trim(adjustl(infname)), t/tunit_out,
     $           TOTALpracticalLUM/Lunit_out,avgt/numcell/tempunit_out
            close(89)
            write(o,*) ""
         end if
c        ****************************************************************




         
c        ****************************************************************
         if(get_true_luminosity) then
            call integrateSpherical
         end if
c        ****************************************************************




         
c        ****************************************************************
         if(get_closest_particles) then
            call init_grid
            write(o,*) "Finding information on the closest particles"
            call getClosestParticles
         end if
c        ****************************************************************




         
c        ****************************************************************
         if(get_particles_at_pos) then
            write(o,*) "Finding all particles at posx, posy"
            call particlesAtPos(posx,posy,pid,nid)
         end if
c        ****************************************************************



         

c        ****************************************************************
         if(get_integration_at_pos) then
            write(o,*) "Finding data at integration steps at posx,posy"
            call integrationAtPos
         end if
c        ****************************************************************



         

c        ****************************************************************
         if(get_integration_at_all_pos) then
            call init_grid
            write(o,*) "Finding data at integration steps at posx,posy"
            call integrationAtAllPos
         end if
c        ****************************************************************




         
c        ****************************************************************
         if(track_particles) then
            write(o,*) "Tracking particles"
            call trackParticles
         end if
c        ****************************************************************



         
         
c        ****************************************************************
         if(get_info_of_particle) then
            write(o,*) "Finding data of particle ",info_particle
            call output(pinfo_file,info_particle,.false.)
         end if
c        ****************************************************************


         

         write(o,*) ""
         write(o,*) "**************************************************"
         write(o,*) ""
      enddo
      write(o,*) "** Complete **"
      write(o,*) ""
 500  format(A20," ",f15.7," ",A4)

c      close(writefile)          ! Remove later
      
      end program
