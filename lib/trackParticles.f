      subroutine trackParticles
      include 'flux_cal.h'

      character*255 trackingfile
      integer m, i
      logical fileexists

      do i=1,nmax            ! Initialize the track array
         track(i) = 0
      end do

      if(.not.track_all) then   ! Use the trackfile
         write(*,*) "Reading trackfile"
         open(17,file=trim(adjustl(trackfile)),status='old')
         do i=1,nmax
            read(17,*,end=120) track(i)
         enddo
 120     close(17)
      else                      ! Track all the particles
         write(*,*) "Tracking all particles"
         do i=1,n
            track(i) = i
         end do
      end if
      
      if(innit.le.9999) then
         write(trackingfile,"('fluxcal_',i4.4,'.track')") innit
      else if(innit.le.99999) then
         write(trackingfile,"('fluxcal_',i5.5,'.track')") innit
      else
         write(trackingfile,"('fluxcal_',i6.6,'.track')") innit
      end if

      call makeOutputFile(trackingfile)
      write(*,*) "Writing to '",trim(adjustl(trackingfile)),"'"
      do i=1,n
         if(track(i).gt.0) then
            call output(trackingfile,i,binary_tracking_file)
         end if
      end do
      
c      if(binary_tracking_file) then
c         open(44,file=trackingfile,form='unformatted',status='unknown')
c         write(*,*) "Writing to binary file '",trim(trackingfile),"'"
c         write(44) t*tunit
c         do i=1,nmax
c            if(track(i).gt.0) then
c               m=track(i)
c               write(44) x(m),y(m),z(m),am(m),hp(m),rho(m),vx(m),vy(m),
c     $              vz(m),a(m),wmeanmolecular(m),localg(m),tempp(m),
c     $              pp(m),track(i)
c            end if
c         enddo
c         close(44)
c      else
c 104     format(14ES22.14,i22)
c         open(44,file=trackingfile,status='unknown')
c         write(*,*) "Writing to text file '",trim(trackingfile),"'"
c         write(44,"(ES22.14)") t*tunit
c         do i=1,nmax
c            if(track(i).gt.0) then
c               m=track(i)
c               write(44,104) x(m),y(m),z(m),am(m),hp(m),rho(m),vx(m),
c     $              vy(m),vz(m),a(m),wmeanmolecular(m),localg(m),
c     $              tempp(m),pp(m),track(i)
c            end if
c         enddo
c         close(44)
c      end if
      
      end subroutine
