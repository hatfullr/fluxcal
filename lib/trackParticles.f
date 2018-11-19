      subroutine trackParticles
      include 'flux_cal.h'

      character*255 trackingfile
      integer m, i

      if(innit.le.9999) then
         write(trackingfile,"('fluxcal_',i4.4,'.track')") innit
      else if(innit.le.99999) then
         write(trackingfile,"('fluxcal_',i5.5,'.track')") innit
      else
         write(trackingfile,"('fluxcal_',i6.6,'.track')") innit
      end if

      if(binary_tracking_file) then
         open(44,file=trackingfile,form='unformatted',status='unknown')
         write(*,*) "Writing to binary file '",trim(trackingfile),"'"
         write(44) t*tunit
         do i=1,nmax
            if(track(i).gt.0) then
               m=track(i)
               write(44) x(m),y(m),z(m),am(m),hp(m),rho(m),vx(m),vy(m),
     $              vz(m),a(m),wmeanmolecular(m),localg(m),tempp(m),
     $              pp(m),track(i)
            end if
         enddo
         close(44)
      else
 104     format(14E15.7,i15)
         open(44,file=trackingfile,status='unknown')
         write(*,*) "Writing to text file '",trim(trackingfile),"'"
         write(44,"(E15.7)") t*tunit
         do i=1,nmax
            if(track(i).gt.0) then
               m=track(i)
               write(44,104) x(m),y(m),z(m),am(m),hp(m),rho(m),vx(m),
     $              vy(m),vz(m),a(m),wmeanmolecular(m),localg(m),
     $              tempp(m),pp(m),track(i)
            end if
         enddo
         close(44)
      end if
      
      end subroutine
