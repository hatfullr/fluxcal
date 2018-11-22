      subroutine createGrid
c     L232 to L292
      include 'optical_depth.h'

c     Find boundaries of the fluid:
      xmin=1d30
      xmax=-1d30
      ymin=1d30
      ymax=-1d30
      do IP=1,N
         if(a(IP).gt.0.d0)then
            xmin=min(xmin,x(IP)-2*HP(IP))
            xmax=max(xmax,x(IP)+2*HP(IP))
            ymin=min(ymin,y(IP)-2*HP(IP))
            ymax=max(ymax,y(IP)+2*HP(IP))
         endif
      enddo
c      print *,'xmin,xmax=',real(xmin),real(xmax)
c      print *,'ymin,ymax=',real(ymin),real(ymax)


c     Set the initial boundaries of the grid near the outer edge of the fluid:
      xminmapnext=1d30
      xmaxmapnext=-1d30
      yminmapnext=1d30
      ymaxmapnext=-1d30
      do IP=1,N
         if(a(IP).gt.0.d0)then
            xminmapnext=min(xminmapnext,x(IP)+2*HP(IP))
            xmaxmapnext=max(xmaxmapnext,x(IP)-2*HP(IP))
            yminmapnext=min(yminmapnext,y(IP)+2*HP(IP))
            ymaxmapnext=max(ymaxmapnext,y(IP)-2*HP(IP))
         endif
      enddo

      lastTOTALLUM=-1d30

c     Use a 3x3 grid for starters:
      NXMAPnext=3
      NYMAPnext=3

c     Initialize the file names for this iteration
 33   format ('dimen',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 34   format ('dimen',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 35   format ('fluxes',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 36   format ('fluxes',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
      if(innit.le.99999) then
         write (fname2,33) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
         write (fname,35) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      else
         write (fname2,34) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
         write (fname,36) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      end if

      end subroutine
