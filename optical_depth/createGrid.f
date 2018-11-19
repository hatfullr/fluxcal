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

      IF(innit.LE.99999) THEN
         WRITE (FNAME,103) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
         WRITE (FNAME2,193) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      ELSE
         WRITE (FNAME,104) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
         WRITE (FNAME2,194) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      endif
c      WRITE (FNAME3,394) innit,nint(anglez*180.d0/pi),
c     $     nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      WRITE (FNAME4,494) innit,nint(anglez*180.d0/pi),
     $     nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
 103  FORMAT ('fluxes',I5.5,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')
 104  FORMAT ('fluxes',I6.6,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')
 193  FORMAT ('dimen',I5.5,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')
 194  FORMAT ('dimen',I6.6,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')
c 394  FORMAT ('spectrum',I6.6,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')
 494  FORMAT ('trangevsrho',I6.6,'_',I3.3,'_',I3.3,'_',I3.3,'.sph')

      end subroutine
