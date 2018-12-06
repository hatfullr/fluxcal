      subroutine init_grid
      include 'optical_depth.h'

      logical fileexists
      integer i,j

      write(*,*) "Establishing the grid"
      
c     Initialize the file names for this iteration
 33   format ('dimen',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 34   format ('dimen',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 35   format ('teffs',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 36   format ('teffs',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
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
      
      inquire(file=trim(adjustl(fname2)),exist=fileexists)

      if(fileexists) then ! Dimen file found
         write(*,*) "Reading pre-existing dimen file '" //
     $        trim(adjustl(fname2)) // "'"
         open (73,file=trim(adjustl(fname2)))
         read(73,*) xminmap,hxmap,nxmap,yminmap,hymap,nymap
         close(73)

         ! Use these to set hxmap and hymap
         xmaxmap=xminmap+hxmap*(nxmap-1)
         ymaxmap=yminmap+hymap*(nymap-1)
         if(nxmap.gt.nxmapmax) then ! Too many cells on x axis
            nxmap=nxmapmax
            hxmap=(xmaxmap-xminmap)/(nxmap-1)
            print *,'Reset nxmap=',nxmap,' hxmap=',hxmap
         endif
         if(nymap.gt.nymapmax) then ! Too many cells on y axis
            nymap=nymapmax
            hymap=(ymaxmap-yminmap)/(nymap-1)
            print *,'Reset nymap=',nymap,' hymap=',hymap
         endif

c         ! Setup the 3D integrating grid
         call prepareIntegration
         return
      else ! No dimen file found
         call createGrid
         call writeDimenFile
      end if
      end subroutine
