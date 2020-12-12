      subroutine init_grid
      include 'optical_depth.h'

      logical fileexists
      integer i,j

      if(isInitGrid) return
      
      write(o,*) "Establishing the grid"
      prepareIntegrationCalled = .false.
      
c     Initialize the file names for this iteration
 33   format ('dimen',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 34   format ('dimen',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 35   format ('teffs',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 36   format ('teffs',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
      if(innit.le.99999) then
         write (fname2,33) innit,nint(anglex*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglez*180.d0/pi)
         write (fname,35) innit,nint(anglex*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglez*180.d0/pi)
      else
         write (fname2,34) innit,nint(anglex*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglez*180.d0/pi)
         write (fname,36) innit,nint(anglex*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglez*180.d0/pi)
      end if

c     Use a 3x3 grid for starters:
      NXMAPnext=min_Nx          ! default = 3
      NYMAPnext=min_Ny          ! default = 3
      
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
      xminmapnext = xmin
      xmaxmapnext = xmax
      yminmapnext = ymin
      ymaxmapnext = ymax
      hxmap = (xmaxmap-xminmap)/dble(nxmap-1)
      hymap = (ymaxmap-yminmap)/dble(nymap-1)

c      xminmapnext=1d30
c      xmaxmapnext=-1d30
c      yminmapnext=1d30
c      ymaxmapnext=-1d30
c      do IP=1,N
c         if(a(IP).gt.0.d0)then
c            xminmapnext=min(xminmapnext,x(IP)+2*HP(IP))
c            xmaxmapnext=max(xmaxmapnext,x(IP)-2*HP(IP))
c            yminmapnext=min(yminmapnext,y(IP)+2*HP(IP))
c            ymaxmapnext=max(ymaxmapnext,y(IP)-2*HP(IP))
c         endif
c      enddo

      lastTOTALLUM=-1d30
      lastTOTALflux = 0.d0

      
      inquire(file=trim(adjustl(fname2)),exist=dimenFileAlreadyExists)

      if(dimenFileAlreadyExists) then ! Dimen file found
         write(o,*) "Reading pre-existing dimen file '" //
     $        trim(adjustl(fname2)) // "'"
         open (73,file=trim(adjustl(fname2)))
         read(73,*) xminmap,hxmap,nxmap,yminmap,hymap,nymap
         close(73)

         xminmap = xminmap
         hxmap = hxmap
         yminmap = yminmap
         hymap = hymap

         ! Use these to set hxmap and hymap
         xmaxmap=xminmap+hxmap*(nxmap-1)
         ymaxmap=yminmap+hymap*(nymap-1)
         if(nxmap.gt.max_Nx) then ! Too many cells on x axis
            nxmap=max_Nx
            hxmap=(xmaxmap-xminmap)/dble(nxmap-1)
            print *,'Reset nxmap=',nxmap,' hxmap=',hxmap
         endif
         if(nymap.gt.max_Ny) then ! Too many cells on y axis
            nymap=max_Ny
            hymap=(ymaxmap-yminmap)/dble(nymap-1)
            print *,'Reset nymap=',nymap,' hymap=',hymap
         endif

         ! Setup the 3D integrating grid
         call prepareIntegration
         goto 1000
      else                      ! No dimen file found
         call createGrid
         call writeDimenFile
      end if


 1000 isInitGrid=.true.

 100  format(A25," = ",ES22.14,", ",ES22.14)
 101  format(A25," = ",I22,", ",I22)
      write(o,*) ""
      write(o,*) "Grid details:"
      write(o,100) "xminmap,xmaxmap      ",xminmap/runit_out,
     $     xmaxmap/runit_out
      write(o,100) "yminmap,ymaxmap      ",yminmap/runit_out,
     $     ymaxmap/runit_out
      write(o,100) "hxmap,hymap          ",hxmap/runit_out,
     $     hymap/runit_out
      write(o,101) "nxmap,nymap          ",nxmap,nymap
c      write(o,100) "minstpsize,maxstpsize",min_step_size,max_step_size
c      write(o,101) "minNstp,maxNstp      ",min_steps_taken,
c     $     max_steps_taken
      write(o,*) ""
      
c      write(o,*) ""
c      write(o,*) "Grid details (cgs):"
c      write(o,*) "   xminmap,xmaxmap = ",xminmap,xmaxmap
c      write(o,*) "   yminmap,ymaxmap = ",yminmap,ymaxmap
c      write(o,*) "   hxmap,hymap     = ",hxmap,hymap
c      write(o,*) "   nxmap,nymap     = ",nxmap,nymap
c      write(o,*) ""
      
      end subroutine
