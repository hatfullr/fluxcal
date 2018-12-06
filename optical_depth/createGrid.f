      subroutine createGrid
      include 'optical_depth.h'

      integer counter
      integer i,j
      
cc     Find boundaries of the fluid:
c      xmin=1d30
c      xmax=-1d30
c      ymin=1d30
c      ymax=-1d30
c      do IP=1,N
c         if(a(IP).gt.0.d0)then
c            xmin=min(xmin,x(IP)-2*HP(IP))
c            xmax=max(xmax,x(IP)+2*HP(IP))
c            ymin=min(ymin,y(IP)-2*HP(IP))
c            ymax=max(ymax,y(IP)+2*HP(IP))
c         endif
c      enddo
cc      print *,'xmin,xmax=',real(xmin),real(xmax)
cc      print *,'ymin,ymax=',real(ymin),real(ymax)
c
c
cc     Set the initial boundaries of the grid near the outer edge of the fluid:
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
c
c      lastTOTALLUM=-1d30
c
cc     Use a 3x3 grid for starters:
c      NXMAPnext=3
c      NYMAPnext=3

 102  format(" ",A5," ",A15  ," ",A15  ," ",A5," ",A15  ," ",A15  ," ",
     $     A5," ",A15  )
 103  format(" ",i5," ",ES15.7," ",ES15.7," ",i5," ",ES15.7," ",
     $     ES15.7," ",i5," ",ES15.7)
      write(*,102) "iter","xmin","hx","Nx","ymin","hy","Ny","L"
      write(*,102) repeat("-",5), repeat("-",15), repeat("-",15),
     $     repeat("-",5), repeat("-",15), repeat("-",15),
     $     repeat("-",5), repeat("-",15)
      
      counter = 1


      
 41   continue
      
c     nxmap and nymap are how many grid cells to use in the X and Y
c     directions, respectively:
      nxmap=nxmapnext
      nymap=nymapnext
      
c     xminmap and xmaxmap are the left & rightmost coordinates of grid
c     yminmap and ymaxmap are the bottom & topmost coordinates of grid
      xminmap=xminmapnext
      xmaxmap=xmaxmapnext
      yminmap=yminmapnext
      ymaxmap=ymaxmapnext
      
c     Compute cell widths:             
      hxmap=(xmaxmap-xminmap)/dble(nxmap-1)
      hymap=(ymaxmap-yminmap)/dble(nymap-1)


      ! These are the maximum and minimum cell positions that still
      ! contribute significantly to the total luminosity. Initializing
      !iminglow = nxmap+1
      !jminglow = nymap+1
      !imaxglow = 0
      !jmaxglow = 0
      
      
      ! Setup the 3D integrating grid
      prepareIntegrationCalled = .false. ! Remake the grid
      call prepareIntegration

      ! Integrate over it to get the temperatures
      call integrateTau

      ! Determine the luminosity from the grid
      call getFlux

      
      write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap,
     $     TOTALpracticalLUM
      
      ! integrateTau determines the glowing cells. When we find an
      ! optical depth greater than tau_thick at a cell, we reset the min
      ! and max glowing cell positions to encompass that grid cell.

      deltai=IMAXGLOW-IMINGLOW
      deltaj=JMAXGLOW-JMINGLOW
      if( (deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy)) .or. 
     $     deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy)) .or.
     $     abs(TOTALpracticalLUM-lastTOTALLUM).gt.
     $     fracaccuracy*TOTALpracticalLUM)
     $     .and.(nxmap.lt.nxmapmax .or. nymap.lt.nymapmax) ) then
!     The code gets inside this if statement if not enough of the grid
!     has "glowing" cells, or if the total luminosity is not in
!     sufficient agreement with the total luminosity calculated from the
!     previous grid size.
c     print *,'no convergence yet...try again'
         
         ! In case there are two well-separated stars, don't start
         ! making smaller grids too quickly:
         if(IMAXGLOW-IMINGLOW.le.2) then
            IMINGLOW=min(IMINGLOW,2)
            IMAXGLOW=max(IMAXGLOW,NXMAP-1)
         endif
         if(JMAXGLOW-JMINGLOW.le.2) then
            JMINGLOW=min(JMINGLOW,2)
            JMAXGLOW=max(JMAXGLOW,NYMAP-1)
         endif
            
         ! Decide upon appropriate boundaries for the next grid:
         xmaxmapnext=min(Imaxglow*hxmap+xminmap + 3*hymap ,xmax)
         xminmapnext=max((Iminglow-2)*hxmap+xminmap - 3*hymap ,xmin)
         ymaxmapnext=min(Jmaxglow*hymap+yminmap + 3*hxmap ,ymax)
         yminmapnext=max((Jminglow-2)*hymap+yminmap - 3*hxmap ,ymin)
         
         ! Determine whether we should resolve the x or y direction
         ! better in the next grid:
         hxtmp=(xmaxmapnext-xminmapnext)/(nxmap-1)
         hytmp=(ymaxmapnext-yminmapnext)/(nymap-1)
         if(hxtmp.gt.hytmp .and. nxmap.lt.nxmapmax) then
            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+7,nxmapmax)
         else if(hytmp.gt.hxtmp .and. nymap.lt.nymapmax) then
            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+7,nymapmax)
         else
            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+6,nxmapmax)
            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+6,nymapmax)
         endif
         
         counter = counter + 1
         !lasttotallum = totallum ! Record for comparison
         lasttotallum = TOTALpracticalLUM
         goto 41
      end if

c 100  format(A19," = ",ES25.17e3,", ",ES25.17e3)
c 101  format(A19," = ",I25,", ",I25)
      write(*,*) ""
      write(*,*) "Grid details:"
      write(*,*) "   xminmap,xmaxmap = ",xminmap/runit_out,
     $     xmaxmap/runit_out
      write(*,*) "   yminmap,ymaxmap = ",yminmap/runit_out,
     $     ymaxmap/runit_out
      write(*,*) "   hxmap,hymap     = ",hxmap/runit_out,hymap/runit_out
      write(*,*) "   nxmap,nymap     = ",nxmap,nymap
      write(*,*) ""


      end subroutine
