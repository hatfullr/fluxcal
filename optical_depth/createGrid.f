      subroutine createGrid
      include 'optical_depth.h'

      integer counter
      integer i,j
      integer lastnumcell
      character*10 reason
      
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

 102  format(
     $     " ",A5,              !iter
     $     " ",A11,             !xmin
     $     " ",A11,             !hx
     $     " ",A5,              !Nx
     $     " ",A11,             !ymin
     $     " ",A11,             !hy
     $     " ",A5,              !Ny
     $     " ",A12,             !|minstpsize|
     $     " ",A12,             !|maxstpsize|
     $     " ",A8,              !minNstp
     $     " ",A8,              !maxNstp
     $     " ",A11,             !<F>
     $     " ",A10              !reason
     $     )
 103  format(
     $     " ",i5,              !iter
     $     " ",ES11.4,          !xmin
     $     " ",ES11.4,          !hx
     $     " ",i5,              !Nx
     $     " ",ES11.4,          !ymin
     $     " ",ES11.4,          !hy
     $     " ",i5,              !Ny
     $     " ",ES12.5,          !|minstpsize|
     $     " ",ES12.5,          !|maxstpsize|
     $     " ",i8,              !minNstp
     $     " ",i8,              !maxNstp
     $     " ",ES11.4,          !<F>
     $     " ",A10              !reason
     $     )
      write(*,102)
     $     "iter",
     $     "xmin",
     $     "hx",
     $     "Nx",
     $     "ymin",
     $     "hy",
     $     "Ny",
     $     "|minstpsize|",
     $     "|maxstpsize|",
     $     "minNstp",
     $     "maxNstp",
     $     "<F>",
     $     "reason"
      write(*,102)
     $     repeat("-",5),       !iter
     $     repeat("-",11),      !xmin
     $     repeat("-",11),      !hx
     $     repeat("-",5),       !Nx
     $     repeat("-",11),      !ymin
     $     repeat("-",11),      !hy
     $     repeat("-",5),       !Ny
     $     repeat("-",12),      !|minstpsize|
     $     repeat("-",12),      !|maxstpsize|
     $     repeat("-",8),       !minNstp
     $     repeat("-",8),       !maxNstp
     $     repeat("-",11),      !<F>
     $     repeat("-",10)       !reason
      
      counter = 1
      lastnumcell = int(1.d-30)
      
      reason = "init"

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

c      write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
c     $     yminmap/runit_out,hymap/runit_out,nymap,
cc     $     TOTALpracticalLUM/Lunit_out
c     $     TOTALflux/numcell,reason
      
      ! integrateTau determines the glowing cells. When we find an
      ! optical depth greater than tau_thick at a cell, we reset the min
      ! and max glowing cell positions to encompass that grid cell.

      ! IMINGLOW is the minimum cell index in the x direction for which
      ! we detect a photosphere, and IMAXGLOW is the maximum. Same for
      ! JMINGLOW and JMAXGLOW but in the y direction.
      ! Thus, deltai and deltaj are the number of cells within which
      ! we have detected photospheres.

      deltai=IMAXGLOW-IMINGLOW
      deltaj=JMAXGLOW-JMINGLOW
      reason = ""
c      reasons(1) = deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy))
c      reasons(2) = deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy))
c      reasons(3) = abs(TOTALflux/numcell-lastTOTALflux/lastnumcell).gt.
c     $     fracaccuracy*TOTALflux/numcell
      if( (deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy)) .or. 
     $     deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy)) .or.
c      if ((deltai.le.0.5*nxmap .or.
c     $     deltaj.le.0.5*nymap .or.
c     $     abs(TOTALpracticalLUM-lastTOTALLUM).gt.
c     $     fracaccuracy*TOTALpracticalLUM)
     $     abs((TOTALflux/numcell)-(lastTOTALflux/lastnumcell)).gt.
     $     fracaccuracy*TOTALflux/numcell) .and.
     $     (nxmap.lt.nxmapmax .or. nymap.lt.nymapmax) ) then
!     The code gets inside this if statement if not enough of the grid
!     has "glowing" cells, or if the total luminosity is not in
!     sufficient agreement with the total luminosity calculated from the
!     previous grid size.
c     print *,'no convergence yet...try again'

c     Default to complaining about F convergence, otherwise complain about
c     not having enough grid cells in X or Y directions
         reason = "F cnvgnce"
         
c         ! In case there are two well-separated stars, don't start
c         ! making smaller grids too quickly:
c         if(IMAXGLOW-IMINGLOW.le.2) then
c            IMINGLOW=min(IMINGLOW,2)
c            IMAXGLOW=max(IMAXGLOW,NXMAP-1)
c         endif
c         if(JMAXGLOW-JMINGLOW.le.2) then
c            JMINGLOW=min(JMINGLOW,2)
c            JMAXGLOW=max(JMAXGLOW,NYMAP-1)
c         endif
            
         ! Decide upon appropriate boundaries for the next grid:
         xmaxmapnext=min(IMAXGLOW*hxmap+xminmap + 3*hymap ,xmax)
         xminmapnext=max((IMINGLOW-2)*hxmap+xminmap - 3*hymap ,xmin)
         ymaxmapnext=min(JMAXGLOW*hymap+yminmap + 3*hxmap ,ymax)
         yminmapnext=max((JMINGLOW-2)*hymap+yminmap - 3*hxmap ,ymin)
         
         ! Determine whether we should resolve the x or y direction
         ! better in the next grid:
         hxtmp=(xmaxmapnext-xminmapnext)/(nxmap-1)
         hytmp=(ymaxmapnext-yminmapnext)/(nymap-1)
         if(hxtmp.gt.hytmp .and. nxmap.lt.nxmapmax) then
            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+7,nxmapmax)
            reason = "Nx cells"
         else if(hytmp.gt.hxtmp .and. nymap.lt.nymapmax) then
            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+7,nymapmax)
            reason = "Ny cells"
         else
            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+6,nxmapmax)
            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+6,nymapmax)
            reason = "Nxy cells"
         endif

         !lasttotallum = totallum ! Record for comparison
         lasttotallum = TOTALpracticalLUM
         lasttotalflux = TOTALflux
         lastnumcell = numcell

         if(TOTALflux.eq.0) then
            write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $           yminmap/runit_out,hymap/runit_out,nymap,min_step_size,
     $           max_step_size,min_steps_taken,max_steps_taken,
     $           0.d0,"Finished  "
            goto 42
         else
            write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $           yminmap/runit_out,hymap/runit_out,nymap,min_step_size,
     $           max_step_size,min_steps_taken,max_steps_taken,
     $           TOTALflux/numcell,reason
         end if

         counter = counter+1
         
         goto 41
      end if

      if(TOTALflux.eq.0) then
         write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $        yminmap/runit_out,hymap/runit_out,nymap,min_step_size,
     $        max_step_size,min_steps_taken,max_steps_taken,
     $        0.d0,"Finished  "
      else
         write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $        yminmap/runit_out,hymap/runit_out,nymap,min_step_size,
     $        max_step_size,min_steps_taken,max_steps_taken,
     $        TOTALflux/numcell,"Finished  "
      end if

 42   end subroutine
