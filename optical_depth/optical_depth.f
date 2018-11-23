      subroutine optical_depth
      include 'optical_depth.h'
      real*8 dummy
      integer counter, i,j

      counter = 0
      codeerror = .false.
      
      nvar=4+numfilters

      do insl=0,nsl
         do ilogwavelength=ilogwavelengthmin,ilogwavelengthmax
            spectrum(ilogwavelength,insl)=0.d0
         enddo
      enddo

c      vunit=sqrt(gravconst*munit/runit)*1d-5 ! convert to km/s ! StarSmasher specific
c      tunit=sqrt(runit**3/(gravconst*munit))/3600/24 ! convert to days ! StarSmasher specific
c      print *,'unit of velocity (in km/s)=',real(vunit)
c      print *,'unit of time (in days)=',real(tunit)

      write(*,*) "Estabilishing the grid"
      call createGrid

 102  format(" ",A5," ",A22  ," ",A22  ," ",A5," ",A22  ," ",A22  ," ",
     $     A5," ",A22  )
 103  format(" ",i5," ",ES22.14," ",ES22.14," ",i5," ",ES22.14," ",
     $     ES22.14," ",i5," ",ES22.14)
      
      inquire(file=trim(adjustl(fname2)),exist=dimenFileAlreadyExists)
      if(dimenFileAlreadyExists) then
         call useDimenFile
         call prepareIntegration
         call integrateTau
      else
         write(*,102) "iter","xmin","hx","Nx","ymin","hy","Ny","L"
         write(*,102) repeat("-",5), repeat("-",22), repeat("-",22),
     $        repeat("-",5), repeat("-",22), repeat("-",22),
     $        repeat("-",5), repeat("-",22)

 41      continue
         counter=counter+1
c        nxmap and nymap are how many grid cells to use in the X and Y
c        directions, respectively:
         nxmap=nxmapnext
         nymap=nymapnext

c        xminmap and xmaxmap are the left & rightmost coordinates of grid
c        yminmap and ymaxmap are the bottom & topmost coordinates of grid
         xminmap=xminmapnext
         xmaxmap=xmaxmapnext
         yminmap=yminmapnext
         ymaxmap=ymaxmapnext

C        Compute cell widths:             
         hxmap=(xmaxmap-xminmap)/dble(nxmap-1)
         hymap=(ymaxmap-yminmap)/dble(nymap-1)

         prepareIntegrationCalled=.false. !Remake the integrating grid
         call prepareIntegration
         call integrateTau
         call getFlux

         write(*,103) counter,xminmap/runit_out,hxmap/runit_out,nxmap,
     $        yminmap/runit_out,hymap/runit_out,nymap,
     $        TOTALpracticalLUM

         if(codeerror) goto 626

         deltai=IMAXGLOW-IMINGLOW
         deltaj=JMAXGLOW-JMINGLOW
         if( (deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy)) .or. 
     $        deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy)) .or.
     $        abs(TOTALLUM-lastTOTALLUM).gt.fracaccuracy*TOTALLUM) .and.
     $        (nxmap.lt.nxmapmax .or. nymap.lt.nymapmax) ) then
            
!     The code gets inside this if statement if not enough of the grid
!     has "glowing" cells, or if the total luminosity is not in
!     sufficient agreement with the total luminosity calculated from the
!     previous grid size.
c     print *,'no convergence yet...try again'
            
!     In case there are two well-separated stars, don't start making
!     smaller grids too quickly:
            if(IMAXGLOW-IMINGLOW.le.2) then
               IMINGLOW=min(IMINGLOW,2)
               IMAXGLOW=max(IMAXGLOW,NXMAP-1)
            endif
            if(JMAXGLOW-JMINGLOW.le.2) then
               JMINGLOW=min(JMINGLOW,2)
               JMAXGLOW=max(JMAXGLOW,NYMAP-1)
            endif
            
C     Decide upon appropriate boundaries for the next grid:
            xmaxmapnext=min(Imaxglow*hxmap+xminmap + 3*hymap ,xmax)
            xminmapnext=max((Iminglow-2)*hxmap+xminmap - 3*hymap ,xmin)
            ymaxmapnext=min(Jmaxglow*hymap+yminmap + 3*hxmap ,ymax)
            yminmapnext=max((Jminglow-2)*hymap+yminmap - 3*hxmap ,ymin)
            
C     Determine whether we should resolve the x or y direction better
C     in the next grid:
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
            
c     Record the total luminosity for a later comparison:
            lastTOTALLUM=TOTALLUM
            goto 41
         end if
         call writeDimenFile
      end if
      
      call writeTempsFile
      call getFlux
      
c      inquire(file=trim(adjustl(fname2)),exist=dimenFileAlreadyExists)
c      if(.not.dimenFileAlreadyExists) then
c         write(*,102) "iter","xmin","hx","Nx","ymin","hy","Ny","L"
c         write(*,102) repeat("-",5), repeat("-",22), repeat("-",22),
c     $        repeat("-",5), repeat("-",22), repeat("-",22),
c     $        repeat("-",5), repeat("-",22)
c      end if
c 41   continue
c 
c      counter = counter + 1
c      call useDimenFile
c      call prepareIntegration
c      call integrateTau
c      call getFlux
c 
c      ! Write current grid information to standard output
c      if(.not.dimenFileAlreadyExists) then
c         write(*,103) counter,xminmap,hxmap,nxmap,yminmap,hymap,nymap,
c     $        TOTALpracticalLUM
c      end if
c 
c      if(codeerror) goto 626
c 
c      deltai=IMAXGLOW-IMINGLOW
c      deltaj=JMAXGLOW-JMINGLOW
c      if( (deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy)) .or. 
c     $     deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy)) .or.
c     $     abs(TOTALLUM-lastTOTALLUM).gt.fracaccuracy*TOTALLUM) .and.
c     $     .not. dimenFileAlreadyExists
c     $     .and. (nxmap.lt.nxmapmax .or. nymap.lt.nymapmax) ) then
c 
c!     The code gets inside this if statement if not enough of the grid
c!     has "glowing" cells, or if the total luminosity is not in
c!     sufficient agreemecnt with the total luminosity calculated from the
c!     previous grid size.
cc         print *,'no convergence yet...try again'
c 
c!     In case there are two well-separated stars, don't start making
c!     smaller grids too quickly:
c         if(IMAXGLOW-IMINGLOW.le.2) then
c            IMINGLOW=min(IMINGLOW,2)
c            IMAXGLOW=max(IMAXGLOW,NXMAP-1)
c         endif
c         if(JMAXGLOW-JMINGLOW.le.2) then
c            JMINGLOW=min(JMINGLOW,2)
c            JMAXGLOW=max(JMAXGLOW,NYMAP-1)
c      endif
c      
cC     Decide upon appropriate boundaries for the next grid:
c         xmaxmapnext=min(Imaxglow*hxmap+xminmap + 3*hymap ,xmax)
c         xminmapnext=max((Iminglow-2)*hxmap+xminmap - 3*hymap ,xmin)
c         ymaxmapnext=min(Jmaxglow*hymap+yminmap + 3*hxmap ,ymax)
c         yminmapnext=max((Jminglow-2)*hymap+yminmap - 3*hxmap ,ymin)
c 
cC     Determine whether we should resolve the x or y direction better
cC     in the next grid:
c         hxtmp=(xmaxmapnext-xminmapnext)/(nxmap-1)
c         hytmp=(ymaxmapnext-yminmapnext)/(nymap-1)
c         if(hxtmp.gt.hytmp .and. nxmap.lt.nxmapmax) then
c            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+7,nxmapmax)
c         else if(hytmp.gt.hxtmp .and. nymap.lt.nymapmax) then
c            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+7,nymapmax)
c         else
c            nxmapnext=min(2*max(IMAXGLOW-IMINGLOW,0)+6,nxmapmax)
c            nymapnext=min(2*max(JMAXGLOW-JMINGLOW,0)+6,nymapmax)
c         endif
c 
cc     Record the total luminosity for a later comparison:
c         lastTOTALLUM=TOTALLUM
c         prepareIntegrationCalled=.false. ! Remake the integrating grid
c         goto 41
c      endif
     
 100  format(A19," = ",ES22.14,", ",ES22.14)
 101  format(A19," = ",I22,", ",I22)
      write(*,*) "Grid details:"
      write(*,100) "xminmap,xmaxmap",xminmap/runit_out,xmaxmap/runit_out
      write(*,100) "yminmap,ymaxmap",yminmap/runit_out,ymaxmap/runit_out
      write(*,100) "hxmap,hymap    ",hxmap/runit_out,hymap/runit_out
      write(*,101) "nxmap,nymap    ",nxmap,nymap
c      call writeDimenFile
      write(*,*) ""
      
c      write (*,*) 'Total luminosity in enclosed area=',real(TOTALLUM)
c      write (*,*) 'Total practical luminosity in enclosed area=',
c     $     real(TOTALpracticalLUM)
      write(*,'(A,ES22.14)')' Total luminosity in enclosed area=',
     $     TOTALpracticalLUM/Lunit_out
      write(*,'(A,ES22.14)')' Maximum Teff=',TMAX/tempunit_out
      write(*,'(A,ES22.14)')' Average Teff=',avgt/numcell/tempunit_out
      
      call peakWavelengths


      if (mlog.eq.1) then
C     Transform to log scale:
         do j=1,nymap
            do i=1,nxmap
               if (TXY(I,J).ge.TMIN) then
c                  if(TXY(I,J).eq.TMAX) print *, 'MAX AT',i,j
                  TXY(I,J)=LOG10(TXY(I,J))
               else
                  TXY(I,J)=0.d0
               end if
            end do
         end do
      end if

 626  return
      end
