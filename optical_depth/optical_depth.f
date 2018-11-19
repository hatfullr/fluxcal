      subroutine optical_depth
c For the light emission, the code assume blackbody.  There is no
c treatment of recombination.  What the light curve routine does is to
c integrate along a grid of parallel lines of sight that approach the
c system from any desired viewing angle.  We integrate
c     d(tau)=kappa*rho*ds
c to get find where the optical depth tau is 1.  Here ds is just an
c infinitesimal length along the line of sight.  The integration starts
c from outside and works into the gas.  The opacity kappa comes from OPAL
c with a supplemental table for low temperatures. Kappa takes temperature
c T and density rho from the SPH results (pretending that they can be
c believed even in the outskirts of the system).  At tau=1 we record the
c temperature T for each grid cell.  In this way, we get a 2-D grid of
c temperatures.  The flux from each of these cells is assumed to be
c sigma*T^4.  The luminosity from that cell is therefore A*sigma*T^4,
c where A is the area of the cell.  The total energy per time being
c emitted toward the observer is the sum of A*sigma*T^4 over all the
c cells.  The "luminosity" is 4 times this value.  I put quotes around
c luminosity because it's not actually the energy per time being emitted
c in all directions; however it is what an observer would calculate as
c the energy per time if he or she assumed the objected was spherically
c symmetric.  The 4 comes from the fact that the surface area of
c a sphere is 4 times its projected area on the plane of the sky. There
c are some lines of sight that pass through gas but never reach tau=1.
c These lines of sight do not contribute to the luminosity.
c 
c The main output is done with these lines
c       write(6,'(g11.3,11g12.4,a)') t*tunit, TOTALLUM,tavg,sigmat,
c      $     t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,sigmavz*vunit,
c      $     vravg*vunit,sigmavr*vunit,' lc'
c       write(42,'(g11.3,99g12.4)') t, TOTALLUM,tavg,sigmat,t4avg**0.25d0,
c      $     rmid,sigmar,area,vzavg*vunit,sigmavz*vunit,vravg*vunit,
c      $     sigmavr*vunit
c The first writes to standard output and can be most easily found by
c storing that standard output in a file and grepping on "lc".  The second
c writes to fort.42.  The columns are time (in days), luminosity (solar
c units), temperature (K), standard deviation of temperature (K), (average
c of T^4)^0.25 (K), average distance from origin, standard deviation of
c average distance from origin, line of sight velocity (km/s), standard
c deviation of line of sight velocity (km/s), radial velocity (km/s)m and
c     the standard deviation of the radial velocity (km/s).

      include 'optical_depth.h'
      real*8 dummy
      integer counter, i

      counter = 0
      codeerror = .false.
      
      NVAR=4+numfilters

      do insl=0,nsl
         do ilogwavelength=ilogwavelengthmin,ilogwavelengthmax
            spectrum(ilogwavelength,insl)=0.d0
         enddo
      enddo

c      vunit=sqrt(gravconst*munit/runit)*1d-5 ! convert to km/s ! StarSmasher specific
c      tunit=sqrt(runit**3/(gravconst*munit))/3600/24 ! convert to days ! StarSmasher specific
c      print *,'unit of velocity (in km/s)=',real(vunit)
c      print *,'unit of time (in days)=',real(tunit)
ccccc L164

      write(*,*) "Estabilishing the grid"
      call createGrid
      
c 102  format(A5," ",A9," ",A9," ",A5," ",A9," ",A9," ",A5)
c     103  format(i5," ",f9.2," ",f9.2," ",i5," ",f9.2," ",f9.2," ",i5)
 102  format(" ",A5," ",A15  ," ",A15  ," ",A5," ",A15  ," ",A15  ," ",
     $     A5," ",A15  )
 103  format(" ",i5," ",E15.7," ",E15.7," ",i5," ",E15.7," ",E15.7," ",
     $     i5," ",E15.7)
      
      inquire(file=FNAME2,exist=dimenFileAlreadyExists)
      if(.not.dimenFileAlreadyExists) then
         write(*,102) "iter","xmin","hx","Nx","ymin","hy","Ny","L"
c         write(*,102) "----","---------","---------","-----",
c     $        "---------","---------","-----"
         write(*,102) repeat("-",5), repeat("-",15), repeat("-",15),
     $        repeat("-",5), repeat("-",15), repeat("-",15),
     $        repeat("-",5), repeat("-",15)
      end if

 41   continue

      counter = counter + 1
      call useDimenFile
c      if(.not.dimenFileAlreadyExists) then
cc         write(*,103) counter,xminmap,xmaxmap,yminmap,ymaxmap,
cc     $        nxmap,nymap
c         write(*,103) counter,xminmap,hxmap,nxmap,yminmap,hymap,nymap
c      end if
c      write(6,*)'*** dimen information:****************'
c      write(6,*)xminmap,hxmap,nxmap,yminmap,hymap,nymap
c     write(6,*)'**************************************'
      call prepareIntegration
      call integrateTau
      call getFlux

      if(.not.dimenFileAlreadyExists) then
         write(*,103) counter,xminmap,hxmap,nxmap,yminmap,hymap,nymap,
     $        TOTALpracticalLUM
      end if

      if(codeerror) goto 626

c      write(*,*) "Luminosity=",TOTALpracticalLUM,"Lsun"
c      write(*,*) "   V1309 L=",5.85562637712,"Lsun"
ccccc L819 to L842
c      if(area.gt.0.d0) then
c         write(6,'(a11,99a12)') 't', 'TOTALLUM','Tavg','sigmaT',
c     $        'T^4avg^0.25','rmid','sigmar','area','vzavg',
c     $        'sigmavz','vravg','sigmavr','area_photo',
c     $        'T^4_ph^0.25','area_therm',
c     $        'T^4_th^0.25','L_ph','L_th','N_photo',
c     $        'Tavg_ph','Tavg_th','Warnings',
c     $        'ccphoto ','gT^4_ph^0.25','L(good)','...'
cc         write(6,'(g11.3,99g12.4)') t*tunit, TOTALpracticalLUM,tavg,
cc     $        sigmat,
cc     $        t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,
cc     $        sigmavz*vunit,vravg*vunit,sigmavr*vunit,
cc     $        areaphoto,tphoto4avg**0.25d0,
cc     $        areatherm,ttherm4avg**0.25d0,
cc     $        konstant*areaphoto*tphoto4avg,
cc     $        konstant*areatherm*ttherm4avg,
cc     $        nphotoavg,tphotoavg,tthermavg,warning,ccphoto,
cc     $        (goodtphoto4avg(insl)**0.25d0,
cc     $        konstant*areaphoto*goodtphoto4avg(insl),insl=nsl,1,-1),
cc     $        tpractical4avg**0.25d0,
c          write(6,'(g11.3,99g12.4)') dummy, TOTALpracticalLUM,tavg,
c     $        sigmat,
c     $        t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,
c     $        sigmavz*vunit,vravg*vunit,sigmavr*vunit,
c     $        areaphoto,tphoto4avg**0.25d0,
c     $        areatherm,ttherm4avg**0.25d0,
c     $        konstant*areaphoto*tphoto4avg,
c     $        konstant*areatherm*ttherm4avg,
c     $        nphotoavg,tphotoavg,tthermavg,warning,ccphoto,
c     $        (goodtphoto4avg(insl)**0.25d0,
c     $        konstant*areaphoto*goodtphoto4avg(insl),insl=nsl,1,-1),
c     $        tpractical4avg**0.25d0,
cc     $        9000d0/(mag(1)-mag(2)+0.93d0), ! temperature in K from eq. (13.36)
cc                                                of Ryden and Peterson's textbook
c     $        (mag(ifilter),ifilter=1,numfilters),vPCygni*vunit
c      endif
ccccc L842
      

ccccc L844 to L906

c     fracaccuracy=0.1d0      ! make faster with some sacrifice in accuracy
c      fracaccuracy=0.01d0  ! A good default value

c      if(numfilters.lt.5)then
c         fracaccuracy=0.0025d0   ! A good default value
c      else
c         fracaccuracy=0.025d0  ! A good default value
c      endif

c      print *, 'Using fracaccuracy=',real(fracaccuracy)
c      fracaccuracy=0.9d0        ! Use this for less strict convergence criteria and debugging

      deltai=IMAXGLOW-IMINGLOW
      deltaj=JMAXGLOW-JMINGLOW
      if( (deltai.le.max(NXMAP/2,nint(20*0.025d0/fracaccuracy)) .or. 
     $     deltaj.le.max(NYMAP/2,nint(20*0.025d0/fracaccuracy)) .or.
     $     abs(TOTALLUM-lastTOTALLUM).gt.fracaccuracy*TOTALLUM) .and.
     $     .not. dimenFileAlreadyExists
     $     .and. (NXMAP.lt.nxmapmax .or. NyMAP.lt.nymapmax) ) then

!     The code gets inside this if statement if not enough of the grid
!     has "glowing" cells, or if the total luminosity is not in
!     sufficient agreement with the total luminosity calculated from the
!     previous grid size.
c         print *,'no convergence yet...try again'

         call flush(6)

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
         xmaxmapnext=min(Imaxglow*HXMAP+XMINMAP + 3*HYMAP ,xmax)
         xminmapnext=max((Iminglow-2)*HXMAP+XMINMAP - 3*HYMAP ,xmin)
         ymaxmapnext=min(Jmaxglow*HyMAP+yMINMAP + 3*Hxmap ,ymax)
         yminmapnext=max((Jminglow-2)*HyMAP+yMINMAP - 3*Hxmap ,ymin)

C     Determine whether we should resolve the x or y direction better
C     in the next grid:
         hxtmp=(XMAXMAPNEXT-XMINMAPNEXT)/(NXMAP-1)
         hytmp=(YMAXMAPNEXT-YMINMAPNEXT)/(NYMAP-1)
         if(HXtmp.gt.HYtmp .and. nxmap.lt.nxmapmax) then
            NXMAPnext=min(2*max(IMAXGLOW-IMINGLOW,0)+7,NXMAPMAX)
         else if(HYtmp.gt.HXtmp .and. nymap.lt.nymapmax) then
            NYMAPnext=min(2*max(JMAXGLOW-JMINGLOW,0)+7,NYMAPMAX)
         else
            NXMAPnext=min(2*max(IMAXGLOW-IMINGLOW,0)+6,NXMAPMAX)
            NYMAPnext=min(2*max(JMAXGLOW-JMINGLOW,0)+6,NYMAPMAX)
         endif

c     Record the total luminosity for a later comparison:
         lastTOTALLUM=TOTALLUM

         goto 41
      endif

 100  format(A19," = ",E11.4,", ",E11.4)
 101  format(A19," = ",I11,", ",I11)
      write(*,*) "Grid details:"
      write(*,100) "xminmap,xmaxmap",xminmap,xmaxmap
      write(*,100) "yminmap,ymaxmap",yminmap,ymaxmap
      write(*,100) "hxmap,hymap    ",hxmap,hymap
      write(*,101) "nxmap,nymap    ",nxmap,nymap
      write(*,*) ""
      
      write (*,*) '    TXY_max=',real(TMAX),real(log10(TMAX))
      write (*,*) 'Total luminosity in enclosed area=',real(TOTALLUM)
      write (*,*) 'Total practical luminosity in enclosed area=',
     $     real(TOTALpracticalLUM)
      write (*,*) 'Average Teff=',real(avgt)/real(numcell)
      
      call peakWavelengths


      
ccccc L952 to L985
c     Record the total luminosity both in the std output
cc      write(6,'(f11.3,53g12.4,a,99g12.4)')t*tunit,TOTALLUM,tavg,sigmat,
cc     $     t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,sigmavz*vunit, ! 10 items so far
cc     $     vravg*vunit,sigmavr*vunit,areaphoto,tphoto4avg**0.25d0,
cc     $     areatherm,ttherm4avg**0.25d0, konstant*areaphoto*tphoto4avg,
cc     $     konstant*areatherm*ttherm4avg, nphotoavg, tphotoavg,
cc     $     tthermavg,warning,ccphoto,  ! 23 items so far
cc     $     (goodtphoto4avg(insl)**0.25d0, 
cc     $     konstant*areaphoto*goodtphoto4avg(insl),insl=nsl,1,-1), ! 43 items so far
cc     $     (Tfit(insl),insl=nsl,0,-1),   ! 54 items so far
cc     $     ' lc',               ! the text 'lc' should be in the 55th column
cc     $     tpractical4avg**0.25d0,
cc     $     (mag(ifilter),ifilter=1,numfilters),vPCygni*vunit

c      write(6,'(f11.3,53g12.4,a,99g12.4)')dummy,TOTALLUM,tavg,sigmat,
c     $     t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,sigmavz*vunit, ! 10 items so far
c     $     vravg*vunit,sigmavr*vunit,areaphoto,tphoto4avg**0.25d0,
c     $     areatherm,ttherm4avg**0.25d0, konstant*areaphoto*tphoto4avg,
c     $     konstant*areatherm*ttherm4avg, nphotoavg, tphotoavg,
c     $     tthermavg,warning,ccphoto,  ! 23 items so far
c     $     (goodtphoto4avg(insl)**0.25d0, 
c     $     konstant*areaphoto*goodtphoto4avg(insl),insl=nsl,1,-1), ! 43 items so far
c     $     (Tfit(insl),insl=nsl,0,-1),   ! 54 items so far
c     $     ' lc',               ! the text 'lc' should be in the 55th column
c     $     tpractical4avg**0.25d0,
c     $     (mag(ifilter),ifilter=1,numfilters),vPCygni*vunit
cc      write(42,'(g11.3,99g12.4)') t*tunit,TOTALLUM,tavg,sigmat,
cc     $     t4avg**0.25d0,rmid,sigmar,area,vzavg*vunit,sigmavz*vunit,
cc     $     vravg*vunit,sigmavr*vunit,areaphoto,tphoto4avg**0.25d0,
cc     $     areatherm,ttherm4avg**0.25d0

      call flush(6)

      IF (MLOG.EQ.1) THEN                                              
C     Transform to log scale:  
         DO J=1,NYMAP
            DO I=1,NXMAP       
               IF (TXY(I,J).GE.TMIN) THEN 
c                  IF(TXY(I,J).eq.TMAX) PRINT *, 'MAX AT',I,J
                  TXY(I,J)=LOG10(TXY(I,J))
               ELSE
                  TXY(I,J)=0.d0
               END IF     
            END DO        
         END DO           
      END IF    
ccccc L985



      call writeDimenFile

 626  return
      end
