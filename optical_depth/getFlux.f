      subroutine getFlux
c     L669 to L817
      include 'optical_depth.h'

      integer i,j
      real*8 Acell, myxpos,myypos
      integer nillum
      
      TMAX=0.d0
      TMIN=1.d+30
      TOTALLUM=0.d0
      TOTALpracticalLUM=0.d0
      TOTALflux=0.d0
c     We might want this later
c      tavg=0.d0
c      t2avg=0.d0
c      t4avg=0.d0
c      tphoto4avg=0.d0
c      tpractical4avg=0.d0
c      ttherm4avg=0.d0
c      tphotoavg=0.d0
c      tthermavg=0.d0
c      cellcount=0
c      cctherm=0
c
c      do ifilter=1,numfilters
c         totalfluxdensity(ifilter)=0.d0
c      enddo
      
      if(codeerror) goto 626

      if(dimenFileAlreadyExists) then
         write(o,*) "Calculating the flux at each grid point"
      end if
      Atot = 0.d0
      DO J=1,NYMAP
         DO I=1,NXMAP
            myxpos=(i-1)*hxmap+xminmap ! x-coordinate of line of sight
            myypos=(j-1)*hymap+yminmap ! y-coordinate of line of sight

            ! Do it flat
c            Acell = fourPointArea(
c     $           myxpos,myypos,0.d0,
c     $           myxpos,myypos+hymap,0.d0,
c     $           myxpos+hxmap,myypos+hymap,0.d0,
c     $           myxpos+hxmap,myypos,0.d0)
            
c            Acell = fourPointArea(
c     $           myxpos,myypos,zintpos(i,j),
c     $           myxpos,myypos+hymap,zintpos(i,j+1),
c     $           myxpos+hxmap,myypos+hymap,zintpos(i+1,j+1),
c     $           myxpos+hxmap,myypos,zintpos(i+1,j))

c     We might want this later
c            TOTALLUM=TOTALLUM+TXY(I,J)**4.d0
            TOTALLUM=TOTALLUM+TOTALTpracticalXY(i,j)**4.d0

            ! Attenuate the luminosity from the thick particles if
            ! there is otpically thin stuff in front of it.
            ! Contributions from optically thin stuff are already
            ! attenuated by the derivs.f function.
c            TOTALpracticalLUM=TOTALpracticalLUM
c     $           + TpracticalXYthick(I,J)**4.d0 * exp(-tauthin(i,j))
c     $           + TpracticalXYthin(I,J)**4.d0
c            TOTALpracticalLUM=TOTALpracticalLUM
c     $           + (TpracticalXYthick(I,J)**4.d0 * exp(-tauthin(i,j))
c     $           + TpracticalXYthin(I,J)**4.d0) * hxmap*hymap

c            ! Count the vertices that have Teff > 0
c            nillum = 0
c            if(TOTALTpracticalXY(i,j) .gt. 0)     nillum = nillum + 1
c            if(TOTALTpracticalXY(i+1,j) .gt. 0)   nillum = nillum + 1
c            if(TOTALTpracticalXY(i,j+1) .gt. 0)   nillum = nillum + 1
c            if(TOTALTpracticalXY(i+1,j+1) .gt. 0) nillum = nillum + 1

            Acell = hxmap*hymap
c            if(nillum.eq.1) then
c               !Acell = Acell*0.125d0
c               Acell = 0.25d0 * Acell
c               Atot = Atot + Acell
c            else if (nillum .eq. 2) then
c               !Acell = Acell*0.5d0
c               Acell = 0.5d0 * Acell
c               Atot = Atot + Acell
c            else if (nillum .eq. 3) then
c               !Acell = Acell*0.875
c               Acell = 0.75d0 * Acell
c               Atot = Atot + Acell
c            else if (nillum.eq.4) then
c               Atot = Atot + Acell
c            end if
            if(TOTALTpracticalXY(i,j).gt.0) Atot = Atot+Acell
c            TOTALpracticalLUM=TOTALpracticalLUM
c     $           + TOTALTpracticalXY(i,j)**4.d0

c     We need to find the convergence in the total flux from the
c     optically thick regions. THEN apply attenuation to each
c     contribution. This way, we can calculate the specific
c     luminosity without worrying about surface areas.

            TOTALflux = TOTALflux + TOTALTpracticalXY(i,j)**4.d0

            if(TOTALflux.ne.TOTALflux) then
               write(o,*) "i,j,xpos,ypos = ",i,j,xpos,ypos
               write(o,*) "TOTALflux = ",TOTALflux
               write(o,*) "TOTALTpracticalXY(i,j) = ",
     $              TOTALTpracticalXY(i,j)
               write(o,*) "Acell = ", Acell
               write(o,*) "ERROR: NaN flux"
               error stop "getFlux.f"
            end if

            
            if(totalpracticallum.ne.totalpracticallum)then
               print *,'impractical practical luminsoity at ',i,j
               write(o,*) tauthin(i,j)
               write(o,*) exp(-tauthin(i,j))
               write(o,*) TOTALpracticalLUM
               write(o,*) TOTALTpracticalXY(i,j)
               write(o,*) thick_part(i,j)
               write(o,*) Avis(thick_part(i,j))
               error stop "getFlux.f"
            endif

            TMAX=MAX(TMAX,TOTALTpracticalXY(i,j))
            if(TOTALTpracticalXY(i,j).gt.0.d0) then
               TMIN=MIN(TMIN,TOTALTpracticalXY(i,j))
            end if
            
c     We might want this later
c            do ifilter=1,numfilters
c               totalfluxdensity(ifilter)=totalfluxdensity(ifilter)+
c     $              fluxdensityXY(I,J,ifilter)*HXMAP*HYMAP
c     $              /distance**2
c            enddo
c
c            TMAX=MAX(TMAX,TXY(I,J))
c            IF(TXY(I,J).GT.0.d0) then
c               TMIN=MIN(TMIN,TXY(I,J))
c               tavg=tavg+TXY(I,J)
c               t2avg=t2avg+TXY(I,J)**2
c               t4avg=t4avg+TXY(I,J)**4
c               cellcount=cellcount+1
c            endif
c
c            IF(TphotoXY(I,J).GT.0.d0) then
c               tphotoavg=tphotoavg+TphotoXY(I,J)
c               tphoto4avg=tphoto4avg+TphotoXY(I,J)**4
c               tpractical4avg=tpractical4avg
c     $              + TOTALTpracticalXY(I,J)**4.d0
c
c            endif
c
c            IF(TthermXY(I,J).GT.0.d0) then
c               tthermavg=tthermavg+TthermXY(I,J)
c               ttherm4avg=ttherm4avg+TthermXY(I,J)**4
c               cctherm=cctherm+1
c            endif
         ENDDO
      ENDDO


c     We might want this later
c      do ifilter=1,numfilters
c         mag(ifilter)=-2.5d0*log10(
c     $        totalfluxdensity(ifilter)/absoluteflux(ifilter))
c      enddo

!     IMINGLOW = leftmost   grid cell that "glows" (has tau>1)
!     IMAXGLOW = rightmost  grid cell that "glows" (has tau>1)
!     JMINGLOW = bottommost grid cell that "glows" (has tau>1)
!     JMAXGLOW = topmost    grid cell that "glows" (has tau>1)
c      print *,'IMINGLOW,IMAXGLOW,JMINGLOW,JMAXGLOW=',IMINGLOW,IMAXGLOW,
c     $     JMINGLOW,JMAXGLOW

c      tavg=tavg/cellcount
c      t2avg=t2avg/cellcount
c      t4avg=t4avg/cellcount
c     We might want this later
c      do insl=1,nsl
c         if (goodccphoto(insl).gt.0) then
c            goodtphoto4avg(insl)=goodtphoto4avg(insl)/goodccphoto(insl)
c         endif
c      enddo
cc      if(warning+goodccphoto(10).ne.ccphoto) then
cc         print *,'counting error',warning,goodccphoto(10),ccphoto
c!     warning = number of cells for which the photosphere is less
c!               than 10 smoothing lengths from surface
c!     goodccphoto(10) = number of cells for which the photosphere is
c!                       at least 10 smoothing lengths from surface
c!     ccphoto = total number of cells spanned by the photosphere
cc         stop
cc      endif
c      nphotoavg=nphotoavg/ccphoto
c      if(cctherm.gt.0)then
c         ttherm4avg=ttherm4avg/cctherm
c         tthermavg=tthermavg/cctherm
c      endif
c      rmid=rmid/cellcount
c      r2mid=r2mid/cellcount
c      call getLocalvz(xcosthetamax,ycosthetamax,zcosthetamax,
c     $     localvx,localvy,vPCygni)
      
c      vzavg=vzavg/tphoto4avg   ! Use this now that weighting by flux
c      vz2avg=vz2avg/tphoto4avg ! Use this now that weighting by flux
c      tphotoavg=tphotoavg/ccphoto
c      tphoto4avg=tphoto4avg/ccphoto
c      tpractical4avg=tpractical4avg/ccphoto
c      vravg=vravg/ccphoto
c      vr2avg=vr2avg/ccphoto
c      area=cellcount*HXMAP*HYMAP
c      areaphoto=ccphoto*HXMAP*HYMAP
c      areatherm=cctherm*HXMAP*HYMAP
c     write (6,*) 'Total luminosity in enclosed area=',TOTALLUM

c     Get luminosity to solar units:
c     The factor of 4 in the following line comes from thinking about the
c     case when the object is spherical, with radius R. We are seeing a
c     projection of area pi*R^2, while the true surface area is 4*pi*R^2.
c     So we multiply by 4 to get what we would presume the full
c     luminosity L should be. The actual luminosity could of course be
c     different than this, depending, for example, on how brightly
c     glowing the back side of the system is in reality. This approach
c     is more or less what one would do in "real life" too, as you can't
c     see the back side of the object you're measuring.
c      TOTALLUM=TOTALLUM*sigma*4d0*HXMAP*HYMAP/3.839d33
      TOTALLUM=TOTALLUM*sigma*4d0*HXMAP*HYMAP
c      TOTALpracticalLUM=
c     $     TOTALpracticalLUM*sigma*4d0*HXMAP*HYMAP/Lunit_out
      TOTALflux = TOTALflux*sigma
      TOTALpracticalLUM=TOTALflux*4.d0*hxmap*hymap

c      write(o,*) Atot/runit_out**2.d0
c      write(o,*) Atot/runit_out**2.d0 / nxmap/nymap,
c     $     hxmap*hymap / runit_out**2.d0

c     We might want this later
c      sigmat=(t2avg-tavg**2)**0.5d0
c      if(cellcount.le.1) then
c         sigmar=0.d0
c      else
c         sigmar=(r2mid-rmid**2)**0.5d0
c      endif
c
c      if(sigmar.ne.sigmar) then
c         print *,'sigmar=',sigmar
c         print *,'r2mid=',r2mid
c         print *,'rmid=',rmid
c         print *,'r2mid-rmid**2=',r2mid-rmid**2
c         print *,'cellcount=',cellcount
c         sToP
c      endif
c
c      sigmavz=(vz2avg-vzavg**2)**0.5d0
c      sigmavr=(vr2avg-vravg**2)**0.5d0
c
c      konstant=4d0*sigma/3.839d33

 626  continue
      end subroutine
