      subroutine getFlux
c     L669 to L817
      include 'optical_depth.h'

      real*8 TOTALTpracticalXY
      integer i
      
      TMAX=0.d0
      TMIN=1.d+30
      TOTALLUM=0.d0
      TOTALpracticalLUM=0.d0
      tavg=0.d0
      t2avg=0.d0
      t4avg=0.d0
      tphoto4avg=0.d0
      tpractical4avg=0.d0
      ttherm4avg=0.d0
      tphotoavg=0.d0
      tthermavg=0.d0
      cellcount=0
      cctherm=0

      natalum = 0.d0
      sphlum = 0.d0
      natalumhi = 0.d0
      natalumlo = 0.d0

      do ifilter=1,numfilters
         totalfluxdensity(ifilter)=0.d0
      enddo
      
      if(codeerror) goto 626

      if(dimenFileAlreadyExists) then
         write(*,*) "Calculating the flux at each grid point"
      end if
      DO J=1,NYMAP
         DO I=1,NXMAP
            TOTALLUM=TOTALLUM+TXY(I,J)**4.d0
c            TOTALTpracticalXY=TpracticalXYthin(i,j)**4.d0+
c     $           TpracticalXYthick(i,j)**4.d0
c            if(thick_part(i,j).gt.0) then
c               TOTALpracticalLUM=TOTALpracticalLUM
c     $              + Avis(thick_part(i,j))
c     $              *TOTALTpracticalXY
c            else
c               TOTALpracticalLUM=TOTALpracticalLUM
c     $              + HXMAP*HYMAP*TOTALTpracticalXY
c            end if

            ! Attenuate the luminosity from the thick particles if
            ! there is otpically thin stuff in front of it.
            ! Contributions from optically thin stuff are already
            ! attenuated by the derivs2.f function.
            TOTALpracticalLUM=TOTALpracticalLUM
     $           + TpracticalXYthick(I,J)**4.d0 * exp(-tauthin(i,j))
     $           + TpracticalXYthin(I,J)**4.d0

            natalum=natalum+nataT(I,J)**4
            sphlum=sphlum+sphT(I,J)**4
            natalumhi=natalumhi+nataThi(I,J)**4
            natalumlo=natalumlo+nataTlo(I,J)**4

            
            if(totalpracticallum.ne.totalpracticallum)then
               print *,'impractical practical luminsoity at ',i,j,
     $              TpracticalXYthick(I,J),TpracticalXYthin(I,J)
               write(*,*) tauthin(i,j)
               write(*,*) exp(-tauthin(i,j))
               write(*,*) TOTALpracticalLUM
               write(*,*) TOTALTpracticalXY
               write(*,*) thick_part(i,j)
               write(*,*) Avis(thick_part(i,j))
               stop
            endif

            do ifilter=1,numfilters
               totalfluxdensity(ifilter)=totalfluxdensity(ifilter)+
     $              fluxdensityXY(I,J,ifilter)*HXMAP*HYMAP
     $              /distance**2
            enddo

            TMAX=MAX(TMAX,TXY(I,J))
            IF(TXY(I,J).GT.0.d0) then
               TMIN=MIN(TMIN,TXY(I,J))
               tavg=tavg+TXY(I,J)
               t2avg=t2avg+TXY(I,J)**2
               t4avg=t4avg+TXY(I,J)**4
               cellcount=cellcount+1
            endif
            IF(TphotoXY(I,J).GT.0.d0) then
               tphotoavg=tphotoavg+TphotoXY(I,J)
               tphoto4avg=tphoto4avg+TphotoXY(I,J)**4
               tpractical4avg=tpractical4avg
     $              + TpracticalXYthick(I,J)**4.d0
     $              + TpracticalXYthin(I,J)**4.d0
            endif
            IF(TthermXY(I,J).GT.0.d0) then
               tthermavg=tthermavg+TthermXY(I,J)
               ttherm4avg=ttherm4avg+TthermXY(I,J)**4
               cctherm=cctherm+1
            endif
         ENDDO
      ENDDO

   
      do ifilter=1,numfilters
         mag(ifilter)=-2.5d0*log10(
     $        totalfluxdensity(ifilter)/absoluteflux(ifilter))
      enddo

!     IMINGLOW = leftmost   grid cell that "glows" (has tau>1)
!     IMAXGLOW = rightmost  grid cell that "glows" (has tau>1)
!     JMINGLOW = bottommost grid cell that "glows" (has tau>1)
!     JMAXGLOW = topmost    grid cell that "glows" (has tau>1)
c      print *,'IMINGLOW,IMAXGLOW,JMINGLOW,JMAXGLOW=',IMINGLOW,IMAXGLOW,
c     $     JMINGLOW,JMAXGLOW

      tavg=tavg/cellcount
      t2avg=t2avg/cellcount
      t4avg=t4avg/cellcount
      do insl=1,nsl
         if (goodccphoto(insl).gt.0) then
            goodtphoto4avg(insl)=goodtphoto4avg(insl)/goodccphoto(insl)
         endif
      enddo
c      if(warning+goodccphoto(10).ne.ccphoto) then
c         print *,'counting error',warning,goodccphoto(10),ccphoto
!     warning = number of cells for which the photosphere is less
!               than 10 smoothing lengths from surface
!     goodccphoto(10) = number of cells for which the photosphere is
!                       at least 10 smoothing lengths from surface
!     ccphoto = total number of cells spanned by the photosphere
c         stop
c      endif
      nphotoavg=nphotoavg/ccphoto
      if(cctherm.gt.0)then
         ttherm4avg=ttherm4avg/cctherm
         tthermavg=tthermavg/cctherm
      endif
      rmid=rmid/cellcount
      r2mid=r2mid/cellcount
      call getLocalvz(xcosthetamax,ycosthetamax,zcosthetamax,
     $     localvx,localvy,vPCygni)
      
      vzavg=vzavg/tphoto4avg   ! Use this now that weighting by flux
      vz2avg=vz2avg/tphoto4avg ! Use this now that weighting by flux
      tphotoavg=tphotoavg/ccphoto
      tphoto4avg=tphoto4avg/ccphoto
      tpractical4avg=tpractical4avg/ccphoto
      vravg=vravg/ccphoto
      vr2avg=vr2avg/ccphoto
      area=cellcount*HXMAP*HYMAP
      areaphoto=ccphoto*HXMAP*HYMAP
      areatherm=cctherm*HXMAP*HYMAP
c     write (6,*) 'Total luminosity in enclosed area=',TOTALLUM

c     Get luminosity to solar units:
c     The factor of 4 in the following line comes from thinking about the
c     case when the object is spherical, with radius R.  We are seeing a
c     projection of area pi*R^2, while the true surface area is 4*pi*R^2.
c     So we multiply by 4 to get what we would presume the full
c     luminosity L should be.  The actal luminosity could of course be
c     different than this, depending, for example, on how brightly
c     glowing the back side of the system is in reality.  This approach
c     is more or less what one would do in "real life" too, as you can't
c     see the back side of the object you're measuring.
c      TOTALLUM=TOTALLUM*sigma*4d0*HXMAP*HYMAP/3.839d33
c      TOTALpracticalLUM=
c     $     TOTALpracticalLUM*sigma*4d0*HXMAP*HYMAP/3.839d33
      TOTALLUM=TOTALLUM*sigma*4d0*HXMAP*HYMAP/Lunit_out
      TOTALpracticalLUM=
     $     TOTALpracticalLUM*sigma*4d0*HXMAP*HYMAP/Lunit_out
c      TOTALpracticalLUM=TOTALpracticalLUM*sigma

      natalum=natalum*sigma*4d0*HXMAP*HYMAP
      sphlum=sphlum*sigma*4d0*HXMAP*HYMAP
      natalumhi=natalumhi*sigma*4d0*HXMAP*HYMAP
      natalumlo=natalumlo*sigma*4d0*HXMAP*HYMAP
      
      
c      WRITE (6,*) '    TXY_max=',real(TMAX),real(log10(TMAX))
c      write (6,*) 'Total luminosity in enclosed area=',real(TOTALLUM)
c      write (6,*) 'Total practical luminosity in enclosed area=',
c     $     real(TOTALpracticalLUM)

      sigmat=(t2avg-tavg**2)**0.5d0
      if(cellcount.le.1) then
         sigmar=0.d0
      else
         sigmar=(r2mid-rmid**2)**0.5d0
      endif

      if(sigmar.ne.sigmar) then
         print *,'sigmar=',sigmar
         print *,'r2mid=',r2mid
         print *,'rmid=',rmid
         print *,'r2mid-rmid**2=',r2mid-rmid**2
         print *,'cellcount=',cellcount
         sToP
      endif

      sigmavz=(vz2avg-vzavg**2)**0.5d0
      sigmavr=(vr2avg-vravg**2)**0.5d0

      konstant=4d0*sigma/3.839d33

 626  continue
      end subroutine
