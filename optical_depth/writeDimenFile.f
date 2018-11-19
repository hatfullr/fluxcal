      subroutine writeDimenFile
c     L987 to L1013
      include 'optical_depth.h'
      integer i

c     write to dimen*.sph file
      if(.not. dimenFileAlreadyExists) then
         write(*,*) "Writing dimen file"
         OPEN (73,FILE=FNAME2)
         write(73,*)xminmap,hxmap,nxmap,yminmap,hymap,nymap
         close(73)
      endif

      write(*,*) "Writing '" // trim(adjustl(fname)) // "'"
      OPEN (72,FILE=FNAME)
      write(72,*)xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap,t/tunit_out
      DO J=1,NYMAP
         WRITE (72,'(2001f8.3)') (log10((TpracticalXYthick(I,J)**4.d0
     $        + TpracticalXYthin(I,J)**4.d0)**0.25d0), I=1,NXMAP)
      ENDDO
      do ifilter=1,numfilters
         write(72,*)real(mag(ifilter)),'= ',filtername(ifilter),
     $        ' magnitude ( ',real(1d4*wavelength(ifilter)),
     $        ' micron) for filter number', ifilter,
     $        'which comes from total flux density=',
     $        real(totalfluxdensity(ifilter)),
     $        '& zero magnitude flux density=',
     $        real(absoluteflux(ifilter)),
     $        '(both in erg cm^-2 sec^-1 Hz^-1)'
         DO J=1,NYMAP
            WRITE (72,'(2001g11.4)') (log10(fluxdensityXY(I,J,ifilter)),
     $           I=1,NXMAP)
         ENDDO
      enddo
      CLOSE(72)

      end subroutine
