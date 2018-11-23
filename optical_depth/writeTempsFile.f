      subroutine writeTempsFile
      include 'optical_depth.h'
      integer i,j

      write(*,*) "Writing '" // trim(adjustl(fname)) // "'"
      open (72,file=fname)
 103  format(" ",ES22.14," ",ES22.14," ",i5," ",ES22.14," ",
     $     ES22.14," ",i5," ",ES22.14)
      write(72,103)xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap,t/tunit_out
      do j=1,nymap
         write (72,'(2001f8.3)') (log10((TpracticalXYthick(i,j)**4.d0
     $        + TpracticalXYthin(i,j)**4.d0)**0.25d0 / tempunit_out),
     $        i=1,nxmap)
      enddo
c      do ifilter=1,numfilters
c         write(72,*)real(mag(ifilter)),'= ',filtername(ifilter),
c     $        ' magnitude ( ',real(1d4*wavelength(ifilter)),
c     $        ' micron) for filter number', ifilter,
c     $        'which comes from total flux density=',
c     $        real(totalfluxdensity(ifilter)),
c     $        '& zero magnitude flux density=',
c     $        real(absoluteflux(ifilter)),
c     $        '(both in erg cm^-2 sec^-1 Hz^-1)'
c         do j=1,nymap
c            write (72,'(2001g11.4)') (log10(fluxdensityXY(i,j,ifilter)),
c     $           i=1,nxmap)
c         enddo
c      enddo
      close(72)
      
      end subroutine
