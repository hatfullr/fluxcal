      subroutine writeFluxesFile
      include 'optical_depth.h'
      integer i,j

      write(*,*) "Writing '" // trim(adjustl(fname)) // "'"
      open (72,file=fname)
      write(72,*)xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap,t/tunit_out
      do j=1,nymap
         write (72,'(2001f8.3)') (log10((TpracticalXYthick(i,j)**4.d0
     $        + TpracticalXYthin(i,j)**4.d0)**0.25d0), i=1,nxmap)
      enddo
      do ifilter=1,numfilters
         write(72,*)real(mag(ifilter)),'= ',filtername(ifilter),
     $        ' magnitude ( ',real(1d4*wavelength(ifilter)),
     $        ' micron) for filter number', ifilter,
     $        'which comes from total flux density=',
     $        real(totalfluxdensity(ifilter)),
     $        '& zero magnitude flux density=',
     $        real(absoluteflux(ifilter)),
     $        '(both in erg cm^-2 sec^-1 Hz^-1)'
         do j=1,nymap
            write (72,'(2001g11.4)') (log10(fluxdensityXY(i,j,ifilter)),
     $           i=1,nxmap)
         enddo
      enddo
      close(72)
      
      end subroutine
