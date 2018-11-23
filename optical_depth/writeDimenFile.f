      subroutine writeDimenFile
      include 'optical_depth.h'
      logical fileexists
      
 103  format(" ",ES22.14," ",ES22.14," ",i5," ",ES22.14," ",
     $     ES22.14," ",i5)
      
c     write to dimen*.dat file
      write(*,*) "Writing dimen file"
      open (73,file=fname2)
      write(73,103)xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap
      close(73)

      end subroutine
