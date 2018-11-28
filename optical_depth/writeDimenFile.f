      subroutine writeDimenFile
      include 'optical_depth.h'
      logical fileexists

c     We need to write the dimen files with a certain amount of
c     precision, so that reading them back in results in the same
c     numbers regardless of runtime environment. We have to use list-
c     directed input formatting to do this. It does not work any other
c     way.
      
c 103  format(" ",ES23.15e3," ",ES23.15e3," ",i5," ",ES23.15e3," ",
c     $     ES23.15e3," ",i5)
      
c     write to dimen*.dat file
      write(*,*) "Writing dimen file"
      open (73,file=trim(adjustl(fname2)))
      write(73,*)xminmap/runit_out,hxmap/runit_out,nxmap,
     $     yminmap/runit_out,hymap/runit_out,nymap
      close(73)

      end subroutine
