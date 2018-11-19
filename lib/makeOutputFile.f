      subroutine makeOutputFile(filename)
      character*255 filename

      write(*,*) "Writing '",trim(adjustl(filename)),"'"
      
 101  format(14A15)
      open(10,file=trim(adjustl(filename)),status='replace')
      write(10,101) "x","y","z","m","h","rho","u","mu","g","T_SPH",
     $     "Teff","P","tau","ID"
      close(10)
      end subroutine
