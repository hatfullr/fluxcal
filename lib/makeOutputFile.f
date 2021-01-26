      subroutine makeOutputFile(filename)
      include 'flux_cal.h'
      character*255 filename

      write(o,*) "Initializing '",trim(adjustl(filename)),"'"
      
 101  format(20A22)
      open(10,file=trim(adjustl(filename)),status='replace')
      write(10,101)
     $     "x",
     $     "y",
     $     "z",
     $     "m",
     $     "h",
     $     "rho",
     $     "u",
     $     "mu",
     $     "g",
     $     "T_SPH",
     $     "Teff",
     $     "P",
     $     "P_env",
     $     "P_surf",
     $     "opacity",
     $     "opacity_surf",
     $     "rho_surf",
     $     "tau",
     $     "entropy",
     $     "ID"
      close(10)
      end subroutine
