c     Taken from https://stackoverflow.com/a/8508757 and edited

      module iso_fortran_env
      
      ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
      ! See Subclause 13.8.2 of the Fortran 2003 standard. 

      implicit none
      public 
      
      integer, parameter :: stderr = 0
      integer, parameter :: stdin = 5 
      integer, parameter :: stdout = 6
      integer, parameter :: e = stderr
      integer, parameter :: o = stdout
      
      end module iso_fortran_env
