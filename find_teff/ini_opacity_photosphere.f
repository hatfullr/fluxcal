      subroutine ini_opacity_photospehre(filename)

      include "../lib/flux_cal.h"
      
      integer il,jl
      character*4 dummy
      character*255 filename
       
c     initialize opacities for a photospheric fit
      write(*,*) "Reading photospheric opacity file"
      open(1,file=trim(adjustl(filename)),status='old')

      numr=19
      read(1,*) dummy
      read(1,*) dummy
      read(1,*) dummy, dummy, (gridR(il),il=1,numr)

      numo=85
      do jl=1,numo
         read(1,*) gridT(jl), (kap(il,jl),il=1,numr)
      end do
      close(1)

      return
      end

      subroutine ini_opacity_photospehre_lowT(filename)

      include "../lib/flux_cal.h"
      
      integer il,jl, k
      character*4 dummy
      character*255 filename
       
c     initialize opacities for a photospheric fit
c     open(1,file="lowT_fa05_gs98_z0.02_x0.7.data ",status='old')
      write(*,*) "Reading photospheric cold opacity file"
      open(1,file=trim(adjustl(filename)),status='old')

      numr=46
      read(1,*) dummy
      read(1,*) dummy
      read(1,*) dummy
      read(1,*) dummy
      read(1,*) (gridR(il),il=1,numr)

      numo=105
      do k=1,numo
         jl=numo+1-k
         read(1,*) gridT(jl), (kap(il,jl),il=1,numr)
      end do
      close(1)
      return
      end


