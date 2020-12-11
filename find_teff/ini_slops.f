      subroutine ini_slops(filename)

      include "../lib/flux_cal.h"
      integer i,j
      character*4 dummy
      character*255 filename

c     initialize opacity file for photosphere fit
      write(o,*) "Reading envfit nabla file"
      open(1,file=trim(adjustl(filename)),status='old')
      num_g=51
      num_t=271
 
      read(1,*) dummy, dummy, (grid_gn(i),i=1,num_g)

      do j=1,num_t
         read(1,*) grid_tn(j), (nabla(i,j),i=1,num_g)
      end do
      close(1)

      table_nabla_Tmin = 10.**grid_tn(1)
      table_nabla_Tmax = 10.**grid_tn(num_t)
      table_nabla_gmin = 10.**grid_gn(1)
      table_nabla_gmax = 10.**grid_gn(num_g)
      
      return
      end
