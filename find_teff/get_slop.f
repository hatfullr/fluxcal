      real*8 function get_slop(ttt, ggg,tmin,slop_types)

      include "../lib/flux_cal.h"

      real*8 temp, gsurf,ttt, ggg, tmin, f_g
      real*8  slop
      integer found,i,j, k,look_this_g,slop_types

c     default slop
      slop=0.3
      slop_type=0

      temp=log10(ttt)
      gsurf=log10(ggg)

      tmin=-1      

      do i=1,num_g-1
         look_this_g=0
         if(gsurf.ge.grid_gn(i).and.gsurf.le.grid_gn(i+1)) look_this_g=1
         if(i.eq.1.and.gsurf.le.grid_gn(i)) look_this_g=1
         f_g=(gsurf-grid_gn(i))/(grid_gn(i+1)-grid_gn(i))
         if(look_this_g.eq.1) then
            slop_type=10        ! for this slop, we have located g
            do j=1,num_t-1
               if(nabla(i,j).gt.-10.and.tmin.lt.0) tmin=grid_tn(j) ! identified minimum temperature
               if(temp.ge.grid_tn(j).and.temp.le.grid_tn(j+1).
     &            and.nabla(i,j).gt.-10.) then
                  slop_type=1 ! slop was calculated
                  slop=nabla(i,j)+f_g*(nabla(i+1,j)-nabla(i,j))
c                  write(*,*) slop,grid_tn(j),nabla(i,j),nabla(i+1,j),f_g
                  exit
               end if
               if(temp.ge.grid_tn(j).and.temp.le.grid_tn(j+1).
     &            and.nabla(i,j).le.-10.) then
                  slop_type=2   ! slop was taken from the closets point
                  do k=j,num_t
                     if(nabla(i,j).gt.-10.) then
                        slop=nabla(i,j)
                        exit
                     end if
                  end do
               end if
            end do
            if(temp.gt.grid_tn(num_t)) then
               slop_type=3     ! slop is extrapolated towards higher T
               slop=nabla(i,num_t)
            end if
         end if
      end do

      
      get_slop=slop
      slop_types=slop_type


      return
      end
