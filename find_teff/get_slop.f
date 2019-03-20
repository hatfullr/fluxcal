      real*8 function get_slop(ttt, ggg,tmin,slop_types)
c     Finds nabla = dlnT/dlnP for values of ttt and ggg from the tabulated nabla values
c     in defaults/table_nabla.dat.
c      
c     ttt = t_sph, the central temperature for this particle
c     ggg = g_sph, the local gravitational acceleration for this particle
c     tmin = identified minimum temperature in the nabla table
c     slop_types = The type of slope found
c         10 = Slope has been found
c          1 = Slope has been calculated
c          2 = Slope was calculated from the closest points
c          3 = Slope was extrapolated towards higher T
      include "../lib/flux_cal.h"

      real*8 temp, gsurf,ttt, ggg, tmin, f_g
      real*8  slop
      integer found,i,j, k,look_this_g,slop_types

c     default slop
      slop=0.3
      slop_type=0

      temp=log10(ttt)           ! Table has log values in it
      gsurf=log10(ggg)

      tmin=-1

      ! grid_gn and grid_tn are the lgg and lgT values from the nabla table
      do i=1,num_g-1            ! Loop over the g values in the nabla table
         look_this_g=0
         ! Check if ggg is in this part of the table
         if(gsurf.ge.grid_gn(i).and.gsurf.le.grid_gn(i+1)) look_this_g=1
         if(i.eq.1.and.gsurf.le.grid_gn(i)) look_this_g=1 ! It's between the first two values
         ! Calculate the fractional position of ggg within the range grid_gn(i) <= ggg <= grid_gn(i+1)
         f_g=(gsurf-grid_gn(i))/(grid_gn(i+1)-grid_gn(i))
         if(look_this_g.eq.1) then
            slop_type=10        ! for this slop, we have located g
            do j=1,num_t-1      ! Loop over the T values in the nabla table
               ! Technically there is one nabla(i,j)=-10 value, but it should probably be -20 (the lower limit)
               if(nabla(i,j).gt.-10.and.tmin.lt.0) tmin=grid_tn(j) ! identified minimum temperature
               if(temp.ge.grid_tn(j).and.temp.le.grid_tn(j+1) ! Check if ttt is in this part of the table
     &            .and.nabla(i,j).gt.-10.) then
                  slop_type=1   ! slop was calculated
                  
                  ! Here, f_g*(nabla(i+1,j)-nabla(i,j)) is the value of nabla that corresponds with
                  ! the fractional position of ggg within grid_gn(i) <= ggg <= grid_gn(i+1).
                  ! We add nabla(i,j) to this to get the actual value of nabla for this ggg.
                  slop=nabla(i,j)+f_g*(nabla(i+1,j)-nabla(i,j))
c                  write(*,*) slop,grid_tn(j),nabla(i,j),nabla(i+1,j),f_g
                  exit
               end if
               ! If the slope value is very small (-20, or -10 or so) for this value of ggg and ttt,
               ! we will simply take the closest nabla above -10 for this ggg value
               if(temp.ge.grid_tn(j).and.temp.le.grid_tn(j+1)
     &            .and.nabla(i,j).le.-10.) then
                  slop_type=2   ! slop was taken from the closets point
                  do k=j,num_t  ! Loop over T values outside of the previously found T value
                     if(nabla(i,j).gt.-10.) then ! Until we find a suitable nabla
                        slop=nabla(i,j)
                        exit
                     end if
                  end do        ! Stop looping over T values outside of the previously found T value
               end if
            end do              ! Stop looping over T values in the nabla table

            ! If ttt is larger than the maximum tabulated T, take nabla to be the value at the
            ! largest T available for this value of ggg.
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
