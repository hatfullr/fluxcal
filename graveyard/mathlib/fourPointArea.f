      real*8 function fourPointArea(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
c     This function calculates the surface area traced out by four
c     points in 3D space by splitting the area into four triangles
c     and using Heron's formula to calculate the area of each one.
c
c     To identify the points, assume all z=0.d0, then
c     "1" is the bottom-left
c     "2" is the top-left
c     "3" is the top-right
c     "4" is the bottom-right
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      real*8 d1,d2,d3,d4
      real*8 xmin,xamx,ymin,ymax,zmin,zmax
      real*8 xm,ym,zm
      real*8 s, A1,A2,A3,A4
      real*8 d1m,d2m,d3m,d4m
      integer nbad(4),i,total
      real*8 A, sfactor(4)

      do i=1,4
         nbad(i)=0
      end do
c     If we have more than 1 point w/ infinite depth, area=0
      if((abs(z1).gt.1.d30).or.(z1.ne.z1)) then
         nbad(1) = 1
      end if
      if((abs(z2).gt.1.d30).or.(z2.ne.z2)) then
         nbad(2) = 1
      end if
      if((abs(z3).gt.1.d30).or.(z3.ne.z3)) then
         nbad(3) = 1
      end if
      if((abs(z4).gt.1.d30).or.(z4.ne.z4)) then
         nbad(4) = 1
      end if

      total = 0
      do i=1,4
         total = total + nbad(i)
      end do
      
      if(total.gt.1) then
         fourPointArea=0.d0
         goto 20
      end if
      
c     If we have 1 point with an infinite depth, find triangle area
      if(total.eq.1) then
         d1 = 0.d0
         d2 = 0.d0
         d3 = 0.d0
         d4 = 0.d0
         s = 0.d0
         do i=1,4
            sfactor(i)=1.d0
         end do
         if((nbad(1).eq.0).and.(nbad(2).eq.0)) then
            d1 = ((x1-x2)**2.d0 + (y1-y2)**2.d0 + (z1-z2)**2.d0)**0.5d0
            s = s + 0.5d0*d1
            sfactor(1) = s-d1
         end if
         if((nbad(2).eq.0).and.(nbad(3).eq.0)) then
            d2 = ((x2-x3)**2.d0 + (y2-y3)**2.d0 + (z2-z3)**2.d0)**0.5d0
            s = s + 0.5d0*d2
            sfactor(2) = s-d2
         end if
         if((nbad(3).eq.0).and.(nbad(4).eq.0)) then
            d3 = ((x3-x4)**2.d0 + (y3-y4)**2.d0 + (z3-z4)**2.d0)**0.5d0
            s = s + 0.5d0*d3
            sfactor(3) = s-d3
         end if
         if((nbad(4).eq.0).and.(nbad(1).eq.0)) then
            d4 = ((x4-x1)**2.d0 + (y4-y1)**2.d0 + (z4-z1)**2.d0)**0.5d0
            s = s + 0.5d0*d4
            sfactor(4) = s-d4
         end if

         fourPointArea=
     $        (s*sfactor(1)*sfactor(2)*sfactor(3)*sfactor(4))**0.5d0
         goto 20
      end if
      
      
         


      
      
c     Calculate outer triangle distances
      d1 = ((x1-x2)**2.d0 + (y1-y2)**2.d0 + (z1-z2)**2.d0)**0.5d0
      d2 = ((x2-x3)**2.d0 + (y2-y3)**2.d0 + (z2-z3)**2.d0)**0.5d0
      d3 = ((x3-x4)**2.d0 + (y3-y4)**2.d0 + (z3-z4)**2.d0)**0.5d0
      d4 = ((x4-x1)**2.d0 + (y4-y1)**2.d0 + (z4-z1)**2.d0)**0.5d0

c     Calculate minimums and maximums of all dimensions
      xmin = min(x1,x2,x3,x4)
      xmax = max(x1,x2,x3,x4)
      ymin = min(y1,y2,y3,y4)
      ymax = max(y1,y2,y3,y4)
      zmin = min(z1,z2,z3,z4)
      zmax = max(z1,z2,z3,z4)

c     Calculate the position of the midpoint
      xm = (xmax+xmin)/2.
      ym = (ymax+ymin)/2.
      zm = (zmax+zmin)/2.

c     Calculate distances to the midpoint
      d1m = ((x1-xm)**2.d0 + (y1-ym)**2.d0 + (z1-zm)**2.d0)**0.5d0
      d2m = ((x2-xm)**2.d0 + (y2-ym)**2.d0 + (z2-zm)**2.d0)**0.5d0
      d3m = ((x3-xm)**2.d0 + (y3-ym)**2.d0 + (z3-zm)**2.d0)**0.5d0
      d4m = ((x4-xm)**2.d0 + (y4-ym)**2.d0 + (z4-zm)**2.d0)**0.5d0
      
c     Calculate area of the triangles
      s = (d1+d1m+d2m)/2.
      A1 = (s*(s-d1m)*(s-d2m)*(s-d1))**0.5d0

      s = (d2+d2m+d3m)/2.
      A2 = (s*(s-d2m)*(s-d3m)*(s-d2))**0.5d0

      s = (d3+d3m+d4m)/2.
      A3 = (s*(s-d3m)*(s-d4m)*(s-d3))**0.5d0

      s = (d4+d4m+d1m)/2.
      A4 = (s*(s-d4m)*(s-d1m)*(s-d4))**0.5d0

c     Return the surface area
      fourPointArea = A1+A2+A3+A4

 20   end function
