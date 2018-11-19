      subroutine particlesAtPos(myposx,myposy,pid,nid)
      include 'optical_depth.h'

      integer pid(n)
      integer nid, i
      real*8 myposx, myposy, rprime
      real*8 zs(n),temp1(n),temp2(n)
      character*255 filename

      do ip=1,N
         pid(ip) = -1
         zs(ip) = -1.d30 ! As far back as possible
      end do

      nid = 0
      do ip=1,N
         ! Find out if the particle overlaps with the myposx,myposy
         rprime = ((myposx-x(ip))**2.d0+(myposy-y(ip))**2.d0)**0.5d0
         if(rprime.le.2.d0*hp(ip)) then
            nid=nid+1
            zs(nid) = z(ip) + (4.d0*hp(ip)**2.d0
     $           -(myposx-x(ip))**2.d0 - (myposy-y(ip))**2.d0)**0.5d0
            pid(nid) = ip
         end if
      end do

      ! Sort the particles in order of surface locations along the zaxis
      ! So the first particle's surface is closest to the observer
      call quicksort2(zs,pid,1,nid)

      ! Reverse the order
      do i=1,nid
         temp1(i) = zs(nid-i+1)
         temp2(i) = pid(nid-i+1)
      end do
      do i=1,nid
         zs(i) = temp1(i)
         pid(i) = temp2(i)
      end do


      if(innit.le.9999) then
         write(filename,"('part_at_pos_',i4.4,'.dat')") innit
      else if(innit.le.99999) then
         write(filename,"('part_at_pos_',i5.5,'.dat')") innit
      else
         write(filename,"('part_at_pos_',i6.6,'.dat')") innit
      end if
      call makeOutputFile(filename)
      do i=1,nid
         call output(filename,pid(i))
      end do
      
      return
      end subroutine
