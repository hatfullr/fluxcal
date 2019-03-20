      subroutine getClosestParticles
      ! Name of the game is to sort the particles by their distance
      ! to the observer and then omit particles that are completely
      ! obscured by the ones in front of them.
      include 'optical_depth.h'
      integer pid(n)
      integer nid,count
      integer cidsize, i,j
      real*8 myhxmap, myhymap
      real*8 xparam, yparam

      integer min_val, max_val
      integer closest1D(nxmapmax*nymapmax)
      real*8 maxdzs1D(nxmapmax*nymapmax)
      integer unique(nxmapmax*nymapmax)
      real*8 surfaces(nxmapmax*nymapmax)
      real*8 maxdzs(nxmapmax,nymapmax)
      integer myn
      character*255 filename
      integer sizeunique

      integer,dimension(:),allocatable :: holder1
      real*8, dimension(:), allocatable :: holder2

      real time1,time2

      real*8 hpip2,dyip2
      
      integer m

      real*8 zmaxx(nxmapmax,nymapmax)
      
c      call createGrid
c      call useDimenFile
c      call prepareIntegration

      count=0
      do j=1,nymap
         do i=1,nxmap
            count = count+1
            closest(i,j)=-1
            zmaxx(i,j)=-1.d30
            maxdzs(i,j) = 0.d0
         enddo
      enddo

      do ip=1,n
         if(a(ip).gt.0.d0) then
            hpip2 = hp(ip)**2.d0
            yposmin=y(ip)-2.d0*hp(ip)
            yposmax=y(ip)+2.d0*hp(ip)
            jmin=max(int((yposmin-yminmap)/hymap+2),1)
            jmax=min(int((yposmax-yminmap)/hymap+1),nymap)
            DO j=jmin,jmax
               YPOS=(J-1)*HYMAP+YMINMAP ! y of line of sight
               dyip2 = (y(ip)-ypos)**2.d0
               maxdx=(4.d0*hpip2-dyip2)**0.5d0
               xposmin=x(ip)-maxdx
               xposmax=x(ip)+maxdx
               imin=max(int((xposmin-xminmap)/hxmap+2),1)
               imax=min(int((xposmax-xminmap)/hxmap+1),nxmap)
               do i=imin,imax
                  XPOS=(i-1)*HXMAP+XMINMAP ! x of line of sight
                  maxdz=(4*hpip2
     $                 -(x(ip)-xpos)**2.d0-dyip2)**0.5d0
                  if(z(ip)+maxdz.gt.zmaxx(i,j)) then
                     closest(i,j) = ip
                     maxdzs(i,j) = maxdz
                     zmaxx(i,j)=z(ip)+maxdz
                  end if
               enddo
            enddo
         endif
      end do

c     Condense the "closest" and "maxdzs" arrays into 1D arrays
      do i=1,nxmapmax*nymapmax
         closest1D(i) = -1
         maxdzs1D(i) = 0.d0
      end do
      
      count=0
      do j=1,nymap
         do i=1,nxmap
            if(closest(i,j).ge.0) then
               count=count+1
               closest1D(count) = closest(i,j)
               maxdzs1D(count) = maxdzs(i,j)
            end if
         end do
      end do
      
c     This is a sorting algorithm I found online to sort an array and
c     only worry about unique values. Thus, we get rid of repeated
c     particle IDs
      m=0
      min_val = minval(closest1D)
      max_val = maxval(closest1D)
      do while (min_val<max_val)
         m = m+1
         min_val = minval(closest1D, mask=closest1D>min_val)
         unique(m) = min_val
      enddo

      allocate(holder1(m))
      allocate(holder2(m))
      do i=1,m
         holder1(i) = unique(i)
         holder2(i) = z(unique(i)) + maxdzs1D(unique(i))
      end do
      
      call quicksort2(holder2,holder1,1,m)

c     Reverse the order
      do i=1,m
         closest1D(i) = holder1(m-i+1)
         surfaces(i) = holder2(m-i+1)
      end do
      

 35   format ('closest',i5.5,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
 36   format ('closest',i6.6,'_',i3.3,'_',i3.3,'_',i3.3,'.dat')
      
      if(innit.le.99999) then
         write (filename,35) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      else
         write (fname2,36) innit,nint(anglez*180.d0/pi),
     $        nint(angley*180.d0/pi),nint(anglex*180.d0/pi)
      end if

      call makeOutputFile(filename)
      do i=1,m
        ! write(*,*) holder1(i)
         if(closest1D(i).ge.0) then
            call output(filename,closest1D(i),.false.)
         end if
      end do

      end subroutine

