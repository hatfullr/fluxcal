      subroutine integrateSpherical
      include 'optical_depth.h'

      integer i,j
      integer nthetamap, nphimap
      real*8 phipos, thetapos,rip
      real*8 dtheta, dphi
      
c     Initialize the polar grid
      anglexdeg = 0.d0
      angleydeg = 0.d0
      anglezdeg = 0.d0

      rmap = 0.d0 ! The maximum particle surface distance from origin
      do i=1,n
         rip=(x(i)**2.d0+y(i)**2.d0+z(i)**2.d0)**0.5d0
         rmap = max(rmap,rip+2.d0*hp(i))
      end do

c     Start with 3x3 grid
      nthetamap = 8             ! Move in increments of 2
      nphimap = 8               ! Move in increments of 2

c     Don't start on phi=0 each time
c      phimin = 0.d0
c      phimax = 360.d0*(1.d0-1.d0/dble(nphimap))
c      thetamin = 0.d0
c      thetamax = 360.d0*(1.d0-1.d0/dble(nthetamap))

      dphi = 360.d0/dble(nphimap)
      dtheta = 360.d0/dble(nthetamap)

c     Start with all the angles around the equator

      anglezdeg = dphi
      anglexdeg = dtheta
      x(2)=1.d0
      y(2)=0.d0 
      z(2)=0.d0
      open(66,file='data.dat',status='unknown')
      do i=1, nthetamap
         anglexdeg = dtheta*mod((i-1),2)
         do j=1, nphimap
            call setViewingAngle ! Rotates by degrees
            write(66,*)anglexdeg,angleydeg,anglezdeg,x(2),y(2),z(2)
         end do
      end do
      close(66)
      end subroutine
