      subroutine tabulinit
c     initialize look-up tables:
      include 'flux_cal.h'
      real*8 u,w,w33,wc4
      real*8 rlow,rhigh
      integer i
      
c     compute normalization constant:
      ctab=dble(ntab-1)/4.d0
c     compute tabulations of w(u):
      do i=1,ntab
         u = 2d0*(dble(i-1)/dble(ntab-1))**0.5d0
         if(u.lt.0 .or. u.gt.2.d0) then
            write(*,*) "u=r/h is out of range"
            write(*,*) "u=",u
            error stop "kernels.f line 14"
         endif

         if (nkernel.eq.0) then
c     cubic spline
            wtab(i)=w(u)
         elseif (nkernel.eq.1) then
c     Wendland 3,3 kernel
            wtab(i)=w33(u)
         elseif (nkernel.eq.2) then
c     Wendland C4 kernel
            wtab(i)=wc4(u)
         else
            write(*,*) 'nkernel must equal 0, 1, or 2'
            stop
         endif
      enddo

      rlow=0.d0
      do i=1,ntab
         rhigh=(i-0.5d0)**0.5d0
         taucoef=taucoef+wtab(i)*(rhigh-rlow)
         rlow=rhigh
      enddo
      taucoef=taucoef/ctab**0.5d0 !  where ctab=dble(ntab-1)/4.d0
      
      return
      end

***********************************************************************
      real*8 function w(u)
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W is the usual cubic spline kernel with compact support 2h
      implicit none
      real*8 u
      if (u < 0) then
         w = 1
      else if (u.lt.1.d0) then
         w=1.d0+u**2.d0*(-1.5d0+0.75d0*u)
      else if (u.lt.2.d0) then
         w=0.25d0*(2.d0-u)**3
      else
         w=0.d0
      endif
      w=w/3.1415926535897932384626d0

      return
      end
************************************************************************
c     Wendland 33 kernel
      real*8 function w33(u)  ! NOTE... function of u=r/h
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W=Wendland 3,3 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =1365.d0/(512.d0*3.1415926535897932384626d0)) ! Note the 512
!     We're using 512 instead of 64 because the output is W*h^3 and not W*(2h)^3
      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      w33 = sigW*(1.d0-q)**8.d0*(32.d0*q**3.d0+25.d0*q**2.d0+8.d0*q+1.d0)

      return
      end
 

************************************************************************
c     Wendland c4 kernel
      real*8 function wc4(u)  ! NOTE... function of u=r/h
c     input: u=r/h from 0 to 2 
c     output: the dimensionless quantity W*h^3, where h=smoothing length
c     and W=Wendland C4 kernel scaled to have a compact support 2h
      implicit none
      real*8 q,sigW,u
      parameter(sigW =495.d0/(256.d0*3.1415926535897932384626d0)) ! Note the 256 instead of 2 (this is because the output is W*h^3 and not W*(2h)^3

      q=u/2.d0                  !  q=r/(2h) is the fractional radius
      wc4 = sigW*(1.d0-q)**6.d0*(35.d0/3d0*q**2.d0+6.d0*q+1.d0)
      return
      end
