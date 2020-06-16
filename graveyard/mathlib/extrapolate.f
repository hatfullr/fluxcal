      subroutine extrapolate(x,y,mx,xarr,my,yarr,zarr,z)
c     Input:
c        x is the x position to extrapolate to
c        y is the y position to extrapolate to
c        mx is the size of xarr
c        xarr is the 1D array containing x positions
c        my is the size of yarr
c        yarr is the 1D array containing y positions
c        zarr is the 2D array containing z values at (xarr(i),yarr(j))
c     Output:
c        z is the z value at (x,y)
      implicit none
      real x,y,z
      integer mx,my,csize
      real zresult_array(mx*my)
      
      integer iopt,kx,ky,nxest,nyest,nx,ny,nmax,lwrk,kwrk,ier
c      parameter(kx=3,ky=3)      ! Polynomial degrees
c      parameter(nxest=mx+kx+1,nyest=my+ky+1) ! Maximum # of knots
c      parameter(nx=10,ny=10)    ! Total # of knots
c      parameter(nmax=max(nx,ny))
c      parameter(lwrk=4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
c     $     my*(ky+1)+max(my,nxest)) ! Size of 'workspace'
c      parameter(kwrk=3+mx+my+nxest+nyest) ! Size of 'workspace'
      real xarr(mx)
      real yarr(my)
      real zarr(mx,my)
      real zarr1D(mx*my)
      real xb,xe,yb,ye
      real s
c      real tx(nmax), ty(nmax) ! Knots of the spline
c      real c((nxest-kx-1)*(nyest-ky-1)) ! Coefficients of the spline
      real fp                   ! Sum of squared residuals
c     real wrk(lwrk)          ! Used as 'workspace'
c     integer iwrk(kwrk)      ! Used as 'workspace'
      
      integer i,j

      iopt = 0                 ! Choose least squares fitting

      kx=3
      ky=3
      
      nx=3
      ny=3
      nmax=max(nx,ny)

      nxest = mx+kx+1
      nyest = my+ky+1

      lwrk=4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
     $     my*(ky+1)+max(my,nxest)
      kwrk=3+mx+my+nxest+nyest

      csize = (nxest-kx-1)*(nyest-ky-1)
c      csize = (nx-kx-1)*(ny-ky-1)
      
      do i=1,mx
         do j=1,my
            zarr1D(my*(i-1)+j) = zarr(i,j)
         end do
      end do
      
      xb = minval(xarr)         ! Boundaries
      xe = maxval(xarr)
      yb = minval(yarr)
      ye = maxval(yarr)
      
      s = 0.                    ! Smoothness factor

      write(*,*) "Trying to get past extrapolate_helper"
      call extrapolate_helper(x,y,z,zarr,
     $     iopt,mx,xarr,my,yarr,zarr1D,xb,xe,yb,ye,
     $     kx,ky,s,nxest,nyest,nx,ny,fp,lwrk,kwrk,ier,
     $     csize,nmax)

      write(*,*) "I made it here"
      end subroutine













      
      
      subroutine extrapolate_helper(x,y,z,zarr,
     $     iopt,mx,xarr,my,yarr,zarr1D,xb,xe,yb,ye,
     $     kx,ky,s,nxest,nyest,nx,ny,fp,lwrk,kwrk,ier,
     $     csize,nmax)
      implicit none
      real x,y,z
      
      integer nmax,csize

      integer iopt,mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk,ier
      real xarr(mx)
      real yarr(my)
      real zarr(mx,my)
      real zarr1D(mx*my)
      real zresult_array(mx*my)
      real xb,xe,yb,ye
      real s
      real fp                   ! Sum of squared residuals
      real tx(nmax), ty(nmax)
      real c(csize)
      real wrk(lwrk)
      real iwrk(kwrk)

      integer i,j

      do i=1,kx+1
         tx(i) = xb
         tx(nx) = xe
      end do

      call regrid(iopt,mx,xarr,my,yarr,zarr1D,xb,xe,yb,ye,kx,ky,s,
     $     nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)

      if ( ier .gt. 0 ) then
         write(*,*) "iopt = ",iopt
         write(*,*) "kx = ",kx
         write(*,*) "ky = ",ky
         write(*,*) "mx = ",mx
         write(*,*) "my = ",my
         write(*,*) "nxest = ",nxest,2*kx+2
         write(*,*) "nyest = ",nyest,2*ky+2
         write(*,*) "kwrk = ",kwrk,3+mx+my+nxest+nyest
         write(*,*) "lwrk = ",lwrk,4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+
     $        mx*(kx+1)+my*(ky+1) + max(my,nxest)
         write(*,*) "xb = ",xb
         write(*,*) "xe = ",xe
         write(*,*) "yb = ",yb
         write(*,*) "ye = ",ye
         write(*,*) "ier = ",ier

         if (ier.eq.10) then
            write(*,*) -1.le.iopt .and. iopt.le.1," -1<=iopt<=1"
            write(*,*) 1.le.kx," 1<=kx"
            write(*,*) ky.le.5," ky<=5"
            write(*,*) mx.gt.kx," mx>kx"
            write(*,*) my.gt.ky," my>ky"
            write(*,*) nxest.ge.2*kx+2," nxest>=2*kx+2"
            write(*,*) nyest.ge.2*ky+2," nyest>=2*ky+2"
            write(*,*) kwrk.ge.3+mx+my+nxest+nyest," kwrk>=3+mx+my"//
     $           "+nxest+nyest"
            write(*,*) lwrk.ge.4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+
     $           mx*(kx+1)+my*(ky+1) + max(my,nxest)," lwrk>=4+nxest*"//
     $           "(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+my*(ky+1)"//
     $           "+max(my,nxest)"
            
            do i=2,mx
               if (.not.(xb.le.xarr(i-1) .and. xarr(i-1).lt.xarr(i).and.
     $              xarr(i).le.xe)) then
                  write(*,*) "F","  xb<=x(i-1)<x(i)<=xe"
                  exit
               end if
            end do
            
            do j=2,my
               if (.not.(yb.le.yarr(j-1) .and. yarr(j-1).lt.yarr(j).and.
     $              yarr(j).le.ye)) then
                  write(*,*) "F","  yb<=y(j-1)<y(j)<=ye"
                  exit
               end if
            end do

            if (iopt.ge.0) then
               write(*,*) s.lt.0," s<0"
               if ( s.ge.0 ) then
                  write(*,*) s.eq.0," s=0"
                  write(*,*) nxest.ge.mx+kx+1," nxest>=mx+kx+1"
                  write(*,*) nyest.ge.my+ky+1," nyest>=my+ky+1"
               end if
            end if
         end if
         
         
         write(*,*) "ier > 0. Check regrid.f for error codes."
         error stop "extrapolate.f"
      end if

      
      call bispev(tx,nx,ty,ny,c,kx,ky,xarr,mx,yarr,my,
     $     zresult_array,wrk,lwrk,iwrk,kwrk,ier)

      if ( ier .gt. 0 ) then
         write(*,*) "ier = ",ier
         write(*,*) "ier returned non-zero. Check bispev.f for error "//
     $        "codes."
         error stop "extrapolate.f"
      end if

      write(*,*) "I am finishing in extrapolate_helper"

      end subroutine



