      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER n,NMAX
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs, rkck
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      double precision errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),
     &SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,
     &ERRCON=1.89d-4)
      real*8 ysav
      common/ysave/ysav
c      integer ierrmax
c      integer ithatgiveserrmax

      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)

      errmax=0.d0

c      do 11 i=1,min(n,2)        ! Using min(n,2) because we don't want
c                                ! to let the i=3 integration (which counts
c                                ! the number of smoothing lengths) to
c                                ! affect the stepsize
c         errmax=max(errmax,abs(yerr(i)/yscal(i)))
c 11   continue

c      print *, yerr(1),yscal(1),yerr(4),yscal(4),
c     $     yerr(1)/yscal(1),yerr(4)/yscal(4)

      errmax=max(errmax,abs(yerr(1)/yscal(1))) ! necessary to get small enough step sizes in some cases
c      ithatgiveserrmax = 1
      do i=3,n
c         if (abs(yerr(i)/yscal(i)).gt.errmax)then
c            ithatgiveserrmax = i
c         endif
         errmax=max(errmax,abs(yerr(i)/yscal(i))) ! necessary to get small enough step sizes in some cases
         
      enddo
c      write(25,*) ithatgiveserrmax,yerr(ithatgiveserrmax),
c     $     yscal(ithatgiveserrmax)
c      do 11 i=1,n
c         errmax=max(errmax,abs(yerr(i)/yscal(i)))
c11    continue

      errmax=errmax/eps

      if(errmax.gt.1.d0)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1d0*abs(h)),h)
        xnew=x+h

        if(xnew.eq.x)then
           write(6,*) 'stepsize underflow in rkqs'
           write(6,*) 'z    = ',x
           write(6,*) 'znew = ',xnew
           write(6,*) 'tau  = ',y(1)
           write(6,*) 'yscal(1)=',yscal(1)
           write(6,*) 'yscal(2)=',yscal(2)
           write(6,*) 'yscal(3)=',yscal(3)
           write(6,*) 'yscal(4)=',yscal(4)
           write(6,*) "Try increasing some of your step variables"
           error stop "rkqs.f"
        endif
        goto 1

      else

        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.d0*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
