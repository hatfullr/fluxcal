      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmaximum,nok,nbad,
     $     derivs,rkqs)
      implicit none
      INTEGER nbad,nok,nvar,KMAXX,NMAX,maxstp
      real*8 eps,h1,hmaximum,x1,x2,ystart(nvar) ! ,TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200) !,TINY=1.d-30)
      INTEGER i,kmax,kount,kount1 !,kount2
      real*8 dxsav,h,hdid,hnext,xint,xsav,dydx(NMAX),xp(KMAXX),
     *yint(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      common /kounts/kount1 !,kount2
      integer numfiltersmax
      parameter (numfiltersmax=15)
      real*8 yscalfactor(numfiltersmax)
      common /yscales/ yscalfactor
      real*8 step1, step2, step3, step4
      common/steps/ step1,step2,step3,step4
      real*8 taulimit,posx,posy
      common/inputstuff/ taulimit,posx,posy
c      logical get_integration_at_pos,get_integration_at_all_pos
c      common/getint/ get_integration_at_pos,get_integration_at_all_pos
      real*8 rhocgs,xhp,gcgs,pcgs,tcgs,ucgs
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
      real*8 opacit,opacit1,opacit2
      real*8 fcondense,nn
      real*8 Rform
      common/dust/ Rform
      integer intout
      common/intoutput/ intout
      real*8 getLocalAngle

      integer mymax
      parameter(mymax=420000)
      real*8 x(mymax),y(mymax),z(mymax),am(mymax),hp(mymax),rho(mymax),
     $      vx(mymax),vy(mymax),vz(mymax),a(mymax),wmeanmolecular(mymax),
     $      localg(mymax),metal(mymax),tempp(mymax),pp(mymax)
      common/part/ x,y,z,am,hp,rho,vx,vy,vz,a,wmeanmolecular,localg,
     $      metal,tempp,pp
      
      integer n,ip,lp
      real*8 hmin,hmax,dr2,dx2,dy2,dz2,bigd,maxbigd
      common/intpar/ hmin,hmax,n

      external usetable,usetabledust
      external usetable2,usetabledust2
      
      integer lastpart
      common/lastparticle/ lastpart

      real*8 tau_thick_integrator
      common/tauthick/ tau_thick_integrator

      real*8 rayout1(10,maxstp)
      integer rayout2(maxstp),nstp
      common/rayout/ rayout1,rayout2,nstp

      logical printIntegrationSteps
      common/intatpos/ printIntegrationSteps
      
      real*8 munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $     runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $     rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $     Lunit_out
      common/units/ munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $     runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $     rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $     Lunit_out

      real*8 xpos,ypos
      common/taugrid/xpos,ypos
      
      

c      logical resolved
c      common /resolvedboolean/ resolved
      xint=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=1
      kount1=0
c      kount2=0
      do 11 i=1,nvar
        yint(i)=ystart(i)
 11   continue

      do 16 nstp=1,MAXSTP
c        print *,'looping',x,y(1)
        call derivs(xint,yint,dydx)
        do 12 i=1,nvar
           yscal(i)=abs(yint(i))+abs(h*dydx(i)) ! +TINY
12      continue

c     yscal(1)=yscal(1)+1d3
        yscal(1)=yscal(1)+step1 ! Step size near tau=1
        yscal(2)=yscal(2)+step2
        yscal(3)=yscal(3)+step3
        yscal(4)=yscal(4)+step4

        do i=5, nvar
           yscal(i)=yscal(i)
     $          +yscalfactor(nvar-4)*exp(-yint(1))
c     $          +abs(yint(i))
        enddo

c     I want to make sure the penultimate step (the step right before the photosphere or thermalization depth) is always saved
        xp(kount)=xint
        do i=1,nvar
           yp(i,kount)=yint(i)
        enddo

        if((xint+h-x2)*(xint+h-x1).gt.0.d0) h=x2-xint
c        write(*,*) xpos/runit,ypos/runit
        call rkqs(yint,dydx,nvar,xint,h,eps,yscal,hdid,hnext,derivs)

        if(printIntegrationSteps) then
           call getLocalQuantities(xpos,ypos,xint)
           call getOpacitySub(dble(t6*1d6),rhocgs,yint(1),opacit)
           write(intout,'(12ES22.14,I22)')xpos/runit_out,ypos/runit_out,
     $          xint/runit_out,xhp/runit_out,rhocgs/rhounit_out,
     $          ucgs/Eunit_out*munit_out,xh/muunit_out,
     $          gcgs/gunit_out,tcgs/tempunit_out,pcgs/punit_out,
     $          opacit,yint(1),lastpart
        end if
        
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif

c     yint(1) = optical depth tau
c     yint(3) = number of smoothing lengths into the system
c     yint(4) Integral[ T^4 exp(-tau'), {tau', 0, tau(1)}], where temperature T=T(tau')

c        if(yint(1).lt.1 .and. yint(4).ge.(10d0**3.70492d0)**4 .and. yint(3).lt.2)
c        if(yint(4).ge.(10d0**3.70492d0)**4 .and. yint(3).lt.4)
c     $       then
c           resolved=.false.
c           return 
c        endif

c        print *,'ok,bad',nok,nbad,h,hdid,yint(1)
c        if((yint(1).ge.1.d0 .and. kount1.eq.0) .or.
c     $       (yint(1)*yint(2).ge.1.d0/3.d0 .and. kount2.eq.0))then
c     write(*,*) yint(1)
        if(yint(1).ge.tau_thick_integrator .and. kount1.eq.0) then
c           print *,'ok,bad',nok,nbad,h,hdid,yint(1)
           kount=kount+1
           if(yint(1).ge.tau_thick_integrator .and. kount1.eq.0)
     $          kount1=kount
c           if(yint(1)*yint(2).ge.1.d0/3.d0 .and. kount2.eq.0) kount2=kount
           xp(kount)=xint
           do i=1,nvar
              yp(i,kount)=yint(i)
           enddo
           kount=kount+1
        endif
        if((xint-x2)*(x2-x1).ge.0.d0 .or.
c     $       yint(1).gt.1000)      ! integrate only to tau=1000
c     Don't waste time probing the tau>1 region that isn't needed
     $       yint(1).ge.taulimit)  ! integrate only to tau=taulimit
     $       then
          do i=1,nvar
             ystart(i)=yint(i)
          enddo
c     kount=kount-1       xxxxx        Maybe this line should still be here????
          return
       endif
       
        if(abs(hnext).gt.abs(hmaximum)) then
           hnext=sign(hmaximum,hnext)
c           print *,'limiting step size to',hnext
        endif

        h=hnext
16    continue
       print *, 'too many steps in odeint'
c     pause 'too many steps in odeint'
      return
      END
