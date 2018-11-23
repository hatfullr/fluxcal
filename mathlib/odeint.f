      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmaximum,nok,nbad,
     $     derivs,rkqs)
      implicit none
      INTEGER nbad,nok,nvar,KMAXX,NMAX,maxstp
      real*8 eps,h1,hmaximum,x1,x2,ystart(nvar) ! ,TINY
      EXTERNAL derivs,rkqs
      PARAMETER (NMAX=50,KMAXX=200,maxstp=10000) !,TINY=1.d-30)
      INTEGER i,kmax,kount,kount1 !,kount2
      real*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      common /kounts/kount1 !,kount2
      integer numfiltersmax
      parameter (numfiltersmax=15)
      real*8 yscalfactor(numfiltersmax)
      common /yscales/ yscalfactor
      real*8 step1, step2, step3, step4
      common/steps/ step1,step2,step3,step4
      logical rossonly,get_closest_taus
      common/opacitytype/ rossonly,get_closest_taus
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

      external usetable,usetabledust
      external usetable2,usetabledust2
      
      integer lastpart
      common/lastparticle/ lastpart

      real*8 intz
      common/intedz/ intz

      real*8 tau_thick
      common/tauthick/ tau_thick

      real*8 rayout1(10,maxstp)
      integer rayout2(maxstp),nstp
      common/rayout/ rayout1,rayout2,nstp

      logical dointatpos,dointatallpos
      common/intatpos/ dointatpos,dointatallpos
      
      real*8 munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $      runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $      rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $      Lunit_out
      common/units/ munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $      runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $      rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $      Lunit_out
      

c      logical resolved
c      common /resolvedboolean/ resolved
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=1
      kount1=0
c      kount2=0
      do 11 i=1,nvar
        y(i)=ystart(i)
 11   continue


      do 16 nstp=1,MAXSTP
c        print *,'looping',x,y(1)
        call derivs(x,y,dydx)
        do 12 i=1,nvar
           yscal(i)=abs(y(i))+abs(h*dydx(i)) ! +TINY
12      continue

c     yscal(1)=yscal(1)+1d3
        yscal(1)=yscal(1)+step1 ! Step size near tau=1
        yscal(2)=yscal(2)+step2
        yscal(3)=yscal(3)+step3
        yscal(4)=yscal(4)+step4

        do i=5, nvar
           yscal(i)=yscal(i)
     $          +yscalfactor(nvar-4)*exp(-y(1))
c     $          +abs(y(i))
        enddo

c     I want to make sure the penultimate step (the step right before the photosphere or thermalization depth) is always saved
        xp(kount)=x
        do i=1,nvar
           yp(i,kount)=y(i)
        enddo

        if((x+h-x2)*(x+h-x1).gt.0.d0) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        intz = x
        if(dointatpos) then
 800       format(10ES22.14,I22)
           call getLocalQuantities(posx,posy,x)
           call getOpacitySub(posx,posy,x,dble(t6*1d6),rhocgs,y(1),
     $          Rform,opacit)
           write(intout,800) x/runit_out,xhp/runit_out,
     $          rhocgs/rhounit_out,ucgs/Eunit_out*munit_out,
     $          gcgs/gunit_out,xh/muunit_out,pcgs/punit_out,
     $          tcgs/tempunit_out,opacit,y(1),lastpart
        end if

        if(dointatallpos) then
           call getLocalQuantities(posx,posy,x)
           call getOpacitySub(posx,posy,x,dble(t6*1d6),rhocgs,y(1),
     $          Rform,opacit)
           rayout1(1,nstp)  = x/runit_out
           rayout1(2,nstp)  = xhp/runit_out
           rayout1(3,nstp)  = rhocgs/rhounit_out
           rayout1(4,nstp)  = ucgs/Eunit_out*munit_out
           rayout1(5,nstp)  = gcgs/gunit_out
           rayout1(6,nstp)  = xh/muunit_out
           rayout1(7,nstp)  = pcgs/punit_out
           rayout1(8,nstp)  = tcgs/tempunit_out
           rayout1(9,nstp)  = opacit
           rayout1(10,nstp) = y(1)
           rayout2(nstp)    = lastpart
        end if
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif

c     y(1) = optical depth tau
c     y(3) = number of smoothing lengths into the system
c     y(4) Integral[ T^4 exp(-tau'), {tau', 0, tau(1)}], where temperature T=T(tau')

c        if(y(1).lt.1 .and. y(4).ge.(10d0**3.70492d0)**4 .and. y(3).lt.2)
c        if(y(4).ge.(10d0**3.70492d0)**4 .and. y(3).lt.4)
c     $       then
c           resolved=.false.
c           return 
c        endif

c        print *,'ok,bad',nok,nbad,h,hdid,y(1)
c        if((y(1).ge.1.d0 .and. kount1.eq.0) .or.
c     $       (y(1)*y(2).ge.1.d0/3.d0 .and. kount2.eq.0))then
c     write(*,*) y(1)
        if(y(1).ge.tau_thick .and. kount1.eq.0) then
c           print *,'ok,bad',nok,nbad,h,hdid,y(1)
           kount=kount+1
           if(y(1).ge.tau_thick .and. kount1.eq.0) kount1=kount
c           if(y(1)*y(2).ge.1.d0/3.d0 .and. kount2.eq.0) kount2=kount
           xp(kount)=x
           do i=1,nvar
              yp(i,kount)=y(i)
           enddo
           kount=kount+1
        endif
        if((x-x2)*(x2-x1).ge.0.d0 .or.
c     $       y(1).gt.1000)      ! integrate only to tau=1000
c     Don't waste time probing the tau>1 region that isn't needed
     $       y(1).ge.taulimit)  ! integrate only to tau=taulimit
     $       then
          do i=1,nvar
             ystart(i)=y(i)
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
