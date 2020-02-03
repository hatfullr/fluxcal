      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      real*8 xpos,ypos
      common/taugrid/xpos,ypos
      real*8 rhocgs,xhp,gcgs,pcgs,tcgs,ucgs
      real xh,t6
      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
      real*8 munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $     runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $     rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $     Lunit_out
      common/units/ munit,runit,tunit,vunit,Eunit,rhounit,muunit,gunit,
     $     runit_out,munit_out,tunit_out,vunit_out,Eunit_out,
     $     rhounit_out,muunit_out,gunit_out,tempunit_out,punit_out,
     $     Lunit_out

      integer j,k
      real*8 zpos,dz,rhostart,factor,opacstart,opacity
      
      INTEGER n,intNMAX
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs, rkck
      PARAMETER (intNMAX=50)
C     U    USES derivs,rkck
      INTEGER i
      double precision errmax,h,htemp,xnew,yerr(intNMAX),ytemp(intNMAX),
     &     SAFETY,PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,
     &     ERRCON=1.89d-4)
      real*8 ysav
      common/ysave/ysav
c     integer ierrmax
c     integer ithatgiveserrmax

      h=htry
 1    call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)

      errmax=0.d0

c     do 11 i=1,min(n,2)        ! Using min(n,2) because we don't want
c     ! to let the i=3 integration (which counts
c     ! the number of smoothing lengths) to
c     ! affect the stepsize
c     errmax=max(errmax,abs(yerr(i)/yscal(i)))
c     11   continue

c     print *, yerr(1),yscal(1),yerr(4),yscal(4),
c     $     yerr(1)/yscal(1),yerr(4)/yscal(4)

      errmax=max(errmax,abs(yerr(1)/yscal(1))) ! necessary to get small enough step sizes in some cases
c     ithatgiveserrmax = 1
      do i=3,n
c     if (abs(yerr(i)/yscal(i)).gt.errmax)then
c     ithatgiveserrmax = i
c     endif
         errmax=max(errmax,abs(yerr(i)/yscal(i))) ! necessary to get small enough step sizes in some cases
         
      enddo
c     write(25,*) ithatgiveserrmax,yerr(ithatgiveserrmax),
c     $     yscal(ithatgiveserrmax)
c     do 11 i=1,n
c     errmax=max(errmax,abs(yerr(i)/yscal(i)))
c11   continue

      errmax=errmax/eps

      if(errmax.gt.1.d0)then
         htemp=SAFETY*h*(errmax**PSHRNK)
         ! This occurs 
         if(x+sign(max(abs(htemp),0.1d0*abs(h)),h).eq.x)then
            write(6,*) 'stepsize underflow in rkqs'
            write(6,*) 'xpos, ypos = ',xpos/runit_out,ypos/runit_out
            write(6,*) 'Before step (z):'
            write(6,*) '   z    = ',x/runit_out
            call getLocalQuantities(xpos,ypos,x)
            write(6,*) '   rho = ',rhocgs/rhounit_out
            write(6,*) '   T = ',dble(t6*1d6)/tempunit_out
            depth = x
            opacstart = getOpacity(dble(t6*1d6),rhocgs,xh)
            write(6,*) '   opacity = ',opacstart/kunit_out
            write(6,*) '   xhp = ',xhp/runit_out
            write(6,*) '   h = ',h/runit_out

            rhostart = rhocgs
            
            write(6,*)
            
            write(6,*) 'After step (z+h):'
            write(6,*) '   z+h = ',(x+h)/runit_out
            call getLocalQuantities(xpos,ypos,x+sign(max(abs(htemp),
     $           0.1d0*abs(h)),h))
            depth = x+sign(max(abs(htemp),0.1d0*abs(h)),h)
            write(6,*) '   rho = ',rhocgs/rhounit_out
            write(6,*) '   T = ',dble(t6*1d6)/tempunit_out
            write(6,*) '   opacity = ',getOpacity(dble(t6*1d6),rhocgs,
     $           xh)/kunit_out
            write(6,*) '   xhp = ',xhp/runit_out
            write(6,*) '   h = sign(max(abs(htemp),0.1d0*abs(h)),h)'
            write(6,*) '   abs(htemp) = ',abs(htemp)/runit_out
            write(6,*) '   0.1d0*abs(h) = ',0.1d0*abs(h)/runit_out
            write(6,*) '   max(abs(htemp),0.1d0*abs(h)) = ',
     $           max(abs(htemp),0.1d0*abs(h))/runit_out
            write(6,*) '   h = ',sign(max(abs(htemp),0.1d0*abs(h)),h)/
     $           runit_out

            write(6,*)

            ! Write out where the nearest material is located
            write(6,*) 'Collecting data on where the nearest material'//
     $           ' is located...'
            write(6,*) ''
            zpos = x
            ! We know we are always going from +z to -z.
            dz = -1.d0*abs(h)
            factor = 2.d0
            findrho : do j=1,1000
               zpos = zpos + dz
               call getLocalQuantities(xpos,ypos,zpos)

               ! Double the step size in the -z direction until we find
               ! a different density than what we started out with.
               if ( rhocgs.eq.rhostart ) then
                  dz = factor * dz
               else             ! We may have overshot
                  ! Search in the region between zpos and zpos-dz/2
                  dz = (dz/factor)/1000.d0
                  do k=1,1000

                     ! If we have found where density started changing
                     if ( rhocgs.eq.rhostart ) then
                        exit findrho
                     end if

                     zpos = zpos - dz
                     call getLocalQuantities(xpos,ypos,zpos)
                     
                  end do
               end if
            end do findrho

            if ( rhocgs.ne.rhostart ) then
               write(6,*) 'Density begins changing at:'
               write(6,*) '   z = ',zpos/runit_out
               write(6,*) '   rho = ',rhocgs/rhounit_out
               write(6,*) '   T = ',dble(t6*1d6)/tempunit_out
               depth = zpos
               write(6,*) '   opacity = ',getOpacity(dble(t6*1d6),
     $              rhocgs,xh)/kunit_out
               write(6,*) '   xhp = ',xhp/runit_out
            else
               write(6,*) 'Could not find where density begins '//
     $              'changing'
            end if

            write(6,*)


            ! Write out where opacity starts changing

            zpos = x
            ! We know we are always going from +z to -z.
            dz = -1.d0*abs(h)
            factor = 2.d0
            findopac : do j=1,1000
               zpos = zpos + dz
               call getLocalQuantities(xpos,ypos,zpos)
               depth = zpos
               opacity = getOpacity(dble(t6*1d6),rhocgs,xh)
               
               ! Double the step size in the -z direction until we find
               ! a different density than what we started out with.
               if ( opacity.eq.opacstart ) then
                  dz = factor * dz
               else             ! We may have overshot
                  ! Search in the region between zpos and zpos-dz/2
                  dz = (dz/factor)/1000.d0
                  do k=1,1000

                     ! If we have found where density started changing
                     if ( opacity.eq.opacstart ) then
                        exit findopac
                     end if

                     zpos = zpos - dz
                     call getLocalQuantities(xpos,ypos,zpos)
                     depth = zpos - dz
                     opacity = getOpacity(dble(t6*1d6),rhocgs,xh)
                     
                  end do
               end if
            end do findopac

            if ( opacity.ne.opacstart ) then
               write(6,*) 'Opacity begins changing at:'
               write(6,*) '   z = ',zpos/runit_out
               write(6,*) '   rho = ',rhocgs/rhounit_out
               write(6,*) '   T = ',dble(t6*1d6)/tempunit_out
               write(6,*) '   opacity = ',opacity/kunit_out
               write(6,*) '   xhp = ',xhp/runit_out
            else
               write(6,*) 'Could not find where opacity begins '//
     $              'changing'
            end if
            
            write(6,*)
            
            write(6,*) 'Other variables in rkqs.f:'
            write(6,*) '   tau  = ',y(1)
            write(6,*) '   step1 = ',yscal(1)
            write(6,*) '   step2 = ',yscal(2)
            write(6,*) '   step3 = ',yscal(3)
            write(6,*) '   step4 = ',yscal(4)
            
            write(6,*)
            
            write(6,*) "This error is a result of the integrator "//
     $           "taking a step forward in the z direction so small"//
     $           " that the new z position equals the old z"//
     $           " position."
            write(6,*) "You can try increasing some of your step "//
     $           "variables"
            error stop "rkqs.f"
         else
            h=sign(max(abs(htemp),0.1d0*abs(h)),h)
            xnew=x+h
         end if

         

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
