      function getOpacity(tem,rhocgs,xhy) ! real*8
      use bilinear_interpolation
      include 'flux_cal.h'
c     INPUT:
c     tem    - temperature
c     rhocgs - density in cgs
c     xhy    - hydrogen fraction
      
      real*8 tem,rhocgs,xhy,logR,logtem,logrho
      real*8 opacity,tempopacity
      integer i,j
      real*8 calculate_table_opacity
      
      external cop
      real*8 eD_R, eG_R, eD_P, eG_P
      DIMENSION eD_R(5,6), eD_P(5,6)
      DIMENSION eG_R(71,71), eG_P(71,71)
      common/opac/eD_R,eG_R,eD_P,eG_P
      logical exist
      
      real*8 XPOS, YPOS, TMAX, TMIN, anglex, angley, anglez
      COMMON/TAUGRID/XPOS,YPOS
      real*8 nnbfrac,nnb,nnbdr2,hp2
      integer k,nnbr,nnbi,nnbitot
      real*8 neighborslist(1000)
      integer writefile
      common/writefilee/writefile
      real*8 dxy,particle_pos,mydr2,minmydr2,dist
      integer closestparticle

      integer lastpart
      common/lastparticle/ lastpart
      integer index
      real*8 r2,wpc,rhoxyz

c      if ( rhocgs.lt.1.d-18 ) then
c         r2=(x(lastpart)-xpos)**2.d0+(y(lastpart)-ypos)**2.d0
c     $        +(z(lastpart)-depth)**2.d0
c         index=int(ctab*r2/hp(lastpart)**2.d0)+1
c         wpc=wtab(index)/hp(lastpart)**3.d0
c         rhoxyz = am(lastpart)*wpc
c         write(*,*) "xpos, ypos, zpos = ",xpos/runit_out,ypos/runit_out,
c     $        depth/runit_out
c         write(*,*) "x(i), y(i), z(i), 2*hp(i) = ",
c     $        x(lastpart)/runit_out,
c     $        y(lastpart)/runit_out,z(lastpart)/runit_out,
c     $        2.d0*hp(lastpart)/runit_out
c         write(*,*) "rho(i) = ",rho(lastpart)
c         write(*,*) "rhoxyz = ",rhoxyz
c         write(*,*) "sqrt(r2), sqrt(r2)/(2*hp(i)) = ",
c     $        (r2**0.5d0)/runit_out,
c     $        (r2**0.5d0)/(2.d0*hp(lastpart))
c         write(*,*) "rhocgs, tem, lastpart = ",rhocgs,tem,lastpart
c         getOpacity = -1.d30
c         return
c      end if
      
      
c      write(*,*) "rhocgs, tem = ",rhocgs,tem
      if ( rhocgs .le. opacity_rho_cutoff ) then
         getOpacity = 0.d0
         return
      end if

      if ( tem.eq.0.d0 ) then
         write(*,*) "(rho,T) = (",rhocgs,tem,")"
         write(*,*) "Cannot compute the opacity of a T = 0 K fluid."
         error stop "getOpacity.f"
      end if

      logtem = log10(tem)
      logrho = log10(rhocgs)
      logR = logrho - 3.d0*logtem + 18.d0

c      write(*,*) "logR, logtem = ",logR,logtem

      ! We do not care about these temperature ranges. There is no
      ! way to calculate the opacity there in a way that makes sense,
      ! so just take the last known opacity as constant throughout the
      ! temperature range.
      !
      ! We are assuming here that we will always have our Planck opacityfile
      if ( logtem.lt.logTmins(1) ) then
         getOpacity = calculate_table_opacity(
     $        10.d0**(0.5d0*(logRmaxs(1)+logRmins(1))),
     $        10.d0**logTmins(1))
c         write(*,*) "log10(opacity) = ",log10(opacity)
      else
      
         opacity = calculate_table_opacity(10.d0**logR,tem)

         if ( opacity.eq.-1.d30 ) then
            call extrapolate_opacity(tem,rhocgs,xhy,opacity)
         end if


      
         getOpacity = opacity
      end if
      
c$$$c     Remove the following later:
c$$$      write(writefile,'(3E15.7,I10)'),log10(rhocgs),log10(tem),
c$$$     $     log10(getOpacity),n
c$$$      
c$$$      if ( n.gt.0 ) then
c$$$
c$$$         ! Find distance to closest contributing particle
c$$$         minmydr2 = 1.d30
c$$$         closestparticle = 0
c$$$         do i=1,n
c$$$            mydr2=(xpos-x(i))**2.d0+(ypos-y(i))**2.d0+(depth-z(i))**2.d0
c$$$            if ( mydr2.lt.minmydr2 ) then
c$$$               minmydr2 = mydr2
c$$$               closestparticle = i
c$$$            end if
c$$$         end do
c$$$
c$$$         dist = minmydr2**0.5d0 / (2.d0*hp(closestparticle))
c$$$
c$$$         write(writefile,'(3E15.7)') rhocgs,tem,dist
         
c$$$c     write(*,*) "xpos, ypos, depth = ",xpos,ypos,depth
c$$$         
c$$$         
c$$$         nnbr = 0
c$$$         nnbi = 0
c$$$         nnbitot = 0
c$$$         do i=1,n
c$$$            hp2 = (2.d0*hp(i))**2.d0
c$$$            mydr2=(xpos-x(i))**2.d0+(ypos-y(i))**2.d0+(depth-z(i))**2.d0
c$$$c     write(*,*) "mydr2.le.hp2, mydr2, hp2",mydr2.le.hp2,mydr2,hp2
c$$$            if(mydr2.le.hp2) then
c$$$               nnbr = nnbr + 1
c$$$               nnbitot = nnbitot + neighbors(i)
c$$$c     if ( neighbors(i).eq.0 ) then
c$$$c     write(*,*) "Problem i,n = ",i,n
c$$$c     error stop "getOpacity.f"
c$$$c     end if
c$$$            end if
c$$$         end do
c$$$         
c$$$         nnbfrac = float(nnbr) / (float(nnbitot)/float(nnbr))
c$$$c     write(*,*) "nnbr,nnbitot,<nnbi>,nnbfrac = ",
c$$$c     $        nnbr,nnbitot,
c$$$c     $        float(nnbr)/(float(nnbitot)/float(nnbr)),nnbfrac
c$$$         if ( nnbr.eq.0 .or. nnbitot.eq.0 ) then
c$$$            write(*,*) "This should never happen."
c$$$            call getLocalQuantities(xpos,ypos,depth)
c$$$            write(*,*) "Called getLocalQuantities"
c$$$            write(*,*) "xpos,ypos,depth = ",xpos/runit_out,
c$$$     $           ypos/runit_out,depth/runit_out
c$$$            write(*,*) "rho, T = ",rhocgs,tem
c$$$            write(*,*) "nnbr,nnbitot,<nnbi>,nnbfrac = ",
c$$$     $           nnbr,nnbitot,
c$$$     $           float(nnbr)/(float(nnbitot)/float(nnbr)),nnbfrac
c$$$c     error stop "getOpacity.f"
c$$$         end if
c$$$         write(writefile,'(3E15.7)') rhocgs,tem,nnbfrac
c$$$c     write(*,*)
c$$$
c$$$      end if
c$$$  ------------------------------------
         
      return
      end function























      subroutine find_closestT(logR,logtem,logRbound_arr,
     $     logTbound_arr,closestT)
c     1 = Lower-left
c     2 = Upper-left
c     3 = Lower-right
c     4 = Upper-right
      include 'flux_cal.h'
      real*8 logR,logtem,closestT,b,m
      real*8 logRbound_arr(4)
      real*8 logTbound_arr(4)

      ! If we are underneath the opacityfile
      if ( logtem.lt.max(logTbound_arr(1),logTbound_arr(3)) ) then
         ! If we are on the right side of the opacityfile
         if ( logR.gt.min(logRbound_arr(3),
     $        logRbound_arr(4)) ) then
            ! If the opacityfile on the right side is not rectangular
            ! and has a positive slope
            if ( logRbound_arr(4).gt.logRbound_arr(3) ) then
               m = (logTbound_arr(3)-logTbound_arr(4))/
     $              (logRbound_arr(3)-logRbound_arr(4))
               b = logTbound_arr(3) - m*logRbound_arr(3)
               closestT = m*logR + b
                     
            ! The opacityfile is either rectangular or has a negative
            ! slope.
            else
               ! If the bottom border is not rectangular
               if ( logTbound_arr(1).ne.logTbound_arr(3) )then
                  m = (logTbound_arr(1)-logTbound_arr(3))/
     $                 (logRbound_arr(1)-logRbound_arr(3))
                  b = logTbound_arr(1) - m*logRbound_arr(1)
                  closestT = m*logR + b
               ! Bottom border is rectangular
               else
                  closestT = logTbound_arr(1)
               end if
            end if

         ! If we are on the left side of the opacityfile
         else if(logR.lt.max(logRbound_arr(1),
     $           logRbound_arr(2)) ) then
            ! If the opacityfile on the left side is not rectangular
            ! and has a negative slope
            if ( logRbound_arr(1).gt.logRbound_arr(2) ) then
               m = (logTbound_arr(1)-logTbound_arr(2))/
     $              (logRbound_arr(1)-logRbound_arr(2))
               b = logTbound_arr(1) - m*logRbound_arr(1)
               closestT = m*logR + b

            ! The opacityfile is either rectangular or has a positive
            ! slope.
            else
               ! If the bottom border is not rectangular
               if ( logTbound_arr(1).ne.logTbound_arr(3) )then
                  m = (logTbound_arr(1)-logTbound_arr(3))/
     $                 (logRbound_arr(1)-logRbound_arr(3))
                  b = logTbound_arr(1) - m*logRbound_arr(1)
                  closestT = m*logR + b
               ! Bottom border is rectangular
               else
                  closestT = logTbound_arr(1)
               end if
            end if
                  
         ! If we are simply underneath the opacityfile
         else
            ! If the bottom border is not rectangular
            if ( logTbound_arr(1).ne.logTbound_arr(3) )then
               m = (logTbound_arr(1)-logTbound_arr(3))/
     $              (logRbound_arr(1)-logRbound_arr(3))
               b = logTbound_arr(1) - m*logRbound_arr(1)
               closestT = m*logR + b
            ! Bottom border is rectangular
            else
               closestT = logTbound_arr(1)
            end if
         end if



      ! If we are above the opacityfile
      else if(logtem.gt.min(logTbound_arr(2),logTbound_arr(4)) ) then

         ! If we are on the right side of the opacityfile
         if ( logR.gt.min(logRbound_arr(3),
     $        logRbound_arr(4)) ) then
            ! If the opacityfile on the right side is not rectangular
            ! and has a negative slope
            if ( logRbound_arr(4).lt.logRbound_arr(3) ) then
               m = (logTbound_arr(3)-logTbound_arr(4))/
     $              (logRbound_arr(3)-logRbound_arr(4))
               b = logTbound_arr(3) - m*logRbound_arr(3)
               closestT = m*logR + b
                     
            ! The opacityfile is either rectangular or has a positive
            ! slope.
            else
               if ( logTbound_arr(2).ne.logTbound_arr(4) )then
                  m = (logTbound_arr(2)-logTbound_arr(4))/
     $                 (logRbound_arr(2)-logRbound_arr(4))
                  b = logTbound_arr(2) - m*logRbound_arr(2)
                  closestT = m*logR + b
               ! Top border is rectangular
               else
                  closestT = logTbound_arr(2)
               end if
            end if

         ! If we are on the left side of the opacityfile
         else if(logR.lt.max(logRbound_arr(1),
     $           logRbound_arr(2)) ) then
            ! If the opacityfile on the left side is not rectangular
            ! and has a positive slope
            if ( logRbound_arr(1).lt.logRbound_arr(2) ) then
               m = (logTbound_arr(1)-logTbound_arr(2))/
     $              (logRbound_arr(1)-logRbound_arr(2))
               b = logTbound_arr(1) - m*logRbound_arr(1)
               closestT = m*logR + b

            ! The opacityfile is either rectangular or has a negative
            ! slope.
            else
               if ( logTbound_arr(2).ne.logTbound_arr(4) )then
                  m = (logTbound_arr(2)-logTbound_arr(4))/
     $                 (logRbound_arr(2)-logRbound_arr(4))
                  b = logTbound_arr(2) - m*logRbound_arr(2)
                  closestT = m*logR + b
               ! Top border is rectangular
               else
                  closestT = logTbound_arr(2)
               end if
            end if
            
         ! If we are simply above the opacityfile
         else
            ! If the top border is not rectangular
            if ( logTbound_arr(2).ne.logTbound_arr(4) )then
               m = (logTbound_arr(2)-logTbound_arr(4))/
     $              (logRbound_arr(2)-logRbound_arr(4))
               b = logTbound_arr(2) - m*logRbound_arr(2)
               closestT = m*logR + b
            ! Top border is rectangular
            else
               closestT = logTbound_arr(2)
            end if
         end if
      end if

      end subroutine



      subroutine find_closestR(logR,logtem,logRbound_arr,
     $     logTbound_arr,closestR)
c     1 = Lower-left
c     2 = Upper-left
c     3 = Lower-right
c     4 = Upper-right
      include 'flux_cal.h'
      real*8 logR,logtem,closestR,b,m
      real*8 logRbound_arr(4)
      real*8 logTbound_arr(4)

      ! If we are to the left of the opacityfile
      if ( logR.lt.max(logRbound_arr(1),logRbound_arr(2)) ) then
         ! If we are above the opacityfile
         if ( logtem.gt.min(logTbound_arr(2),
     $        logTbound_arr(4)) ) then
            ! If the opacityfile above is not rectangular
            ! and has a positive slope
            if ( logTbound_arr(4).gt.logRbound_arr(2) ) then
               m = (logRbound_arr(2)-logRbound_arr(4))/
     $              (logTbound_arr(2)-logTbound_arr(4))
               b = logRbound_arr(2) - m*logTbound_arr(2)
               closestR = m*logtem + b
                     
            ! The opacityfile is either rectangular or has a negative
            ! slope.
            else
               ! If the left border is not rectangular
               if ( logRbound_arr(1).ne.logRbound_arr(2) )then
                  m = (logRbound_arr(1)-logRbound_arr(2))/
     $                 (logTbound_arr(1)-logTbound_arr(2))
                  b = logRbound_arr(1) - m*logTbound_arr(1)
                  closestR = m*logtem + b
               ! Left border is rectangular
               else
                  closestR = logRbound_arr(1)
               end if
            end if

         ! If we are below the opacityfile
         else if(logtem.lt.max(logTbound_arr(1),
     $           logTbound_arr(3)) ) then
            ! If the opacityfile below is not rectangular
            ! and has a negative slope
            if ( logTbound_arr(1).gt.logTbound_arr(3) ) then
               m = (logRbound_arr(1)-logRbound_arr(3))/
     $              (logTbound_arr(1)-logTbound_arr(3))
               b = logRbound_arr(1) - m*logTbound_arr(1)
               closestR = m*logtem + b

            ! The opacityfile is either rectangular or has a positive
            ! slope.
            else
               ! If the left border is not rectangular
               if ( logRbound_arr(1).ne.logRbound_arr(2) )then
                  m = (logRbound_arr(1)-logRbound_arr(2))/
     $                 (logTbound_arr(1)-logTbound_arr(2))
                  b = logRbound_arr(2) - m*logTbound_arr(2)
                  closestR = m*logtem + b
               ! Left border is rectangular
               else
                  closestR = logRbound_arr(2)
               end if
            end if
                  
         ! If we are simply to the left of the opacityfile
         else
            ! If the left border is not rectangular
            if ( logRbound_arr(1).ne.logRbound_arr(2) )then
               m = (logRbound_arr(1)-logRbound_arr(2))/
     $              (logTbound_arr(1)-logTbound_arr(2))
               b = logRbound_arr(2) - m*logTbound_arr(2)
               closestR = m*logtem + b
            ! Left border is rectangular
            else
               closestR = logRbound_arr(2)
            end if
         end if






      ! If we are to the right of the opacityfile
      else if ( logR.gt.min(logRbound_arr(3),logRbound_arr(4)) ) then
         ! If we are above the opacityfile
         if ( logtem.gt.min(logTbound_arr(2),
     $        logTbound_arr(4)) ) then
            ! If the opacityfile above is not rectangular
            ! and has a positive slope
            if ( logTbound_arr(4).gt.logRbound_arr(2) ) then
               m = (logRbound_arr(2)-logRbound_arr(4))/
     $              (logTbound_arr(2)-logTbound_arr(4))
               b = logRbound_arr(2) - m*logTbound_arr(2)
               closestR = m*logtem + b
                     
            ! The opacityfile is either rectangular or has a negative
            ! slope.
            else
               ! If the right border is not rectangular
               if ( logRbound_arr(3).ne.logRbound_arr(4) )then
                  m = (logRbound_arr(3)-logRbound_arr(4))/
     $                 (logTbound_arr(3)-logTbound_arr(4))
                  b = logRbound_arr(4) - m*logTbound_arr(4)
                  closestR = m*logtem + b
               ! Right border is rectangular
               else
                  closestR = logRbound_arr(4)
               end if
            end if

         ! If we are below the opacityfile
         else if(logtem.lt.max(logTbound_arr(1),
     $           logTbound_arr(3)) ) then
            ! If the opacityfile below is not rectangular
            ! and has a negative slope
            if ( logTbound_arr(1).gt.logTbound_arr(3) ) then
               m = (logRbound_arr(1)-logRbound_arr(3))/
     $              (logTbound_arr(1)-logTbound_arr(3))
               b = logRbound_arr(1) - m*logTbound_arr(1)
               closestR = m*logtem + b

            ! The opacityfile is either rectangular or has a positive
            ! slope.
            else
               ! If the right border is not rectangular
               if ( logRbound_arr(3).ne.logRbound_arr(4) )then
                  m = (logRbound_arr(3)-logRbound_arr(4))/
     $                 (logTbound_arr(3)-logTbound_arr(4))
                  b = logRbound_arr(4) - m*logTbound_arr(4)
                  closestR = m*logtem + b
               ! Right border is rectangular
               else
                  closestR = logRbound_arr(4)
               end if
            end if
                  
         ! If we are simply to the right of the opacityfile
         else
            ! If the right border is not rectangular
            if ( logRbound_arr(3).ne.logRbound_arr(4) )then
               m = (logRbound_arr(3)-logRbound_arr(4))/
     $              (logTbound_arr(3)-logTbound_arr(4))
               b = logRbound_arr(3) - m*logTbound_arr(3)
               closestR = m*logtem + b
            ! Right border is rectangular
            else
               closestR = logRbound_arr(3)
            end if
         end if
      end if

      end subroutine

      


      

      



      subroutine extrapolate_opacity(tem,rhocgs,xhy,opacity)
      include 'flux_cal.h'
      real*8 tem,rhocgs,xhy,opacity
      real*8 closestR,closestT,closestopac
      integer closest_opacityfilerho,closest_opacityfileT,
     $     closest_opacityfileR,closest_opacityfile
      real*8 ext_Teq,ext_rhoeq,ext_Req,calculate_table_opacity
      real*8 mind,d00,d01,d10,d11,tempmind,cornerR,cornerT
      logical containsR,containsT
      real*8 logRbound_arr(4,numopacityfiles)
      real*8 logTbound_arr(4,numopacityfiles)
      real*8 logTbound1_arr(numopacityfiles),
     $     logTbound2_arr(numopacityfiles),
     $     logRbound1_arr(numopacityfiles),
     $     logRbound2_arr(numopacityfiles)
      real*8 dT, dR, weight, logR, logtem,cornerR1,cornerR2,cornerT1,
     $     cornerT2
      real*8 dlogrho,mindlogrho,closestrho,dlogT,mindlogT,mindlogR,dlogR
      integer i,j,k
      real*8 mindT,mindR
      real*8 logTcompare, logRcompare, minTbound, maxTbound,
     $     minRbound, maxRbound
      real*8 m,b,d,closestR_temp,closestT_temp
      real*8 mindlogR2,mindlogT2

      
      logtem = log10(tem)
      logR = log10(rhocgs) - 3.d0*logtem + 18.d0
      
      ! We have not found a domain for (rho,T) within which we can
      ! calculate opacity. We need to use analytic approximations in
      ! a last-ditch effort to calculate opacity.
c     write(*,*) "We are not using a table"
c     write(*,*) "rho,logR,logT = ",rhocgs,logR,logtem
      closestR = -1.d30
      closestT = -1.d30
      closestopac = -1.d30
      closest_opacityfilerho = 0
      closest_opacityfileT = 0
      closest_opacityfileR = 0
      closest_opacityfile = 0

      mindlogR = 1.d30
      mindlogrho = 1.d30
      mindlogT = 1.d30
      mind = 1.d30
      cornerR = 1.d30
      cornerT = 1.d30
      mindT = 1.d30
      mindR = 1.d30

c     Do not edit these values. If you do, also edit opacityTables.f
      logTbound_arr(1,1) = logTmins(1)
      logTbound_arr(2,1) = logTmaxs(1)
      logTbound_arr(3,1) = logTmins(1)
      logTbound_arr(4,1) = logTmaxs(1)
      logTbound_arr(1,2) = logTmins(2)
      logTbound_arr(2,2) = logTmaxs(2)
      logTbound_arr(3,2) = logTmins(2)
      logTbound_arr(4,2) = logTmaxs(2)
      logRbound_arr(1,1) = (log10(rhoPlanckRossmin)+1.d-8)
     $     -3.d0*logTmins(1) + 18.d0
      logRbound_arr(2,1) = (log10(rhoPlanckRossmin)+1.d-8)
     $     -3.d0*logTmaxs(1) + 18.d0
      logRbound_arr(3,1) = (log10(rhoPlanckRossmax)-1.d-8)
     $     -3.d0*logTmins(1) + 18.d0
      logRbound_arr(4,1) = (log10(rhoPlanckRossmax)-1.d-8)
     $     -3.d0*logTmaxs(1) + 18.d0
      logRbound_arr(1,2) = (log10(rhoPlanckRossmin)+1.d-8)
     $     -3.d0*logTmins(2) + 18.d0
      logRbound_arr(2,2) = (log10(rhoPlanckRossmin)+1.d-8)
     $     -3.d0*logTmaxs(2) + 18.d0
      logRbound_arr(3,2) = (log10(rhoPlanckRossmax)-1.d-8)
     $     -3.d0*logTmins(2) + 18.d0
      logRbound_arr(4,2) = (log10(rhoPlanckRossmax)-1.d-8)
     $     -3.d0*logTmaxs(2) + 18.d0
      

      do i=3,numopacityfiles
         ! Lower-left
         logRbound_arr(1,i) = logRmins(i)
         logTbound_arr(1,i) = logTmins(i)

         ! Upper-left
         logRbound_arr(2,i) = logRmins(i)
         logTbound_arr(2,i) = logTmaxs(i)

         ! Lower-right
         logRbound_arr(3,i) = logRmaxs(i)
         logTbound_arr(3,i) = logTmins(i)

         ! Upper-right
         logRbound_arr(4,i) = logRmaxs(i)
         logTbound_arr(4,i) = logTmaxs(i)
      end do
      

      do i=1,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).le.0) cycle
         
         logTbound1_arr(i) = logTmins(i)
         logTbound2_arr(i) = logTmaxs(i)
         logRbound1_arr(i) = logRmins(i)
         logRbound2_arr(i) = logRmaxs(i)

         minRbound = minval(logRbound_arr(:,i))
         maxRbound = maxval(logRbound_arr(:,i))
         minTbound = minval(logTbound_arr(:,i))
         maxTbound = maxval(logTbound_arr(:,i))

         containsR = .false.
         containsT = .false.

         ! See if there is an opacityfile that contains logR AND has 
         ! a logT bound closest to our logtem

         ! closest_opacityfileT = 0 no opacityfile contains our logtem
         ! closest_opacityfileT =/= 0 an opacityfile contains our logtem 
         !
         ! closest_opacityfileR = 0 no opacityfile contains our logR
         ! closest_opacityfileR =/= 0 an opacityfile contains our logR
         !
         ! closest_opacityfile = 0 our (R,T) is within an opacityfile's range
         ! closest_opacityfile =/= 0 this opacityfile is the nearest to (R,T)
         if ( minRbound.le.logR .and.
     $        logR.le.maxRbound ) then
            containsR = .true.
            call find_closestT(logR,logtem,logRbound_arr(:,i),
     $           logTbound_arr(:,i),closestT_temp)
            dlogT = abs(closestT_temp-logtem)
            if ( dlogT.lt.mindlogT ) then
               mindlogT = dlogT
               closestT = closestT_temp
               closest_opacityfileT = i
               
               ! Update which opacityfile is the closest to be this one
               closest_opacityfile = i
            end if
         end if

         if ( minTbound.le.logtem .and.
     $        logtem.le.maxTbound ) then
            containsT = .true.
            call find_closestR(logR,logtem,logRbound_arr(:,i),
     $           logTbound_arr(:,i),closestR_temp)
            dlogR = abs(closestR_temp-logR)
            if ( dlogR.lt.mindlogR ) then
               mindlogR = dlogR
               closestR = closestR_temp
               closest_opacityfileR = i
               
               ! Update which opacityfile is the closest to be this one
               closest_opacityfile = i
            end if
         end if



         ! We need to find the opacityfile with the minimum distance 
         ! to our (R,T), but we care mostly about the closest opacityfile
         ! in the T direction, not the R direction. 
         if ((containsR.eqv..false.).and.(containsT.eqv..false.)) then
            mindlogT2 = 1.d30
            mindlogR2 = 1.d30
            do j=1,4
               dlogT = abs(logtem-logTbound_arr(j,i))
               if ( dlogT.lt.mindlogT2 ) then
                  mindlogT2 = dlogT
                  cornerR = logRbound_arr(j,i)
                  cornerT = logTbound_arr(j,i)
                  closest_opacityfile = i

               ! Check to see if cornerR should change
               else if (dlogT.eq.mindlogT2) then
                  dlogR = abs(logR-logRbound_arr(j,i))
                  if ( dlogR.lt.abs(logR-cornerR) ) then
                     cornerR = logRbound_arr(j,i)
                  end if
               end if
            end do
         end if
      end do

      ! If we have selected the Planck opacityfile as the closest,
      ! but in fact the Rosseland opacity is the closest in the T direction,
      ! reset the value of closest_opacityfile to the Rosseland opacityfile
      if ( closest_opacityfile.eq.1 .and.
     $     logtem.gt.min(logTbound_arr(2,2),logTbound_arr(4,2)) ) then
         closest_opacityfile = 2
      end if

      ! If our (R,T) is not between two opacityfiles and closest_opacityfile
      ! is not rectangular
      if ( logRbound_arr(1,closest_opacityfile).ne.
     $     logRbound_arr(2,closest_opacityfile) .or.
     $     logRbound_arr(3,closest_opacityfile).ne.
     $     logRbound_arr(4,closest_opacityfile) .or.
     $     logTbound_arr(1,closest_opacityfile).ne.
     $     logTbound_arr(3,closest_opacityfile) .or.
     $     logTbound_arr(2,closest_opacityfile).ne.
     $     logTbound_arr(4,closest_opacityfile) ) then
         ! Set the corner values
         mindlogT2 = 1.d30
         mindlogR2 = 1.d30
         do i=1,4
            dlogT = abs(logtem-logTbound_arr(i,closest_opacityfile))
            if ( dlogT.lt.mindlogT2 ) then
               mindlogT2 = dlogT
               cornerR = logRbound_arr(i,closest_opacityfile)
               cornerT = logTbound_arr(i,closest_opacityfile)

            else if (dlogT.eq.mindlogT2) then
               dlogR = abs(logR-logRbound_arr(i,closest_opacityfile))
               if ( dlogR.lt.abs(logR-cornerR) ) then
                  cornerR = logRbound_arr(i,closest_opacityfile)
               end if
            end if
         end do

         ! Force this to be a corner case
         if ( logtem.gt.cornerT .and. logR.gt.cornerR .and.
     $        logtem.gt.logTmaxs(2) ) then
            closest_opacityfileR = 0
            closest_opacityfileT = 0
         end if
      end if

      ! In the "between case", we need to do our temperature smoothing
      ! relative to the nearest (in temperature) opacityfile corner. We
      ! send analytic_extrap a cornerT value to accomplish this.
      if ( closest_opacityfileR.ne.0.and.closest_opacityfileT.ne.0 )then
         if ( closest_opacityfileR.le.2 .or.
     $        closest_opacityfileT.le.2 ) then
            cornerT = logTbound_arr(2,2)
         else
            mindlogT = 1.d30
            do i=1,4
               dlogT = abs(logtem-logTbound_arr(i,closest_opacityfileR))
               if ( dlogT.lt.mindlogT ) then
                  mindlogT = dlogT
                  cornerT = logTbound_arr(i,closest_opacityfileR)
               end if
               
               dlogT = abs(logtem-logTbound_arr(i,closest_opacityfileT))
               if ( dlogT.lt.mindlogT ) then
                  mindlogT = dlogT
                  cornerT = logTbound_arr(i,closest_opacityfileT)
               end if
            end do
         end if
      end if

      call analytic_extrap(logR,logtem,closest_opacityfileR,
     $     closest_opacityfileT,closest_opacityfile,closestR,closestT,
     $     cornerR,cornerT,xhy,opacity)

      end subroutine





      

      subroutine analytic_extrap(logR,logT,closest_opacityfileR,
     $     closest_opacityfileT,closest_opacityfile,closestR,closestT,
     $     cornerR,cornerT,xhy,opacity)
c     If this routine fails to compute the analytic extrapolation of
c     the opacity from the opacityfiles, it returns an opacity equal
c     to -1.d30.
      include 'flux_cal.h'
      real*8 logR,logT,closestR,closestT,cornerR,cornerT,xhy
      integer closest_opacityfileR,closest_opacityfileT,
     $     closest_opacityfile
      real*8 opacity,ext_Teq,ext_Req,calculate_table_opacity,
     $     closestopac,rhocgs
      real*8 opac_test
      real*8 smooth_opacity
      real*8 Rsmooth,Tsmooth
      ! Smoothing distance for R and T directions
      parameter(Rsmooth=0.d0,Tsmooth=0.d0)
      real*8 opacityReq,opacityTeq,weightReq,weightTeq,
     $     closestopacReq,closestopacTeq
      real*8 R1,R2,T1,T2,weightR,weightT,weight
      real*8 smooth

c      write(*,*) "closest_opacityfileR, closest_opacityfileT, "//
c     $     "closest_opacityfile = ",closest_opacityfileR,
c     $     closest_opacityfileT,closest_opacityfile
      
      opacity = -1.d30
      
      rhocgs = 10.d0**(logR + 3.d0*logT - 18.d0)


      ! closest_opacityfileT = 0 no opacityfile contains our logtem
      ! closest_opacityfileT =/= 0 an opacityfile contains our logtem 
      !
      ! closest_opacityfileR = 0 no opacityfile contains our logR
      ! closest_opacityfileR =/= 0 an opacityfile contains our logR
      !
      ! closest_opacityfile = 0 our (R,T) is within an opacityfile's range
      ! closest_opacityfile =/= 0 this opacityfile is the nearest to (R,T)



      ! Corner case
      ! In this case, none of our opacityfiles have logtem nor logR
      ! in their range. We can still try to extrapolate to our
      ! (rho,T), however. First we extrapolate with constant T and
      ! then we extrapolate with constant R.
      !
      ! The second condition in this if statement (in parentheses) is for
      ! when there is an opacityfile in the T direction that does not contain
      ! logR but is closer to our (R,T) than an opacityfile in the T direction
      ! that does contain our logR. We want to favor closeness in the T
      ! direction regardless of R.
      if ((closest_opacityfileR.eq.0 .and.
     $     closest_opacityfileT.eq.0 .and.
     $     closest_opacityfile.ne.0) .or.
     $    (closest_opacityfileR.eq.0 .and. ! Favor closeness in T direction
     $     closest_opacityfileT.ne.0 .and. 
     $     closest_opacityfile.ne.closest_opacityfileT .and.
     $     closest_opacityfile.ne.0)) then
c         write(*,*) "Corner case"
         closestopac = calculate_table_opacity(10.d0**cornerR,
     $        10.d0**cornerT)
         R1 = cornerR
         T1 = cornerT
         opacity = ext_Teq(10.d0**cornerR,10.d0**logR,10.d0**cornerT,
     $        closestopac,xhy,metallicity,.false.)

c         write(*,*) "cornerR = ",cornerR
c         write(*,*) "cornerT = ",cornerT
c         write(*,*) "closestopac = ",closestopac
c         write(*,*) "opacity = ",opacity

         if ( opacity.ne.-1.d30 ) then
            ! Test the T direction (constant R)
            opac_test = ext_Req(10.d0**logR,10.d0**cornerT,10.d0**logT,
     $           opacity,xhy,metallicity,.false.)
            if ( opac_test .ne. -1.d30 ) then
               opacity = opac_test
            end if
         else
            ! Try the T direction first (constant R)
            opacity = ext_Req(10.d0**logR,10.d0**cornerT,10.d0**logT,
     $           closestopac,xhy,metallicity,.false.)
            if ( opacity.ne.-1.d30 ) then
               ! Test the R direction (constant T)
               opac_test = ext_Teq(10.d0**cornerR,10.d0**logR,
     $              10.d0**cornerT,opacity,xhy,metallicity,.false.)
c               write(*,*) "opac_test = ",opac_test
c               write(*,*) "opacity, opac_test = ",opacity,opac_test
               if ( opac_test .ne. -1.d30 ) then
                  opacity = opac_test
               end if
            end if
         end if

      ! Teq case
      ! We have found an opacityfile that has logtem in its range but no
      ! opacityfile containing logR
      else if((closest_opacityfileR.ne.0 .and.
     $        closest_opacityfileT.eq.0) .or.
     $        (closest_opacityfileR.le.2 .and.
     $        closest_opacityfileR.gt.0)) then

c         write(*,*) "Teq case"
         closestopac = calculate_table_opacity(10.d0**closestR,
     $        10.d0**logT)
         if ( cornerR.ne.1.d30 .and. cornerT.ne.1.d30 ) then
            R1 = cornerR
            T1 = cornerT
         else
            R1 = closestR
            T1 = logT
         end if
c         write(*,*) "cornerR,cornerT = ",cornerR,cornerT
         opacity = ext_Teq(10.d0**closestR,10.d0**logR,10.d0**logT,
     $        closestopac,xhy,metallicity,.false.) ! If error, opacity=-1.d30
         if ( opacity.eq.-1.d30 ) then
            ! This is the last thing we want. It's non-physical, but
            ! so is returning opacity = -1.d30. This is the lesser of two evils.
            opacity = closestopac
         end if
         

      ! Req case
      ! We have found an opacityfile that has logR in its logR range
      else if(closest_opacityfileR.eq.0 .and.
     $        closest_opacityfileT.ne.0)then
c         write(*,*) "Req case"

         closestopac = calculate_table_opacity(10.d0**logR,
     $        10.d0**closestT)

         ! For the sake of smoothing, we need to smooth from the corners
         ! of opacityfiles when the opacityfile is not rectangular
         if ( cornerR.ne.1.d30 .and. cornerT.ne.1.d30 ) then
            R1 = cornerR
            T1 = cornerT
         else
            R1 = logR
            T1 = closestT
         end if

c         write(*,*) "cornerR = ",cornerR
c         write(*,*) "cornerT = ",cornerT

c         write(*,*) "closestopac = ",closestopac
         opacity = ext_Req(10.d0**logR,10.d0**closestT,10.d0**logT,
     $        closestopac,xhy,metallicity,.false.) ! If error, opacity=-1.d30
         if ( opacity.eq.-1.d30 ) then
            ! This is the last thing we want. It's non-physical, but
            ! so is returning opacity = -1.d30. This is the lesser of two evils.
            opacity = closestopac
         end if

      ! Between case
      ! We have found 2 opacityfiles: one with logtem, and one with logR
      else if((closest_opacityfileR.ne.0 .and.
     $         closest_opacityfileT.ne.0) .and.
     $        (closest_opacityfileR.gt.2)) then
c         write(*,*) "Between case"
         ! Try extrapolating with constant T
         closestopacTeq = calculate_table_opacity(10.d0**closestR,
     $        10.d0**logT)
         opacityTeq = ext_Teq(10.d0**closestR,10.d0**logR,10.d0**logT,
     $        closestopacTeq,xhy,metallicity,.false.) ! If error, opacity=-1.d30

         ! Try extrapolating with constant R
         closestopacReq = calculate_table_opacity(10.d0**logR,
     $        10.d0**closestT)
         opacityReq = ext_Req(10.d0**logR,10.d0**closestT,10.d0**logT,
     $        closestopacReq,xhy,metallicity,.false.) ! If error, opacity=-1.d30
         
         ! By default, use opacityTeq because this is the behavior defined
         ! in the extrapolation routines.
         if ( opacityTeq.ne.-1.d30 ) then
            opacity = opacityTeq
            closestopac = closestopacTeq
            if ( cornerT.ne.1.d30 ) then
               T1 = cornerT
            else
               T1 = closestT
            end if
            R1 = closestR
c            write(*,*) "closestR = ",closestR
c            write(*,*) "cornerT = ",cornerT
         else if( opacityReq.ne.-1.d30 ) then ! If Teq extrapolation failed
            opacity = opacityReq ! Take the Req extrapolation result
            closestopac = closestopacReq
            R1 = logR
            if ( cornerT.ne.1.d30 ) then
               T1 = cornerT
            else
               T1 = closestT
            end if
         else                   ! Rely on closest opacityfile in R direction
            opacity = closestopacTeq
            closestopac = closestopacTeq
            R1 = closestR
            if ( cornerT.ne.1.d30 ) then
               T1 = cornerT
            else
               T1 = closestT
            end if
         end if
         
      else if(closest_opacityfileR.eq.0 .and.
     $        closest_opacityfileT.eq.0 .and.
     $        closest_opacityfile.eq.0)then
         write(*,*) "(rho,T) = (",rhocgs,10.d0**logT,")"
         write(*,*) "Could not find the closest opacityfile to"//
     $        " (rho,T). Check your input file for missing "//
     $        "opacityfiles or incorrect logTmins, logTmaxs, "//
     $        "logRmins, and logRmaxs values."
         error stop "getOpacity.f"
      end if


      if ( opacity.eq.-1.d30 ) then ! Error
         if ( opacity_oob_warning .or. opacity_oob_error ) then
            write(*,*) "(rho,T) = (",rhocgs,10.d0**logT,")"
            write(*,*) "WARNING (getOpacity.f): (rho,T) is out "//
     $           "of bounds of your opacityfiles and FluxCal is"//
     $           " unable to extrapolate an opacity value from "//
     $           "your provided opacityfiles."
            if ( opacity_oob_error ) error stop "getOpacity.f"
         end if

      ! No error. Smooth between closestopac and opacity
      else
         weightR = 1.d0
         weightT = 1.d0
         if ( logR .ge. R1 ) then
            R2 = R1+Rsmooth
            if ( logR.ge.R1 .and. logR.le.R2 ) then
               weightR = (R2-logR)/(R2-R1)
            end if
         else
            R2 = R1-Rsmooth
            if ( logR.ge.R2 .and. logR.le.R1 ) then
               weightR = (logR-R2)/(R1-R2)
            end if
         end if

         if ( logT .ge. T1 ) then
            T2 = T1+Tsmooth
            if ( logT.ge.T1 .and. logT.le.T2 ) then
               weightT = (T2-logT)/(T2-T1)
            end if
         else
            
            if ( closest_opacityfileR.ge.2 ) then
               T2 = T1-Tsmooth
               if ( logT.ge.T2 .and. logT.le.T1 ) then
                  weightT = (logT-T2)/(T1-T2)
               end if
            end if
         end if

         if ( logT.lt.max(T1,T2) .and.
     $        logT.gt.min(T1,T2) ) then
            weight = weightT
         else
            weight = 0.d0
         end if
         
c         write(*,*) "T1,logT,T2 = ",T1,logT,T2
c         write(*,*) "weight = ",weight
c         write(*,*) "weightR, weightT = ",weightR,weightT
         
         ! weight = 1 when logR = R1 or logT = T1
         ! weight = 0 when logR = R2 or logT = T2
c         opacity = weight*closestopac + (1.d0-weight)*opacity

         opacity = smooth(10.d0**logT,10.d0**T1,10.d0**T2,
     $        closestopac,opacity)
         
      end if

      return
      end subroutine





      real*8 function calculate_table_opacity(myR,myT)
c     We need to be able to calculate the opacity from the opacityfiles
c     only, so that we can extrapolate outside the files later if need be.
c     Hence, we must include all blending considerations, stepping through
c     all opacityfiles to find the opacity. Also, we should not blend
c     opacities unless there is actual overlap between two opacityfiles.
c
c     Returns the opacity, not log opacity.
      use bilinear_interpolation
      include "flux_cal.h"
      integer i
      real*8 myR,myT
      real*8 opacity_P,opacity_R,myrho,logmyT,logmyR,opacity
      real*8 opacities(numopacityfiles-1),weights(numopacityfiles-1)
      real*8 weightT, weightR, weight,weight1,weight2
      external cop
      real*8 eD_R, eG_R, eD_P, eG_P
      DIMENSION eD_R(5,6), eD_P(5,6)
      DIMENSION eG_R(71,71), eG_P(71,71)
      common/opac/eD_R,eG_R,eD_P,eG_P
      integer numtabused
      real*8 smooth
      
      ! If the function returns this value, it indicates an error (move
      ! to attempts at extrapolation)
      calculate_table_opacity = -1.d30
      
      logmyT = log10(myT)
      logmyR = log10(myR)
      myrho = 10.d0**(logmyR + 3.d0*logmyT - 18.d0)

      do i=1,numopacityfiles-1
         opacities(i) = 0.d0
         weights(i) = 0.d0
      end do


      numtabused = 0
      
      ! No matter what, start with the cop subroutine first. Check if (rho,T) is
      ! within the bounds of the cop subroutine and, if so, calculate those 
      ! opacities. We treat these first two opacityfiles as one.
      if ((myrho.gt.rhoPlanckRossmin.and.myrho.lt.rhoPlanckRossmax).and.
     $    (myT.gt.10.d0**logTmins(1).and.myT.lt.10.d0**logTmaxs(2)))then
         opacity_P = 0.d0
         opacity_R = 0.d0
         
         if(len(trim(adjustl(opacityfiles(1)))).gt.0) then
            call cop(eD_P,eG_P,myrho,myT,opacity_P,0) ! Planck
            numtabused = numtabused + 1
         end if
         if(len(trim(adjustl(opacityfiles(2)))).gt.0) then
            call cop(eD_R,eG_R,myrho,myT,opacity_R,0) ! Rosseland
            numtabused = numtabused + 1
         end if

         
         
         opacities(1) = smooth(myT,10.d0**logT_blend2(1),
     $        10.d0**logT_blend1(2),opacity_P,opacity_R)

c         calculate_table_opacity = opacities(1)
c         return

c         if(myT.lt.10.d0**(2)) then
c            write(*,*) "This is happening", log10(myrho),log10(myT)
c            write(*,*) opacities(1),opacity_P,opacity_R
c         end if

c         ! We need to stitch the Planck and Rosseland files together
c         ! prior to considering any blending
c         
c         ! Within only the Planck opacities
c         if ( logTmins(1).le.logmyT .and.
c     $        logmyT.le.logT_blend2(1) ) then
c            opacities(1) = opacity_P
c            
c         ! Within only the Rosseland opacities
c         else if(logT_blend1(2).le.logmyT .and.
c     $           logmyT.le.logTmaxs(2)) then
c            opacities(1) = opacity_R
c            
c         ! Between Planck and Rosseland
c         else if(logmyT.le.logT_blend1(2) .and.
c     $           logmyT.ge.logT_blend2(1)) then
c            weight1 = (logmyT-logT_blend2(1))/ ! Planck weight
c     $           (logTmaxs(1)-logT_blend2(1))
c            weight2 = (logT_blend1(2)-logmyT)/ ! Rosseland weight
c     $           (logT_blend1(2)-logTmins(2))
c            opacities(1) = ((1.d0-weight1)*opacity_P +
c     $           (1.d0-weight2)*opacity_R)/(2.d0-weight1-weight2)
c         end if

         ! Record any blending that could potentially occur. If
         ! opacities(1) = 0.d0, keep the weight associated with
         ! it equal to zero so it doesn't influence the weighted
         ! sum.
         if ( opacities(1).ne.0.d0 ) then
            weightT = 1.d0
            if ( logT_blend1(1).ne.logTmins(1) .and.
     $           logmyT.le.logT_blend1(1) .and.
     $           logmyT.ge.logTmins(1) )then
               weightT = weightT * (1.d0 - ((logT_blend1(1)-logmyT)/
     $              (logT_blend1(1)-logTmins(1))))
            end if
            if ( logT_blend2(2).ne.logTmaxs(2) .and.
     $           logmyT.ge.logT_blend2(2) .and.
     $           logmyT.le.logTmaxs(2) )then
               weightT = weightT * (1.d0 - ((logmyT-logT_blend2(2))/
     $              (logTmaxs(2)-logT_blend2(2))))
            end if
            
            ! Blending doesn't occur in R for the first 2 files
            
            weights(1) = weightT
            
         end if
      end if

c      write(*,*) "logmyT, logmyR = ",logmyT,logmyR

      
      ! Search through the rest of the opacityfiles for
      ! anything to blend with.
      do i=3,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).eq.0) cycle
         
         ! Check to see if we are in bounds
         if ( logmyT.le.logTmaxs(i) .and. logmyT.ge.logTmins(i) .and.
     $        logmyR.le.logRmaxs(i) .and. logmyR.ge.logRmins(i) ) then

            ! Calculate the opacity
            opacities(i-1) = bilinear_interpolate(
     $           nrows_tables(i),
     $           logTs(i,1:nrows_tables(i)),
     $           ncols_tables(i),
     $           logRs(i,1:ncols_tables(i)),
     $           logopacitytables(i,
     $           1:nrows_tables(i),1:ncols_tables(i)),
     $           logmyT,logmyR)

            weightT = 1.d0
            ! Lower T blending region
            if ( logT_blend1(i).ne.logTmins(i) .and.
     $           logmyT.le.logT_blend1(i) ) then
               weightT = weightT * (1.d0 - ((logT_blend1(i)-logmyT)/
     $              (logT_blend1(i)-logTmins(i))))
            end if

            ! Upper T blending region
            if ( logT_blend2(i).ne.logTmaxs(i) .and.
     $           logmyT.ge.logT_blend2(i) ) then
               weightT = weightT * (1.d0 - ((logmyT-logT_blend2(i))/
     $              (logTmaxs(i)-logT_blend2(i))))
            end if

            weightR = 1.d0
            ! Lower R blending region
            if ( logR_blend1(i).ne.logRmins(i) .and.
     $           logmyR.le.logR_blend1(i) ) then
               weightR = weightR * (1.d0 - ((logR_blend1(i)-logmyR)/
     $              (logR_blend1(i)-logRmins(i))))
            end if

            ! Upper R blending region
            if ( logR_blend2(i).ne.logRmaxs(i) .and.
     $           logmyR.ge.logR_blend2(i) ) then
               weightR = weightR * (1.d0 - ((logmyR-logR_blend2(i))/
     $              (logRmaxs(i)-logR_blend2(i))))
            end if

            weights(i-1) = weightT*weightR
            
         end if
         numtabused = numtabused + 1
      end do

      ! Just in case we only use the default low T opacities,
      ! make sure we return the correct opacity. The correct
      ! opacity will include no blending at all in this case
      ! with any other opacity files.
      if ( numtabused .eq. 2 .and. opacities(1).ne.0 ) then
         weights(1) = 1.d0
      end if
      
      ! cop routine returns non-logarithmic opacities
      calculate_table_opacity = weights(1)*opacities(1)
      do i=2,numopacityfiles-1
         calculate_table_opacity =
     $        calculate_table_opacity + weights(i)*10.d0**opacities(i)
      end do

      if ( sum(weights).gt.0.d0 ) then
         calculate_table_opacity = calculate_table_opacity/
     $        sum(weights)
      else if ( calculate_table_opacity.eq.0.d0 ) then
         ! This implies we are outside all opacityfiles. Return -1.d30 so
         ! that the main opacity routine can try extrapolating
         calculate_table_opacity = -1.d30
      else
         write(*,*) "sum(opacities), sum(weights) = ",sum(opacities),
     $        sum(weights)
         write(*,*) "When trying to calculate the opacity from the "//
     $        "opacityfiles, the sum of the weights was <= 0 and"//
     $        " the sum of the opacities was =/= 0."
         error stop "getOpacity.f"
      end if
      
      return
      end function








c     All of the following functions are to be used ONLY if the given
c     (rho,T) is outside all provided opacityfiles.
c      
c     We need to extrapolate from some point on an opacityfile,
c     (rho_o,T_o) to a new point (rho,T), where (rho_T) lies outside
c     the boundaries of the opacityfile. If rho_o = rho, the only
c     difference is in temperature, and if T_o=T, the only difference
c     is in density.
c
c        *-----------------*
c        |                 |
c        |                 |
c        |                 |
c        |                 | 
c        |     (rho_o,T_o) X--X (rho,T)
c     T  |                 |
c     ^  *-----------------*
c     |
c     o -- > rho
c     
c     When T_o=T, we can easily extrapolate outside the opacityfile
c     depending on the value of T. We calculate the opacity under
c     three main effects, (1) Kramer opacity, (2) e- scattering, and
c     (3) H- ion absorption. Most of the H- ion absorption region is
c     covered by the cop subroutine, but 10^-7 < rho < 10^-5 is not.

      ! A version that scales by R instead
      real*8 function ext_Teq(R_o,R,T,opac_o,X,Z,error)
      implicit none
      real*8 R_o,R,T,opac_o,X,Z
      real*8 Rbound1,Rbound2
      parameter(Rbound1=1.d0/270.d0,Rbound2=1250.d0/27.d0)
      logical error
      real*8 weighte, weightK
      logical isHion, isescattering, isKramer
      real*8 opacities(5)
      integer i
      real*8 Rsmooth
      parameter(Rsmooth=1.d2)
      real*8 smooth
      real*8 smoothing_window_T,smoothing_window_rho
      common/smoothingwindow/ smoothing_window_T,smoothing_window_rho
      real*8 rho, rho_o

      ext_Teq = -1.d30
      
c      write(*,*) "Extrapolating Teq"

      ! R_o H- ion, R H- ion
      opacities(1) = (R/R_o)**(0.5d0) * opac_o
      ! T Kramer
      opacities(2) = R/R_o * opac_o
      ! R_o H- ion, R < Rbound1
      opacities(3) = (Rbound1/R_o)**(0.5d0) * opac_o
      ! R_o H- ion, R > Rbound2
      opacities(4) = (Rbound2/R_o)**(0.5d0) * opac_o
      ! T e- scattering
      opacities(5) = opac_o
      
      
      ! For H- ion absorption, we need to make sure we are
      ! extrapolating from H- region to H- region (constrains rho
      ! AND rho_o)
      if ( isHion(R_o,T,Z) .and. isHion(R,T,Z) ) then
         ext_Teq = opacities(1)
      else if(isKramer(T))then
         ext_Teq = opacities(2) + opacities(5)
      else if(isescattering(T))then
         ext_Teq = opacities(5)
         
c      else if (isescattering(T) .or. isKramer(T)) then
c         ext_Teq = smooth(T,1.d4,1.d6,opacities(2),
c     $        opacities(5))

      ! If analytic extrapolation failed, we take the boundary case
      ! as the result instead (constant opacity)
      else if(isHion(R_o,T,Z) .and.
     $        (Rbound1.gt.R .and. R.le.Rbound2)) then
         ! Here we are below Rbound1 only
         ext_Teq = opacities(3)
      else if(isHion(R_o,T,Z) .and.
     $        (Rbound1.le.R .and. R.gt.Rbound2)) then
         ! Here we are above Rbound2 only
         ext_Teq = opacities(4)
      else if(T.lt.3.d3 .or.
     $        T.gt.1.d8 ) then  ! We have no information on this region
         ext_Teq = opacities(5) ! Trust the opacityfiles
      end if


      ! We can only get here if there has been an error
      if (ext_Teq.eq.-1.d30 .and. error) then
         write(*,*) "R_o,R,T,opac_o,Z = ",R_o,R,T,opac_o,Z
         write(*,*) "ext_Teq can only be called for 3.d3<T<6.d3"//
     $        " with ",Rbound1,"<R<",Rbound2," and 0.001d0<Z<"//
     $        "0.02d0 or for 1.d4<T<1.d8."
         write(*,*) "ext_Teq has been called for a (rho,T) that is"//
     $        " out of its valid bounds. Check to make sure your "//
     $        "opacityfiles and their ranges are correct in your "//
     $        "input file."
         error stop "getOpacity.f"
      end if

      rho = 10.d0**(log10(R)+3.d0*log10(T)-18.d0)
      rho_o = 10.d0**(log10(R_o)+3.d0*log10(T)-18.d0)
      if ( rho .gt. rho_o ) then
         ext_Teq = smooth(rho,rho_o,rho_o-smoothing_window_rho,opac_o,
     $        ext_Teq)
      else if ( rho .lt. rho_o ) then
         ext_Teq = smooth(rho,rho_o,rho_o+smoothing_window_rho,opac_o,
     $        ext_Teq)
      else
         ext_Teq = opac_o
      end if
     
      return
      end function
      




      logical function isKramer(T)
      implicit none
      real*8 T
      if ( 1.d4.le.T   .and.   T.le.1.d6 ) then
         isKramer = .true.
      else
         isKramer = .false.
      end if
      return
      end function
      logical function isHion(R,T,Z)
      implicit none
      real*8 R,T,Z
      real*8 Rbound1,Rbound2
      parameter(Rbound1=1.d0/270.d0,Rbound2=1250.d0/27.d0)
      if ( 3.d3.le.T    .and. T.le.6.d3 .and.
     $     Rbound1.le.R .and. R.le.Rbound2 .and.
     $     0.001d0.le.Z .and. Z.le.0.02d0 ) then
         isHion = .true.
      else
         isHion = .false.
      end if
      return
      end function
      logical function isescattering(T)
      implicit none
      real*8 T
      if ( 1.d4.le.T .and. T.le.1.d8 ) then
         isescattering = .true.
      else
         isescattering = .false.
      end if
      return
      end function



      real*8 function smooth2(T,N,Ts,ks)
      implicit none
      integer N, mink,maxk,i
      real*8 T, minT, maxT
      real*8 Ts(N), ks(N), B(N), C(N), D(N)
      real*8 SEVAL

      minT = 1.d30
      maxT = -1.d30
      do i=1,N
         if ( Ts(i).lt.minT ) then
            minT = Ts(i)
            mink = ks(i)
         end if
         if ( Ts(i).gt.maxT ) then
            maxT = Ts(i)
            maxk = ks(i)
         end if
      end do

      if ( T.le.minT ) then
         smooth2 = mink
      else if (T.ge.maxT)then
         smooth2 = maxk
      else
         do i=1,N
            Ts(i) = log10(Ts(i))
            ks(i) = log10(ks(i))
         end do
         call SPLINE(N,Ts,ks,B,C,D)
         smooth2 = 10.d0**(SEVAL(N,log10(T),Ts,ks,B,C,D))
      end if
      
      return
      end function

      

      real*8 function smooth(x,x1,x2,y1,y2)
c     Do smoothing in log space for greater effect. Otherwise,
c     opacities change too sharply!
      implicit none
      real*8 x,y1,y2,x1,x2
      real*8 weight
      real*8 SEVAL
      real*8 xarr(2), yarr(2), B(2), C(2), D(2)

      if ( x.le.x1 ) then
         smooth = y1
      else if ( x.ge.x2 ) then
         smooth = y2
      else
         ! Do a cubic spline interpolation
         xarr(1) = log10(x1)
         xarr(2) = log10(x2)
         yarr(1) = log10(y1)
         yarr(2) = log10(y2)
         call SPLINE(2,xarr,yarr,B,C,D)
         smooth = 10.d0**(SEVAL(2,log10(x),xarr,yarr,B,C,D))
c         weight = (T-T1)/(T-T2)
c         smooth = (k1*exp(weight) + k2*exp(1.d0/weight))/
c     $        (exp(weight)+exp(1.d0/weight))
      end if
      
      return
      end function

      


      
      ! A version that uses R instead
      real*8 function ext_Req(R,T_o,T,opac_o,X,Z,error)
      implicit none
      real*8 R,T_o,T,opac_o,X,Z
      real*8 Rbound1,Rbound2
      parameter(Rbound1=1.d0/270.d0,Rbound2=1250.d0/27.d0)
      logical error
      real*8 opacities(7)
      real*8 weights(7)
      integer i,nused
      logical isKramer, isHion, isescattering
      real*8 sumweights
      real*8 smooth

      ext_Req = -1.d30

c      write(*,*) "Extrapolating Req"

      ! T_o Kramer, T Kramer
      opacities(1) = (T/T_o)**(-0.5d0) * opac_o
      ! T_o H-, T H-
      opacities(2) = (T/T_o)**(21.d0/2.d0) * opac_o
      ! T e-, T_o anything
      opacities(3) = 0.2d0*(1.d0+X)
      ! H- ion absorption
      opacities(4) = 2.5d-31*(Z/0.02d0)*R**(0.5d0)*T**(21.d0/2.d0) *
     $     10.d0**(-18.d0/2.d0)
      ! Bound-free
      opacities(5) = 4.3d25*(1.d0+X)*Z*R*T**(-0.5d0) * 10.d0**(-18.d0)
      ! Free-free
      opacities(6) = 3.8d22*(1.d0+X)*R*T**(-0.5d0) * 10.d0**(-18.d0)
      ! Original opacity
      opacities(7) = opac_o


c     For us, there are 5 empty regions, or "gaps" in our analytic
c     function range,
c     
c     Gap 1:
c          3.d3 <= T <= 6.d3
c        Rbound1 > R > Rbound2
c       0.001d0 <= Z <= 0.02d0
c     Gap 2:
c          3.d3 <= T <= 6.d3
c        Rbound1 > R > Rbound2
c        0.001d0 > Z > 0.02d0
c     Gap 3:
c           6.d3 < T < 1.d4
c     Gap 4:
c                  T < 3.d3
c     Gap 5:
c           1.d8 < T
c     
c     Across these gaps, we need to apply some kind of blending solution
c     in the case where the user doesn't provide a large enough opacityfile
c     range (which is almost certain).

      
      if ( isKramer(T) ) then
         ! If it's Kramer, then it is also e- scattering. Need to weight.
         if ( isKramer(T_o) ) then
            ext_Req = smooth(T,1.d4,1.d6,opacities(1),
     $           opacities(3))
         else if ( Z.ge.1.d-3) then
            ext_Req = smooth(T,1.d4,1.d6,opacities(5),
     $           opacities(3))
         else
            ext_Req = smooth(T,1.d4,1.d6,opacities(6),
     $           opacities(3))
         end if

      else if(isHion(R,T,Z))then
         if ( isHion(R,T_o,Z) ) then
            ext_Req = opacities(2)
         else
            ext_Req = opacities(4)
         end if
      else if(isescattering(T)) then
         ext_Req = opacities(3)
         
         
      else                      ! We are in one of the gaps
         ! Gap 1
         if ( 3.d3.le.T .and. T.le.6.d3 .and.
     $        (Rbound1.gt.R .or. R.gt.Rbound2) .and.
     $        0.001d0.le.Z .and. Z.le.0.02d0 ) then
            if ( 3.d3.le.T_o .and. T_o.le. 6.d3 ) then ! Trust the opacityfiles
c               write(*,*) "Gap 1"
               ext_Req = opacities(7)
            end if

            
         ! Gap 2
         else if(3.d3.le.T .and. T.le.6.d3 .and.
     $           (Rbound1.gt.R .or. R.gt.Rbound2) .and.
     $           (0.001d0.gt.Z .or. Z.gt.0.02d0) ) then
            if ( 3.d3.le.T_o .and. T_o.le.6.d3 ) then ! Trust the opacityfiles
c               write(*,*) "Gap 2"
               ext_Req = opacities(7)
            end if

            
         ! Gap 3
         else if(6.d3.lt.T .and. T.lt.1.d4)then
            if ( 6.d3.lt.T_o .and. T_o.lt.1.d4 ) then ! Trust the opacityfiles
c               write(*,*) "Gap 3"
               ext_Req = opacities(7)
            end if
            
            
         ! Gap 4
         else if(T.lt.3.d3) then
            if ( T_o .lt. 3.d3 ) then ! Trust the opacityfiles
c               write(*,*) "Gap 4"
               ext_Req = opacities(7)
            end if

            
         ! Gap 5
         else if(1.d8.lt.T) then
c            write(*,*) "Gap 5"
            if ( 1.d8.lt.T_o ) then ! Trust the opacityfiles
               ext_Req = opacities(7)
            else                ! Just use analytic e- scattering
               ext_Req = opacities(3)
            end if
            
         end if
         
      end if

      
      ! We can only get here if there has been an error
      if (ext_Req.eq.-1.d30 .and. error) then
         write(*,*) "R,T_o,T,opac_o,X,Z = ",R,T_o,T,opac_o,X,Z
         write(*,*) "ext_Req has been called for a (R,T) that is "//
     $        "out of its valid bounds. Check to make sure your "//
     $        "opacityfiles and their ranges are correct in your "//
     $        "input file."
         error stop "getOpacity.f"
      end if
      
      
      return
      end function















      


      



         
