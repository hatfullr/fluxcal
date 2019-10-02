      function getOpacity(tem,rhocgs,xhy) ! real*8
      use bilinear_interpolation
      include 'flux_cal.h'
c     INPUT:
c     tem    - temperature
c     rhocgs - density in cgs
c     Xh     - hydrogen fraction
      
      real*8 tem,rhocgs,logR,logtem,logrho
      real*8 logT0,logT1,logR0,logR1,
     $     logopac00,logopac01,logopac10,logopac11,
     $     logrho00,logrho01,logrho10,logrho11
      real*8 opacity,tempopacity
      integer i,j,k,numtabused
      real*8 weight, weightT, weightR, sumweights
      real*8 xhy

      external cop
      real*8 eD_R, eG_R, eD_P, eG_P
      DIMENSION eD_R(5,6), eD_P(5,6)
      DIMENSION eG_R(71,71), eG_P(71,71)
      common/opac/eD_R,eG_R,eD_P,eG_P

      

      logtem = log10(tem)
      logrho = log10(rhocgs)
      logR = logrho - 3.d0*logtem + 18.d0

      numtabused = 0
      opacity = 0.d0
      sumweights = 0.d0

      ! No matter what, start with the cop subroutine first. Check if (rho,T) is
      ! within the bounds of the cop subroutine and, if so, calculate those 
      ! opacities.
      if ( (rhocgs.gt.1.d-19 .and. rhocgs.lt.1.d-7) .and.
     $     (tem.gt.500.d0 .and. tem.lt.10000.d0) ) then
         
         ! Planck opacity
         if ( len(trim(adjustl(opacityfiles(1)))).gt.0 ) then
            if ( logtem.le.logTmaxs(1) ) then
               call cop(eD_P,eG_P,rhocgs,tem,tempopacity,1)
               weight = 1.d0
               if ( logtem.ge.logT_blend2(1) ) then
                  weight = (logtem-logT_blend2(1))/
     $                 (logTmaxs(1)-logT_blend2(1))
               end if
               opacity = opacity + weight*tempopacity
               sumweights = sumweights + weight
               
               numtabused = numtabused + 1
            end if
         end if

         ! Rosseland opacity
         if ( len(trim(adjustl(opacityfiles(2)))).gt.0 ) then
            if ( logtem .ge. logTmins(2) ) then
               call cop(eD_R,eG_R,rhocgs,tem,tempopacity,1)
               weight = 1.d0
               if ( logtem.le.logT_blend1(2) ) then
                  weight = (logT_blend1(2)-logtem)/
     $                 (logT_blend1(2)-logTmins(2))
               end if
               opacity = opacity + weight*tempopacity
               sumweights = sumweights + weight
               
               numtabused = numtabused + 1
            end if
         end if

         if ( opacity .le. 0.d0 ) then
            write(*,*) "opacity = ",opacity
            write(*,*) "Fatal error in getOpacity.f. The cop"//
     $           " subroutine failed to calculate a sensible opacity."
            error stop "getOpacity.f"
         end if

      end if





      
      ! Now check all the other opacityfiles for any opacities we want to 
      ! blend in to the cop subroutine, or to calculate opacities outside the
      ! cop subroutine.
      do i=3,numopacityfiles
         if(len(trim(adjustl(opacityfiles(i)))).le.0) cycle
         
         ! The following if statement only handles the boundary conditions
         ! defined in the user's input file
         if ( logtem.le.logTmins(i) .or.
     $        logtem.ge.logTmaxs(i) .or.
     $        logR.le.logRmins(i) .or.
     $        logR.ge.logRmaxs(i) ) cycle

         weightT = 1.d0
         weightR = 1.d0

         if ( logtem.le.logT_blend1(i) ) then
            weightT = (logtem-logTmins(i))/
     $           (logT_blend1(i)-logTmins(i))
         else if ( logtem.ge.logT_blend2(i) ) then
            weightT = (logTmaxs(i)-logtem)/
     $           (logTmaxs(i)-logT_blend2(i))
         end if
         
         if ( logR.le.logR_blend1(i) ) then
            weightR = (logR-logRmins(i))/
     $           (logR_blend1(i)-logRmins(i))
         else if ( logR.ge.logR_blend2(i) ) then
            weightR = (logRmaxs(i)-logR)/
     $           (logRmaxs(i)-logR_blend2(i))
         end if

         weight = weightT*weightR
         sumweights = sumweights + weight

         opacity = opacity+ weight*10.d0**bilinear_interpolate(
     $        nrows_tables(i),
     $        logTs(i,1:nrows_tables(i)),
     $        ncols_tables(i),
     $        logRs(i,1:ncols_tables(i)),
     $        logopacity_highT_tables(i,
     $        1:nrows_tables(i),1:ncols_tables(i)),
     $        logtem,logR)

c         write(*,*) "opacity = ",opacity,bilinear_interpolate(
c     $        nrows_tables(i),
c     $        logTs(i,1:nrows_tables(i)),
c     $        ncols_tables(i),
c     $        logRs(i,1:ncols_tables(i)),
c     $        logopacity_highT_tables(i,
c     $        1:nrows_tables(i),1:ncols_tables(i)),
c     $        logtem,logR)
         
         numtabused = numtabused + 1
      end do



      ! Now that we have checked all the opacityfiles as well as the cop
      ! subroutine, if we still have not found a domain for (rho,T) within which
      ! we can calculate opacity, we need to use analytic approximations in a
      ! last-ditch effort to calculate opacity.
      if ( numtabused.eq.0 ) then
         ! Bound-free and electron scattering opacity
         if ( 1.d4.le.tem .and. tem.le.1.d6 .and.
     $        metallicity .ge. 1.d-3 ) then
            weight = (1.d6-tem)/(1.d6-1.d4)
            opacity = weight*boundfree(rhocgs,tem,xhy,metallicity,
     $           opacity_analytic_warning)
     $           + (1.d0-weight)*escattering(rhocgs,tem,xhy,
     $           opacity_analytic_warning)
         ! Free-free and electron scattering opacity
         else if ( 1.d4.le.tem .and. tem.le.1.d6 .and.
     $           metallicity .lt. 1.d-3 ) then
            weight = (1.d6-tem)/(1.d6-1.d4)
            opacity = weight*freefree(rhocgs,tem,xhy,metallicity,
     $           opacity_analytic_warning)
     $           + (1.d0-weight)*escattering(rhocgs,tem,xhy,
     $           opacity_analytic_warning)
         ! Purely electron scattering opacity
         else if ( tem.gt.1.d6 .and. tem.le.1.d8 ) then
            opacity = escattering(rhocgs,tem,xhy,
     $           opacity_analytic_warning)
         ! Negative H ion opacity
         else if ( 3.d3.le.tem .and. tem.le.1.d4 .and.
     $           1.d-10.le.rhocgs .and. rhocgs.le.1.d-5 .and.
     $           0.001.le.metallicity .and. metallicity.le.0.02 ) then
            opacity = negativeHion(rhocgs,tem,metallicity,
     $           opacity_analytic_warning)
         end if

         ! There's no more tricks to try. We have to admit defeat and notify the
         ! user that we failed to calculate the opacity.
         if ( opacity.eq.0.d0 ) then
            if ( opacity_oob_warning .or. opacity_oob_error) then
               write(*,'(A,ES11.4,", ",ES11.4,A)')
     $              "WARNING (getOpacity.f): (rho,T) = (",
     $              rhocgs,tem,") is in an undefined "//
     $              "region for calculating opacities. Setting "//
     $              "opacity = 0.d0."
               if ( opacity_oob_error ) error stop "getOpacity.f"
            end if
         end if

         ! We succeeded in calculating analytic approximation of opacity
         getOpacity = opacity
         
         return
      end if


      if ( opacity .le. 0.d0 ) then
         write(*,*) "logrho,logR,logtem = ",logrho,logR,logtem
         write(*,*) "numtabused = ",numtabused
         write(*,*) "opacity = ",opacity
         write(*,*) "WHY?!"
         stop
      end if
      
      ! We calculated opacity using the opacityfiles and no analytic
      ! approximations.
      getOpacity = opacity/sumweights

      return
      end function

















      


c     The following functions are approximations and should ONLY
c     be used in emergency situations. The user will be vehemently
c     warned in the output that these functions have been used.
c     These functions have a small degree of accuracy compared to
c     using opacity tables, which is why we are so adament to warn
c     the user.
      
      real*8 function boundfree(rho,T,X,Z,warn)
c     Taken from Onno Pols notes on Stellar Evolution
c     Should be used for Z >= 10^(-3), 1.d4 < T < 1.d6
      implicit none
      real*8 rho, T, X, Z
      logical warn
      if ( warn ) then
         write(*,'(A,ES11.4,",",ES11.4,A)')
     $        " WARNING (getOpacity.f): Using analytic boundfree "//
     $        "opacities for (rho,T) = (",rho,T,")"
      end if
      if ( T.lt.1.d4 .or. T.gt.1.d6 .or.
     $     Z.lt.1.d-3 ) then
         write(*,*) "rho,T,X,Z = ",rho,T,X,Z
         write(*,*) "boundfree opacity can only be used in the "//
     $        "temperature range 1.d4 < T < 1.d6."
         error stop "getOpacity.f"
      end if
      boundfree = 4.3d25*(1.d0 + X)*Z*rho*T**(-7.d0/2.d0)
      return
      end function
      real*8 function freefree(rho,T,X,Z,warn)
c     Taken from Onno Pols notes on Stellar Evolution
c     Should be used for Z < 10^(-3), 1.d4 < T < 1.d6
      implicit none
      real*8 rho, T, X, Z
      logical warn
      if ( warn ) then
         write(*,'(A,ES11.4,",",ES11.4,A)')
     $        " WARNING (getOpacity.f): Using analytic freefree "//
     $        "opacities for (rho,T) = (",rho,T,")"
      end if
      if ( T.lt.1.d4 .or. T.gt.1.d6 .or.
     $     Z.gt.1.d-3 ) then
         write(*,*) "rho,T,X = ",rho,T,X
         write(*,*) "freefree opacity can only be used in the "//
     $        "temperature range 1.d4 < T < 1.d6."
         error stop "getOpacity.f"
      end if
      freefree = 3.8d22*(1.d0+X)*rho*T**(-7.d0/2.d0)
      return
      end function
      real*8 function escattering(rho,T,X,warn)
c     Taken from Onno Pols notes on Stellar Evolution
c     Should only be used for 1.d4 < T < 1.d8
      implicit none
      real*8 rho, T, X
      logical warn
      if ( warn ) then
         write(*,'(A,ES11.4,",",ES11.4,A)')
     $        " WARNING (getOpacity.f): Using analytic electron "//
     $        "scattering opacities for (rho,T) = (",rho,T,")"
      end if
      if ( T.lt.1.d4 .or. T.gt.1.d8 ) then
         write(*,*) "T,X = ",T,X
         write(*,*) "electron scattering opacity can only be used "//
     $        "in the temperature range 1.d4 < T < 1.d8."
         error stop "getOpacity.f"
      end if
      escattering = 0.2d0*(1.d0+X)
      return
      end function
      real*8 function negativeHion(rho,T,Z,warn)
c     Taken from Onno Pols notes on Stellar Evolution
c     Should only be used for 3.d3 < T < 1.d4,
c     1.d-10 < rho < 1.d-5, and 0.001 < Z < 0.02
      implicit none
      real*8 rho,T,Z
      logical warn
      if ( warn ) then
         write(*,'(A,ES11.4,",",ES11.4,A)')
     $        " WARNING (getOpacity.f): Using analytic negative H ion"//
     $        " opacities for (rho,T) = (",rho,T,")"
      end if
      if ( T.lt.3.d3 .or. T.gt.6.d3 .or.
     $     rho.lt.1.d-10 .or. rho.gt.1.d-5 .or.
     $     Z.lt.0.001d0 .or. Z.gt.0.02d0 ) then
         write(*,*) "rho,T,Z = ",rho,T,Z
         write(*,*) "negative H ion opacity can only be used "//
     $        "in the temperature range 3.d3 < T < 6.d3, density"//
     $        " range 1.d-10 < rho < 1.d-5, and metallicity range"//
     $        " 0.001 < Z < 0.02."
         error stop "getOpacity.f"
      end if
      negativeHion = 2.5d-31*(Z/0.02d0)*rho**(0.5d0)*T**9.d0
      return
      end function
         
