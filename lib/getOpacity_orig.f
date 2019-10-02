      function getOpacity(tem,rhocgs,tau1) ! real*8
      use bilinear_interpolation
      include 'flux_cal.h'
c     INPUT:
c     tem    - temperature
c     rhocgs - density in cgs
c     tau1   - some variable that is known as tau(1) in derivs2.f
      
      real*8 tem,rhocgs,tau1,logR,logtem,logrho

      logtem = log10(tem)
      logrho = log10(rhocgs)
      logR = logrho - 3.d0*logtem + 18.d0

      ! The following if statement only handles the boundary conditions
      ! defined in the user's input file
      if ( logtem.lt.minval(logT_opac_array) .or.
     $     logtem.gt.maxval(logT_opac_array) .or.
     $     logR.lt.minval(logR_opac_array) .or.
     $     logR.gt.maxval(logR_opac_array) ) then
         if ( opacity_oob_warning ) then
            write(*,*) "WARNING: (rho,T) is out of range of your "//
     $           "opacityfiles."
            write(*,*) "rho, R, T = ",rhocgs,10.d0**(logR),T
            write(*,*) "minT,maxT = ",10.d0**minval(logT_opac_array),
     $           10.d0**maxval(logT_opac_array)
            write(*,*) "minR,maxR = ",10.d0**minval(logR_opac_array),
     $           10.d0**maxval(logR_opac_array)
         end if
         if ( opacity_oob_error ) then
            write(*,*) "rho, R, T = ",rhocgs,10.d0**(logR),T
            write(*,*) "minT,maxT = ",10.d0**minval(logT_opac_array),
     $           10.d0**maxval(logT_opac_array)
            write(*,*) "minR,maxR = ",10.d0**minval(logR_opac_array),
     $           10.d0**maxval(logR_opac_array)
            write(*,*) "(rho,T) is out of range of your opacityfiles."
            error stop "getOpacity.f"
         end if
      end if
      
      ! T is too low, rho is in bounds
      if ( tem.lt.opacity_low_T .and.
     $     rhocgs.gt.opacity_low_rho .and.
     $     rhocgs.lt.opacity_high_rho ) then
         getOpacity=opacity_at_low_T
      ! T is too low, rho is too low
      else if(tem.lt.opacity_low_T .and.
     $        rhocgs.lt.opacity_low_rho ) then
         if(opacity_oob_prioritize_rho) then
            getOpacity = opacity_at_low_rho
         else
            getOpacity = opacity_at_low_T
         end if
      ! T is too low, rho is too high
      else if(tem.lt.opacity_low_T .and.
     $        rhocgs.gt.opacity_low_rho)then
         if(opacity_oob_prioritize_rho) then
            getOpacity = opacity_at_high_rho
         else
            getOpacity = opacity_at_low_T
         end if
      ! T is too high, rho is in bounds
      else if(tem.gt.opacity_high_T .and.
     $        rhocgs.ge.opacity_low_rho .and.
     $        rhocgs.le.opacity_high_rho) then
         getOpacity = opacity_at_high_T
      ! T is too high, rho is too low
      else if(tem.gt.opacity_high_T .and.
     $        rhocgs.lt.opacity_low_rho) then
         if(opacity_oob_prioritize_rho) then
            getOpacity = opacity_at_low_rho
         else
            getOpacity = opacity_at_high_T
         end if
      ! T is too high, rho is too high
      else if(tem.gt.opacity_high_T .and.
     $        rhocgs.gt.opacity_high_rho) then
         if(opacity_oob_prioritize_rho) then
            getOpacity = opacity_at_high_rho
         else
            getOpacity = opacity_at_high_T
         end if
      ! T is in bounds, rho is too low
      else if(tem.ge.opacity_low_T .and.
     $        tem.le.opacity_high_T .and.
     $        rhocgs.lt.opacity_low_rho)then
         getOpacity = opacity_at_low_rho
      ! T is in bounds, rho is too high
      else if(tem.ge.opacity_low_T .and.
     $        tem.le.opacity_high_T .and.
     $        rhocgs.gt.opacity_high_rho)then
         getOpacity = opacity_at_high_rho
      else
         getOpacity = bilinear_interpolate(
     $        nrows_opac_table,
     $        logT_opac_array,
     $        ncols_opac_table,
     $        logR_opac_array,
     $        logopacity_highT_table,
     $        logtem,logR)
      end if

      return
      end function
