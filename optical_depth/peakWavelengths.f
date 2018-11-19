      subroutine peakWavelengths
c     L908 to L950
      include 'optical_depth.h'
      
      do insl=0,nsl
         peak=0.d0
         do ilogwavelength=ilogwavelengthmin,ilogwavelengthmax
            if(spectrum(ilogwavelength,insl).gt.peak) then
               peak=spectrum(ilogwavelength,insl)
               ilambdamax=ilogwavelength
            endif
         enddo

         if(ilambdamax.gt.ilogwavelengthmin .and.
     $        ilambdamax.lt.ilogwavelengthmax)then
            ratio=(spectrum(ilambdamax,insl)
     $            -spectrum(ilambdamax-1,insl))/
     $           (spectrum(ilambdamax,insl)
     $            -spectrum(ilambdamax+1,insl))

            if(ratio.lt.0) then
               print *,'ratio problem...'
               stop
            endif

            print *, insl,ilambdamax,ratio,0.5d0*(ratio-1)/(ratio+1)
            loglambdamax=ilambdamax+0.5d0*(ratio-1)/(ratio+1) ! This
!     makes a simple estimate as to where the peak actually occurs
!     If ratio is very small, then the peak is likely halfway toward the
!     next lower wavelength.
!     If the ratio is 1, then the peak is likely right of the wavelength
!     we found.
!     If ratio is very large, then the peak is likely halfway toward the
!     next higher wavelength.
         else
            loglambdamax=ilambdamax
         endif

         logwavelength=0.01d0*loglambdamax ! logwavelength has wavelength
!                                            in micrometers
         lambdamax=10**logwavelength *1d-4 ! lambdamax is in cm
         if(peak.ne.0) then
            Tfit(insl)=0.2897755d0/lambdamax
         else
            Tfit(insl)=0.d0
         endif
      enddo

      end subroutine
