      real*8 function getOpacity(ncooling,tem,rhocgs,tau1,
     $     table,table2)
c     This is a subroutine taken directly from derivs2.f.
c     The point is to make StarSmasher-specific code easier to separate
c     from the main code. In this way, the code will be more modular and
c     we will be able to more easily substitute routines (such as this
c     one) for another routine from a different code.

c     INPUT:
c     ncooling
c     tem    - temperature
c     rhocgs - density in cgs
c     tau1   - some variable that is known as tau(1) in derivs2.f
c     table  - subroutine found in usekappatable.f (usetable).
c              Takes rho,tem, and localross as parameters.
c              To be used when ncooling <= 1
c     table2 - subroutine found in usekappatable.f (usetable2).
c              Takes rho,tem,localross, and localplanck as parameters.
c              To be used when ncooling > 1

      integer ncooling
      real*8 tem,rhocgs,tau1,opacitross,opacitplanck,opacit
      external table,table2,cop

      real*8 eD, eG
      DIMENSION eD(5,6)
      DIMENSION eG(71,71)
      common/opac_dust/eD,eG
      
c     silly test
c         if(tem.le.8000.d0)  then
c            do i=1,10
c               tem=1000.*i
c               call cop(eD,eG,rhocgs,tem,opacit)
c               write(*,*) "COP:", tem,rhocgs,opacit,ncooling               
c               call table(rhocgs,tem,opacit)
c               write(*,*) "TAB:", tem,rhocgs,opacit,ncooling                 
c            end do
c         end if

      
c     We should use both rosseland and planck opacities, always.

      if(tem.le.1000.d0) then
         call cop(eD,eG,rhocgs,tem,opacit)
         getOpacity=opacit
         return
      end if
      
      if(ncooling.le.1)then
c     opacit is local Rosseland opacity
         call table(rhocgs,tem,opacit)
         getOpacity=opacit
      else
         write(*,*) "ncooling = ", ncooling
         write(*,*) "getOpacity.f: ncooling > 1?"
         error stop
c     opacitross is local Rosseland opacity
c     opacitplanck is local Planck opacity
         call table2(rhocgs,tem,opacitross,opacitplanck)
         if(tau1.gt.0d0) then
            getOpacity=exp(-2*tau1)*opacitplanck
     $           +(1d0-exp(-2*tau1))*opacitross
         else
            getOpacity=opacitplanck
         endif
        
      endif

      return
      end function
