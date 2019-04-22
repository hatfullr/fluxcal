      real*8 function getOpacity(tem,rhocgs,tau1,table,table2)
c     This is a subroutine taken directly from derivs2.f.
c     The point is to make StarSmasher-specific code easier to separate
c     from the main code. In this way, the code will be more modular and
c     we will be able to more easily substitute routines (such as this
c     one) for another routine from a different code.

c     INPUT:
c     tem    - temperature
c     rhocgs - density in cgs
c     tau1   - some variable that is known as tau(1) in derivs2.f
c     table  - subroutine found in usekappatable.f (usetable).
c              Takes rho,tem, and localross as parameters.
c              To be used when ncooling <= 1
c     table2 - subroutine found in usekappatable.f (usetable2).
c              Takes rho,tem,localross, and localplanck as parameters.
c              To be used when ncooling > 1

      real*8 tem,rhocgs,tau1,opacitross,opacitplanck,opacit
      external table,table2,cop

      real*8 eD, eG
      DIMENSION eD(5,6)
      DIMENSION eG(71,71)
      common/opac_dust/eD,eG
      real*8 Topac_Planck
      common/Topacity_Planck/ Topac_Planck
      
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

      
      if(tem.le.Topac_Planck) then
         call cop(eD,eG,rhocgs,tem,opacit)
      else
         call table(rhocgs,tem,opacit)
      end if
      
      getOpacity=opacit

c     This is old code, but I want to keep it here to have easy
c     access to the smoothing function Jamie wrote
c      if(rossonly) then         ! User-set input file variable
c         ! Rosseland opacities only
c         call table(rhocgs,tem,opacit)
c         getOpacity=opacit
c      else
c         ! Rosseland and Planck opacities transfer function
c         ! opacitross is local Rosseland opacity
c         ! opacitplanck is local Planck opacity
c         call table2(rhocgs,tem,opacitross,opacitplanck)
c         if(tau1.gt.0d0) then
c            getOpacity=exp(-2*tau1)*opacitplanck
c     $           +(1d0-exp(-2*tau1))*opacitross
c         else
c            getOpacity=opacitplanck
c         endif
c        
c      endif
      
      return
      end function
