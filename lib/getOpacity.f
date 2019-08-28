      real*8 function getOpacity(tem,rhocgs,tau1)
c     This is a subroutine taken directly from derivs2.f.
c     The point is to make StarSmasher-specific code easier to separate
c     from the main code. In this way, the code will be more modular and
c     we will be able to more easily substitute routines (such as this
c     one) for another routine from a different code.

c     INPUT:
c     tem    - temperature
c     rhocgs - density in cgs
c     tau1   - some variable that is known as tau(1) in derivs2.f
      
      real*8 tem,rhocgs,tau1,opacit
      external cop

      real*8 eD_R, eG_R, eD_P, eG_P
      DIMENSION eD_R(5,6), eD_P(5,6)
      DIMENSION eG_R(71,71), eG_P(71,71)
      common/opac/eD_R,eG_R,eD_P,eG_P
      real*8 Tplanck
      common/Topacity_Planck/ Tplanck

      
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

c      inquire(file="testing.dat",exist=exist)
c      if(exist)then
c         open(55,file="testing.dat",status='old',position='append',
c     $        action='write')
c      else
c         open(55,file="testing.dat",status='new',action='write')
c      end if
      
      if(rhocgs.lt.1.d-30) then
         ! There is no material in this case
         getOpacity = 0.d0
         return
      end if


      opacit = 0.d0
      if(tem.le.Tplanck) then
         ! Use Planck opacities
         call cop(eD_P,eG_P,rhocgs,tem,opacit)
      else
         ! Use Rosseland opacities
         call cop(eD_R,eG_R,rhocgs,tem,opacit)
      end if
      
      getOpacity = opacit

c      close(55)
      
c      if(tem.le.Topac_Planck) then
c         call cop(eD,eG,rhocgs,tem,opacit)
c      else
c      end if
      
c      getOpacity=opacit

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
