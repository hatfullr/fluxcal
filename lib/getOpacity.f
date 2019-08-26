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
c     table  - subroutine found in usekappatable.f (usetable).
c              Takes rho,tem, and localross as parameters.
c              To be used when ncooling <= 1
c     table2 - subroutine found in usekappatable.f (usetable2).
c              Takes rho,tem,localross, and localplanck as parameters.
c              To be used when ncooling > 1
      use bilinear_interpolation
      real*8 tem,rhocgs,tau1,opacitross,opacitplanck,opacit
      external table,table2,cop

      real*8 eD, eG
      DIMENSION eD(5,6)
      DIMENSION eG(71,71)
      common/opac_dust/eD,eG
      real*8 Topac_Planck
      common/Topacity_Planck/ Topac_Planck

      integer Nrho_highT,NT_highT
      parameter(Nrho_highT=37,NT_highT=138)
      real*8 logR_highT(Nrho_highT), logT_highT(NT_highT)
      real*8 opacit_highT(Nrho_highT,NT_highT)
      common/opac_highT/opacit_highT,logR_highT,logT_highT
      integer Nrho_lowT,NT_lowT
      parameter(Nrho_lowT=46,NT_lowT=105)
      real*8 logR_lowT(Nrho_lowT), logT_lowT(NT_lowT)
      real*8 opacit_lowT(Nrho_lowT,NT_lowT)
      common/opac_lowT/opacit_lowT,logR_lowT,logT_lowT

c      real*8 opac_high,opac_low

      real*8 logR,logT,dlogT
      real*8 kappa_highT,kappa_lowT

      logical exist
      
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
      
      logR = log10(rhocgs)-3.d0*log10(tem)+18.d0

      if(logR.gt.max(logR_highT(Nrho_highT),logR_lowT(Nrho_lowT)) .or.
     $  (logR.lt.min(logR_highT(1),logR_lowT(1)))) then
         write(*,*) "Warning: out of range of highT opacity file for"
         write(*,*) "(logR,logT)=(",logR,",",log10(tem),")"
         write(*,*) "rhocgs = ",rhocgs
         write(*,*) "setting opacit=0.d0. You should double check"
         write(*,*) "all your particles and integrations for any"//
     $        " suspicious opacities."
         getOpacity = 0.d0
         return
      end if

      logT = log10(tem)
c      opac_low = nan
c      opac_high = nan
      
      if((logT.le.logT_lowT(NT_lowT) .and. logT.ge.logT_lowT(1)) ! We are in range of the lowT file
     $     .and.(logT.lt.logT_highT(1))) then ! We are exclusively inside the lowT file
         
         getOpacity = 10.d0**bilinear_interpolate(Nrho_lowT,
     $        logR_lowT,NT_lowT,logT_lowT,opacit_lowT,logR,logT)

c         opac_low = 10.d0**bilinear_interpolate(Nrho_lowT,
c     $     logR_lowT,NT_lowT,logT_lowT,opacit_lowT,logR,logT)
         
      else if((logT.le.logT_highT(NT_highT).and.(logT.ge.logT_highT(1))) ! We are in range of the highT file
     $        .and.(logT.gt.logT_lowT(NT_lowT))) then ! We are exclusively inside the highT file

         getOpacity = 10.d0**bilinear_interpolate(Nrho_highT,
     $        logR_highT,NT_highT,logT_highT,opacit_highT,logR,logT)

c         opac_high = 10.d0**bilinear_interpolate(Nrho_highT,
c     $     logR_highT,NT_highT,logT_highT,opacit_highT,logR,logT)

      else if((logT.lt.logT_lowT(NT_lowT))
     $        .and.(logT.gt.logT_highT(1))) then ! We are in-between the lowT and highT files

         ! Let's smooth between the two tables by doing a weighted average by how close logT
         ! is to the boundaries of the tables. The closer logT is to a boundary, the more
         ! even the weighting becomes.
         dlogT = logT_lowT(NT_lowT) - logT_highT(1)
         
         kappa_lowT = 10.d0**bilinear_interpolate(Nrho_lowT,
     $        logR_lowT,NT_lowT,logT_lowT,opacit_lowT,logR,logT)
         kappa_highT = 10.d0**bilinear_interpolate(Nrho_highT,
     $        logR_highT,NT_highT,logT_highT,opacit_highT,logR,logT)

         getOpacity = kappa_lowT * (logT_lowT(NT_lowT)-logT)/dlogT +
     $        kappa_highT * (logT-logT_highT(1))/dlogT

c         write(55,*) logT, (logT_lowT(NT_lowT)-logT)/dlogT,
c     $        (logT-logT_highT(1))/dlogT
         
      else ! Out of range of both files!
         write(*,*) "Warning: out of range of highT opacity file for"
         write(*,*) "(logR,logT)=(",logR,",",log10(tem),")"
         write(*,*) "rhocgs = ",rhocgs
         write(*,*) "setting opacit=0.d0. You should double check"
         write(*,*) "all your particles and integrations for any"//
     $        " suspicious opacities."
         getOpacity = 0.d0
      end if

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
