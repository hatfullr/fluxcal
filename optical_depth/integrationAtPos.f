      subroutine integrationAtPos
      include 'optical_depth.h'

      character*255 myfname
      real*8 mytau,mt1,mt2
      integer i
      real*8 xminmap_temp,hxmap_temp,yminmap_temp,hymap_temp
      integer nxmap_temp,nymap_temp
      real*8 xpos_temp,ypos_temp
      
c      call createGrid

      xpos_temp = xpos
      ypos_temp = ypos
      xpos = posx
      ypos = posy
      
c     hxmap and hymap need to be exactly equal to 1 in whatever units
c     are being used. This is so that prepareIntegration correctly sets
c     xpos=posx and ypos=posy in its calculations.
      xminmap_temp = xminmap
      hxmap_temp = hxmap
      nxmap_temp = nxmap
      yminmap_temp = yminmap
      hymap_temp = hymap
      nymap_temp = nymap
      
      xminmap=xpos
      hxmap=runit
      nxmap=1
      yminmap=ypos
      hymap=runit
      nymap=1

      prepareIntegrationCalled=.false. ! Pretend we haven't called yet
      call prepareIntegration
      ! Reset this boolean because we only prepared the integration for
      ! one single point and not the whole driving grid.
      prepareIntegrationCalled=.false.

      if(innit.le.9999) then
         write(myfname,"('int_at_pos_',i4.4,'.dat')") innit
      else if(innit.le.99999) then
         write(myfname,"('int_at_pos_',i5.5,'.dat')") innit
      else
         write(myfname,"('int_at_pos_',i6.6,'.dat')") innit
      end if

      intout=51
      open(intout,file=trim(adjustl(myfname)),status='unknown')
 801  format(2ES22.14)
 802  format(11A22)
      write(*,*) "Writing to ", trim(adjustl(myfname))
c      write(intout,801) posx/runit_out, posy/runit_out
c      write(intout,802)"z","h","rho","u","g","mu","P","T","kappa","tau",
c     $     "particle"
 101  format(13A22)
      write(intout,101) "posx","posy","posz","<h>","<rho>",
     $     "<u>","<mu>","<g>","<T_SPH>","<P>","opacity","sumtau",
     $     "lastpart"
      
      printIntegrationSteps=.true.
      call integrateTau
      printIntegrationSteps=.false.
      
      close(intout)

      xminmap = xminmap_temp
      hxmap = hxmap_temp
      nxmap = nxmap_temp
      yminmap = yminmap_temp
      hymap = hymap_temp
      nymap = nymap_temp

      xpos = xpos_temp
      ypos = ypos_temp
      
      end subroutine
