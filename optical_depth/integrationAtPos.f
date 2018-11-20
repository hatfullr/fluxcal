      subroutine integrationAtPos
      include 'optical_depth.h'

      character*255 myfname
      real*8 mytau,mt1,mt2
      integer i
      
c      call createGrid

c     hxmap and hymap need to be exactly equal to 1 in whatever units
c     are being used. This is so that prepareIntegration correctly sets
c     xpos=posx and ypos=posy in its calculations.
      xminmap=posx
      hxmap=runit
      nxmap=1
      yminmap=posy
      hymap=runit
      nymap=1

      call prepareIntegration

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
      write(intout,801) posx/runit_out, posy/runit_out
      write(intout,802) "z","h","rho","u","g","mu","P","T","kappa","tau",
     $     "particle"


      call integrateTau

      close(intout)
      
      end subroutine
