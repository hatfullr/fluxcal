      subroutine integrationAtPos
      include 'optical_depth.h'

      character*255 myfname
      real*8 mytau,mt1,mt2
      integer i
      
c      call createGrid

      xminmap=posx
      hxmap=1.d0
      nxmap=1
      yminmap=posy
      hymap=1.d0
      nymap=1

      xmaxmap=xminmap+hxmap*(nxmap-1)
      ymaxmap=yminmap+hymap*(nymap-1)

      call prepareIntegration

      if(innit.le.9999) then
         write(myfname,"('int_at_pos_',i4.4,'.dat')") innit
      else if(innit.le.99999) then
         write(myfname,"('int_at_pos_',i5.5,'.dat')") innit
      else
         write(myfname,"('int_at_pos_',i6.6,'.dat')") innit
      end if

      i=nint((posx-xminmap)/hxmap + 1)
      j=nint((posy-yminmap)/hymap + 1)
c      xpos=(i-1)*hxmap+xminmap  ! x of line of sight should be close to posx
c      ypos=(j-1)*hymap+yminmap  ! y of line of sight should be close to posy
c      hzmap=(zmax(i,j)-zmin(i,j))/(nzmap-1)

      intout=51
      open(intout,file=myfname,status='unknown')
 801  format(2E15.7)
 802  format(11A15)
      write(*,*) "Writing to ", trim(adjustl(myfname))
      write(intout,801) posx, posy
c      write(intout,801) posy,"posy"
c      write(intout,801) zmax(i,j), "zmax [cm]"
c      write(intout,801) zmin(i,j), "zmin [cm]"
c      write(intout,801) taulimit, "taulimit"
c      write(intout,801) tau_thick, "tau_thick"
      write(intout,802) "z","h","rho","u","g","mu","P","T","kappa","tau",
     $     "particle"


      call integrateTau

      close(intout)
      
      end subroutine
