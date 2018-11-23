      subroutine integrationAtAllPos
      include 'optical_depth.h'
      character*255 myfname
      integer i
      real*8 posx_temp,posy_temp

      call useDimenFile
      call prepareIntegration

 700  format('int_at_all_pos_',i4.4,'.dat')
 701  format('int_at_all_pos_',i5.5,'.dat')
 702  format('int_at_all_pos_',i6.6,'.dat')
      if(innit.le.9999) then
         write(myfname,700) innit
      else if(innit.le.99999) then
         write(myfname,701) innit
      else
         write(myfname,702) innit
      end if

      posx_temp=posx
      posy_temp=posy
      
      intout=655
      open(intout,file=trim(adjustl(myfname)),status='unknown')
      write(*,*) "Writing to '",trim(adjustl(myfname)),"'"
      
      dointatallpos=.true.
      call integrateTau
      dointatallpos=.false.

      close(intout)

      posx=posx_temp
      posy=posy_temp
      
      end subroutine
      
