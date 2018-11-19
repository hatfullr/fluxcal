      subroutine progressBar(pgress)
      implicit none
      integer(kind=4)::pgress,k
      character(len=17)::bar="???% |          |"
      write(unit=bar(1:3),fmt="(i3)") pgress
      do k=1, pgress/10
         bar(6+k:6+k)="*"
      enddo
! print the progress bar.
      write(*,fmt="(a1,a17)",advance="no") char(13), bar
      if (pgress/=100) then
         flush(6)
      else
         write(*,*)
      endif
      return
      
      end subroutine
