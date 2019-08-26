      subroutine ini_highT_opacities(filename)
      implicit none
      character*255 filename
      real*8 logT_min, logT_max
      integer Nrho_highT,NT_highT
      parameter(Nrho_highT=37,NT_highT=138)
      integer form, version,logRs,logTs
      real*8 X,Z,logR_min,logR_max
      real*8 logR_highT(Nrho_highT), logT_highT(NT_highT)
      real*8 opacit_highT(Nrho_highT,NT_highT)
      integer i,j
      character*255 fmt
      common/opac_highT/opacit_highT,logR_highT,logT_highT
      
      open(42,file=trim(adjustl(filename)),status='old')
c     At time of writing this code, expected to read in
c     gs98_z0.02_x0.7.data.
      read(42,*)
      read(42,*)

      read(42,*) form,version,X,Z,logRs,logR_min,logR_max,logTs,
     $     logT_min,logT_max

      if(Nrho_highT.ne.logRs .or. NT_highT.ne.logTs)then
         write(*,*) "Nrho_highT=/=logRs or NT_highT=/=logTs in"//
     $        " ini_highT_opacities.f!"
         error stop
      end if
      
      read(42,*)
      read(42,*)

      read(42,*) (logR_highT(i),i=1,logRs)
      read(42,*)
      do i=1, logTs
         read(42,*) logT_highT(i), (opacit_highT(j,i),j=1,logRs)
      end do
      
      close(42)

      return
      end subroutine
