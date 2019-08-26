      subroutine ini_lowT_opacities(filename)
      implicit none
      character*255 filename
      real*8 logT_min, logT_max
      integer Nrho_lowT,NT_lowT
      parameter(Nrho_lowT=46,NT_lowT=105)
      integer form, version,logRs,logTs
      real*8 X,Z,logR_min,logR_max
      real*8 logR_lowT(Nrho_lowT), logT_lowT(NT_lowT)
      real*8 opacit_lowT(Nrho_lowT,NT_lowT)
      integer i,j
      character*255 fmt
      common/opac_lowT/opacit_lowT,logR_lowT,logT_lowT
      
      open(42,file=trim(adjustl(filename)),status='old')
c     At time of writing this code, expected to read in
c     lowT_fa05_gs98_z0.02_x0.7.data
      read(42,*)
      read(42,*)

      read(42,*) form,version,X,Z,logRs,logR_min,logR_max,logTs,
     $     logT_min,logT_max

      if(Nrho_lowT.ne.logRs .or. NT_lowT.ne.logTs)then
         write(*,*) "Nrho_lowT, logRs = ",Nrho_lowT,logRs
         write(*,*) "NT_lowT, logTs   = ",NT_lowT,logTs
         write(*,*) "Nrho_lowT=/=logRs or NT_lowT=/=logTs in"//
     $        " ini_lowT_opacities.f!"
         error stop
      end if
      
      read(42,*)
      read(42,*)

      read(42,*) (logR_lowT(i),i=1,logRs)
      read(42,*)
      do i=1, logTs
         read(42,*) logT_lowT(i), (opacit_lowT(j,i),j=1,logRs)
      end do
      
      close(42)

      return
      end subroutine
