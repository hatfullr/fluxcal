      subroutine readineostable(eosfile)
      implicit none
      integer n_lower,n_upper,myrank,nprocs,qthreads,q
      common/nlimits/n_lower,n_upper,nprocs,myrank,qthreads,q
      integer numrho,numu,iu,irho
      integer maxtablesize
      parameter(maxtablesize=1000)
      real*8 rhotable(maxtablesize),
     $     utable(maxtablesize)

      real*8 eostable(maxtablesize,maxtablesize,6)
      real*8 steprho,stepu,rhotable1,utable1,
     $     rhotablelast,utablelast
      common/eoscom/ numrho,numu,rhotable1,utable1,
     $     steprho,stepu,eostable
      integer i
      real*8 steprhotest,steputest,dummy
      real*8 xxx,yyy,zzz,starmu
      common/abundances/XXX,YYY
      logical fileInPath

      character*255 eosfile

      write(*,*) "Reading eos file"
      open (43,file=trim(adjustl(eosfile)))
      read(43,*) xxx
      read(43,*) yyy
      read(43,*) zzz
c     See equation (2-15) of Clayton:
      starmu=1.67262158d-24/(2*xxx+0.75d0*yyy+0.5d0*zzz)
      if(abs(xxx+yyy+zzz-1.d0).gt.1d-16)then
         write(*,*) "ERROR: readineostable.f (L37)"
         write(*,*) "X+Y+Z-1=",xxx+yyy+zzz-1.d0
         stop
      endif
      do i=1,2
         read(43,*)
      enddo
      read(43,*) numrho,rhotable1,rhotablelast,steprho
      read(43,*) numu,utable1,utablelast,stepu
      do i=1,2
         read(43,*)
      enddo

      do irho=1,numrho
         do iu=1,numu
            read(43,*) rhotable(irho),utable(iu),
     $             (eostable(iu,irho,i),i=1,4),dummy,dummy,dummy,
     $              eostable(iu,irho,5),eostable(iu,irho,6)
            eostable(iu,irho,1)=10d0**eostable(iu,irho,1) ! record temperature instead of log temperature
            eostable(iu,irho,2)=eostable(iu,irho,2)*1.67262158d-24 ! record mean molecular mass in cgs units
         enddo
      enddo

      close(43)

      if(rhotable1.ne.rhotable(1))then
         print *,'rhotable(1) mismatch',rhotable1,rhotable(1),
     $        rhotable1-rhotable(1)
         stop
      endif
      if(utable1.ne.utable(1))then
         print *,'utable(1) mismatch'
         stop
      endif
      if(rhotablelast.ne.rhotable(numrho))then
         print *,'rhotable(numrho) mismatch'
         stop
      endif
      if(utablelast.ne.utable(numu))then
         print *,'utable(numu) mismatch'
         stop
      endif

      steprhotest=(rhotable(numrho)-rhotable(1))/(numrho-1)
      steputest=(utable(numu)-utable(1))/(numu-1)

      if(abs(steprhotest-steprho).gt.1d-9) then
         print *,'steprho problem'
         stop
      endif
      if(abs(steputest-stepu).gt.1d-9) then
         print *,'stepu problem'
         stop
      endif


      end
