      subroutine readinkappatable(filename)
      include 'flux_cal.h'
      integer numrho,numtem,item,irho
      integer maxtablesize
      parameter(maxtablesize=500)
      real*8 rhotable(maxtablesize),temtable(maxtablesize)
      real dummy
      real localrosstable(maxtablesize,maxtablesize)
      real localplancktable(maxtablesize,maxtablesize)
      real*8 steprho,steptem,rhotable1,temtable1
      common/kappacom/ numrho,numtem,rhotable1,temtable1,
     $     steprho,steptem,localrosstable,
     $     localplancktable
      logical fileexists
      character*255 filename
      
c      inquire(file=opacityfile,exist=fileexists)
c      if(.not.fileexists) then
c         write(*,*) "Opacity file not found"
c         error stop "usekappatable.f"
c      end if
      write(*,*) "Reading opacity file"
      open (43,file=trim(adjustl(filename)))
      read(43,*)
      read(43,*)
      read(43,*) numrho,numtem

      if(rossonly) then
         do irho=1,numrho
            do item=1,numtem
               read(43,*) rhotable(irho),temtable(item),
     $              localrosstable(irho,item)
            enddo
         enddo
      else
         do irho=1,numrho
            do item=1,numtem
               read(43,*) rhotable(irho),temtable(item),
     $              localrosstable(irho,item),
     $              dummy,localplancktable(irho,item)
            enddo
         enddo
      endif


      rhotable1=rhotable(1)
      temtable1=temtable(1)

      steprho=(rhotable(numrho)-rhotable(1))/(numrho-1)
      steptem=(temtable(numtem)-temtable(1))/(numtem-1)

      end

      subroutine readinkappatabledust(filename)
      include 'flux_cal.h'
      integer numrhodust,numtemdust,itemdust,irhodust
      integer maxtablesize
      parameter(maxtablesize=500)
      real*8 rhotabledust(maxtablesize),temtabledust(maxtablesize)
      real localrosstabledust(maxtablesize,maxtablesize)
      real localplancktabledust(maxtablesize,maxtablesize)
      real*8 steprhodust,steptemdust,rhotabledust1,temtabledust1
      common/kappacomdust/numrhodust,numtemdust,rhotabledust1,
     $     temtabledust1,steprhodust,steptemdust,
     $     localrosstabledust,localplancktabledust
      real dummy
      logical fileexists
      character*255 filename
      
c      inquire(file=opacitydustfile,exist=fileexists)
c      if(.not.fileexists) then
c         write(*,*) "Opacity file for dust not found"
c         error stop "usekappatable.f"
c      end if
      write(*,*) "Reading opacity file for dust"
      open (44,file=trim(adjustl(filename)))
      read(44,*)
      read(44,*)
      read(44,*) numrhodust,numtemdust


      if(rossonly) then
         do irhodust=1,numrhodust
            do itemdust=1,numtemdust
               read(44,*) rhotabledust(irhodust),temtabledust(itemdust),
     $              localrosstabledust(irhodust,itemdust)
            enddo
         enddo
      else
         do irhodust=1,numrhodust
            do itemdust=1,numtemdust
               read(44,*) rhotabledust(irhodust),temtabledust(itemdust),
     $              localrosstabledust(irhodust,itemdust),
     $              dummy,localplancktabledust(irhodust,itemdust)
            enddo
         enddo
      endif

      rhotabledust1=rhotabledust(1)
      temtabledust1=temtabledust(1)

      steprhodust=(rhotabledust(numrhodust)-rhotabledust(1))
     $     /(numrhodust-1)
      steptemdust=(temtabledust(numtemdust)-temtabledust(1))
     $     /(numtemdust-1)

      end
      
      subroutine usetable(rho,tem,localross)
c     Use this routine for the rossonly=.true. case
      implicit none
      real*8 rho,tem
      real*8 localross
      integer item,irho
      integer numrho,numtem
      real*8 temtable1,rhotable1
      integer maxtablesize
      parameter(maxtablesize=500)
      real localrosstable(maxtablesize,maxtablesize)
      real localplancktable(maxtablesize,maxtablesize)

      real*8 steprho,steptem
      common/kappacom/ numrho,numtem,rhotable1,temtable1,
     $     steprho,steptem,localrosstable
     $     ,localplancktable

      real*8 f00,f01,f10,f11,log10rho,log10tem,
     $     rholow,rhohigh,temlow,temhigh

c     temtable(item)=temtable(1)+(item-1)*steptem
c     so, item= (temtable(item)-temtable(1))/steptem + 1

      log10rho=log10(rho)
      log10tem=log10(tem)

      irho = min(max(1,int((log10rho-rhotable1)/steprho + 1)),numrho-1)

      item = min(int((log10tem-temtable1)/steptem + 1) ,numtem-1)

      rholow=log10rho-(rhotable1+(irho-1)*steprho)
      rhohigh=rhotable1+irho*steprho-log10rho
      
c      write(90,*)rho,irho,rholow,rhohigh
c      write(90,*)tem,item

      if(item.ge.1) then

         temlow=log10tem-(temtable1+(item-1)*steptem)
         temhigh=temtable1+item*steptem-log10tem
         
c     Use bi-linear interpolation among the four cartesian
c     grid point (irho,item), (irho+1,item), (irho,item+1), and (irho+1,item+1)
         f00=rholow*temlow
         f10=rhohigh*temlow
         f01=rholow*temhigh
         f11=rhohigh*temhigh
         
         localross=(f00*localrosstable(irho+1,item+1)
     $        + f10*localrosstable(irho,  item+1)
     $        + f01*localrosstable(irho+1,item)
     $        + f11*localrosstable(irho,  item))/(steprho*steptem)

      else
         localross=(rholow*localrosstable(irho+1,1)
     $        + rhohigh*localrosstable(irho,  1))/steprho
     $        *(tem/10d0**temtable1)**2
      endif

      end

      subroutine usetable2(rho,tem,localross,localplanck)
c     Use this routine for the rossonly=.false. case
      implicit none
      real*8 rho,tem
      real*8 localross,localplanck
      integer item,irho
      integer numrho,numtem
      real*8 temtable1,rhotable1
      integer maxtablesize
      parameter(maxtablesize=500)
      real localrosstable(maxtablesize,maxtablesize)
      real localplancktable(maxtablesize,maxtablesize)

      real*8 steprho,steptem
      common/kappacom/ numrho,numtem,rhotable1,temtable1,
     $     steprho,steptem,localrosstable,
     $     localplancktable

      real*8 f00,f01,f10,f11,log10rho,log10tem,
     $     rholow,rhohigh,temlow,temhigh

c     temtable(item)=temtable(1)+(item-1)*steptem
c     so, item= (temtable(item)-temtable(1))/steptem + 1

      log10rho=log10(rho)
      log10tem=log10(tem)

      irho = min(max(1,int((log10rho-rhotable1)/steprho + 1)),numrho-1)

      item = min(int((log10tem-temtable1)/steptem + 1) ,numtem-1)

c      print *, numrho,numtem,rhotable1,temtable1,steprho,steptem,
c     $     localrosstable(1,1),localplancktable(1,1)

      rholow=log10rho-(rhotable1+(irho-1)*steprho)
      rhohigh=rhotable1+irho*steprho-log10rho
      
c      write(90,*)rho,irho,rholow,rhohigh
c      write(90,*)tem,item

      if(item.ge.1) then

         temlow=log10tem-(temtable1+(item-1)*steptem)
         temhigh=temtable1+item*steptem-log10tem
         
c     Use bi-linear interpolation among the four cartesian
c     grid point (irho,item), (irho+1,item), (irho,item+1), and (irho+1,item+1)
         f00=rholow*temlow
         f10=rhohigh*temlow
         f01=rholow*temhigh
         f11=rhohigh*temhigh
         
         localross=(f00*localrosstable(irho+1,item+1)
     $        + f10*localrosstable(irho,  item+1)
     $        + f01*localrosstable(irho+1,item)
     $        + f11*localrosstable(irho,  item))/(steprho*steptem)
         localplanck=(f00*localplancktable(irho+1,item+1)
     $        + f10*localplancktable(irho,  item+1)
     $        + f01*localplancktable(irho+1,item)
     $        + f11*localplancktable(irho,  item))/(steprho*steptem)

      else
         localross=(rholow*localrosstable(irho+1,1)
     $        + rhohigh*localrosstable(irho,  1))/steprho
     $        *(tem/10d0**temtable1)**2
         localplanck=(rholow*localplancktable(irho+1,1)
     $        + rhohigh*localplancktable(irho,  1))/steprho
     $        *(tem/10d0**temtable1)**2
      endif

c      print *, item, localross,localplanck
c      stop

      end

      subroutine usetabledust(rho,tem,localross)
c     Use this routine for the rossonly=.true. case
      implicit none
      real*8 rho,tem
      real*8 localross
      integer item,irho
      integer numrhodust,numtemdust
      real*8 temtabledust1,rhotabledust1
      integer maxtablesize
      parameter(maxtablesize=500)
      real localrosstabledust(maxtablesize,maxtablesize)
      real localplancktabledust(maxtablesize,maxtablesize)

      real*8 steprhodust,steptemdust
      common/kappacomdust/numrhodust,numtemdust,rhotabledust1,
     $     temtabledust1,steprhodust,steptemdust,
     $     localrosstabledust,localplancktabledust

      real*8 f00,f01,f10,f11,log10rho,log10tem,
     $     rholow,rhohigh,temlow,temhigh

c     temtable(item)=temtable(1)+(item-1)*steptemdust
c     so, item= (temtable(item)-temtable(1))/steptemdust + 1

      log10rho=log10(rho)
      log10tem=log10(tem)

      irho=min(max(1,int((log10rho-rhotabledust1)/steprhodust)+1),
     $     numrhodust-1)

      item=min(int((log10tem-temtabledust1)/steptemdust)+1,numtemdust-1)

      rholow=log10rho-(rhotabledust1+(irho-1)*steprhodust)
      rhohigh=rhotabledust1+irho*steprhodust-log10rho
      
c      write(90,*)rho,irho,rholow,rhohigh
c      write(90,*)tem,item

      if(item.ge.1) then

         temlow=log10tem-(temtabledust1+(item-1)*steptemdust)
         temhigh=temtabledust1+item*steptemdust-log10tem
         
c     Use bi-linear interpolation among the four cartesian
c     grid point (irho,item), (irho+1,item), (irho,item+1), and (irho+1,item+1)
         f00=rholow*temlow
         f10=rhohigh*temlow
         f01=rholow*temhigh
         f11=rhohigh*temhigh
         
         localross=(f00*localrosstabledust(irho+1,item+1)
     $        + f10*localrosstabledust(irho,item+1)
     $        + f01*localrosstabledust(irho+1,item)
     $        + f11*localrosstabledust(irho,item))/
     $        (steprhodust*steptemdust)

      else
         localross=(rholow*localrosstabledust(irho+1,1)
     $        + rhohigh*localrosstabledust(irho,  1))/steprhodust
     $        *(tem/10d0**temtabledust1)**2

      endif

      end

      subroutine usetabledust2(rho,tem,localross,localplanck)
c     Use this routine for the rossonly=.false. case
      implicit none
      real*8 rho,tem
      real*8 localross,localplanck
      integer item,irho
      integer numrhodust,numtemdust
      real*8 temtabledust1,rhotabledust1
      integer maxtablesize
      parameter(maxtablesize=500)
      real localrosstabledust(maxtablesize,maxtablesize)
      real localplancktabledust(maxtablesize,maxtablesize)

      real*8 steprhodust,steptemdust
      common/kappacomdust/numrhodust,numtemdust,rhotabledust1,
     $     temtabledust1,steprhodust,steptemdust,
     $     localrosstabledust,localplancktabledust

      real*8 f00,f01,f10,f11,log10rho,log10tem,
     $     rholow,rhohigh,temlow,temhigh

c     temtable(item)=temtable(1)+(item-1)*steptemdust
c     so, item= (temtable(item)-temtable(1))/steptemdust + 1

      log10rho=log10(rho)
      log10tem=log10(tem)

      irho=min(max(1,int((log10rho-rhotabledust1)/steprhodust)+1),
     $     numrhodust-1)

      item=min(int((log10tem-temtabledust1)/steptemdust)+1,numtemdust-1)

      rholow=log10rho-(rhotabledust1+(irho-1)*steprhodust)
      rhohigh=rhotabledust1+irho*steprhodust-log10rho
      
c      write(90,*)rho,irho,rholow,rhohigh
c      write(90,*)tem,item

      if(item.ge.1) then

         temlow=log10tem-(temtabledust1+(item-1)*steptemdust)
         temhigh=temtabledust1+item*steptemdust-log10tem
         
c     Use bi-linear interpolation among the four cartesian
c     grid point (irho,item), (irho+1,item), (irho,item+1), and (irho+1,item+1)
         f00=rholow*temlow
         f10=rhohigh*temlow
         f01=rholow*temhigh
         f11=rhohigh*temhigh
         
         localross=(f00*localrosstabledust(irho+1,item+1)
     $        + f10*localrosstabledust(irho,item+1)
     $        + f01*localrosstabledust(irho+1,item)
     $        + f11*localrosstabledust(irho,item))/
     $        (steprhodust*steptemdust)
         localplanck=(f00*localplancktabledust(irho+1,item+1)
     $        + f10*localplancktabledust(irho,item+1)
     $        + f01*localplancktabledust(irho+1,item)
     $        + f11*localplancktabledust(irho,item))/
     $        (steprhodust*steptemdust)

      else
         localross=(rholow*localrosstabledust(irho+1,1)
     $        + rhohigh*localrosstabledust(irho,  1))/steprhodust
     $        *(tem/10d0**temtabledust1)**2
         localplanck=(rholow*localplancktabledust(irho+1,1)
     $        + rhohigh*localplancktabledust(irho,  1))/steprhodust
     $        *(tem/10d0**temtabledust1)**2

      endif

      end
