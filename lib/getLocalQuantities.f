      subroutine getLocalQuantities(xpos,ypos,zpos)
      include 'flux_cal.h'
      real*8 R2, R2MAX, WPC, XPOS, YPOS, ZPOS, dens, rhocgs,xhp,
     $     gcgs,pcgs,tcgs
      integer INDEX
      real xh,t6, zz, uu, gg,ppp,tt
      real*8 ucgs,temperatur
      common/localQuantities/ rhocgs,xh,t6, xhp,ucgs,gcgs,pcgs,tcgs
      real xhi
c      real*8 useeostable
      integer NXMAPMAX,NYMAPMAX,NZMAP
c      PARAMETER (NXMAPMAX=999,NYMAPMAX=999,NZMAP=202)
c      real*8 rhoxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      PARAMETER (NXMAPMAX=499,NYMAPMAX=499,NZMAP=416)
      real rhoxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real uxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real hpxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real xhxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real gxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real pxyz(NXMAPMAX,NYMAPMAX,NZMAP)
      real txyz(NXMAPMAX,NYMAPMAX,NZMAP)
      integer last_part(NXMAPMAX,NYMAPMAX,NZMAP)
      integer lastpart
      common/lastparticle/ lastpart,last_part
      real*8 xminmap,yminmap,hxmap,hymap,hzmap
      real*8 zmin(NXMAPMAX,NYMAPMAX),zmax(NXMAPMAX,NYMAPMAX)
      common/densitygrid/xminmap,yminmap,zmin,zmax,
     $     hxmap,hymap,rhoxyz,uxyz,hpxyz,xhxyz,gxyz,pxyz,txyz
      integer ix,iy,iz
      real*8 realiz,zzz
      common/metals/zzz
      
      zz=metallicity

      ix=(xpos-xminmap)/hxmap+1.5d0
      iy=(ypos-yminmap)/hymap+1.5d0
      hzmap=(zmax(ix,iy)-zmin(ix,iy))/(nzmap-1)

      realiz=(zpos-zmin(ix,iy))/hzmap+1.d0
      iz=realiz
      if(iz.eq.nzmap) then
         iz=nzmap-1
      else if(iz.eq.0) then
         iz=1
      end if
c      write(*,*) xpos,ypos,zpos !,zmin(ix,iy),zmax(ix,iy),hzmap
c     write(*,*) ix,iy,iz,SIZE(rhoxyz)
c      if((zpos.eq.zmin(ix,iy)).and.(zmin(ix,iy).eq.zmax(ix,iy))) then
c         write(*,*) "ERROR!",ix,iy,zpos,zmin(ix,iy),zmax(ix,iy),hzmap
c         stop
c      end if
c      write(*,*) realiz,zpos,zmin(ix,iy),hzmap
c     Take a weighted average between the two nearest grid cells
      dens=rhoxyz(ix,iy,iz)*(iz+1-realiz)
     $     +rhoxyz(ix,iy,iz+1)*(realiz-iz)
c      write(*,*) "#### ix,iy,iz,realiz = ", ix,iy,iz,realiz
c      write(*,*) "#### rhoxyz(ix,iy,iz) = ", rhoxyz(ix,iy,iz)
c      write(*,*) "#### rhoxyz(ix,iy,iz+1) = ", rhoxyz(ix,iy,iz+1)
c      write(*,*) "#### dens = ", dens

      if(dens.le.0.d0) then
         rhocgs=0.d0
         xh=0.d0
         t6=0.d0
         xhp=0.d0
         ucgs=0.d0
         gcgs=0.d0
         pcgs=0.d0
         tcgs=0.d0
         return
      end if
c      if(dens.le.0.d0) then
c         t6=0.d0
c         return
c      end if

      xhp=hpxyz(ix,iy,iz)*(iz+1-realiz)
     $   +hpxyz(ix,iy,iz+1)*(realiz-iz)

      xh=xhxyz(ix,iy,iz)*(iz+1-realiz)
     $  +xhxyz(ix,iy,iz+1)*(realiz-iz)
         
      rhocgs=dens
      if(rhocgs.eq.0.d0) then
         t6=0.
c         write(*,*) "getLocalQuantities returned because rhocgs=0"
         return
      endif
      uu=uxyz(ix,iy,iz)*(iz+1-realiz)
     $    +uxyz(ix,iy,iz+1)*(realiz-iz)
      if(uu.lt.0d0) uu=0d0

      gg=gxyz(ix,iy,iz)*(iz+1-realiz)
     $    +gxyz(ix,iy,iz+1)*(realiz-iz)

      ppp=pxyz(ix,iy,iz)*(iz+1-realiz)
     $    +pxyz(ix,iy,iz+1)*(realiz-iz)

      tt=txyz(ix,iy,iz)*(iz+1-realiz)
     $    +txyz(ix,iy,iz+1)*(realiz-iz)

      tt=tt/dens
      ppp=ppp/dens
      gg=gg/dens
      xh=xh/dens
      uu=uu/dens
      xhp=xhp/dens              ! now xhp is the average smoothing length
                                ! at this position
      wmeanmu=0.8d0*1.67262158d-24/(xh+0.6d0-0.2d0*zz)

      ucgs=uu
      gcgs=gg
      pcgs=ppp
      tcgs=tt
c      if(iz.eq.0) then
c         write(*,*) ix,iy,iz,zpos,zmin(ix,iy),hzmap,realiz
c         write(*,*) uxyz(ix,iy,1),last_part(ix,iy,1)
c      end if
      lastpart=last_part(ix,iy,iz)
      t6=tcgs*1.e-6

c      if(t6.eq.0.d0) then
c         write(*,*) "getLocalQuantities: t6 = 0"
c         write(*,*) "dens         = ", dens
c         write(*,*) "t6, tt, tcgs = ", t6, tt, tcgs
c      end if

      if(t6.ne.t6) then
         write(*,*) "getLocalQuantities.f: t6 = ",t6
         write(*,*) "dens = ",dens
         write(*,*) "tt = ", tt
         write(*,*) "ix,iy,iz,realiz = ",ix,iy,iz,realiz
         write(*,*) "txyz(ix,iy,iz) = ", txyz(ix,iy,iz)
         write(*,*) "txyz(ix,iy,iz+1) = ",txyz(ix,iy,iz+1)
c         stop
      end if

      RETURN
      END
