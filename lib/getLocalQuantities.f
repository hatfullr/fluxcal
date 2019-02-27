      subroutine getLocalQuantities(myx,myy,myz)
      include '../optical_depth/optical_depth.h'

      real*8 myx,myy,myz
      real*8 realiz
      integer ix,iy

      zz=metallicity

      ix=(myx-xminmap)/hxmap+1.5d0
      iy=(myy-yminmap)/hymap+1.5d0
      hzmap=(zmax(ix,iy)-zmin(ix,iy))/(nzmap-1)

      realiz=(myz-zmin(ix,iy))/hzmap+1.d0
      iz=realiz
c      write(*,*) "hzmap = ",hzmap
      if(iz.eq.nzmap) then
         iz=nzmap-1
      else if(iz.eq.0) then
         iz=1
      end if

c      if(iy.gt.60 .and. innit.eq.5) then
c         write(*,*) ix,iy,iz, realiz
c         write(*,*) hzmap, nzmap, zmax(ix,iy),zmin(ix,iy)
c      end if
c      write(*,*) "I am here"
      if((rhoxyz(ix,iy,iz+1).le.0.d0).or.
     $     (rhoxyz(ix,iy,iz).le.0.d0)) then ! One or both cells are empty
c        Find a more precise value for boundary cases

         rhocgs = 0.d0
         ucgs = 0.d0
         xhp = 0.d0
         xh = 0.d0
         gcgs = 0.d0
         pcgs = 0.d0
         tcgs = 0.d0
         lastpart = 0
         if((rhoxyz(ix,iy,iz+1).le.0.d0).and.
     $        (rhoxyz(ix,iy,iz).le.0.d0)) then ! Both are empty
            ! There exists still the extremely rare case where part of
            ! a particle's smoothing length is just barely visible
            ! along the line of sight in between two empty grid cells.
            ! This case is so negligible that it realy isn't worth
            ! coding in. Our grid spacing is small enough to ignore cases
            ! like this.
c            write(*,*) "I finished in spot"
            return
         else if((rhoxyz(ix,iy,iz).le.0.d0).and.
     $           (rhoxyz(ix,iy,iz+1).gt.0.d0)) then !Only behind is empty
            ip = last_part(ix,iy,iz+1)
         else if((rhoxyz(ix,iy,iz+1).le.0.d0).and.
     $           (rhoxyz(ix,iy,iz).gt.0.d0)) then   !Only front is empty
            ip = last_part(ix,iy,iz)
         end if

         if(a(ip).gt.0.d0) then ! Not a sink particle
            r2 = (x(ip)-myx)**2.d0+(y(ip)-myy)**2.d0
     $           +(z(ip)-myz)**2.d0
            if(r2 .lt. 4.d0*hp(ip)**2.d0) then ! We're inside
               index = int(ctab*r2/hp(ip)**2.d0)+1
               wpc = wtab(index)/hp(ip)**3.d0
               rhocgs = am(ip)*wpc
               ucgs = am(ip)*wpc*a(ip)
               xhp = am(ip)*wpc*hp(ip)
               xhi =  -0.6d0+0.2d0*metallicity+0.8d0*1.67262158d-24/
     $              wmeanmolecular(ip)
               xh = am(ip)*wpc*xhi
               gcgs = am(ip)*wpc*localg(ip)
               pcgs = am(ip)*wpc*pp(ip)
               tcgs = am(ip)*wpc*tempp(ip)
               lastpart = ip
            end if
         end if
      else
c     Take a weighted average between the two nearest grid cells
         rhocgs=rhoxyz(ix,iy,iz)*(iz+1-realiz)
     $        +rhoxyz(ix,iy,iz+1)*(realiz-iz)
         xhp=hpxyz(ix,iy,iz)*(iz+1-realiz)
     $        +hpxyz(ix,iy,iz+1)*(realiz-iz)
         xh=xhxyz(ix,iy,iz)*(iz+1-realiz)
     $        +xhxyz(ix,iy,iz+1)*(realiz-iz)
         ucgs=uxyz(ix,iy,iz)*(iz+1-realiz)
     $        +uxyz(ix,iy,iz+1)*(realiz-iz)
         gcgs=gxyz(ix,iy,iz)*(iz+1-realiz)
     $        +gxyz(ix,iy,iz+1)*(realiz-iz)
         pcgs=pxyz(ix,iy,iz)*(iz+1-realiz)
     $        +pxyz(ix,iy,iz+1)*(realiz-iz)
         tcgs=txyz(ix,iy,iz)*(iz+1-realiz)
     $        +txyz(ix,iy,iz+1)*(realiz-iz)
         lastpart=last_part(ix,iy,iz)
      end if

      if(rhocgs.le.0.d0) then
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

      if(ucgs.lt.0d0) ucgs=0d0
      
      tcgs=tcgs/rhocgs
      pcgs=pcgs/rhocgs
      gcgs=gcgs/rhocgs
      xh=xh/rhocgs
      ucgs=ucgs/rhocgs
      xhp=xhp/rhocgs            ! now xhp is the average smoothing length
                                ! at this position
      
      wmeanmu=0.8d0*1.67262158d-24/(xh+0.6d0-0.2d0*zz)
      
      t6=tcgs*1.e-6

      if(t6.ne.t6) then
         write(*,*) "getLocalQuantities.f: t6 = ",t6
         write(*,*) "rhocgs = ",rhocgs
         write(*,*) "tcgs = ", tcgs
         write(*,*) "ix,iy,iz,realiz = ",ix,iy,iz,realiz
         write(*,*) "txyz(ix,iy,iz) = ", txyz(ix,iy,iz)
         write(*,*) "txyz(ix,iy,iz+1) = ",txyz(ix,iy,iz+1)
         write(*,*) ""
c         stop
      end if

c      write(*,*) "I finished"
      return
      end
