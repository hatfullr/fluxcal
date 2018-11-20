      subroutine prepareIntegration
c     L346 to L468
      include 'optical_depth.h'

      real*8 zmmax
      real*8 opacit

      integer count, i

      eps=1.0d-6                ! parameter for desired accuracy in
                                ! integration
      kmax=200
      IMINGLOW=NXMAP+1
      JMINGLOW=NYMAP+1
      IMAXGLOW=0
      JMAXGLOW=0
      rmid=0.d0
      r2mid=0.d0
      costhetamax=-1d30
      vzavg=0.d0
      vz2avg=0.d0
      vravg=0.d0
      vr2avg=0.d0
      dssav=1d30
      nphotoavg=0.d0
      warning=0

      do insl=1,nsl
         goodtphoto4avg(insl)=0
         goodccphoto(insl)=0
      enddo

      ccphoto=0

      if(dimenFileAlreadyExists) then
         write(*,*) "Creating the driving grid"
      end if

c     Find limits of integration:
      count=0
      DO J=1,NYMAP
         DO I=1,NXMAP
            count = count+1
            zmin(i,j)=1d30
            zmax(i,j)=-1d30
            zmax_thick(i,j)=-1d30
            thick_part(i,j)=0
         enddo
      enddo

      count=0
      do IP=1,N
         if(a(IP).gt.0.d0) then
            yposmin=y(ip)-2*hp(ip)
            yposmax=y(ip)+2*hp(ip)
            jmin=max(int((yposmin-yminmap)/hymap+2),1)
            jmax=min(int((yposmax-yminmap)/hymap+1),nymap)
            DO J=jmin,jmax
               YPOS=(J-1)*HYMAP+YMINMAP ! y of line of sight
               maxdx=(4*hp(ip)**2-(y(ip)-ypos)**2)**0.5d0
               xposmin=x(ip)-maxdx
               xposmax=x(ip)+maxdx
               imin=max(int((xposmin-xminmap)/hxmap+2),1)
               imax=min(int((xposmax-xminmap)/hxmap+1),nxmap)
               DO I=imin,imax
                  XPOS=(I-1)*HXMAP+XMINMAP ! x of line of sight
                  maxdz=(4*hp(ip)**2
     $                 -(x(ip)-xpos)**2-(y(ip)-ypos)**2)**0.5d0
                  ! Protect against infinitessimally small integration
                  if(maxdz.gt.0) then
                     if(Z(IP)-maxdz.lt.zmin(i,j)) then
                        zmin(i,j)=Z(IP)-maxdz
                        h1(i,j)=hp(ip)
                     endif
c                     if(Z(IP)+maxdz.gt.zmax(i,j)) then
c                        count=count+1
c                        closest(count) = IP
c                     end if
                     zmax(i,j)=max(zmax(i,j),Z(IP)+maxdz)
                     if(tauA(ip) .gt. tau_thick) then
                        if((Z(IP)+maxdz).gt.zmax_thick(i,j)) then
                           zmax_thick(i,j) = Z(IP)+maxdz
                           thick_part(i,j) = ip
                        end if
                     end if
                  endif
               enddo
            enddo
         endif
      enddo

      if(dimenFileAlreadyExists) then
         write(*,*) "Creating the integrating grid"
      end if
c     Find 3d density grid for quick look-up later
      DO J=1,NYMAP
         DO I=1,NXMAP
            do k=1,nzmap
               rhoxyz(i,j,k)=0. ! density
               uxyz(i,j,k)=0.   ! specific internal energy
               hpxyz(i,j,k)=0.  ! smoothing length
               xhxyz(i,j,k)=0.  ! hydrogen abundance
               gxyz(i,j,k)=0.   ! local gravitational acceleration
               pxyz(i,j,k)=0.   ! pressure
               txyz(i,j,k)=0.   ! temperature
               last_part(i,j,k)=0 ! Last particle touched
            enddo
         enddo
      enddo

      do IP=1,N
         if(a(IP).gt.0.d0) then
            yposmin=y(ip)-2*hp(ip)
            yposmax=y(ip)+2*hp(ip)
            jmin=max(int((yposmin-yminmap)/hymap+2),1)
            jmax=min(int((yposmax-yminmap)/hymap+1),nymap)
            
            DO J=jmin,jmax
               YPOS=(J-1)*HYMAP+YMINMAP ! y of line of sight
               maxdx=(4*hp(ip)**2-(y(ip)-ypos)**2)**0.5d0
               xposmin=x(ip)-maxdx
               xposmax=x(ip)+maxdx
               imin=max(int((xposmin-xminmap)/hxmap+2),1)
               imax=min(int((xposmax-xminmap)/hxmap+1),nxmap)
               DO I=imin,imax
                  XPOS=(I-1)*HXMAP+XMINMAP ! x of line of sight
                  maxdz=(4*hp(ip)**2
     $                 -(x(ip)-xpos)**2-(y(ip)-ypos)**2)**0.5d0
                  zposmin=z(ip)-maxdz
                  zposmax=z(ip)+maxdz
                  hzmap=(zmax(i,j)-zmin(i,j))/(nzmap-1)
                  izmin=max(int((zposmin-zmin(i,j))/hzmap+2),1)
                  izmax=min(int((zposmax-zmin(i,j))/hzmap+1),nzmap)
                  zmmax=-1.d30
                  do iz=izmin,izmax
                     zPOS=(iz-1)*HzMAP+zMIN(i,j) ! z along line of sight
                     R2=(X(ip)-XPos)**2+(Y(ip)-YPos)**2+(Z(ip)-ZPos)**2
                     index=INT(CTAB*R2/HP(Ip)**2)+1
                     WPC=WTAB(index)/HP(IP)**3
                     rhoxyz(i,j,iz)=rhoxyz(i,j,iz)+am(ip)*wpc
                     uxyz(i,j,iz)=uxyz(i,j,iz)+am(ip)*wpc*a(ip)
                     hpxyz(i,j,iz)=hpxyz(i,j,iz)+am(ip)*wpc*hp(ip)
                     xhi=-0.6d0+0.2d0*metallicity
     $                    +0.8d0*1.67262158d-24/wmeanmolecular(Ip)
                     xhxyz(i,j,iz)=xhxyz(i,j,iz)+am(ip)*wpc*xhi
                     gxyz(i,j,iz)=gxyz(i,j,iz)+am(ip)*wpc*
     $                    localg(ip)
                     pxyz(i,j,iz)=pxyz(i,j,iz)+am(ip)*wpc*pp(ip)
                     txyz(i,j,iz)=txyz(i,j,iz)+am(ip)*wpc*tempp(ip)

                     if(zpos.gt.zmmax) then
                        last_part(i,j,iz)=ip
                        zmmax=zpos
                     end if
                  enddo
                  
               enddo
            enddo
         endif
      enddo

      end subroutine
