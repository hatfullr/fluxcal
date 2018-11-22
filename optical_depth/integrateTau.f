      subroutine integrateTau
      include 'optical_depth.h'
      real*8 pressure,g_sph,rpos,m_sph
      real*8 avgrpos,avgxhp
      integer i,j, cnt,ii,q
      common/analysis2/ avgrpos,avgxhp
      integer counter1, counter2, counter3
      real*8 ta,tb,rprimex,rprimey,deltax,deltay,
     $     arcx,arcy

      real*8 opacit

      real*8 zstop
      real*8 TOTALTpracticalXY
      character*255 myfname


      ! Set the visible surface areas of each particle to zero
      do ip=1,nmax
         Avis(ip)=0.d0
      end do
      
      do j=1,NYMAP
         do i=1,NXMAP
            nataT(I,J) = 0.d0
            sphT(I,J) = 0.d0
            tauthin(i,j) = 0.d0
         end do
      end do

      numcell = 0
      avgt = 0.d0
      counter1=0
      counter2=0
      counter3=0
      avgrpos=0.d0
      avgxhp=0.d0
      if(dimenFileAlreadyExists) then
         write(*,*) "Integrating through each point on the driving grid"
      end if

      DO J=1,NYMAP
         DO I=1,NXMAP
            XPOS=(I-1)*HXMAP+XMINMAP ! x-coordinate of line of sight
            YPOS=(J-1)*HYMAP+YMINMAP ! y-coordinate of line of sight
            TpracticalXYthin(I,J)=0.d0
            TpracticalXYthick(I,J)=0.d0
            
            if(zmin(i,j).lt.1d30)then

               if(dointatallpos) then
                  ! Initialize the integration results array
                  do ii=1,maxstp
                     rayout1(1,ii)  = 1d30
                     rayout1(2,ii)  = 0.d0
                     rayout1(3,ii)  = 0.d0
                     rayout1(4,ii)  = 0.d0
                     rayout1(5,ii)  = 0.d0
                     rayout1(6,ii)  = 0.d0
                     rayout1(7,ii)  = 0.d0
                     rayout1(8,ii)  = 0.d0
                     rayout1(9,ii)  = 0.d0
                     rayout1(10,ii) = 0.d0
                     rayout2(ii)    = 0
                  end do
                  posx=xpos
                  posy=ypos
               end if
               
               call getTpractical(zmin(i,j),zmax(i,j),
     $              zmax_thick(i,j),thick_part(i,j),h1(i,j),
     $              TpracticalXYthin(i,j),TpracticalXYthick(i,j),
     $              tauthin(i,j))
               TOTALTpracticalXY=TpracticalXYthin(i,j)+
     $              TpracticalXYthick(i,j)

               if(dointatallpos) then
                  write(intout,'(2ES22.14,2I22)') xpos/runit_out,
     $                 ypos/runit_out,nstp,size(rayout1(:,1))+1
                  do ii=1,nstp
                     write(intout,'(10ES22.14,I22)')
     $                    (rayout1(q,ii),q=1,size(rayout1(:,1))),
     $                    rayout2(ii)
                  end do
               end if
c               if(tauthin(i,j).gt.0.d0) then
c                  write(*,*) "Optically thin material at ",i,j,xpos,ypos
c               end if
               
               if(thick_part(i,j).gt.0) then
c     This idea didn't work, but it was a good idea, so save it in case we want it in the future.
c                  ! Add to the total visible surface area for this
c                  ! particle. We want to add the "true" visible
c                  ! surface area: the area traced out by the cell on
c                  ! top of the surface of the spherical particle.
c                  ! This will be equal to the arc length traced out
c                  ! in the x direction times the arc length traced out
c                  ! in the y direction.
c                  ! arc_x = 2r*arcsin(hxmap/2r)
c                  ! arc_y = 2r*arcsin(hymap/2r)
c                  ! A = (arc_x) * (arc_y)
c                  rprimex = sqrt((2.d0*hp(thick_part(i,j)))**2.d0 -
c     $                 (YPOS - y(thick_part(i,j)))**2.d0)
c                  rprimey = sqrt((2.d0*hp(thick_part(i,j)))**2.d0 -
c     $                 (XPOS - x(thick_part(i,j)))**2.d0)
c                  if(abs(XPOS-x(thick_part(i,j))) .gt. rprimex) then ! LHS of particle
c                     deltax = rprimex - (XPOS+hxmap-x(thick_part(i,j)))
c                  else          ! RHS of particle
c                     deltax = rprimex - (XPOS - x(thick_part(i,j)))
c                  end if
c
c                  if(abs(YPOS-y(thick_part(i,j))) .gt. rprimey) then ! Bottom of particle
c                     deltay = rprimey - (YPOS+hymap-y(thick_part(i,j)))
c                  else          ! Top of particle
c                     deltay = rprimey - (YPOS - y(thick_part(i,j)))
c                  end if
c
c                  arcx = 2.d0 * rprimex * ASIN(deltax / (2.d0*rprimex))
c                  arcy = 2.d0 * rprimey * ASIN(deltay / (2.d0*rprimey))
c
c                  Avis(thick_part(i,j))=Avis(thick_part(i,j))+
c     $                 arcx*arcy
                  Avis(thick_part(i,j))=Avis(thick_part(i,j))+
     $                 hxmap*hymap
               end if
               
               if(envfit) then
                  if(thick_part(i,j).gt.0) then
                     if(I.lt.IMINGLOW) IMINGLOW=I
                     if(I.gt.IMAXGLOW) IMAXGLOW=I
                     if(J.lt.JMINGLOW) JMINGLOW=J
                     if(J.gt.JMAXGLOW) JMAXGLOW=J
                  end if
               else
                  if(kount1.gt.1) then
                     if(I.lt.IMINGLOW) IMINGLOW=I
                     if(I.gt.IMAXGLOW) IMAXGLOW=I
                     if(J.lt.JMINGLOW) JMINGLOW=J
                     if(J.gt.JMAXGLOW) JMAXGLOW=J
                  end if
               end if
               
               if(TOTALTpracticalXY.gt.0.d0) then ! There is detectable gas here
                  avgt = avgt + TOTALTpracticalXY
                  numcell = numcell + 1
c     Flux density here will be the thin stuff + the thick stuff
c                  
c     taustart(4+ifilter) has the same dimensions as the Planck function
c     Bnu=2*planck*crad/wavelengthcm(ifilter)**3d0/
c     $                 (exp(exponent)-1.d0)
c     derivs2 is using cgs units to get taustart(4+ifilter), and this is
c     erg/s/cm^2/Hz/sr  (where sr=steradian)
c                  
c     Use Planck blackbody spectrum for presumed reasonable temperature for
c     optically thin material
                  do ifilter=1,numfilters
                     exponent=coeff/(wavelength(ifilter)*t6)
                     eexponent=exp(exponent)
                     denominator=eexponent-1
                     fluxdensityXY(I,J,ifilter)=taustart(4+ifilter)+
     $                    (2*planck*crad/
     $                    wavelength(ifilter)**3d0/
     $                    denominator)
                  end do

                  ! If there is a photosphere:
                  ! If thick_part(i,j) > 0,
                  ! OR if tauthin >= tau_thick
c     if(TpracticalXYthick(I,J).gt.0.d0) then !"detectable photopshere"
                  if((thick_part(i,j) .gt. 0) .or.
     $                 (tauthin(i,j) .ge. tau_thick)) then ! "detectable photosphere"
c                     write(*,*) thick_part(i,j), tauthin, tau_thick
                     
                     if(I.lt.IMINGLOW) IMINGLOW=I
                     if(I.gt.IMAXGLOW) IMAXGLOW=I
                     if(J.lt.JMINGLOW) JMINGLOW=J
                     if(J.gt.JMAXGLOW) JMAXGLOW=J
                    
                     
                     ccphoto=ccphoto+1
                     TphotoXY(I,J) = TOTALTpracticalXY

                     ! This part of the code may not work because the main part of the code
                     ! that drove this has been removed.
                     ! nphoto is part of the legacy code
                     na=tau(3,kount1-1) ! tau(3) = # of smoothing lengths into the surface
                     nb=tau(3,kount1)
                     nphoto=(na*(tau(1,kount1)-1.d0)
     $                    +nb*(1.d0-tau(1,kount1-1)))/
     $                    (tau(1,kount1)-tau(1,kount1-1))
                     
                     do insl=1,nsl
                        if(nphoto.ge.insl) then
                           goodtphoto4avg(insl)=goodtphoto4avg(insl)
     $                          +TphotoXY(I,J)**4
                           goodccphoto(insl)=goodccphoto(insl)+1
                        endif
                     enddo
                     do insl=0,nsl
                        if(nphoto.ge.insl) then
                           do ilogwavelength=ilogwavelengthmin,
     $                          ilogwavelengthmax
                              
c     log of wavelength in micrometers:
                              logwavelength=0.01d0*ilogwavelength
c     wavelength is in centimeters
                              wavelengthcm=10.d0**logwavelength * 1d-4 
                              
                              exponent=planck*crad/
     $                             (wavelengthcm*boltz*TphotoXY(I,J))
                              spectrum(ilogwavelength,insl)=
     $                             spectrum(ilogwavelength,insl)
     $                             +wavelengthcm**(-5d0)/(exp(exponent)
     $                             -1.d0)
                           enddo
                        endif
                     enddo

                     call getLocalvz(XPOS,YPOS,sphoto,localvx,localvy,
     $                    localvz)
               
                     costheta=
     $                   -localvz/sqrt(localvx**2+localvy**2+localvz**2)
                     if(costheta.gt.costhetamax)then
                        costhetamax=costheta
                        xcosthetamax=XPOS
                        ycosthetamax=YPOS
                        zcosthetamax=sphoto
                     endif
                     
                     vzavg=vzavg+localvz*TphotoXY(I,J)**4 ! will now weight vz average based on cbrightness
                     vz2avg=vz2avg+localvz**2*TphotoXY(I,J)**4 ! will now weight vz average basecd on brightness
                     localvr=(localvx*XPOS+localvy*YPOS+localvz*sphoto)/
     $                    (XPOS**2+YPOS**2+sphoto**2)**0.5d0
                     vravg=vravg+localvr
                     vr2avg=vr2avg+localvr**2
                     
                  else
                     TphotoXY(I,J) = 0.d0
                  end if
               end if

                  
            else
               TpracticalXYthick(I,J)=0.d0
               TpracticalXYthin(I,J)=0.d0
               do ifilter=1,numfilters
                  fluxdensityXY(I,J,ifilter)=0d0
               enddo
               TphotoXY(I,J) = 0.d0
            end if
            TthermXY(I,J)=0.d0
            TXY(I,J)=TphotoXY(I,J)

         end do
      end do

      
      end subroutine
