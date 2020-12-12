      subroutine integrateTau
      USE OMP_LIB
      include 'optical_depth.h'
      real*8 pressure,g_sph,rpos,m_sph
      real*8 avgrpos,avgxhp
      integer i,j, cnt,ii,q
      common/analysis2/ avgrpos,avgxhp
      integer counter1, counter2, counter3
      real*8 ta,tb,rprimex,rprimey,deltax,deltay,
     $     arcx,arcy

      real*8 zstop
      character*255 myfname
      real*8 TpracticalXYthin,TpracticalXYthick
      integer nan
      real start_time,finish_time
      
      call cpu_time(start_time)
      
      nan = 0
      ! Set the visible surface areas of each particle to zero
      do ip=1,nmax
         Avis(ip)=0.d0
      end do


      do j=1,NYMAP
         do i=1,NXMAP
            nataT(I,J) = 0.d0
            sphT(I,J) = 0.d0
            tauthin(i,j) = 0.d0
            TOTALTpracticalXY(i,j) = 0.d0
            zintpos(i,j) = -1.d30
         end do
      end do

      numcell = 0
      avgt = 0.d0
      counter1=0
      counter2=0
      counter3=0
      avgrpos=0.d0
      avgxhp=0.d0
c      if(dimenFileAlreadyExists) then
c         write(o,*) "Integrating through each point on the driving grid"
c     end if
      min_step_size = 0.d0
      max_step_size = 0.d0
      min_steps_taken = 2147483647
      max_steps_taken = -1

      DO J=1,NYMAP
         DO I=1,NXMAP
            XPOS=(I-1)*HXMAP+XMINMAP ! x-coordinate of line of sight
            YPOS=(J-1)*HYMAP+YMINMAP ! y-coordinate of line of sight
            TpracticalXYthin=0.d0
            TpracticalXYthick=0.d0

            if(zmin(i,j).lt.1d30)then
               nstp = 0
c               write(o,*) "zmin(i,j), zmax(i,j) = ",zmin(i,j)/runit_out,
c     $              zmax(i,j)/runit_out
               if ( debug ) then
                  write(o,*) ""
                  write(o,*) "i,j = ",i,j,"fluid"
                  write(o,*) "zmin(i,j),zmax(i,j)=",
     $                 zmin(i,j)/runit_out,zmax(i,j)/runit_out
                  write(o,*) "thick_part(i,j)=",thick_part(i,j)
                  write(o,*) "zmax_thick(i,j)=",zmax_thick(i,j),
     $                 zmax_thick(i,j)/runit_out
               end if
               call getTpractical(zmin(i,j),zmax(i,j),
     $              zmax_thick(i,j),thick_part(i,j),h1(i,j),
     $              TOTALTpracticalXY(i,j),
     $              tauthin(i,j),zintpos(i,j))
               min_steps_taken = min(min_steps_taken,nstp)
               max_steps_taken = max(max_steps_taken,nstp)
               
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
               
               if(TOTALTpracticalXY(i,j).gt.0.d0) then ! There is detectable gas here
                  avgt = avgt + TOTALTpracticalXY(i,j)
                  numcell = numcell + 1
c     Flux density here will be the thin stuff + the thick stuff
c                  
c     taustart(4+ifilter) has the same dimensions as the Planck function
c     Bnu=2*planck*crad/wavelengthcm(ifilter)**3d0/
c     $                 (exp(exponent)-1.d0)
c     derivs is using cgs units to get taustart(4+ifilter), and this is
c     erg/s/cm^2/Hz/sr  (where sr=steradian)
c                  
c     Use Planck blackbody spectrum for presumed reasonable temperature for
c     optically thin material

c     We might want this later
c                  do ifilter=1,numfilters
c                     exponent=coeff/(wavelength(ifilter)*t6)
c                     eexponent=exp(exponent)
c                     denominator=eexponent-1
c                     fluxdensityXY(I,J,ifilter)=taustart(4+ifilter)+
c     $                    (2*planck*crad/
c     $                    wavelength(ifilter)**3d0/
c     $                    denominator)
c                  end do

                  ! If there is a photosphere:
                  ! If thick_part(i,j) > 0,
                  ! OR if tauthin >= tau_thick
                     
                  if(I.lt.IMINGLOW) IMINGLOW=I
                  if(I.gt.IMAXGLOW) IMAXGLOW=I
                  if(J.lt.JMINGLOW) JMINGLOW=J
                  if(J.gt.JMAXGLOW) JMAXGLOW=J

c     We might want this later
c                     ccphoto=ccphoto+1
c                     TphotoXY(I,J) = TOTALTpracticalXY(i,j)
c
c                     ! This part of the code may not work because the main part of the code
c                     ! that drove this has been removed.
c                     ! nphoto is part of the legacy code
c                     na=tau(3,kount1-1) ! tau(3) = # of smoothing lengths into the surface
c                     nb=tau(3,kount1)
c                     nphoto=(na*(tau(1,kount1)-1.d0)
c     $                    +nb*(1.d0-tau(1,kount1-1)))/
c     $                    (tau(1,kount1)-tau(1,kount1-1))
c                     
c                     do insl=1,nsl
c                        if(nphoto.ge.insl) then
c                           goodtphoto4avg(insl)=goodtphoto4avg(insl)
c     $                          +TphotoXY(I,J)**4
c                           goodccphoto(insl)=goodccphoto(insl)+1
c                        endif
c                     enddo
c                     do insl=0,nsl
c                        if(nphoto.ge.insl) then
c                           do ilogwavelength=ilogwavelengthmin,
c     $                          ilogwavelengthmax
c                              
cc     log of wavelength in micrometers:
c                              logwavelength=0.01d0*ilogwavelength
cc     wavelength is in centimeters
c                              wavelengthcm=10.d0**logwavelength * 1d-4 
c                              
c                              exponent=planck*crad/
c     $                             (wavelengthcm*boltz*TphotoXY(I,J))
c                              spectrum(ilogwavelength,insl)=
c     $                             spectrum(ilogwavelength,insl)
c     $                             +wavelengthcm**(-5d0)/(exp(exponent)
c     $                             -1.d0)
c                           enddo
c                        endif
c                     enddo
c
c                     call getLocalvz(XPOS,YPOS,sphoto,localvx,localvy,
c     $                    localvz)
c               
c                     costheta=
c     $                   -localvz/sqrt(localvx**2+localvy**2+localvz**2)
c                     if(costheta.gt.costhetamax)then
c                        costhetamax=costheta
c                        xcosthetamax=XPOS
c                        ycosthetamax=YPOS
c                        zcosthetamax=sphoto
c                     endif
c                     
c                     vzavg=vzavg+localvz*TphotoXY(I,J)**4 ! will now weight vz average based on cbrightness
c                     vz2avg=vz2avg+localvz**2*TphotoXY(I,J)**4 ! will now weight vz average basecd on brightness
c                     localvr=(localvx*XPOS+localvy*YPOS+localvz*sphoto)/
c     $                    (XPOS**2+YPOS**2+sphoto**2)**0.5d0
c                     vravg=vravg+localvr
c                     vr2avg=vr2avg+localvr**2
c                     
c                  else
c                     TphotoXY(I,J) = 0.d0
               end if
            else
               if ( debug ) then
                  write(o,*) "i,j=",i,j,"no fluid"
               end if
               TOTALTpracticalXY(i,j)=0.d0
c     We might want this later
c               do ifilter=1,numfilters
c                  fluxdensityXY(I,J,ifilter)=0d0
c               enddo
c               TphotoXY(I,J) = 0.d0
            end if
c     We might want this later
c            TthermXY(I,J)=0.d0
c            TXY(I,J)=TphotoXY(I,J)

         end do
      end do
      
      call cpu_time(finish_time)
c      write(o,*) "TOTALTpracticalXY(i/2,j/2) = ",i/2,j/2,
c     $     TOTALTpracticalXY(i/2,j/2)
c      write(o,*) "integrateTau took ",finish_time-start_time

      end subroutine
