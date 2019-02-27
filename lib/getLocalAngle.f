      function getLocalAngle(xin,yin,zin)
c     This function calculates the local angle cos(psi) given a
c     coordinate within the simulated fluid. It considers a local area
c     element dA, which is oriented in exactly the same direction as an
c     element at the surface dS. Both dA and dS are along the same ray
c     that travels from the surface of the simulation to the origin. dS
c     is located on the surface of the farthest kernel along r and dA is
c     located somewhere closer in (given by user input).
c
c     Input:
c        xin = x coordinate
c        yin = y coordinate   
c        zin = z coordinate
c
c     Output:
c        cos(psi) = the local angle to the observer from dA

      include 'flux_cal.h'
      real*8 xin,yin,zin
      real*8 r,theta,phi,ri,thetai,phii,dtheta,dphi,rn
      real*8 myx,myy,myz
      integer i,maxid
      real*8 myzmax,maxdz
      real*8 maxdz_temp
      real*8 alpha,beta
      
c     Let's find the components of r
      r = sqrt(xin**2.d0 + yin**2.d0 + zin**2.d0)
      theta = acos(zin/r)
      phi = atan2(yin,xin)
      
c     Let's pretend we're at a different viewing angle looking down r
c     Need to rotate ourselves to (theta,phi)
c     The following rotation is done the same way as in SPLASH
c     (and setViewingAngle.f)
c     open(51,file="angle_testing.dat")
      myzmax = -1.d30
      do i=1,n
         ri = sqrt(x(i)**2.d0 + y(i)**2.d0 + z(i)**2.d0)
         thetai = acos(z(i)/ri)
         phii = atan2(y(i),x(i))
         dtheta = theta-thetai
         dphi = phi-phii
         myx = ri*sin(dtheta)*cos(dphi)
         myy = ri*sin(dtheta)*sin(dphi)
         myz = ri*cos(dtheta)
         maxdz = sqrt((4.d0*hp(i)**2.d0
     $        -myx**2.d0 - myy**2.d0))
         maxdz_temp = myz+maxdz
         if(maxdz_temp .gt. myzmax) then
            myzmax = maxdz_temp
            ! Find the surface normal
c            rn = sqrt(myx**2.d0+myy**2.d0+(myz-myzmax)**2.d0)
            maxid = i
         end if
c         write(51,*) x(i),y(i),z(i),myx(i),myy(i),myz(i)
      end do
c      close(51)

      ri = sqrt(x(maxid)**2.d0 + y(maxid)**2.d0 + z(maxid)**2.d0)
      thetai = acos(z(maxid)/ri)
      dtheta = theta-thetai
      myz = ri*cos(dtheta)

c      ! The surface normal
c      rn = sqrt(myx**2.d0+myy**2.d0+(myz-myzmax)**2.d0)

      ! The angle our viewing angle makes with the surface normal
      alpha = acos((myzmax-myz)/(2.d0*hp(maxid)))

      ! Orient to the surface normal
      
      
      write(*,*) alpha,theta,cos(alpha+theta)
! The angle our viewing angle makes with the observer
c      write(*,*) xin,yin,zin,cos(alpha+theta),thetai,theta
      getLocalAngle = cos(alpha+theta)
      
      
c     The particle farthest away from the center has been found (maxid),
c     as well as the position on its surface that is the farthest along
c     this rotated viewing angle. We also found the surface normal rn.


      end function
