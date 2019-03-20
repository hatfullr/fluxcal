      subroutine setViewingAngle
      include 'optical_depth.h'
c     L166 to L229
c     Set viewing angle in radians:
      anglex=anglexdeg/180*pi
      angley=angleydeg/180*pi
      anglez=anglezdeg/180*pi

c     Rather than use different viewing angles, it's easier just
c     to rotate the particles in the system and continue to look
c     along the z-axis.  The effect is the same.

c     The following rotation is done the same way as in SPLASH
      do ip=1,N
c     Rotate position vector of particle i:
         if (anglez.ne.0d0) then ! rotate about z
            rold = sqrt(x(ip)**2 + y(ip)**2)
            phi = ATAN2(y(ip),x(ip))
            phi = phi - anglez
            x(ip) = rold*COS(phi)
            y(ip) = rold*SIN(phi)
         endif
         if (angley.ne.0d0) then ! rotate about y
            rold = sqrt(z(ip)**2 + x(ip)**2)
            phi = ATAN2(z(ip),x(ip))
            phi = phi - angley
            z(ip) = rold*SIN(phi)
            x(ip) = rold*COS(phi)
         endif
         if (anglex.ne.0d0) then ! rotate about x
            rold = sqrt(y(ip)**2 + z(ip)**2)
            phi = ATAN2(z(ip),y(ip))
            phi = phi - anglex
            y(ip) = rold*COS(phi)
            z(ip) = rold*SIN(phi)
         endif

c     Rotate velocity vector of particle i:
         if(anglez.ne.0) then   ! rotate about z
            rold = sqrt(vx(ip)**2 + vy(ip)**2)
            phi = ATAN2(vy(ip),vx(ip))
            phi = phi - anglez
            vx(ip) = rold*COS(phi)
            vy(ip) = rold*SIN(phi)
         endif
         if (angley.ne.0d0) then ! rotate about y
            rold = sqrt(vz(ip)**2 + vx(ip)**2)
            phi = ATAN2(vz(ip),vx(ip))
            phi = phi - angley
            vz(ip) = rold*SIN(phi)
            vx(ip) = rold*COS(phi)
         endif
         if (anglex.ne.0d0) then ! rotate about x
            rold = sqrt(vy(ip)**2 + vz(ip)**2)
            phi = ATAN2(vz(ip),vy(ip))
            phi = phi - anglex
            vy(ip) = rold*COS(phi)
            vz(ip) = rold*SIN(phi)
         endif

      enddo

      end subroutine
