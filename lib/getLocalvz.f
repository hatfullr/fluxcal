      subroutine getLocalvz(xpos,ypos,zpos,localvx,localvy,localvz)
      include 'flux_cal.h'
      real*8 R2, R2MAX, WPC, XPOS, YPOS, ZPOS, dens
      integer INDEX, i
      real*8 localvx,localvy,localvz

      dens=0.d0
      localvx=0.d0
      localvy=0.d0
      localvz=0.d0
      DO I=1,N
         if(a(i).gt.0.d0)then
            R2=(Xpos-X(I))**2+(Ypos-Y(I))**2+(Zpos-Z(I))**2
            R2MAX=4.d0*HP(I)**2
            IF (R2.LT.R2MAX) THEN
               INDEX=INT(CTAB*R2/HP(I)**2)+1
               WPC=WTAB(INDEX)/HP(I)**3
               dens=dens+AM(I)*WPC
               localvx=localvx+AM(I)*WPC*vx(i)
               localvy=localvy+AM(I)*WPC*vy(i)
               localvz=localvz+AM(I)*WPC*vz(i)
            ENDIF
c     print *,'dens',i,dens
         endif
      ENDDO

      if(dens.ne.0.d0) then
         localvx=localvx/dens
         localvy=localvy/dens
         localvz=localvz/dens
      else
c         print *,'density=0 at',xpos,ypos,zpos
c         stop
      endif


      RETURN
      END
