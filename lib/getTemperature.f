      subroutine getTemperature(q,r,x3)

      implicit none
c      include 'sphu.h'

      double precision q,r,k,B,piece1,piece2
      double precision y1,y2,yy,AA,b2,c2,x3,kh
c      double precision pgas,prad,beta1,gam

C     subroutine to solve 4th order equations to determine the temperature x3
C     for an equation of state with both ideal gas and radiation pressure.
c     Written by Scott Fleming 10/04/02 and James Lombardi 2002-2003

C     the fourth order equation comes from u_gas+ u_rad = u, with
C     u_gas proportional to T and u_rad proportional to T^4

C     In general, we can transform a 4th order equation to x^4+px^2+qx+r=0
C     (see pages 57-58 of Stillwell's "Mathematics and its history" text)
C     but we fortunately don't even have an x^2 term (that is, p=0).

C     Follow Stillwell, we can transform this into a cubic equation:
C     First solve for y by using the fact that B^2-4AC=0
C     equation is then:  y^3=ry+q^2/8
C     using the solution of cubic equations found in Stillwell page 55:

      k = 0.125d0*q**2
      kh=0.5d0*k
      if(kh**2-(r/3.d0)**3.le.0.d0)then
         print *, k,r,kh**2-(r/3.d0)**3
         stop 'bad input: imaginary results?'
      endif

      piece1 = kh+(kh**2-(r/3.d0)**3)**0.5d0
      piece2 = (r/3.d0)**3.d0/piece1
c      piece2 = kh-(kh**2-(r/3.d0)**3)**0.5d0

      y1 = piece1**(1.d0/3.d0)
      
c     Fortran can't handle cube roots of neg. #'s
      y2 = -dabs(piece2)**(1.d0/3.d0)

      yy = y1+y2
      

C     equation to solve: (x^2+p+y)^2=Ax^2+Bx+C
C     Now take square root of both sides with:

      AA = 2.d0*yy
      B = -q
C      C = -r+y**2

C     re-writing Ax^2+Bx+C as a square then solving the equation we
C     obtain 2 results:
C     x^2 + (-(A^(1/2)))x + (-B/(2(A)^(1/2))+p+y) = 0 (1)
C     or
C     x^2 + (A^(1/2))x + (B/(2(A)^(1/2))+p+y) = 0     (2)

C     Our solution we're interested in:

      b2 = AA**0.5d0
      c2 = 0.5d0*B/b2 + yy

C     therefore, we once again have x^2+bx+c=0, and our answer we want is:

c     x3 = 0.5d0*(-b2 + (b2**2-4.d0*c2)**0.5d0)
      x3 = -2.d0*c2/(b2 + (b2**2-4.d0*c2)**0.5d0)

      if(piece1.lt.0.d0)then
         print *, 'piece 1 lt 0',
     $        k,r,piece1,piece2
      endif
      if(b2.eq.0.d0)then
         x3=(-r -q*( -r-q*(-r)**0.25d0 )**0.25d0)**0.25d0
      endif
      if(piece2.ge.0.d0) x3=-(r+(r/q)**4)/q
c      if(piece1.eq.-piece2)then
c         print *, 'piece 1 eq -piece 2 (rad pressure dominates)',
c     $        k,r,piece1,piece2,b2,c2,x3,(-r)**0.25d0
c         x3=(-r -q*( -r-q*(-r)**0.25d0 )**0.25d0)**0.25d0
c         print *, x3
c      endif
c      if(piece2.ge.0.d0)then
c         print *, 'piece 2 ge 0 (gas pressure dominates)',
c     $        k,q,r,piece1,piece2,b2,c2,x3,-r/q

c         x3=-(r+(r/q)**4)/q

c      endif
      
c      if(piece2.gt.0.d0) then
c         pgas=q/qconst*boltz*x3
c         prad=arad*x3**4/3.d0
c         beta1=pgas/(pgas+prad)
c         gam=(32.d0-24.d0*beta1-3.d0*beta1**2) /
c     $        (24.d0-21.d0*beta1)
c         print *, x3, pgas,prad,beta1,gam
c      endif

      end
      
