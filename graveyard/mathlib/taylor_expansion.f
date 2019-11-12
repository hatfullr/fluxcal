      real*8 function taylor_expansion(x,xarr,farr,N,a,maxorder,
     $     leftward)
c     This function calculates the value f(x) at location x by expanding
c     the function defined by farr and xarr around location a.
c     Input: 
c        'x' is the array index to evaluate f(x) at upon expansion
c        'xarr' is an array of size 'maxorder' that holds values of x
c        'farr' is an array of size 'maxorder' that holds values of f(x)
c        'N' is the size of xarr and farr (must be the same size)
c        'a' (integer) is the array index to expand f(x) around
c        'maxorder'-1 is the order to take the expansion out to
c        'leftward' = .false. takes the derivates moving rightward along x
c     Output:
c        f(x) as the taylor expansion about a

      implicit none
      integer order,maxorder,i,j
      integer N,a, x1,x2
      real*8 xarr(N), farr(N)
      real*8 derivs(maxorder-1,maxorder-1)
      real*8 x,result, factorial, derivative
      logical left
      logical,optional :: leftward

      if(present(leftward)) then
         left = leftward
      else
         left = .false.
      end if
      
      if ( ((left.eqv..false.) .and. (N-a .lt. maxorder)) .or.
     $     ((left.eqv..true.) .and. (a .lt. maxorder)) ) then
         write(*,*) "a = ",a
         write(*,*) "N = ",N
         write(*,*) "left = ",left
         write(*,*) "maxorder = ",maxorder
         write(*,*) "Cannot compute the finite difference taylor "//
     $        "expansion to the order specified: not enough data in "//
     $        "the input arrays, or the location 'a' at which to "//
     $        "expand the series is too close to the edge of the "//
     $        "input arrays."
         error stop "taylor_expansion.f"
      end if
      
      ! Initialize the derivs array
      do i=1,maxorder-1
         do j=1,maxorder-1
            derivs(i,j) = 0.d0
         end do
      end do
      
      if (left) then            ! Take derivatives moving leftward
         do i=1,maxorder-1
            x1 = a - i
            x2 = a - (i-1)
c            derivs(1,i)=(farr(x2)-farr(x1))/(xarr(x2)-xarr(x1))
            derivs = derivative(xarr(x1),xarr(x2),farr(x1),farr(x2))
c            write(*,*) "derivs(1,i) = ",derivs(1,i)
         end do
         do j=2,maxorder-1
            do i=1,maxorder-1-j   ! Find higher order derivatives
               x1 = a - i
               x2 = a - (i-1)
c               derivs(j,i)=(derivs(j-1,i)-derivs(j-1,i+1))/
c     $              (xarr(x2)-xarr(x1))
               derivs(j,i) = derivative(xarr(x1),xarr(x2),
     $              derivs(j-1,i+1),derivs(j-1,i))
c               write(*,*) "j,i,derivs(j,i) = ",j,i,derivs(j,i)
            end do
         end do
      else                      ! Take derivatives moving rightward
         do i=1, maxorder-1
            x1 = a + (i-1)
            x2 = a + i
c            derivs(1,i)=(farr(x2)-farr(x1))/(xarr(x2)-xarr(x1))
            derivs(1,i)=derivative(xarr(x1),xarr(x2),farr(x1),farr(x2))
c            write(*,*) "derivs(1,i) = ",derivs(1,i)
         end do
         do j=2,maxorder-1
            do i=1,maxorder-1-j   ! Find higher order derivatives
               x1 = a + (i-1)
               x2 = a + i
c               derivs(j,i)=(derivs(j-1,i+1)-derivs(j-1,i))/
c     $              (xarr(x2)-xarr(x1))
               derivs(j,i) = derivative(xarr(x1),xarr(x2),
     $              derivs(j-1,i),derivs(j-1,i+1))
c               write(*,*) "j,i,derivs(j,i) = ",j,i,derivs(j,i)
            end do
         end do
      end if

      result = farr(a)          ! Zero'th order
      
      do order=1,maxorder-1

         ! Calculate the factorial
         factorial = 1.d0
         do j=1,order
            factorial = factorial * j
         end do

         ! Compute the Taylor expansion
         result = result + (derivs(order,a)*(x-xarr(a))**order /
     $        factorial)
c         write(*,*) "result = ",result
      end do
      
      taylor_expansion = result
      return
      end


      real*8 function derivative(x1,x2,f1,f2)
      implicit none
      real*8 x1,x2,f1,f2

      derivative = (f2-f1)/(x2-x1)
      
      return
      end function
