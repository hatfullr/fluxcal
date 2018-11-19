      recursive subroutine quicksort(a, first, last)
c     quicksort.f -*-f90-*-
c     Author: t-nissie
c     License: GPLv3
c     Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
      implicit none
      real*8  a(*), x, t
      integer first, last
      integer i, j
      
      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         i=i+1
         j=j-1
      end do
      if (first < i-1) call quicksort(a, first, i-1)
      if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort
      
      
      recursive subroutine quicksort2(a, b, first, last)
c     Just a modified version of the original.
c     Sort arrays a and b using information from array a
c     For example, if a is distance on the z-axis and b is the
c     particle IDs, this will sort a while also sorting b.
      implicit none
      real*8  a(*), x, t,k
      integer b(*)
      integer first, last
      integer i, j
      
      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         k = b(i);  b(i) = b(j);  b(j) = k
         i=i+1
         j=j-1
      end do
      if (first < i-1) call quicksort2(a,b, first, i-1)
      if (j+1 < last)  call quicksort2(a,b, j+1, last)
      end subroutine quicksort2
      
