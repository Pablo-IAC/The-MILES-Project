c	---------------------------------------------------------------
c   subrutina que calcula la inversa de una matriz
      SUBROUTINE matinx(a)
c
c  This is a fast explicit routine for calculating the inverse of
c  a 4x4 real matrix. It is generalized i.e. it does not take 
c  account of the symmetry/anti-symmetry of the absorption matrix.
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a(4,4),b(4,4)
      INTEGER i, j
c      double precision a(4,4), b(4,4)
c      double precision absmax, fabsmx, det, fdeta
      absmax = 0.d0
c  Find the largest value in the matrix - assumed to be a
c  suitable factorization factor for the whole matrix.
      do 5010 i = 1, 4
        do 5009 j = 1, 4
	  if ( abs( a(i,j) ) .gt. absmax ) absmax = a(i,j)
 5009   continue
 5010 continue
c
c  If matrix is zero, stop. No other information is available.
c
      if ( absmax .eq. 0.d0 ) stop 'Matrix is zero.'
c
c  One division for best speed.
c
      fabsmx = 1.d00 / absmax
c
c  and then factor out absmax.
c
      do 5020 i = 1, 4
        do 5019 j = 1, 4
	  a(i,j) = a(i,j) * fabsmx
 5019   continue
 5020 continue
c
c  Calculate the 16 elements of the inverse. Determinant is excluded.
c
      b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)
     1 + a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)
     2 - a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
      b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)
     1 + a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)
     2 - a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
      b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)
     1 + a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)
     2 - a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
      b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)
     1 + a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)
     2 - a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
      b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)
     1 + a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)
     2 - a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
      b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)
     1 + a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)
     2 - a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
      b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)
     1 + a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)
     2 - a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
      b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)
     1 + a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)
     2 - a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
      b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)
     1 + a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)
     2 - a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
      b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)
     1 + a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)
     2 - a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
      b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)
     1 + a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)
     2 - a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
      b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)
     1 + a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)
     2 - a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
      b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)
     1 + a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)
     2 - a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
      b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)
     1 + a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)
     2 - a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
      b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)
     1 + a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)
     2 - a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
      b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)
     1 + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)
     2 - a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)
c
      det = a(1,1) * b(1,1) + a(1,2) * b(2,1)
     1 + a (1,3) * b(3,1) + a(1,4) * b(4,1)
c
      fdeta = fabsmx / det
c
c  Finally, prepare matrix inverse and place in original matrix.
c
      do 5001 i = 1, 4
        do 5000 j = 1, 4
          a(i,j) = b(i,j) * fdeta
 5000   continue
 5001 continue
      return
      end
