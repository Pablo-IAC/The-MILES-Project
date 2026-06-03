c*******SUBRUTINA respv************************************************
	subroutine respk(xl,s)
c Subroutine for calculating the response curve of the K filter
c It calls the "hunt" and "polint" routines of the numerical recipes 
c
c WARNING: the provided wavelength should be double-precision
c
       	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	PARAMETER (nf=23)
	PARAMETER (nf=17)
	dimension XV(nf),SV(nf)
c        data XV/19400.,19600.,19800.,20000.,20200.,20800.,21000.,
c     &21200.,21400.,21600.,21800.,22000.,22200.,22400.,22600.,
c     &22800., 23000., 23200., 23400., 23800., 24000., 24400.,
c     &24800./
c	data SV/.0,.119999997,.200000003,.300000012,.550000012,
c     &.769999981,.850000024,.899999976,.939999998,.939999998, 
c     &.949999988,.939999998,.959999979,.980000019,.970000029, 
c     &.959999979,.910000026,.879999995,.839999974,.819999993,
c     &.639999986,.100000001,.0/
c       xli=19400.0d0
c	xlf=24800.0d0
        data XV/18000.,18500.,19000.,19500.,20000.,20500.,21000.,
     &21500.,22000.,22500.,23000.,23500.,24000.,24500.,25000.,
     &25500.,26000./
	data SV/0.0,0.100000001,0.479999989,0.949999988,1.,
     &0.980000019,0.959999979,0.949999988,0.970000029,0.959999979, 
     &0.939999998,0.949999988,0.949999988,0.839999974,0.460000008, 
     &0.0799999982,0.0/
 	call hunt(XV,nf,xl,k1)
        xli=18500.0d0
	xlf=25500.0d0
	if(xl.gt.xli.and.xl.lt.xlf)then
	   call polint(XV(k1-1),SV(k1-1),3,xl,s,dy)
	   if(s.gt.1.0d0) s=1.0d0
	else
	   s=0.0d0
	endif
	end
  
