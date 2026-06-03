      SUBROUTINE respu(xl,s)
c Subroutine to calculate the response of the U Johnson filter
c It calls the "hunt" and "polint" routines of the numerical recipes 
c WARNING: the provided wavelength should be double-precision
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nf=24)
      DIMENSION XV(nf),SV(nf)
c Buser 1978
      data XV/3050.,3100.,3150.,3200.,3250.,3300.,3350.,3400.,3450.,
     &3500.,3550.,3600.,3650.,3700.,3750.,3800.,3850.,3900.,3950.,
     &4000.,4050.,4100.,4150.,4200./
      data SV/0.,0.0199999996,0.0769999996,0.135000005,0.203999996,
     &0.282000005,0.38499999,0.493000001,0.600000024,0.704999983,
     &0.819999993,0.899999976,0.958999991,0.992999971,1.000000000,
     &0.975000024,0.850000024,0.644999981,0.400000006,0.223000005,
     &0.125,0.057,0.00499999989,0./ 
        xli=XV(1)
	xlf=XV(nf)
 	call hunt(XV,nf,xl,k1)
	if(xl.gt.xli.and.xl.lt.xlf)then
	   call polint(XV(k1-1),SV(k1-1),3,xl,s,dy)
	   if(s.gt.1.0d0) s=1.0d0
	   if(s.lt.0.0d0) s=0.0d0
	else
	   s=0.0d0
	endif
	RETURN
	END
