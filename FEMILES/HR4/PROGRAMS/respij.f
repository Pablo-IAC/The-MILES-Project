	subroutine respij(xl,s)
c Subroutine to calculate the response of the I Johnson filter
c It calls the "hunt" and "polint" routines of the numerical recipes 
       	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (nf=27)
	dimension XV(nf),SV(nf)
c Johnson, 1965, ApJ, 141,923 (indicated in Alonso+95,297,197
        data XV/6800.,7000.,7200.,7400.,7600.,7800.,8000.,8200.,8400.
     &,8600.,8800.,9000.,9200.,9400.,9600.,9800.,10000.,10200.,10400.
     &,10600.,10800.,11000.,11200.,11400.,11600.,11800.,12000./
	data SV/0.,0.01,0.17,0.36,0.56,0.76,0.96,0.98,0.99,1.00,0.98,
     &0.93,0.84,0.71,0.58,0.47,0.36,0.28,0.20,0.15,0.10,0.08,0.05,
     &0.03,0.02,0.,0./
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
	end
