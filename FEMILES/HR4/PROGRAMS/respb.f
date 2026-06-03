
      SUBROUTINE respb(xl,s)
c Subroutine for calculating the response curve of the B filter
c It calls the "hunt" and "polint" routines of the numerical recipes 
c Filter in use: Buser & Kurucz 1978 (A&A,70,555) (Filter B3)
c WARNING: the provided wavelength should be double-precision
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nf=42)
      DIMENSION XV(nf),SV(nf)
      data XV/3550.,3600.,3650.,3700.,3750.,3800.,3850.,3900.,3950.,
     &4000.,4050.,4100.,4150.,4200.,4250.,4300.,4350.,4400.,4450.,
     &4500.,4550.,4600.,4650.,4700.,4750.,4800.,4850.,4900.,4950.,
     &5000.,5050.,5100.,5150.,5200.,5250.,5300.,5350.,5400.,5450.,
     &5500.,5550.,5600./
c B3 Buser
      data SV/.0,.0,.006,.030,.060,.134,.302,.567,.841,.959,.983,.996,
     &1.00,.996,.987,.974,.957,.931,.897,.849,.800,.748,.698,.648,.597,
     &.545,.497,.447,.397,.345,.297,.252,.207,.166,.129,.095,.069,.043,
     &.024,.009,.0,.0/
c B2 Buser
c	data SV/.0,.0,.006,.023,.045,.106,.254,.492,.752,.881,.923,.955,
c     &.977,.990,1.0,1.0,.997,.984,.958,.916,.871,.820,.775,.723,.672,.617,
c     &.569,.511,.457,.402,.347,.299,.244,.199,.154,.113,.084,.051,.029,.010,
c     &.0,.0/
c       xli=XV(1)
c	xlf=XV(nf)
        xli=3600.0d0
	xlf=5550.0d0
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
