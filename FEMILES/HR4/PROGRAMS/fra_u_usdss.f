c Subrutina que calcula fraccion: u SDSS / U_Johnson
c Se pasan tiss,giss y devuelve fraccion
c Se uso Pickles library 
c subrutina respusdss.f, respu.f 
c Llamada desde hrsl.f
      SUBROUTINE fra_u_usdss(tiss,giss,fUm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      tissl=dlog10(tiss)
	tissl=tiss
cGIANTS1_MILES
      ug00m=-13.8845
      ug10m=7.97029
      ug20m=-1.07311
cGIANTS2_MILES
      ug01m=-14.8468
      ug11m=7.15258
      ug21m=-0.801627
cDWARFS1_MILES
      ud00m=-22.1331
      ud10m=12.3401
      ud20m=-1.64796
cDWARFS2_MILES
      ud01m=-10.8327
      ud11m=5.20346
      ud21m=-0.566797
ccccccccccccccccccccccccccc    
      tlimd=3.96d0
      tlimg=3.915d0
      tMINI=3.40d0
      tMAXI=4.50d0
c     
      IF(giss.lt.3.5d0)THEN
        if(tissl.lt.tMINI)then
         fUm=ug00m+ug10m*tMINI+ug20m*(tMINI*tMINI)
	elseif(tissl.gt.tMINI.and.tissl.lt.tlimg)then
         fUm=ug00m+ug10m*tissl+ug20m*(tissl*tissl)
	elseif(tissl.gt.tMAXI)then
         fUm=ug01m+ug11m*tMAXI+ug21m*(tMAXI*tMAXI)
	else
         fUm=ug01m+ug11m*tissl+ug21m*(tissl*tissl)
        endif
      ELSE
        if(tissl.lt.tMINI)then
         fUm=ud00m+ud10m*tMINI+ud20m*(tMINI*tMINI)
	elseif(tissl.gt.tMINI.and.tissl.lt.tlimd)then
         fUm=ud00m+ud10m*tissl+ud20m*(tissl*tissl)
	elseif(tissl.gt.tMAXI)then
         fUm=ud01m+ud11m*tMAXI+ud21m*(tMAXI*tMAXI)
	else
         fUm=ud01m+ud11m*tissl+ud21m*(tissl*tissl)
	endif
      ENDIF
      RETURN
      END
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      dtiss=1.2d0/100.d0
c      n=100
c      do j=1,2
c       tiss=3.35d0-dtiss
c       if(j.eq.1) giss=2.d0
c       if(j.eq.2) giss=4.d0
c       do i=1,n
c	tiss=tiss+dtiss
c	call sfra_u_usdss(tiss,giss,fUm)
c	write(93,*)giss,tiss,fUm
c       enddo
c       enddo
c      end
