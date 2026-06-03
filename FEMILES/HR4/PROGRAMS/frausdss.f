c Subrutina que calcula fraccion perdida en u SDSS por MILES Y MIUSC
c Se pasan tiss,giss y devuelve fracciones
c Se uso Pickles library: convolucion para U y filtro U Johnson
c Filtro u SDSS: subrutina respusdss.f 
c Llamada desde hrsl.f
      SUBROUTINE frausdss(tiss,giss,fUm,fUi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tissl=dlog10(tiss)
c	tissl=tiss
cGIANTS1_MILES
      ug00m=-20.1458d0
      ug10m=10.9555d0
      ug20m=-1.46282d0
cGIANTS2_MILES
      ug01m=-12.3592d0
      ug11m=5.82078d0
      ug21m=-0.657489d0
cDWARFS1_MILES
      ud00m=-8.73430d0
      ud10m=4.79230d0
      ud20m=-0.628229d0
cDWARFS2_MILES
      ud01m=-7.77464d0
      ud11m=3.66036d0
      ud21m=-0.403701d0    
ccccccccccccccccccccccccccc    
cGIANTS1_MIS2_MILES
      ug00i=-15.6037d0
      ug10i=8.39083d0
      ug20i=-1.11022d0
cGIANTS2_MIS2_MILES
      ug01i=-9.99142d0
      ug11i=4.67950d0
      ug21i=-0.527290d0
cDWARFS1_MIS2_MILES
      ud00i=-7.35468d0
      ud10i=3.92412d0
      ud20i=-0.503812d0
cDWARFS2_MIS2_MILES
      ud01i=-6.34996d0
      ud11i=2.97132d0
      ud21i=-0.327592d0
ccccccccccccccccccccccccccc    
      tlimd=3.955d0
      tlimg=3.895d0
      tMINI=3.40d0
      tMAXI=4.50d0
      IF(giss.lt.3.5d0)THEN
	if(tissl.lt.tMINI)then
         fUi=ug00i+ug10i*tMINI+ug20i*(tMINI*tMINI)
        elseif(tissl.ge.tMINI.and.tissl.lt.tlimg)then
         fUi=ug00i+ug10i*tissl+ug20i*(tissl*tissl)
	elseif(tissl.gt.tMAXI)then
         fUi=ug01i+ug11i*tMAXI+ug21i*(tMAXI*tMAXI)
        else
         fUi=ug01i+ug11i*tissl+ug21i*(tissl*tissl)
        endif
      ELSE
        if(tissl.lt.tMINI)then
         fUi=ud00i+ud10i*tMINI+ud20i*(tMINI*tMINI)
	elseif(tissl.gt.tMINI.and.tissl.lt.tlimd)then
         fUi=ud00i+ud10i*tissl+ud20i*(tissl*tissl)
	elseif(tissl.gt.tMAXI)then
         fUi=ud01i+ud11i*tMAXI+ud21i*(tMAXI*tMAXI)
        else
         fUi=ud01i+ud11i*tissl+ud21i*(tissl*tissl)
        endif
      ENDIF
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
c	tiss=3.35d0-dtiss
c	if(j.eq.1) giss=2.d0
c	if(j.eq.2) giss=4.d0
c	do i=1,n
c	 tiss=tiss+dtiss
c	 call fraudsss(tiss,giss,fUm,fUi)
c	 write(93,*)giss,tiss,fUm,fUi
c	enddo
c       enddo
c      end
