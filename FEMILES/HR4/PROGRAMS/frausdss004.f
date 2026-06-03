c Calcula fraccion perdida en u SDSS por MIUSC a z=0.04
c Se pasan tiss,giss y devuelve fracciones
c Se uso Pickles library: convolucion para U y filtro U Johnson
c Filtro u SDSS: subrutina respusdss.f 
c Llamada desde hrsl.f
      SUBROUTINE frausdss004(tiss,giss,fUi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tissl=dlog10(tiss)
c	tissl=tiss
cGIANTS1_MIS2_MILES
      ug00i=-21.0043d0
      ug10i=11.4781d0
      ug20i=-1.53484d0
cGIANTS2_MIS2_MILES
      ug01i=-13.7803d0
      ug11i=6.51841d0
      ug21i=-0.737586d0
cDWARFS1_MIS2_MILES
      ud00i=-8.01475d0
      ud10i=4.49463d0
      ud20i=-0.593126d0
cDWARFS2_MIS2_MILES
      ud01i=-8.28794d0
      ud11i=3.92668d0
      ud21i=-0.432653d0
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
