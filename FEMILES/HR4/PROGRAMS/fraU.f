c Subrutina que calcula fraccion perdida en U por MILES Y MIUSC
c Se pasan tiss,giss y devuelve fracciones
c Se uso Pickles library: convolucion para U y filtro U Johnson
c Filtro Johnson U (Buser 1978): subrutina respu.f 
c Llamada desde hrsl.f
      SUBROUTINE fraU(tiss,giss,fUm,fUi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tissl=dlog10(tiss)
c	tissl=tiss
cGIANTS1_MILES
      ug00m=-13.9274d0
      ug10m=7.58171d0
      ug20m=-1.01602d0
cGIANTS2_MILES
      ug01m=-9.61709d0
      ug11m=4.43489d0
      ug21m=-0.492780d0
cDWARFS1_MILES
      ud00m=-10.5861d0
      ud10m=5.76763d0
      ud20m=-0.767486d0
cDWARFS2_MILES
      ud01m=-8.27303d0
      ud11m=3.80499d0
      ud21m=-0.419382d0
ccccccccccccccccccccccccccc    
cGIANTS1_MILES
      ug00i=-9.93241d0
      ug10i= 5.35284d0
      ug20i=-0.711383d0
cGIANTS2_MIS2_MILES
      ug01i=-7.01574d0
      ug11i= 3.22153d0
      ug21i=-0.356937d0
cDWARFS1_MIS1_MILES
      ud00i=-7.45403d0
      ud10i= 4.00057d0
      ud20i=-0.525438d0
cDWARFS2_MIS2_MILES
      ud01i=-6.13479d0
      ud11i= 2.81588d0
      ud21i=-0.310566d0
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
c        tiss=3.35d0-dtiss
c        if(j.eq.1) giss=2.d0
c        if(j.eq.2) giss=4.d0
c        do i=1,n
c         tiss=tiss+dtiss
c         call fraU(tiss,giss,fUm,fUi)
c         write(93,*)giss,tiss,fUm,fUi
c        enddo
c       enddo
c      end

