c Subrutina que calcula el espectros de alta resolucion
c integrados. Hace uso de las subrutinas "lbusc" y "sigmac",
c las cuales calculan los espectros estelares rangos Jones, 
c CaT y MILES para una terna de parametros atmosfericos. 
c*******SUBRUTINA HRSL************************************************
      SUBROUTINE hrsl_alphaCO(npxm,k,wflujb,wflujv,wflujc,wflujm,al99,fjal99,
     &al99c,fja99c,al99m,fja99m,Wfcatt,qum,qu10,qu05,qu15,
     &fluxUm,fluxUi,fluxU,fluUms,fluUis,fluUs,fluUis004)
c      SUBROUTINE hrsl_alphaCO(npxm,k,wflujm,al99m,fja99m,qum,qu10,qu05
c     &,qu15,fluxUm,fluxU,fluUms,fluUs)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nstm=1999) !miles
      DIMENSION aam(nstm,5),starm(nstm) !miles
      DIMENSION ssppm(npxm,2),wflujm(npxm,2) !miles
      DIMENSION alm(3),al99m(3) !dimension 4?
      DIMENSION volp0(5) !calidad
      DIMENSION sigtg(3) !calidad
      DIMENSION wmejor(5,5) !evaluación soluciones
      CHARACTER*20 starm
      character*1 tlabn,glabn,flabn
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      COMMON/qsighr/xnbxst,volpam,volp10,volp05,volp15
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk !!!!!!!!!!!!!
      COMMON/MAMA/R140(12,99000),imam(12)
      COMMON/hrm1/aam,nstarm !miles
      COMMON/hrm2/starm !miles ***Se transmite nombre para ss o aa***
      COMMON/mic/mico
      dfmic4=0.3d0
      t50=5040.0d0
c limites lambdas
      B1=3600.d0
      B2=5550.d0
      DB=1010.0d0
      V1=4750.d0
      V2=7400.d0
      DV=870.0d0
      fiss=fez(k)
      qum=0.0d0 !parametro calidad (MILES) normalizado sigma_min
      qu10=0.0d0 !parametro calidad (MILES) aceptable 1.0 sigma
      qu05=0.0d0 !parametro calidad (MILES) aceptable 0.5 sigma
      qu15=0.0d0 !parametro calidad (MILES) aceptable 1.5 sigma
      starmin=3.0d0 !minimo aceptable estrellas: 3 (1 x subcaja)
      fja99m=0.0d0
c U and missing U filter fluxes in MILES and MIUSC
      fluxUm=0.0d0 !Missing U MILES flux
      fluxU=0.0d0 !Total U flux
      fluUms=0.0d0 !Missing u sdss MILES flux
      fluUs=0.0d0 !Total u sdss flux
c
      do l=1,3
	al99m(l)=0.d0 !miles
      enddo
      do l=1,npxm !miles
	wflujm(l,2)=0.d0
      enddo
      zeroi=0.0d0
      DO m=1,jotai(k)
	tiss=10.**R(k,m,3)
	giss=R(k,m,4)
	biss=(-2.5d0)*dlog10(R(k,m,6))
	viss=(-2.5d0)*dlog10(R(k,m,7))
	uvbus=(-2.5d0)*dlog10(R(k,m,5)/R(k,m,7))
c	call xflabs(uvbus,biss,viss,vikc,flujb,flujv,flujcj,fluju)
      call xflabs(uvbus,biss,viss,vikc,vk,flujb,flujv,flujcj,fluju,
     &flujk)
        call fracat(tiss,giss,frc)
        call fraU(tiss,giss,fUm,fUi)
        call fra_u_usdss(tiss,giss,fusdss)
        call frausdss(tiss,giss,fUms,fUis)
        if(giss.gt.4.25d0)then
         ggiss=giss+0.1d0
         dgiss=-0.1d0
        elseif(giss.ge.3.25d0.and.giss.le.4.25d0)then
         ggiss=giss-0.1d0
         dgiss=+0.1d0
        elseif(giss.ge.1.25d0.and.giss.lt.3.25d0)then
         ggiss=giss+0.1d0
         dgiss=-0.1d0
        else
         ggiss=giss-0.1d0
        dgiss=+0.1d0
        endif
c para elegir mejor set de parametros en el bucle de afinamiento:
c par. optimos segun dif metalicidad (se calcularan estos si no hay solucion)     
        autiss=0.0d0 !temperatura f
        augiss=0.0d0 !gravedad f
c 
        wmej=-999999.
        DO k1=1,5
         ggiss=ggiss+dgiss
         ttiss=tiss-1.0d0
         DO k2=1,5
          nk1=0
          ttiss=ttiss+1.0d0
	  IF(mico.eq.4)THEN
	   CALL sigmam_alphaCO(aam,starm,nstarm,ttiss,ggiss,
     &                    fiss-dfmic4,ssppm,alm,sigtg,volp0)
          ELSE
	   CALL sigmam_alphaCO(aam,starm,nstarm,ttiss,ggiss,
     &                    fiss,ssppm,alm,sigtg,volp0)
          ENDIF	
          IF(mico.eq.4)THEN
c           write(55,'(3(A10,F8.1,2(1X,F5.2,1X)))')
c     &     'REQUESTED=',tiss,giss,fiss-dfmic4,
c     &     'TRIED =   ',ttiss,ggiss,fiss-dfmic4,
c     &     'COMPUTED =',(alm(ni),ni=1,3)
	   if(abs(alm(1)-tiss).le.2.5d0.and.
     &      abs(alm(2)-giss).le.0.101d0.and.
     &      abs(alm(3)-(fiss-dfmic4)).le.0.01d0)then
            nk1=1
            goto 7777
	   endif
	  ELSE
c           write(55,'(3(A10,F8.1,2(1X,F5.2,1X)))')
c     &     'REQUESTED=',tiss,giss,fiss,
c     &     'TRIED =   ',ttiss,ggiss,fiss,
c     &     'COMPUTED =',(alm(ni),ni=1,3)
	   if(abs(alm(1)-tiss).le.2.5d0.and.
     &      abs(alm(2)-giss).le.0.101d0.and.
     &      abs(alm(3)-fiss).le.0.01d0)then
            nk1=1
            goto 7777
	   endif
	  ENDIF
	  IF(mico.eq.4)THEN
           wmejor(k1,k2)
     &     =exp(-((((t50/alm(1))-(t50/tiss))/sigtg(1))**2.d0))
     &     *exp(-(((alm(2)-giss)/(3.0d0*sigtg(2)))**2.d0))     
     &     *exp(-(((alm(3)-(fiss-dfmic4))/(0.75d0*sigtg(3)))**2.d0))
          ELSE
           wmejor(k1,k2)
     &     =exp(-((((t50/alm(1))-(t50/tiss))/sigtg(1))**2.d0))
     &     *exp(-(((alm(2)-giss)/(3.0d0*sigtg(2)))**2.d0))     
     &     *exp(-(((alm(3)-fiss)/(0.75d0*sigtg(3)))**2.d0))
          ENDIF
          if(wmejor(k1,k2).gt.wmej)then
	     wmej=wmejor(k1,k2)
	     autiss=ttiss
	     augiss=ggiss
	  endif
       ENDDO
      ENDDO
cccccccccccccc
7777  continue
      if(nk1.eq.0)then
       if(mico.eq.4)then
        CALL sigmam_alphaCO(aam,starm,nstarm,autiss,augiss,
     &                      fiss-dfmic4,ssppm,alm,sigtg,volp0)
       else
         CALL sigmam_alphaCO(aam,starm,nstarm,autiss,augiss,
     &                       fiss,ssppm,alm,sigtg,volp0)
       endif
      endif
             if(abs(alm(1)-tiss).gt.10.1d0)then
	     tlabn='t'
	     labn=1
	    else
	     tlabn=' '
	     labn=0
	    endif
            if(abs(alm(2)-giss).gt.0.11d0)then
	     glabn='g'
	     labn=1
	    else
	     glabn=' '
	     labn=0
	    endif
      IF(mico.eq.4)THEN
            if(abs(alm(3)-(fiss-dfmic4)).gt.0.011d0)then
	     flabn='f'
	     labn=1
	    else
	     flabn=' '
	     labn=0
	    endif
c          write(*,'(2(A10,F8.1,2(1X,F5.2,1X)),3(A1))')
c     &    'REQUESTED=',tiss,giss,fiss-dfmic4,
c     &    'COMPUTED =',(alm(ni),ni=1,3)
c     &     ,tlabn,glabn,flabn
          if(labn.eq.1)write(57,'(2(A10,F8.1,2(1X,F5.2,1X)),3(A1))')
     &    'REQUESTED=',tiss,giss,fiss-dfmic4,
     &    'COMPUTED =',(alm(ni),ni=1,3)
     &     ,tlabn,glabn,flabn
      ELSE
            if(abs(alm(3)-fiss).gt.0.011d0)then
	     flabn='f'
	     labn=1
	    else
	     flabn=' '
	     labn=0
	    endif
c         write(*,'(2(A10,F8.1,2(1X,F5.2,1X)),3(A1))')
c     &    'REQUESTED=',tiss,giss,fiss,
c     &    'COMPUTED =',(alm(ni),ni=1,3)
c     &     ,tlabn,glabn,flabn
           if(labn.eq.1)write(57,'(2(A10,F8.1,2(1X,F5.2,1X)),3(A1))')
     &    'REQUESTED=',tiss,giss,fiss,
     &    'COMPUTED =',(alm(ni),ni=1,3)
     &     ,tlabn,glabn,flabn
      ENDIF
c Normalizamos espectros miles
        fbmil=0.0d0
        fvmil=0.0d0
	do l=1,npxm
	 rfB=0.0d0
	 rfV=0.0d0
	 if(l.eq.npxm)then
	  dl=abs(ssppm(l-1,1)-ssppm(l,1))
	 else
	  dl=abs(ssppm(l+1,1)-ssppm(l,1))
	 endif
	 if(ssppm(l,1).ge.B1.and.ssppm(l,1).le.B2)then
	  call respb(ssppm(l,1),rfB)
	 endif
	 if(ssppm(l,1).ge.V1.and.ssppm(l,1).le.V2)then
	  call respv(ssppm(l,1),rfV)
	 endif
         fbmil=fbmil+rfB*ssppm(l,2)*dl
         fvmil=fvmil+rfV*ssppm(l,2)*dl
	enddo
	DO l=1,npxm !miles
	 ssppm(l,2)=ssppm(l,2)/fvmil
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R140(k,mma)*flujv
	   enddo
	  endif
	  wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,1,14)*flujv
	 else
	  wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,m,14)*flujv
	 endif
	ENDDO
        if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    fja99m=fja99m+R140(k,mma)*flujv
	   enddo
	  endif
	  fja99m=fja99m+R(k,1,14)*flujv
	else
	  fja99m=fja99m+R(k,m,14)*flujv
	endif
	do l=1,3 !miles
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    al99m(l)=al99m(l)+alm(l)*R140(k,mma)*flujv
	   enddo
	  endif
	  al99m(l)=al99m(l)+alm(l)*R(k,1,14)*flujv
	 else
	  al99m(l)=al99m(l)+alm(l)*R(k,m,14)*flujv
	 endif
	enddo
c Missing MILES and MIUSC U fluxes:
        IF(m.eq.1)THEN
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    fluxUm=fluxUm+R140(k,mma)*fluju*fUm
	    fluxUi=fluxUi+R140(k,mma)*fluju*fUi
	    fluxU=fluxU+R140(k,mma)*fluju
	    fluUms=fluUms+R140(k,mma)*fluju*fusdss*fUms
	    fluUis=fluUis+R140(k,mma)*fluju*fusdss*fUis
	    fluUis004=fluUis004+R140(k,mma)*fluju*fusdss*fUis004
	    fluUs=fluUs+R140(k,mma)*fluju*fusdss
	   enddo
	  endif
	  fluxUm=fluxUm+R(k,1,14)*fluju*fUm
	  fluxUi=fluxUi+R(k,1,14)*fluju*fUi
	  fluxU=fluxU+R(k,1,14)*fluju
	  fluUms=fluUms+R(k,1,14)*fluju*fusdss*fUms
	  fluUis=fluUis+R(k,1,14)*fluju*fusdss*fUis
	  fluUis004=fluUis004+R(k,1,14)*fluju*fusdss*fUis004
	  fluUs=fluUs+R(k,1,14)*fluju*fusdss
	ELSE
	  fluxUm=fluxUm+R(k,m,14)*fluju*fUm
	  fluxUi=fluxUi+R(k,m,14)*fluju*fUi
	  fluxU=fluxU+R(k,m,14)*fluju
	  fluUms=fluUms+R(k,m,14)*fluju*fusdss*fUms
	  fluUis=fluUis+R(k,m,14)*fluju*fusdss*fUis
	  fluUis004=fluUis004+R(k,m,14)*fluju*fusdss*fUis004
	  fluUs=fluUs+R(k,m,14)*fluju*fusdss
	ENDIF
	qum=qum+R(k,m,14)*flujv*(xnbxst/volpam)
	qu10=qu10+R(k,m,14)*flujv*(starmin/volp10)
	qu05=qu05+R(k,m,14)*flujv*(starmin/volp05)
	qu15=qu15+R(k,m,14)*flujv*(starmin/volp15)
      ENDDO
      RETURN
      END
