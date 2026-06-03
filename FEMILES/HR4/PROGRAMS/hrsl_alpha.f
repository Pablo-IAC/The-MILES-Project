c Subrutina que calcula espectros integrados. Usa subrutinas "lbusc" y "sigmac",
c que calculan espectros estelares rangos Jones, CaT y MILES para una terna de parametros atmosfericos. 
c*******SUBRUTINA HRSL************************************************
      SUBROUTINE hrsl_alpha(k,wflujb,wflujv,wflujc,wflujm,wflujcs,
     &wflujms,al99,fjal99,al99c,fja99c,al99m,fja99m,Wfcatt,qum,qu10,
     &qu05,qu15)
cc      SUBROUTINE hrsl_alpha(k,wflujb,wflujv,wflujc,wflujm,wflujcs,
cc     &wflujms,al99,fjal99,al99c,fja99c,al99m,fja99m,Wfcatt,qum,qu10,
cc     &qu05,qu15,fluxUm,fluxUi,fluxU,fluUms,fluUis,fluUs,fluUis004)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nstm=1999) !miles
      PARAMETER (npxm=4300) !miles
      PARAMETER (npxj=1107) !jones
      PARAMETER (npxc=710) !cat
      DIMENSION ssppm(npxm,2),ssppml(npxm,2),ssppmh(npxm,2) !miles
      DIMENSION wflujm(npxm,2),wflujms(npxm,2) !miles
      DIMENSION aam(nstm,6),aamlo(nstm,6),aamhi(nstm,6)
      DIMENSION starm(nstm) !miles: added Mg/Fe column
      DIMENSION ssppb(1120,2),wflujb(1120,2) !jones b
      DIMENSION ssppv(1120,2),wflujv(1120,2) !jones v
      DIMENSION ssppc(710,2),wflujc(710,2),wflujcs(710,2) !cat
      DIMENSION al(3),al99(3),alc(3),al99c(3)
      DIMENSION alm(5),alml(5),almh(5),al99m(5),altest(5)
      DIMENSION sigtg(3),sigtgl(3),sigtgh(3)
      DIMENSION volp0(5),volpl(5),volph(5)
      CHARACTER*80 starm
      character*1 tlabn,glabn,flabn,alabn
      COMMON/enhanc/amgfe,DeltFeH ![Mg/Fe],[Fe/H] correc. value
      COMMON/enhand/damgfe !overlapping Mg/Fe value for 4D interpolation
      COMMON/mimgfe/aam,aamlo,aamhi !miles interpolador 4D
      COMMON/hrm1/nstarm !miles
      COMMON/hrm2/starm !miles
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      common/cfnc/fnc11(2999),fnmg(2999)
      common/cd4000/fra4(2999),fbj(2999)
c      COMMON/qsighr/xnbxst,volpam,volp10,volp05,volp15 !ahora en volp0(5)
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,vkvega,
     &zerovk
      COMMON/MAMA/R140(12,99000),imam(12)
      COMMON/ab/iaba,iabaj
      COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
      COMMON/part0/ipartial
      COMMON/part1/pmassL,pmassH
c      COMMON/FNOR/C2Code,FNORMM,FNORMC
c limites lambdas
      B1=3600.d0
      B2=5550.d0
      DB=1010.0d0
c      BJ1=3856.d0
c      BJ2=4476.d0
      V1=4750.d0
      V2=7400.d0
      DV=870.0d0
c      VJ1=4795.d0
c      VJ2=5465.d0
      XIJ1=6800.0d0 
      XIJ2=12000.0d0 
      C1=8475.0d0
      C2=8807.0d0
      DIC=1226.0d0
c
      fbuser=1.0d0/10.**(0.71d0/2.5d0) !FB=FB3*0.52
ccccccccccccccccccccccccc
c [Fe/H] (si [Mg/Fe].ne.0) a partir de metalicidad isocrona [Z/H]=fez(k)
c solo para sigmam_alpha.f
      fiss=fez(k)
      fiss=fez(k)-DeltFeH
ccccccccccccccccccccccccc
      qum=0.0d0 !parametro calidad (MILES) normalizado sigma_min
      qu10=0.0d0 !parametro calidad (MILES) aceptable 1.0 sigma
      qu05=0.0d0 !parametro calidad (MILES) aceptable 0.5 sigma
      qu15=0.0d0 !parametro calidad (MILES) aceptable 1.5 sigma
      starmin=3.0d0 !minimo aceptable estrellas: 3 (1 x subcaja)
c	if(jotai(k).eq.0)goto 40123
      fjal99=0.0d0
      fja99c=0.0d0
      fja99m=0.0d0
      Wfcatt=0.0d0
c U and missing U filter fluxes in MILES and MIUSC
cc      fluxUm=0.0d0 !Missing U MILES flux
cc      fluxUi=0.0d0 !Missing U MIUSC flux
cc      fluxU=0.0d0 !Total U flux
cc      fluUms=0.0d0 !Missing u sdss MILES flux
cc      fluUis=0.0d0 !Missing u sdss MIUSC flux
cc      fluUis004=0.0d0 !Missing u sdss MIUSC flux at z=0.04
cc      fluUs=0.0d0 !Total u sdss flux
c
      do l=1,3
	al99(l)=0.d0 !jones
	al99c(l)=0.d0 !cat
      enddo
      do l=1,5
	al99m(l)=0.d0 !miles
      enddo
      do l=1,npxj !jones
	wflujb(l,2)=0.d0
	wflujv(l,2)=0.d0
      enddo
      do l=1,npxc !cat
	wflujc(l,2)=0.d0
	wflujcs(l,2)=0.d0
      enddo
      do l=1,npxm !miles
	wflujm(l,2)=0.d0
	wflujms(l,2)=0.d0
      enddo
      zeroi=0.0d0
      DO m=1,jotai(k)
       tiss=10.**R(k,m,3)
       giss=R(k,m,4)
       biss=(-2.5d0)*dlog10(R(k,m,6))
       viss=(-2.5d0)*dlog10(R(k,m,7))
       rij=zeroi
       vrkc=zeroi
       vrj=zeroi
       vij=zeroi
       ckc=zeroi
       dkc=zeroi
       bkc=zeroi
       akc=zeroi
       vrkc=(-2.5d0)*dlog10(R(k,m,7)/R(k,m,8))
       vikc=(-2.5d0)*dlog10(R(k,m,7)/R(k,m,9))
       uvbus=(-2.5d0)*dlog10(R(k,m,5)/R(k,m,7))
c      call cotojo(giss,vrkc,vikc,vij)
       CALL xflabs(uvbus,biss,viss,vikc,vk,flujb,flujv,flujcj,fluju,
     &flujk)
       IF(ipartial.eq.1)THEN
	IF(Rm0(k,m,1).lt.pmassL.or.Rm0(k,m,1).gt.pmassH)THEN
	   flujb=0.0d0
	   flujv=0.0d0
	   flujcj=0.0d0
	   fluju=0.0d0
	   flujk=0.0d0
	ENDIF
       ENDIF
       CALL fracat(tiss,giss,frc)
       CALL fraU(tiss,giss,fUm,fUi)
       CALL fra_u_usdss(tiss,giss,fusdss)
       CALL frausdss(tiss,giss,fUms,fUis)
       CALL frausdss004(tiss,giss,fUis004)
       flujc=frc*flujcj
       fjal=(flujv/DV-flujb/DB)*(224.d0/993.d0)+flujb/DB
c Preparando la dimension Mg/Fe
       wlow=0.5d0
       whig=0.5d0
       IF(iaba.eq.9)THEN
        CALL sigmam_alpha
     &(aamlo,tiss,giss,fiss,ssppml,alml,sigtgl,volpl)
        CALL sigmam_alpha
     &(aamhi,tiss,giss,fiss,ssppmh,almh,sigtgh,volph)
        CALL sigmam_alpha
     &(aam,tiss,giss,fiss,ssppm,alm,sigtg,volp0)
        nwarl=0
        nwarh=0
cc	 do isig=1,3
cc	  sigtg(isig)=MIN(sigtgl(isig),sigtgh(isig))
cc	 enddo
cc	 do imgfe=1,5
cc	   volp0(imgfe)=volpl(imgfe)!to update volp0(1),volp0(2) if 2* combined
cc	 enddo
cc	 do imgfe=1,npxm
cc	  ssppm(imgfe,1)=ssppml(imgfe,1)
cc	 enddo
cccc assessing 3-D solution reliability:
        wtesth=(amgfe-alml(4))/(almh(4)-alml(4))
        wtestl=1.d0-wtesth
        DO imgfe=1,5
	 altest(imgfe)=999.9999
	 altest(imgfe)=wtestl*alml(imgfe)+wtesth*almh(imgfe)
        ENDDO
        IF(abs(tiss-altest(1)).gt.sigtg(1).or.
     &abs(giss-altest(2)).gt.sigtg(2)*5.0d0.or.
     &abs(fiss-altest(3)).gt.sigtg(3)*0.25d0)THEN
	   if(abs(tiss-alml(1)).gt.sigtg(1).or.
     &abs(giss-alml(2)).gt.sigtg(2)*5.0d0.or.
     &abs(fiss-alml(3)).gt.sigtg(3)*0.25d0.or.
     &abs(fiss-alml(3)).ge.abs(amgfe-alml(4)))then
	         nwarl=1
	   else
	         nwarl=0
	   endif
	   if(abs(tiss-almh(1)).gt.sigtg(1).or.
     &abs(giss-almh(2)).gt.sigtg(2)*5.0d0.or.
     &abs(fiss-almh(3)).gt.sigtg(3)*0.25d0.or.
     &abs(fiss-almh(3)).ge.abs(amgfe-almh(4)))then
	         nwarh=1
	   else
	         nwarh=0
	   endif
        ELSE
	   nwarl=0
	   nwarh=0
        ENDIF
c
	IF(nwarl.eq.1.and.nwarh.eq.1)THEN
	  goto 2551
	ELSEIF(nwarl.eq.0.and.nwarh.ne.0)THEN
	  wlow=1.0d0
	  whig=0.0d0
	  do imgfe=1,5
	   alm(imgfe)=alml(imgfe)
	  enddo
	  do imgfe=1,5
	   volp0(imgfe)=volpl(imgfe)
	  enddo
	  do imgfe=1,npxm
	   ssppm(imgfe,2)=ssppml(imgfe,2)
	  enddo
	ELSEIF(nwarl.ne.0.and.nwarh.eq.0)THEN
	  wlow=0.0d0
	  whig=1.0d0
	  do imgfe=1,5
	   alm(imgfe)=almh(imgfe)
	  enddo
	  do imgfe=1,5
	   volp0(imgfe)=volph(imgfe)
	  enddo
	  do imgfe=1,npxm
	   ssppm(imgfe,2)=ssppmh(imgfe,2)
	  enddo
	ELSE
	 IF(abs(almh(4)-alml(4)).lt.(damgfe/10.d0))THEN
	  wlow=0.5d0
	  whig=0.5d0
	  do imgfe=1,5
	   alm(imgfe)=wlow*alml(imgfe)+whig*almh(imgfe)
	  enddo
c	  volp0(1)=MAX(volph(1),volpl(1))
c	  volp0(2)=MIN(volph(2),volpl(2))
	  volp0(1)=wlow*volpl(1)+whig*volph(1)
	  volp0(2)=wlow*volpl(2)+whig*volph(2)
	  do imgfe=1,npxm
	   ssppm(imgfe,2)=wlow*ssppml(imgfe,2)+whig*ssppmh(imgfe,2)
	  enddo
	 ELSE
	  IF(almh(4).lt.amgfe.and.alml(4).lt.amgfe)THEN
	   if(abs(almh(4)-amgfe).lt.abs(alml(4)-amgfe))then
	    wlow=0.0d0
	    whig=1.0d0
	    do imgfe=1,5
	     alm(imgfe)=almh(imgfe)
	    enddo
	    do imgfe=1,5
	     volp0(imgfe)=volph(imgfe)
	    enddo
	    do imgfe=1,npxm
	     ssppm(imgfe,2)=ssppmh(imgfe,2)
	    enddo
	   else
	    wlow=1.0d0
	    whig=0.0d0
	    do imgfe=1,5
	     alm(imgfe)=alml(imgfe)
	    enddo
	    do imgfe=1,5
	     volp0(imgfe)=volpl(imgfe)
	    enddo
	    do imgfe=1,npxm
	     ssppm(imgfe,2)=ssppml(imgfe,2)
	    enddo
	   endif
	  ELSEIF(almh(4).gt.amgfe.and.alml(4).gt.amgfe)THEN
	   if(abs(almh(4)-amgfe).lt.abs(alml(4)-amgfe))then
	    wlow=0.0d0
	    whig=1.0d0
	    do imgfe=1,5
	     alm(imgfe)=almh(imgfe)
	    enddo
	    do imgfe=1,5
	     volp0(imgfe)=volph(imgfe)
	    enddo
	    do imgfe=1,npxm
	     ssppm(imgfe,2)=ssppmh(imgfe,2)
	    enddo
	   else
	    wlow=1.0d0
	    whig=0.0d0
	    do imgfe=1,5
	     alm(imgfe)=alml(imgfe)
	    enddo
	    do imgfe=1,5
	     volp0(imgfe)=volpl(imgfe)
	    enddo
	    do imgfe=1,npxm
	     ssppm(imgfe,2)=ssppml(imgfe,2)
	    enddo
	   endif
	  ELSE
c	  whig=exp(-(((aloeci(1,lp)-tisst)/sigmat)**2.))
c     &  *exp(-(((aloeci(2,lp)-giss)/sigmag)**2.))
c     &  *exp(-(((aloeci(3,lp)-fiss)/sigmaz)**2.))
c          whig=exp(-(((almh(4)-amgfe)/damgfe)**2.))
c          wlow=exp(-(((alml(4)-amgfe)/damgfe)**2.))
c         whig=whig/(whig+wlow)
c	  wlow=wlow/(whig+wlow)
	   whig=(amgfe-alml(4))/(almh(4)-alml(4))
	   wlow=1.d0-whig
	   do imgfe=1,5
	    alm(imgfe)=wlow*alml(imgfe)+whig*almh(imgfe)
	   enddo
	   volp0(1)=wlow*volpl(1)+whig*volph(1)
	   volp0(2)=wlow*volpl(2)+whig*volph(2)
c	   volp0(1)=MAX(volph(1),volpl(1))
c	   volp0(2)=MIN(volph(2),volpl(2))
	   do imgfe=1,npxm
	    ssppm(imgfe,2)=wlow*ssppml(imgfe,2)+whig*ssppmh(imgfe,2)
	   enddo
	  ENDIF
	 ENDIF
        ENDIF
       ELSE
        CALL sigmam_alpha(aam,tiss,giss,fiss,ssppm,alm,sigtg,volp0)
       ENDIF
2551   continue
c NOTICE NOT WELL REPRODUCED STELLAR PARAMETERS:
       if(abs(alm(1)-tiss).gt.25.1d0)then
	     tlabn='t'
	     labn=1
       else
	     tlabn=' '
	     labn=0
       endif
       if(abs(alm(2)-giss).gt.0.15d0)then
	     glabn='g'
	     labn=1
       else
	     glabn=' '
	     labn=0
       endif
       if(abs(alm(3)-fiss).gt.0.025d0)then
	     flabn='f'
	     labn=1
       else
	     flabn=' '
	     labn=0
       endif
       if(abs(alm(4)-amgfe).gt.0.025d0)then
	     alabn='a'
	     labn=1
       else
	     alabn=' '
	     labn=0
       endif
1319   FORMAT(4A1,4(1X,A3,1X,F7.0,4(1X,F6.2)),2(2X,F10.6),3(1X,F9.3))
c      IF(labn.eq.1)then
c        WRITE(57,1319)tlabn,glabn,flabn,alabn,
c     &'REQ',tiss,giss,fiss,amgfe,fiss+0.75d0*amgfe,
c     &'SYN',(alm(ikol),ikol=1,5),
c     &'LOW',(alml(nl),nl=1,5),
c     &'HIG',(almh(nl),nl=1,5),
c     &wlow,whig,(sigtg(nl),nl=1,3)
c      ENDIF
c      if(m.eq.32)then
c        write(72,'(A3,X,F6.0,5(X,F6.3),2(X,F10.6))')
c     &  'alm',(alm(nl),nl=1,4),alml(4),almh(4),wlow,whig
c       do nl=1,npxm
c        write(72,*)ssppm(nl,1),ssppm(nl,2)
c       enddo
c      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       IF(iaba.eq.0.and.iabaj.eq.1)THEN
	 CALL lbusc(tiss,giss,fiss,viss,ssppb,ssppv,al,fra4000)
	 fra4(m)=fra4000
	 fbj(m)=flujb
       ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       CALL sigmac(tiss,giss,fiss,ssppc,alc,Wcatt,fnc1,fnmgi)
       fnc11(m)=fnc1
       fnmg(m)=fnmgi
       IF(iaba.eq.0.and.iabaj.eq.1)THEN
        BJ1=ssppb(1,1)
c        BJ1=4000.
        BJ2=ssppb(npxj,1)
        VJ1=ssppv(1,1)
        VJ2=ssppv(npxj,1)
c Normalizamos espectros jones
        fbjon=0.0d0
        fvjon=0.0d0
        do l=1,npxj
         if(l.eq.npxj)then
	  dlb=abs(ssppb(l-1,1)-ssppb(l,1))
	  dlv=abs(ssppv(l-1,1)-ssppv(l,1))
         else
	  dlb=abs(ssppb(l+1,1)-ssppb(l,1))
	  dlv=abs(ssppv(l+1,1)-ssppv(l,1))
         endif
         CALL respb(ssppb(l,1),rfB) !respob factor entre 0 y 1
         CALL respv(ssppv(l,1),rfV) !respob factor entre 0 y 1
         if(ssppb(l,1).ge.BJ1.and.ssppb(l,1).le.BJ2)then
	  fbjon=fbjon+rfB*ssppb(l,2)*dlb
         endif
         fvjon=fvjon+rfV*ssppv(l,2)*dlv
        enddo
        do l=1,npxj
         ssppb(l,2)=ssppb(l,2)/(fbjon*fbuser)
         ssppv(l,2)=ssppv(l,2)/fvjon
        enddo
       ENDIF
c cat
       fcat=0.0d0
       rfC=0.0d0
       do l=1,npxc
	 if(ssppc(l,1).ge.C1.and.ssppc(l,1).le.C2)then
	  if(l.eq.npxc)then
	   dl=abs(ssppc(l-1,1)-ssppc(l,1))
	  else
	   dl=abs(ssppc(l+1,1)-ssppc(l,1))
	  endif
	  call respic(ssppc(l,1),rfC)
	  if(rfC.gt.0.0d0) fcat=fcat+rfC*ssppc(l,2)*dl
	 endif
       enddo
c miles
       fbjon=0.0d0
       fvjon=0.0d0
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
          IF(iaba.eq.0.and.iabaj.eq.1)THEN
	   if(ssppm(l,1).ge.BJ1.and.ssppm(l,1).le.BJ2)then
	   if(rfB.gt.0.0d0) fbjon=fbjon+rfB*ssppm(l,2)*dl
	   endif
          ENDIF
	 endif
	 if(ssppm(l,1).ge.V1.and.ssppm(l,1).le.V2)then
	  call respv(ssppm(l,1),rfV)
          IF(iaba.eq.0.and.iabaj.eq.1)THEN
	   if(ssppm(l,1).ge.VJ1.and.ssppm(l,1).le.VJ2)then
	    if(rfV.gt.0.0d0) fvjon=fvjon+rfV*ssppm(l,2)*dl
	   endif
	  ENDIF
	 endif
         if(rfB.gt.0.0d0) fbmil=fbmil+rfB*ssppm(l,2)*dl
         if(rfV.gt.0.0d0) fvmil=fvmil+rfV*ssppm(l,2)*dl
       enddo
       IF(iaba.eq.0.and.iabaj.eq.1)THEN
        frb=(fbjon*fbuser)/(fbmil*fbuser)
        frv=fvjon/fvmil
	DO l=1,npxj !jones
c	 wflujb(l,2)=wflujb(l,2)+ssppb(l,2)*R(k,m,14)*flujb*frb*FNORM
c	 wflujv(l,2)=wflujv(l,2)+ssppv(l,2)*R(k,m,14)*flujv*frv*FNORM
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	  wflujb(l,2)=wflujb(l,2)+ssppb(l,2)*R140(k,mma)*flujb*frb
	  wflujv(l,2)=wflujv(l,2)+ssppv(l,2)*R140(k,mma)*flujv*frv
	   enddo
	  endif
	  wflujb(l,2)=wflujb(l,2)+ssppb(l,2)*R(k,1,14)*flujb*frb
	  wflujv(l,2)=wflujv(l,2)+ssppv(l,2)*R(k,1,14)*flujv*frv
	 else
	  wflujb(l,2)=wflujb(l,2)+ssppb(l,2)*R(k,m,14)*flujb*frb
	  wflujv(l,2)=wflujv(l,2)+ssppv(l,2)*R(k,m,14)*flujv*frv
	 endif
	ENDDO
       ENDIF
       IF(fcat.gt.0.0d0)THEN
        DO l=1,npxc !cat
c	 wflujc(l,2)=wflujc(l,2)+ssppc(l,2)*R(k,m,14)*flujc*FNORMC
         ssppc(l,2)=(ssppc(l,2)/fcat)
	 spfl=ssppc(l,2)*flujc
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wflujc(l,2)=wflujc(l,2)+R140(k,mma)*spfl
	    wflujcs(l,2)=wflujcs(l,2)+R140(k,mma)*spfl*spfl
	   enddo
	  endif
	  wflujc(l,2)=wflujc(l,2)+R(k,1,14)*spfl
	  wflujcs(l,2)=wflujcs(l,2)+R(k,1,14)*spfl*spfl
	 else
	  wflujc(l,2)=wflujc(l,2)+R(k,m,14)*spfl
	  wflujcs(l,2)=wflujcs(l,2)+R(k,m,14)*spfl*spfl
	 endif
        ENDDO
       ENDIF
       IF(fvmil.gt.0.0d0)THEN
        DO l=1,npxm !miles
c	 ssppm(l,2)=(ssppm(l,2)/fvmil)*FVNOR
c	 wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,m,14)*flujv*FNORMM
c	 wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,m,14)*R(k,m,7)
	 ssppm(l,2)=ssppm(l,2)/fvmil
	 spfl=ssppm(l,2)*flujv
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wflujm(l,2)=wflujm(l,2)+R140(k,mma)*spfl
	    wflujms(l,2)=wflujms(l,2)+R140(k,mma)*spfl*spfl
	   enddo
	  endif
	  wflujm(l,2)=wflujm(l,2)+R(k,1,14)*spfl
	  wflujms(l,2)=wflujms(l,2)+R(k,1,14)*spfl*spfl
	 else
	  wflujm(l,2)=wflujm(l,2)+R(k,m,14)*spfl
	  wflujms(l,2)=wflujms(l,2)+R(k,m,14)*spfl*spfl
	 endif
        ENDDO
       ENDIF
       if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    Wfcatt=Wfcatt+Wcatt*R140(k,mma)*flujc
	    fjal99=fjal99+R140(k,mma)*fjal
	    fja99c=fja99c+R140(k,mma)*flujc
	    fja99m=fja99m+R140(k,mma)*flujv
	   enddo
	  endif
	  Wfcatt=Wfcatt+Wcatt*R(k,1,14)*flujc
          IF(iaba.eq.0.and.iabaj.eq.1) fjal99=fjal99+R(k,1,14)*fjal
	  fja99c=fja99c+R(k,1,14)*flujc
	  fja99m=fja99m+R(k,1,14)*flujv
       else
	  Wfcatt=Wfcatt+Wcatt*R(k,m,14)*flujc
          IF(iaba.eq.0.and.iabaj.eq.1) fjal99=fjal99+R(k,m,14)*fjal
	  fja99c=fja99c+R(k,m,14)*flujc
	  fja99m=fja99m+R(k,m,14)*flujv
       endif
       IF(iaba.eq.0.and.iabaj.eq.1)THEN
	do l=1,3 !jones
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    al99(l)=al99(l)+al(l)*R140(k,mma)*fjal
	   enddo
	  endif
	  al99(l)=al99(l)+al(l)*R(k,1,14)*fjal
	 else
	  al99(l)=al99(l)+al(l)*R(k,m,14)*fjal
	 endif
	enddo
       ENDIF
       DO l=1,3 !cat
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    al99c(l)=al99c(l)+alc(l)*R140(k,mma)*flujc
	   enddo
	  endif
	  al99c(l)=al99c(l)+alc(l)*R(k,1,14)*flujc
	 else
	  al99c(l)=al99c(l)+alc(l)*R(k,m,14)*flujc
	 endif
       ENDDO
      IF(fvmil.gt.0.0d0)THEN
       DO l=1,5 !miles
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
       ENDDO
c Missing MILES and MIUSC U fluxes:
cc       IF(m.eq.1)THEN
cc	  if(imam(k).gt.0)then
cc	   do mma=1,imam(k)
cc	    fluxUm=fluxUm+R140(k,mma)*fluju*fUm
cc	    fluxUi=fluxUi+R140(k,mma)*fluju*fUi
cc	    fluxU=fluxU+R140(k,mma)*fluju
cc	    fluUms=fluUms+R140(k,mma)*fluju*fusdss*fUms
cc	    fluUis=fluUis+R140(k,mma)*fluju*fusdss*fUis
cc	    fluUis004=fluUis004+R140(k,mma)*fluju*fusdss*fUis004
cc	    fluUs=fluUs+R140(k,mma)*fluju*fusdss
cc	   enddo
cc	  endif
cc	  fluxUm=fluxUm+R(k,1,14)*fluju*fUm
cc	  fluxUi=fluxUi+R(k,1,14)*fluju*fUi
cc	  fluxU=fluxU+R(k,1,14)*fluju
cc	  fluUms=fluUms+R(k,1,14)*fluju*fusdss*fUms
cc	  fluUis=fluUis+R(k,1,14)*fluju*fusdss*fUis
cc	  fluUis004=fluUis004+R(k,1,14)*fluju*fusdss*fUis004
cc	  fluUs=fluUs+R(k,1,14)*fluju*fusdss
cc       ELSE
cc	  fluxUm=fluxUm+R(k,m,14)*fluju*fUm
cc	  fluxUi=fluxUi+R(k,m,14)*fluju*fUi
cc	  fluxU=fluxU+R(k,m,14)*fluju
cc	  fluUms=fluUms+R(k,m,14)*fluju*fusdss*fUms
cc	  fluUis=fluUis+R(k,m,14)*fluju*fusdss*fUis
cc	  fluUis004=fluUis004+R(k,m,14)*fluju*fusdss*fUis004
cc	  fluUs=fluUs+R(k,m,14)*fluju*fusdss
cc       ENDIF
       qum=qum+R(k,m,14)*flujv*(volp0(1)/volp0(2))
       qu10=qu10+R(k,m,14)*flujv*(starmin/volp0(3))
       qu05=qu05+R(k,m,14)*flujv*(starmin/volp0(4))
       qu15=qu15+R(k,m,14)*flujv*(starmin/volp0(5))
      ENDIF
      ENDDO
      RETURN
      END
