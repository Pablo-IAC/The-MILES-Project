c Subrutina que calcula el espectros de alta resolucion
c integrados. Hace uso de las subrutinas "lbusc" y "sigmac",
c las cuales calculan los espectros estelares rangos Jones, 
c CaT y MILES para una terna de parametros atmosfericos. 
c*******SUBRUTINA HRSL************************************************
c      SUBROUTINE hrsl(npxm,nstm,k,wflujb,wflujv,wflujc,wflujm,al99,
c     & fjal99,al99c,fja99c,al99m,fja99m,Wfcatt)
      SUBROUTINE hrsl(npxm,k,wflujb,wflujv,wflujc,wflujm,al99,fjal99,
     &al99c,fja99c,al99m,fja99m,Wfcatt,qum,qu10,qu05,qu15,
     &fluxUm,fluxUi,fluxU,fluUms,fluUis,fluUs,fluUis004)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (npxj=1107) !jones
      PARAMETER (npxc=710) !cat
      DIMENSION ssppb(1120,2),wflujb(1120,2) !jones b
      DIMENSION ssppv(1120,2),wflujv(1120,2) !jones v
      DIMENSION ssppc(710,2),wflujc(710,2) !cat
      DIMENSION ssppm(npxm,2),wflujm(npxm,2) !miles
      DIMENSION al(3),al99(3),alc(3),al99c(3),alm(3),al99m(3)
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      common/cfnc/fnc11(2999),fnmg(2999)
      common/cd4000/fra4(2999),fbj(2999)
      COMMON/qsighr/xnbxst,volpam,volp10,volp05,volp15
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu
      COMMON/MAMA/R140(12,99000),imam(12)
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
c
      fiss=fez(k)
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
      fluxUm=0.0d0 !Missing U MILES flux
      fluxUi=0.0d0 !Missing U MIUSC flux
      fluxU=0.0d0 !Total U flux
      fluUms=0.0d0 !Missing u sdss MILES flux
      fluUis=0.0d0 !Missing u sdss MIUSC flux
      fluUis004=0.0d0 !Missing u sdss MIUSC flux at z=0.04
      fluUs=0.0d0 !Total u sdss flux
c
      do l=1,3
	al99(l)=0.d0 !jones
	al99c(l)=0.d0 !cat
	al99m(l)=0.d0 !miles
      enddo
      do l=1,npxj !jones
	wflujb(l,2)=0.d0
	wflujv(l,2)=0.d0
      enddo
      do l=1,npxc !cat
	wflujc(l,2)=0.d0
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
c	call cotojo(giss,vrkc,vikc,vij)
c	call xflabs(biss,viss,vij,flujb,flujv,flujcj)
	call xflabs(uvbus,biss,viss,vikc,flujb,flujv,flujcj,fluju)
        call fracat(tiss,giss,frc)
        call fraU(tiss,giss,fUm,fUi)
        call fra_u_usdss(tiss,giss,fusdss)
        call frausdss(tiss,giss,fUms,fUis)
        call frausdss004(tiss,giss,fUis004)
c	flujc=(frc*flujcj)/(C2-C1)
	flujc=frc*flujcj
c	write(36,*)tiss,vikc,vij,flujb,flujv,flujcj,frc
c	write(*,*)tiss,vikc,vij,flujb,flujv,flujcj,frc
c        write(*,'(6F10.4)')tiss,giss,frc,flujv,flujcj,flujc
c        write(36,'(6F10.4)')tiss,giss,frc,flujv,flujcj,flujc
	fjal=(flujv/DV-flujb/DB)*(224.d0/993.d0)+flujb/DB 
	call sigmam(tiss,giss,fiss,ssppm,alm)
	call lbusc(tiss,giss,fiss,viss,ssppb,ssppv,al,fra4000)
	fra4(m)=fra4000
	fbj(m)=flujb
	call sigmac(tiss,giss,fiss,ssppc,alc,Wcatt,fnc1,fnmgi)
	fnc11(m)=fnc1
	fnmg(m)=fnmgi
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
	call respb(ssppb(l,1),rfB) !respob factor entre 0 y 1
	call respv(ssppv(l,1),rfV) !respob factor entre 0 y 1
	if(ssppb(l,1).ge.BJ1.and.ssppb(l,1).le.BJ2)then
	  fbjon=fbjon+rfB*ssppb(l,2)*dlb
	endif
	  fvjon=fvjon+rfV*ssppv(l,2)*dlv
       enddo
       do l=1,npxj
	ssppb(l,2)=ssppb(l,2)/(fbjon*fbuser)
	ssppv(l,2)=ssppv(l,2)/fvjon
       enddo
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
	  fcat=fcat+rfC*ssppc(l,2)*dl
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
	  if(ssppm(l,1).ge.BJ1.and.ssppm(l,1).le.BJ2)then
	   fbjon=fbjon+rfB*ssppm(l,2)*dl
	  endif
	 endif
	 if(ssppm(l,1).ge.V1.and.ssppm(l,1).le.V2)then
	  call respv(ssppm(l,1),rfV)
	  if(ssppm(l,1).ge.VJ1.and.ssppm(l,1).le.VJ2)then
	   fvjon=fvjon+rfV*ssppm(l,2)*dl
	  endif
	 endif
         fbmil=fbmil+rfB*ssppm(l,2)*dl
         fvmil=fvmil+rfV*ssppm(l,2)*dl
	enddo
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
	DO l=1,npxc !cat
c	 wflujc(l,2)=wflujc(l,2)+ssppc(l,2)*R(k,m,14)*flujc*FNORMC
         ssppc(l,2)=(ssppc(l,2)/fcat)
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wflujc(l,2)=wflujc(l,2)+ssppc(l,2)*R140(k,mma)*flujc
	   enddo
	  endif
	  wflujc(l,2)=wflujc(l,2)+ssppc(l,2)*R(k,1,14)*flujc
	 else
	  wflujc(l,2)=wflujc(l,2)+ssppc(l,2)*R(k,m,14)*flujc
	 endif
	ENDDO
	DO l=1,npxm !miles
c	 ssppm(l,2)=(ssppm(l,2)/fvmil)*FVNOR
c	 wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,m,14)*flujv*FNORMM
c	 wflujm(l,2)=wflujm(l,2)+ssppm(l,2)*R(k,m,14)*R(k,m,7)
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
	    Wfcatt=Wfcatt+Wcatt*R140(k,mma)*flujc
	    fjal99=fjal99+R140(k,mma)*fjal
	    fja99c=fja99c+R140(k,mma)*flujc
	    fja99m=fja99m+R140(k,mma)*flujv
	   enddo
	  endif
	  Wfcatt=Wfcatt+Wcatt*R(k,1,14)*flujc
	  fjal99=fjal99+R(k,1,14)*fjal
	  fja99c=fja99c+R(k,1,14)*flujc
	  fja99m=fja99m+R(k,1,14)*flujv
	else
	  Wfcatt=Wfcatt+Wcatt*R(k,m,14)*flujc
	  fjal99=fjal99+R(k,m,14)*fjal
	  fja99c=fja99c+R(k,m,14)*flujc
	  fja99m=fja99m+R(k,m,14)*flujv
	endif
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
	do l=1,3 !cat
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
	enddo
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
