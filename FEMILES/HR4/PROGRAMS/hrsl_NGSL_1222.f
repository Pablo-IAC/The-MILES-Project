c Subrutina que calcula el espectros de alta resolucion
c integrados. Hace uso de las subrutinas "lbusc" y "sigmac",
c las cuales calculan los espectros estelares rangos Jones, 
c CaT y MILES para una terna de parametros atmosfericos. 
      SUBROUTINE hrsl_NGSL(k,wfluju,wflujus,al99u,fja99u,qum,qu10,qu05
     &,qu15)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c NGSL
      PARAMETER (npxu=6192) !NGSL
      DIMENSION ssppu(npxu,2),wfluju(npxu,2),wflujus(npxu,2) !NGSL
      DIMENSION alu(3),al99u(3)
      COMMON/qsighu/xnbxsu,volpau,volu10,volu05,volu15
c
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk
      COMMON/MAMA/R140(12,99000),imam(12)
      COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
	COMMON/part0/ipartial
	COMMON/part1/pmassL,pmassH
      COMMON/ab/iaba
	COMMON/enhanc/amgfe,DeltFeH ![Mg/Fe],[Fe/H] correc. value
c limites lambdas
      V1=4750.d0
      V2=7400.d0
ccccccccccccccccccccccccc
c [Fe/H] (si [Mg/Fe].ne.0) a partir de metalicidad isocrona [Z/H]=fez(k)
c solo para sigmam.f
      fiss=fez(k)
      fiss=fez(k)-DeltFeH
ccccccccccccccccccccccccc
      qum=0.0d0 !parametro calidad (NGSL) normalizado sigma_min
      qu10=0.0d0 !parametro calidad (NGSL) aceptable 1.0 sigma
      qu05=0.0d0 !parametro calidad (NGSL) aceptable 0.5 sigma
      qu15=0.0d0 !parametro calidad (NGSL) aceptable 1.5 sigma
      starmin=3.0d0 !minimo aceptable estrellas: 3 (1 x subcaja)
      fja99u=0.0d0
      do l=1,3
	 al99u(l)=0.d0 !NGSL
      enddo
      do l=1,npxu !NGSL
	 wfluju(l,2)=0.d0
	 wflujus(l,2)=0.d0
      enddo
      DO m=1,jotai(k)
       tiss=10.**R(k,m,3)
       giss=R(k,m,4)
	 biss=(-2.5d0)*dlog10(R(k,m,6))
	 vk=(-2.5d0)*dlog10(R(k,m,7)/R(k,m,12))
	 viss=(-2.5d0)*dlog10(R(k,m,7))
	 vikc=(-2.5d0)*dlog10(R(k,m,7)/R(k,m,9))
	 uvbus=(-2.5d0)*dlog10(R(k,m,5)/R(k,m,7))
       call xflabs(uvbus,biss,viss,vikc,vk,flujb,flujv,flujcj,fluju,
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
	call sigmam_NGSL(tiss,giss,fiss,ssppu,alu) !NGSL
      fvmil=0.0d0
	do l=1,npxu
	 rfV=0.0d0
	 if(l.eq.npxu)then
	  dl=abs(ssppu(l-1,1)-ssppu(l,1))
	 else
	  dl=abs(ssppu(l+1,1)-ssppu(l,1))
	 endif
	 if(ssppu(l,1).ge.V1.and.ssppu(l,1).le.V2)then
	  call respv(ssppu(l,1),rfV)
        fvmil=fvmil+rfV*ssppu(l,2)*dl
	 endif
	enddo
	DO l=1,npxu !NGSL
	 ssppu(l,2)=ssppu(l,2)/fvmil
	 spfl=ssppu(l,2)*flujv
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wfluju(l,2)=wfluju(l,2)+R140(k,mma)*spfl
	    wflujus(l,2)=wflujus(l,2)+R140(k,mma)*spfl*spfl
	   enddo
	  endif
	  wfluju(l,2)=wfluju(l,2)+R(k,1,14)*spfl
	  wflujus(l,2)=wflujus(l,2)+R(k,1,14)*spfl*spfl
	 else
	  wfluju(l,2)=wfluju(l,2)+R(k,m,14)*spfl
	  wflujus(l,2)=wflujus(l,2)+R(k,m,14)*spfl*spfl
	 endif
	ENDDO
      if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    fja99u=fja99u+R140(k,mma)*flujv
	   enddo
	  endif
	  fja99u=fja99u+R(k,1,14)*flujv
	else
	  fja99u=fja99u+R(k,m,14)*flujv
	endif
	do l=1,3 !NGSL
       if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    al99u(l)=al99u(l)+alu(l)*R140(k,mma)*flujv
	   enddo
	  endif
	  al99u(l)=al99u(l)+alu(l)*R(k,1,14)*flujv
	 else
	  al99u(l)=al99u(l)+alu(l)*R(k,m,14)*flujv
	 endif
	enddo
	qum=qum+R(k,m,14)*fluju*(xnbxsu/volpau) !fluju, pues solo usamos UV
	qu10=qu10+R(k,m,14)*fluju*(starmin/volu10) !fluju, pues solo usamos UV
	qu05=qu05+R(k,m,14)*fluju*(starmin/volu05) !fluju, pues solo usamos UV
	qu15=qu15+R(k,m,14)*fluju*(starmin/volu15) !fluju, pues solo usamos UV
      ENDDO
c      DO l=1,npxu !SBF
c       wflujus(l,2)=wflujus(l,2)/wfluju(l,2)
c      ENDDO
      RETURN
      END
