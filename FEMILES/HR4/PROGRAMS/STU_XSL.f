      SUBROUTINE STU_XSL(m,ageSSP,BETA,reale,ZMISO,xm9,xm9s,a9m,f9m,
     &qum,qu10,qu05,qu15)
cc     SUBROUTINE STU_XSL(m,ageSSP,BETA,reale,ZMISO,xm9,xm9s,a9m,f9m,
cc     &qum,qu10,qu05,qu15,fxUm,fxUi,fxU,fxUms,fxUis,fxUs,fxUis004)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c XSL
      PARAMETER (npxi=106301) !IRTF #pixels in the spectra
      DIMENSION a9m(3),xm9(npxi,2),xm9s(npxi,2)
      DIMENSION z(12,150000,16)
      CHARACTER*1 shape
      COMMON/lmu0/bicl,bicp,bich
      COMMON/lmu1/ZMU
      COMMON/lshape/shape
      COMMON/lson/sson
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      COMMON/stus/ZMASA,tiso(75),num(15),numtis
      COMMON/zsolo/z
      COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/KEYS/keypT(18)
      COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
      COMMON/NC/ncha
      COMMON/masas/xmast,xmaso,xmasa,fmO,fmA,fmR !common subrut e.f
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk
      COMMON/MAMA/R140(12,99000),imam(12)
      fmO=0.0d0 
      fmA=0.0d0
      fmR=0.0d0
      xmast=0.0d0
      xmaso=0.0d0
      xmasa=0.0d0
	barm=0.0d0 !esta es la masa promedio segun la IMF
	BETA=0.d0
	zmax=0.d0
	ZMISO=0.d0
	ZMASOS=ZMASA
	if(ncha.le.0)then
	 edx=ageSSP/1000.0d0 !pues en "a.f" habia multiplicado x 1000
	else
	 edx=dlog10(ageSSP*1.e6)
	endif
	i=0
	call loca(tiso,numtis,edx,n)
	if(ncha.le.0)then
	 reale=tiso(n)
	else
	 reale=10.**tiso(n)/1.e9
	endif
	kkeypT=0 !solo para el caso que ZML no es tan pequeña para TERAMO
	DO 1502 kk=1,num(m) !puede subir hasta 150000
      IF(z(m,kk,1).gt.tiso(n)-.0001d0.and.z(m,kk,1).lt.tiso(n)+
     &.0001d0)THEN
       if(ncha.le.0)then
        kkeypT=kkeypT+1
	 endif
	if((z(m,kk,13).ge.ZML).and.(z(m,kk,13).le.ZMUS))then
       i=i+1
c ISOC TERAMO: SI MASA MINIMA NO TAN < =>keypT(se modifican)
c       if(i.eq.1.and.ncha.le.0)then
c!actualizo los keypT 
c        do iteram=1,18
c         keypT(iteram)=keypT(iteram)-kkeypT-1!-1 porque pto.1 incluido
c        enddo
         kkeypT=0
c        endif
	  jotai(m)=i
	  Rm0(m,i,1)=z(m,kk,13)!Masa_inicial para fases evolutivas PADOVA
c	  R(m,i,15)=z(m,kk,16)!Mbol
	  R(m,i,15)=10.**(-0.4d0*(z(m,kk,16)-bolsol)) !Mbol to L/Lbolsol
	  R(m,i,1)=z(m,kk,15)
	  R(m,i,2)=(-2.5d0)*dlog10(z(m,kk,6)/z(m,kk,7))
	  do jk=3,12
	   R(m,i,jk)=z(m,kk,jk)
	  enddo
	  if(z(m,kk+1,1).ne.tiso(n))then
	   dm=abs(z(m,kk,13)-z(m,kk-1,13))
	  elseif(z(m,kk-1,1).ne.tiso(n))then
	   dm=abs(z(m,kk+1,13)-z(m,kk,13))
	  else
	   dm=abs(z(m,kk,13)-z(m,kk-1,13))
	  endif
	  call limf(z(m,kk,13),fim)
c	  zmax tomara al final la mayor masa que sobrevive
	  zmax=z(m,kk,13)+dm !integra el trozo dm a partir de zmax
c	  estoy integrando en masa para obtener la BETA con las masas vivas
	  BETA=BETA+fim*dm
	  barm=barm+fim*z(m,kk,13)*dm
C	  PARA SUMAR POBLAC. EN MASA ESTE ES EL MOMENTO: ZMASOS,
	  R(m,i,14)=fim*ZMASOS*dm
	  ZMISO=ZMISO+R(m,i,14)*z(m,kk,2)
	  zrem=z(m,kk,2)
c-----------------criterio de clase de luminosidad-------------------------
	  IF(ncha.le.0)THEN
	   if(i.le.keypT(3))then
	    R(m,i,13)=5.0d0
	   else
	    R(m,i,13)=3.0d0
	   endif
	  ELSE
	   if(R(m,i,3).lt.dlog10(5500.d0))then
	    if(R(m,i,4).ge.3.50d0)then
		R(m,i,13)=5.0d0
	    elseif(R(m,i,4).lt.3.5d0)then
		R(m,i,13)=3.0d0
	    endif
	   elseif(R(m,i,3).ge.dlog10(5500.d0))then
		R(m,i,13)=5.0d0
	   else
		R(m,i,13)=9.0d0
	   endif
	  ENDIF
c       para la zona del CaT segun I_Johnson
	  vikc1=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,9))
	  vrkc1=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,8))
	  biss1=-2.5d0*dlog10(z(m,kk,6))
	  viss1=-2.5d0*dlog10(z(m,kk,7))
	  uviss1=-2.5d0*dlog10(z(m,kk,5)/z(m,kk,7))
	  vk1=-2.5d0*dlog10(z(m,kk,7)/z(m,kk,12))
        call xflabs(uviss1,biss1,viss1,vikc1,vk1,flujb1,
     &flujv1,flujj1,fluju1,flujk1)
	 endif
	ENDIF
1502	ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      imam(m)=0
c masas mas pequenyas que la primera masa luminosa
	IF(ZMLow.lt.ZML)THEN
	 dmw=0.000001d0
c	 dmw=0.1d0 !PARA ACELERAR CALCULO
       nimkk=int(abs((ZML-ZMLow)/dmw))+1
	 zmin=ZMLow-dmw
	 imami=0
	 do iil=1,nimkk
	  zmin=zmin+dmw
	  if(iil.eq.nimkk) dmw=abs(ZML-zmin)
        call limf(zmin,fim)
	  if(zmin.gt.xminlu.and.zmin.lt.ZML)then
         imami=imami+1
	   imam(m)=imami
	   R140(m,imami)=fim*ZMASOS*dmw
	   ZMISO=ZMISO+R140(m,imami)*zmin
c	   do lnk=5,12
c          ZZZ(lnk-4)=ZZZ(lnk-4)+R140(m,imami)*R(m,1,lnk)
c	   enddo
        endif
	  BETA=BETA+fim*dmw
	  barm=barm+fim*zmin*dmw
	 enddo
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ZMISA=ZMISO !masa brillante sin incluir los remnants
c masas que ya han muerto hasta la max. inicial (ZMUS)
	IF(ZMUS.gt.zmax)THEN
	 dmw=0.01d0
c zmax ultima masa inicial isocrona + dm    !atencion a dm
       nimkk=int(abs((ZMUS-zmax)/dmw))+1 !int implica entero menor,+1 para resto
	 zmax=zmax-dmw
	 do iil=1,nimkk
	  zmax=zmax+dmw
	  if(iil.eq.nimkk) dmw=abs(ZMUS-zmax) !integracion resto hasta ZMUS
	  call limf(zmax,fim)
c Sumo solo los remnants para las masas mayores (receta: Renzini & Ciotti1993)
	  call remn(m,zmax,zrem)
	  ZMISO=ZMISO+(fim*ZMASOS*dmw)*zrem
	  BETA=BETA+fim*dmw
	  barm=barm+fim*zmax*dmw
	 enddo
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	BETA=1.0d0/BETA
	barm=1.0d0/barm
	BETA=barm
	ZMISO=ZMISO*BETA
	ZMISA=ZMISA*BETA
	fmO=ZMISO/ZMASOS !fraccion masa galaxia en * con remnants 
	fmA=ZMISA/ZMASOS !fraccion masa galaxia en * sin remnants
	fmR=fmO-fmA !fraccion remnants
	xmast=ZMASOS
	xmaso=ZMISO
	xmasa=ZMISA
cc      CALL hrsl_XSL(m,xm9,xm9s,a9m,f9m,qum,qu10,qu05,qu15,
cc     &fxUm,fxUi,fxU,fxUms,fxUis,fxUs,fxUis004)
      CALL hrsl_XSL(m,xm9,xm9s,a9m,f9m,qum,qu10,qu05,qu15)
      RETURN
      END
