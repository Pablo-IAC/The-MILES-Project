c Subrutina que calcula el espectros de alta resolucion
c integrados. Hace uso de las subrutinas "lbusc" y "sigmac",
c las cuales calculan los espectros estelares rangos Jones, 
c CaT y MILES para una terna de parametros atmosfericos. 
      SUBROUTINE hrsl_MIUSC(k,wfluji,al99i,fja99i,qum,qu10,qu05,qu15,
     &fluxUm,fluxUi,fluxU,fluUms,fluUis,fluUs,fluUis004)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c MIUSC
      PARAMETER (npxi=6672) !MIUSC
      DIMENSION ssppi(npxi,2),wfluji(npxi,2) !MIUSC
      DIMENSION ali(3),al99i(3)
      COMMON/qsighi/xnbxsi,volpai,voli10,voli05,voli15
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
      qum=0.0d0 !parametro calidad (MIUSC) normalizado sigma_min
      qu10=0.0d0 !parametro calidad (MIUSC) aceptable 1.0 sigma
      qu05=0.0d0 !parametro calidad (MIUSC) aceptable 0.5 sigma
      qu15=0.0d0 !parametro calidad (MIUSC) aceptable 1.5 sigma
      starmin=3.0d0 !minimo aceptable estrellas: 3 (1 x subcaja)
      fja99i=0.0d0
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
	 al99i(l)=0.d0 !MIUSC
      enddo
      do l=1,npxi !MIUSC
	 wfluji(l,2)=0.d0
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
        call fracat(tiss,giss,frc)
        call fraU(tiss,giss,fUm,fUi)
        call fra_u_usdss(tiss,giss,fusdss)
        call frausdss(tiss,giss,fUms,fUis)
        call frausdss004(tiss,giss,fUis004)
	call sigmam_MIUSC(tiss,giss,fiss,ssppi,ali) !MIUSC
      fvmil=0.0d0
	do l=1,npxi
	 rfV=0.0d0
	 if(l.eq.npxi)then
	  dl=abs(ssppi(l-1,1)-ssppi(l,1))
	 else
	  dl=abs(ssppi(l+1,1)-ssppi(l,1))
	 endif
	 if(ssppi(l,1).ge.V1.and.ssppi(l,1).le.V2)then
	  call respv(ssppi(l,1),rfV)
        fvmil=fvmil+rfV*ssppi(l,2)*dl
	 endif
	enddo
	DO l=1,npxi !MIUSC
	 ssppi(l,2)=ssppi(l,2)/fvmil
         if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    wfluji(l,2)=wfluji(l,2)+ssppi(l,2)*R140(k,mma)*flujv
	   enddo
	  endif
	  wfluji(l,2)=wfluji(l,2)+ssppi(l,2)*R(k,1,14)*flujv
	 else
	  wfluji(l,2)=wfluji(l,2)+ssppi(l,2)*R(k,m,14)*flujv
	 endif
	ENDDO
      if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    fja99i=fja99i+R140(k,mma)*flujv
	   enddo
	  endif
	  fja99i=fja99i+R(k,1,14)*flujv
	else
	  fja99i=fja99i+R(k,m,14)*flujv
	endif
	do l=1,3 !MIUSC
       if(m.eq.1)then
	  if(imam(k).gt.0)then
	   do mma=1,imam(k)
	    al99i(l)=al99i(l)+ali(l)*R140(k,mma)*flujv
	   enddo
	  endif
	  al99i(l)=al99i(l)+ali(l)*R(k,1,14)*flujv
	 else
	  al99i(l)=al99i(l)+ali(l)*R(k,m,14)*flujv
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
	qum=qum+R(k,m,14)*flujv*(xnbxsi/volpai)
	qu10=qu10+R(k,m,14)*flujv*(starmin/voli10)
	qu05=qu05+R(k,m,14)*flujv*(starmin/voli05)
	qu15=qu15+R(k,m,14)*flujv*(starmin/voli15)
      ENDDO
      RETURN
      END
