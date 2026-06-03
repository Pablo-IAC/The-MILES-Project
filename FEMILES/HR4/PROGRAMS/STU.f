c**********************************************************
c SUBRUTINA QUE CALCULA LA TODOS LOS OBSERVABLES DE UNA SSP
c**********************************************************
      SUBROUTINE STU(npxm,m,ageSSP,BETA,reale,ZZZ,ZMISO,vsin,stm,b9
     &  ,v9,c9,xm9,a9,f9,a9c,f9c,a9m,f9m,ZISOL,fZISOL,VM,dwgi,ZZZsbf
     &  ,fhst,fhsbf,t22,g22,f22,Wfcatt,h,fluxx0,qum,qu10,qu05,qu15
     &  ,fxUm,fxUi,fxU,fxUms,fxUis,fxUs,fxUis004)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ZZZ(8),dwgi(2,8),vsin(50),stm(50),a9(3),VM(8),
     &  ZISOL(8),ZZZsbf(50),fhst(4),fhsbf(4),a9c(3),xLsol(8),xmsol(8)
     &  ,fZISOL(8)
      DIMENSION b9(1120,2),v9(1120,2),c9(710,2),xm9(npxm,2)
      DIMENSION z(12,150000,16)
      DIMENSION h(5,9)
      CHARACTER*4 gc
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      COMMON/stus/ZMASA,tiso(75),num(15),numtis
      COMMON/zsolo/z
      COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/rat/ra(3,50)
      COMMON/KEYS/keypT(18)
      COMMON/TOAGB/g(15,2000,7),ng(15)
      COMMON/TOAGBC/gc(15,2000)
      COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
      COMMON/NC/ncha
      COMMON/masas/xmast,xmaso,xmasa,fmO,fmA,fmR !common subrut e.f
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu
      COMMON/MAMA/R140(12,99000),imam(12)
c      COMMON/FNOR/C2Code,FNORMM,FNORMC
      DATA xmsol/5.6d0,5.441d0,4.82d0,4.459d0,4.148d0,
     &  3.711d0,3.392d0,3.334d0/
      xminlu=0.09d0!minima masa luminosa (solo para computacion M/L)
      do i=1,8
       xLsol(i)=10.**((-0.4d0)*xmsol(i))
      enddo
      fmO=0.0d0 
      fmA=0.0d0
      fmR=0.0d0
      xmast=0.0d0
      xmaso=0.0d0
      xmasa=0.0d0
      do nmnm=1,4
       fhst(nmnm)=0.d0
       fhsbf(nmnm)=0.d0
      enddo
      do nmnm=1,8
       ZZZ(nmnm)=0.d0
       ZZZsbf(nmnm)=0.d0
       do nnnn=1,2
        dwgi(nnnn,nmnm)=0.d0
       enddo
      enddo
      do l0=1,5
       do nnnn=1,2
        gsi(l0,nnnn)=0.d0
       enddo
      enddo
c	fluxx0=0.0d0 !lo he puesto un poco mas abajo
	barm=0.0d0 !esta es la masa promedio segun la IMF
	BETA=0.d0
	zmax=0.d0
	ZMISO=0.d0
	ZMASOS=ZMASA
c	if(ncha.eq.2)then
c	 fez(m)=dlog10(Z00(m))-dlog10(0.020d0)
c	else	
c	 fez(m)=dlog10(Z00(m))-dlog10(0.019d0)
c	endif	
	if(ncha.eq.0)then
	 edx=ageSSP/1000.0d0 !pues en "a.f" habia multiplicado x 1000
	else
	 edx=dlog10(ageSSP*1.e6)
	endif
	i=0
	call loca(tiso,numtis,edx,n)
	if(ncha.eq.0)then
	 reale=tiso(n)
	else
	 reale=10.**tiso(n)/1.e9
	endif
	IF(ncha.gt.0)THEN
	do nmnm=1,ng(m)
       if(g(m,nmnm,1).lt.edx+.001d0.and.g(m,nmnm,1).gt.edx-.001d0)then
	  if(gc(m,nmnm).eq.'TOff')then
	   gsi(1,1)=g(m,nmnm,2)                        !masa inicial
	   gsi(1,2)=g(m,nmnm,5)                        !log Teff
	  elseif(gc(m,nmnm).eq.'RGBb')then
	   gsi(2,1)=g(m,nmnm,2)
	   gsi(2,2)=g(m,nmnm,5)                         
	  elseif(gc(m,nmnm).eq.'RGBt')then
	   gsi(3,1)=g(m,nmnm,2)
	   gsi(3,2)=g(m,nmnm,5)                         
	  elseif(gc(m,nmnm).eq.'BHeb')then
	   gsi(4,1)=g(m,nmnm,2)
	   gsi(4,2)=g(m,nmnm,5)                         
	  elseif(gc(m,nmnm).eq.'EHeb')then
	   gsi(5,1)=g(m,nmnm,2)
	   gsi(5,2)=g(m,nmnm,5)
	   goto 33551                         
	  endif
       endif
	enddo
33551   continue
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccc
	fluxx0=0.d0
	do l0=1,5
	  do nnnn=1,9
	   h(l0,nnnn)=0.d0
	  enddo
	enddo
cccccccccccccccccccccccccccccccccccccccccccc
	kkeypT=0!solo para el caso que ZML no es tan pequeña para TERAMO
	DO 1502 kk=1,num(m) !puede subir hasta 150000
c	IF (z(m,kk,1).eq.tiso(n)) THENRm0
      IF(z(m,kk,1).gt.tiso(n)-0.0001d0.and.z(m,kk,1).lt.tiso(n)+.0001d0)THEN
        if(ncha.eq.0)then
		kkeypT=kkeypT+1
	endif
	if((z(m,kk,13).ge.ZML).and.(z(m,kk,13).le.ZMUS))then
		i=i+1
c ISOC TERAMO: SI MASA MINIMA NO TAN < =>keypT(se modifican)
		if(i.eq.1.and.ncha.eq.0)then
!actualizo los keypT 
		  do iteram=1,18
		   keypT(iteram)=keypT(iteram)-kkeypT-1!-1 porque pto.1 incluido
		  enddo
		  kkeypT=0
		endif
		jotai(m)=i
		Rm0(m,i,1)=z(m,kk,13)!Masa_inicial para fases evolutivas PADOVA
		R(m,i,15)=z(m,kk,16)!Mbol
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
c	zmax tomara al final la mayor masa que sobrevive
		zmax=z(m,kk,13)+dm !integra el trozo dm a partir de zmax
c	estoy integrando en masa para obtener la BETA con las masas vivas
		BETA=BETA+fim*dm
c		write(44,*)BETA
		barm=barm+fim*z(m,kk,13)*dm
C	PARA SUMAR POBLAC. EN MASA ESTE ES EL MOMENTO: ZMASOS,
C	ES DECIR ZMASOS1+ZMASOS2=ZMASOS
c		R(m,i,14)=fim*ZMASOS*dm/z(m,kk,13)
		R(m,i,14)=fim*ZMASOS*dm
		ZMISO=ZMISO+R(m,i,14)*z(m,kk,2)
c		write(80,*)z(m,kk,13),z(m,kk,2),R(m,i,14)*z(m,kk,2),ZMISO
c		write(*,*)ZMISO,z(m,kk,2),R(m,i,14)*z(m,kk,2)
		zrem=z(m,kk,2)
c-----------------criterio de clase de luminosidad-------------------------
		IF(ncha.eq.0)THEN
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
c-----------------criterio de clase de luminosidad-------------------------
c		if((R(m,i,3).lt.4.63).and.(R(m,i,1).gt.-1.0))then
c			if(R(m,i,4).ge.4.0)then
c				R(m,i,13)=5.0
c			elseif(R(m,i,4).le.3.5)then
c				R(m,i,13)=3.0
c			elseif((R(m,i,4).gt.3.5).and.(R(m,i,4).lt.4.0))then
c				if(R(m,i,1).le.(2.0*R(m,i,4)-6.0))then
c					R(m,i,13)=5.0
c				elseif(R(m,i,1).gt.(2.0*R(m,i,4)-6.0))then
c					R(m,i,13)=3.0
c				endif
c			endif
c		else
c			R(m,i,13)=9.0
c		endif
c--------------------------------------------------------------------------
c       para la zona del CaT segun I_Johnson
	  vikc1=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,9))
	  vrkc1=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,8))
	  biss1=-2.5d0*dlog10(z(m,kk,6))
	  viss1=-2.5d0*dlog10(z(m,kk,7))
	  uviss1=-2.5d0*dlog10(z(m,kk,5)/z(m,kk,7))
c          call cotojo(R(m,i,4),vrkc1,vikc1,vij1)
      call xflabs(uviss1,biss1,viss1,vikc1,flujb1,flujv1,flujj1,fluju1)
          call fracat(R(m,i,3),R(m,i,4),frx1)
	  fluxx1=(flujj1*frx1)/(8807.0d0-8475.0d0)
	 do lnk=5,12
c	para out_sbf
	  ZZZ(lnk-4)=ZZZ(lnk-4)+R(m,i,14)*R(m,i,lnk)
          ZZZsbf(lnk-4)=ZZZsbf(lnk-4)+R(m,i,14)*R(m,i,lnk)*R(m,i,lnk)
c	para out_contr
	  IF(ncha.eq.0)THEN
           if(i.le.keypT(3))then
            h(1,lnk-4)=h(1,lnk-4)+R(m,i,14)*R(m,i,lnk)     !MS 
	   elseif(i.gt.keypT(3).and.i.lt.keypT(5))then
            h(2,lnk-4)=h(2,lnk-4)+R(m,i,14)*R(m,i,lnk)     !SGB
	   elseif(i.ge.keypT(5).and.i.lt.keypT(9))then
            h(3,lnk-4)=h(3,lnk-4)+R(m,i,14)*R(m,i,lnk)    !RGB
	   elseif(i.ge.keypT(9).and.i.le.keypT(15))then
            h(4,lnk-4)=h(4,lnk-4)+R(m,i,14)*R(m,i,lnk)    !CHeb
	   else
            h(5,lnk-4)=h(5,lnk-4)+R(m,i,14)*R(m,i,lnk)     !AGB or Cb
           endif
	  ELSE
           if(z(m,kk,13).le.gsi(1,1))then
            h(1,lnk-4)=h(1,lnk-4)+R(m,i,14)*R(m,i,lnk)     !MS 
	   elseif(z(m,kk,13).gt.gsi(1,1).and.z(m,kk,13).lt.gsi(2,1))then
            h(2,lnk-4)=h(2,lnk-4)+R(m,i,14)*R(m,i,lnk)     !SGB
	   elseif(z(m,kk,13).ge.gsi(2,1).and.z(m,kk,13).le.gsi(3,1))then
	    if(R(m,i,3).lt.gsi(4,2))then  !eliminamos problema redondeo
             h(3,lnk-4)=h(3,lnk-4)+R(m,i,14)*R(m,i,lnk)    !RGB
            else
             h(4,lnk-4)=h(4,lnk-4)+R(m,i,14)*R(m,i,lnk)    !CHeb
            endif
	   elseif(z(m,kk,13).gt.gsi(3,1).and.z(m,kk,13).le.gsi(5,1))then
	    if(R(m,i,3).gt.gsi(3,2))then  !eliminamos problema redondeo
             h(4,lnk-4)=h(4,lnk-4)+R(m,i,14)*R(m,i,lnk)    !CHeb
            else
             h(3,lnk-4)=h(3,lnk-4)+R(m,i,14)*R(m,i,lnk)    !RGB
            endif
	   else
            h(5,lnk-4)=h(5,lnk-4)+R(m,i,14)*R(m,i,lnk)     !AGB or Cb
           endif
	  ENDIF
c CALCULO CONTRIBUCIONES ESTADIOS EVOLUTIVOS EN ZONA CAT
	  IF(ncha.eq.0)THEN
           if(i.le.keypT(3))then
            h(1,9)=h(1,9)+R(m,i,14)*fluxx1     !MS 
	   elseif(i.gt.keypT(3).and.i.lt.keypT(5))then
            h(2,9)=h(2,9)+R(m,i,14)*fluxx1     !SGB
	   elseif(i.ge.keypT(5).and.i.lt.keypT(9))then
            h(3,9)=h(3,9)+R(m,i,14)*fluxx1     !RGB
	   elseif(i.ge.keypT(9).and.i.le.keypT(15))then
            h(4,9)=h(4,9)+R(m,i,14)*fluxx1     !CHeb
	   else
            h(5,9)=h(5,9)+R(m,i,14)*fluxx1     !AGB or Cb
           endif
	  ELSE
           if(z(m,kk,13).le.gsi(1,1))then
            h(1,9)=h(1,9)+R(m,i,14)*fluxx1     !MS 
	   elseif(z(m,kk,13).gt.gsi(1,1).and.z(m,kk,13).lt.gsi(2,1))then
            h(2,9)=h(2,9)+R(m,i,14)*fluxx1     !SGB
	   elseif(z(m,kk,13).ge.gsi(2,1).and.z(m,kk,13).le.gsi(3,1))then
	    if(R(m,i,3).lt.gsi(4,2))then  !eliminamos problema redondeo
             h(3,9)=h(3,9)+R(m,i,14)*fluxx1    !RGB
            else
             h(4,9)=h(4,9)+R(m,i,14)*fluxx1    !CHeb
            endif
	   elseif(z(m,kk,13).gt.gsi(3,1).and.z(m,kk,13).le.gsi(5,1))then
	    if(R(m,i,3).gt.gsi(3,2))then  !eliminamos problema redondeo
             h(4,9)=h(4,9)+R(m,i,14)*fluxx1    !CHeb
            else
             h(3,9)=h(3,9)+R(m,i,14)*fluxx1    !RGB
            endif
	   else
            h(5,9)=h(5,9)+R(m,i,14)*fluxx1     !AGB or Cb
           endif
	  ENDIF
          fluxx0=fluxx0+R(m,i,14)*fluxx1
c	para out_dg
	  IF(ncha.eq.0)THEN
           if(i.le.keypT(3))then
	    dwgi(1,lnk-4)=dwgi(1,lnk-4)+R(m,i,14)*R(m,i,lnk)
	   else
	    dwgi(2,lnk-4)=dwgi(2,lnk-4)+R(m,i,14)*R(m,i,lnk)
	   endif
	  ELSE
	   if(R(m,i,4).ge.3.50d0)then	          !enanas a grosso
	    dwgi(1,lnk-4)=dwgi(1,lnk-4)+R(m,i,14)*R(m,i,lnk)
	   elseif(R(m,i,4).lt.3.50d0)then            !gigantes a grosso
	    dwgi(2,lnk-4)=dwgi(2,lnk-4)+R(m,i,14)*R(m,i,lnk)
	   endif
	  ENDIF
	 enddo
	 bvhst=(-2.5d0)*dlog10(z(m,kk,6)/z(m,kk,7))
	 vrhst=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,8))
	 vihst=(-2.5d0)*dlog10(z(m,kk,7)/z(m,kk,9))
	 F439WB=(-0.003d0)*bvhst+0.088d0*bvhst*bvhst
	 if(vihst.gt.1.3087d0)then
  	  F555WV=+0.051d0*vihst-0.009d0*vihst*vihst-0.029d0
  	  F675WR=(-0.182d0)*vrhst+0.097d0*vrhst*vrhst-0.040d0
  	  F814WI=+0.124d0*vihst-0.028d0*vihst*vihst-0.076d0
         else
    	  F555WV=+0.052d0*vihst-0.027d0*vihst*vihst
    	  F675WR=(-0.253d0)*vrhst+0.125d0*vrhst*vrhst
    	  F814WI=+0.062d0*vihst-0.025d0*vihst*vihst
    	 endif
         F439WB=R(m,i,6)*10.**((-0.4d0)*F439WB)
         F555WV=R(m,i,7)*10.**((-0.4d0)*F555WV)
         F675WR=R(m,i,8)*10.**((-0.4d0)*F675WR)
         F814WI=R(m,i,9)*10.**((-0.4d0)*F814WI)
         fhst(1)=fhst(1)+R(m,i,14)*F439WB
         fhst(2)=fhst(2)+R(m,i,14)*F555WV
         fhst(3)=fhst(3)+R(m,i,14)*F675WR
         fhst(4)=fhst(4)+R(m,i,14)*F814WI
         fhsbf(1)=fhsbf(1)+R(m,i,14)*F439WB*F439WB
         fhsbf(2)=fhsbf(2)+R(m,i,14)*F555WV*F555WV
         fhsbf(3)=fhsbf(3)+R(m,i,14)*F675WR*F675WR
         fhsbf(4)=fhsbf(4)+R(m,i,14)*F814WI*F814WI
	endif
	ENDIF
1502	ENDDO
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c		write(44,*)'MASAS BAJAS'
        imam(m)=0
c masas mas pequenyas que la primera masa luminosa
	IF(ZMLow.lt.ZML)THEN
ccccccccccc para hrsl.f:
c         if(ZMLow.ge.xminlu.and.ZMLow.lt.ZML)then
c	  mami=1
c	  dmami=abs(z(m,2,13)-ZMlow)
c	  xmami=0.5d0*(ZMLow+z(m,2,13))
c	  call limf(xmami,fmami)
c	  R140(m)=fmami*ZMASOS*dmami
c	   do lnk=5,12
c      ZZZ(lnk-4)=ZZZ(lnk-4)+R140(m)*R(m,1,lnk)
c      ZZZsbf(lnk-4)=ZZZsbf(lnk-4)+R140(m)*R(m,1,lnk)*R(m,1,lnk)
c	   enddo
c	 else
c	  mami=0
c	  dmami=0.d0
c	  xmami=0.0d0
c	  fmami=0.0d0
c	  R140(m)=0.0d0
c	 endif
ccccccccccc
c	 dmw=0.00001d0
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
c	   ZMISO=ZMISO+(fim*ZMASOS*dmw)*zmin
           imami=imami+1
	   imam(m)=imami
	   R140(m,imami)=fim*ZMASOS*dmw
	   ZMISO=ZMISO+R140(m,imami)*zmin
	   do lnk=5,12
      ZZZ(lnk-4)=ZZZ(lnk-4)+R140(m,imami)*R(m,1,lnk)
      ZZZsbf(lnk-4)=ZZZsbf(lnk-4)+R140(m,imami)*R(m,1,lnk)*R(m,1,lnk)
	   enddo
          endif
	  BETA=BETA+fim*dmw
c		write(44,*)BETA
	  barm=barm+fim*zmin*dmw
	 enddo
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ZMISA=ZMISO !masa brillante sin incluir los remnants
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c		write(44,*)'MASAS ALTAS'
c masas que ya han muerto hasta la max. inicial (ZMUS)
	IF(ZMUS.gt.zmax)THEN
	 dmw=0.01d0
c	 dmw=0.0001d0
c	 dmw=0.5d0 !PARA ACELERAR CALCULO 
c zmax ultima masa inicial isocrona + dm    !atencion a dm
         nimkk=int(abs((ZMUS-zmax)/dmw))+1 !int implica entero menor,+1 para resto
c	 nimkk=int(abs(ZMUS-zmax)/dmw)
	 zmax=zmax-dmw
	 do iil=1,nimkk
	  zmax=zmax+dmw
	  if(iil.eq.nimkk) dmw=abs(ZMUS-zmax) !integracion resto hasta ZMUS
	  call limf(zmax,fim)
c Sumo solo los remnants para las masas mayores (receta: Renzini & Ciotti1993)
	  call remn(m,zmax,zrem)
	  ZMISO=ZMISO+(fim*ZMASOS*dmw)*zrem
c	write(80,*)zmax,zrem,(fim*ZMASOS*dmw)*zrem,ZMISO,'rem'
	  BETA=BETA+fim*dmw
c		write(44,*)BETA
	  barm=barm+fim*zmax*dmw
	 enddo
	ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	BETA=1.0d0/BETA
c	write(*,*)'BETA=',BETA
	barm=1.0d0/barm
	BETA=barm
c	write(*,*)'barm=',barm
c	barmf=BETA*barm !esto es para ver la masa promedio
cccc	if(shape.eq.'u'.or.shape.eq.'b')then
cccc	BETA=BETA/barm
cccc	endif
c10000	continue
	continue
	ZMISO=ZMISO*BETA
	ZMISA=ZMISA*BETA
	fmO=ZMISO/ZMASOS !fraccion masa galaxia en * con remnants 
	fmA=ZMISA/ZMASOS !fraccion masa galaxia en * sin remnants
	fmR=fmO-fmA !fraccion remnants
	xmast=ZMASOS
	xmaso=ZMISO
	xmasa=ZMISA
c      write(*,'(4(F10.4,1X))')ZMASOS,fmO,fmA,fmR
	do lnk=1,8
	   VM(lnk)=(-2.5d0)*dlog10(BETA*ZZZ(lnk))
c	   ZISOL(lnk)=ZMISO*xLsol(lnk)/(BETA*ZZZ(lnk))
	   fZISOL(lnk)=(BETA*ZZZ(lnk))/xLsol(lnk)
	   ZISOL(lnk)=ZMISO/fZISOL(lnk)
	enddo
c	write(*,*)'Flux',((BETA*ZZZ(lnk)),lnk=1,8)
	write(*,*)'V',((ZISOL(lnk)),lnk=1,8)
c LLAMADA A LAS SUBRUTINAS ESPECTROSCOPICAS
      CALL hrsl(npxm,m,b9,v9,c9,xm9,a9,f9,a9c,f9c,a9m,f9m,Wfcatt,
     &qum,qu10,qu05,qu15,fxUm,fxUi,fxU,fxUms,fxUis,fxUs,fxUis004)
c LLAMADA A LAS SUBRUTINAS DE LAS FF's
	CALL lcn1(m,vsin(8),stm(8))
	CALL lcn2(m,vsin(9),stm(9))
	CALL ca4227(m,vsin(10),stm(10))
	call lgband(m,vsin(11),stm(11))
	CALL fe4383(m,vsin(12),stm(12))
	CALL ca4455(m,vsin(13),stm(13))
	CALL fe4531(m,vsin(14),stm(14))
	CALL fe4668(m,vsin(15),stm(15))
	CALL hbeta(m,vsin(16),stm(16))
	CALL fe5015(m,vsin(17),stm(17))
	CALL mg1(m,vsin(18),stm(18))
	CALL mg2(m,vsin(19),stm(19))
	CALL mgb(m,vsin(20),stm(20))
	CALL fe5270(m,vsin(21),stm(21))
	CALL fe5335(m,vsin(22),stm(22))
	CALL fe5406(m,vsin(23),stm(23))
	CALL fe5709(m,vsin(24),stm(24))
	CALL fe5782(m,vsin(25),stm(25))
	CALL lnad(m,vsin(26),stm(26))
	CALL ltio1(m,vsin(27),stm(27))
	CALL ltio2(m,vsin(28),stm(28))
	CALL hdelta(m,vsin(29),stm(29))
	CALL hgamaa(m,vsin(30),stm(30))
	CALL hdeltf(m,vsin(31),stm(31))
	CALL hgamaf(m,vsin(32),stm(32))
	CALL d4000(m,vsin(33),stm(33))
c transfiero el coef corrector fnc1 (cont1) de subrout. hrsl
        CALL cattri(m,vsin(34),vsin(35),vsin(36),vsin(37),
     & vsin(38),stm(34),stm(35),stm(36),stm(37),stm(38),t22,g22,f22)
	CALL CO(m,vsin(39),stm(39))
	return
	end
