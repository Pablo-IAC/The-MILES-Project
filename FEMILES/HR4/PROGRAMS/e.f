c*******SUBRUTINA QUE ESCRIBE LOS RESULTADOS******************************
	SUBROUTINE e(alpha,izo,realeo,izy,realey)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	CHARACTER*1 shape
	CHARACTER*2 shapa
	DIMENSION vsin(50),ZZZ(8),ZISOL(8),fZISOL(8),SBF(50),vSBF(50)
     &  ,DWGI(2,8),a90(3),a90c(3),fhst(4),fhsbf(4),cfhst(4),cfhsbf(4)
        dimension tilin(9)
      COMMON/esc/vsin,ZZZ,SBF,DWGI,a90,a90c,ZISOL,fZISOL,ncols,nobs,
     & fhst,fhsbf
	COMMON/tgf/t222,g222,f222
	COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
	COMMON/lmu1/ZMU
	COMMON/lshape/shape
        COMMON/TOAGB2/hh(5,9)
        COMMON/CTFIN/CATFIN
c        COMMON/NC/ncha
	COMMON/fezsol/fez(15)
	COMMON/qua/quaqn,quaq10,quaq05,quaq15
        COMMON/masas/xmast,xmaso,xmasa,fmO,fmA,fmR !common subrut STU.f
	COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !common a.f
	if(shape.eq.'u')then
	 shapa='UN'
	elseif(shape.eq.'b')then
	 shapa='BI'
	elseif(shape.eq.'k')then
	 shapa='Ku'
	elseif(shape.eq.'r')then
	 shapa='Kb'
	endif
	if(alpha.eq.1.0)then
	  zetay=9.99d0
	  reay=9.99d0
	else
c	if(ncha.eq.2)then
c	  zetay=dlog10(Z00(izy)/0.020d0)
c	else	
c	  zetay=dlog10(Z00(izy)/0.019d0)
c	endif
	  zetay=fez(izy)	
	  reay=realey
	endif	
c	if(ncha.eq.2)then
c	  zetao=dlog10(Z00(izo)/0.020d0)
c	else	
c	  zetao=dlog10(Z00(izo)/0.019d0)
c	endif	
	zetao=fez(izo)
c	DO N3=1,npol
c		do J=1,8
c			AZZZ(J)=ZZZ(J)+AD(J,N3)
c		enddo
        vsin(1)=ZZZ(1)-ZZZ(3)
        vsin(2)=ZZZ(2)-ZZZ(3)
	vsin(3)=ZZZ(3)-ZZZ(4)	      
	vsin(4)=ZZZ(3)-ZZZ(5)
	vsin(5)=ZZZ(3)-ZZZ(6)	      
	vsin(6)=ZZZ(3)-ZZZ(7)
	vsin(7)=ZZZ(3)-ZZZ(8)
        vSBF(1)=SBF(1)-SBF(3)
        vSBF(2)=SBF(2)-SBF(3)
	vSBF(3)=SBF(3)-SBF(4)   
	vSBF(4)=SBF(3)-SBF(5)
	vSBF(5)=SBF(3)-SBF(6)   
	vSBF(6)=SBF(3)-SBF(7)
	vSBF(7)=SBF(3)-SBF(8)
	cfhst(1)=fhst(1)-fhst(2)
	cfhst(2)=fhst(2)-fhst(3)
	cfhst(3)=fhst(2)-fhst(4)
	cfhsbf(1)=fhsbf(1)-fhsbf(2)
	cfhsbf(2)=fhsbf(2)-fhsbf(3)
	cfhsbf(3)=fhsbf(2)-fhsbf(4)
c
	write(*,'(A2,1x,6(F7.4,1x))')
     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha
c     	write(*,'(f6.0,2(1x,f6.3))')(a90(nt),nt=1,3)
c     	write(*,'(f6.0,2(1x,f6.3))')(a90c(nt),nt=1,3)
c     	write(*,'(f6.0,2(1x,f6.3))')t222,g222,f222
c
c	escribe en out_li
	write(30,'(A2,1x,(F5.2,1x),2(F7.4,1x),32(F7.3))')
c	write(30,'(A2,1x,6(F5.2,1x),36(F7.3),8(1x,F8.3))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &	(vsin(n31),n31=8,nobs)
c     &	(vsin(n31),n31=1,nobs),(ZISOL(mlum),mlum=1,8)
	write(*,'(A2,1x,(F5.2,1x),2(F7.4,1x),39(F7.3))')
     &  shapa,ZMU-1.,zetao,realeo,(vsin(n31),n31=1,nobs)
c
c	escribe en out_qua
	write(39,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F12.8))')shapa,ZMU-1.,zetao,realeo
     & ,quaqn/quaq10,quaqn/quaq05,quaqn/quaq15,quaqn,quaq10,quaq05,quaq15
	write(*,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F12.8))')shapa,ZMU-1.,zetao,realeo
     & ,quaqn/quaq10,quaqn/quaq05,quaqn/quaq15,quaqn,quaq10,quaq05,quaq15
c
c	escribe en out_ef
     	write(31,'(A2,1x,(F5.2,1x),2(F7.4,1x),3(F6.0,2(1x,F6.3)))')
c     &  shapa,ZMU-1.,zetao,realeo,
     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &	(a90(nt),nt=1,3),(a90c(nt),nt=1,3)
     &  ,t222,g222,f222
c
c	escribe en out_dg
        write(32,'(A2,1x,(F5.2,1x),2(F7.4,1x),2(8(F5.2,1x),2x))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &  (DWGI(1,nt),nt=1,8),(DWGI(2,nt),nt=1,8)
c
c	escribe en out_u
        write(29,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F14.8,1x))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &  fxUtm/fxUt,fxUti/fxUt,fxUt,
     &  fxUtms/fxUts,fxUtis/fxUts,fxUtis004/fxUts,fxUts
c
c	escribe en out_sbf
	write(33,'(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),1x,4(F5.2,1x),3(F7.3))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &	(SBF(nsb),nsb=1,8),(vSBF(nsb),nsb=1,7),
     &  (fhsbf(nt),nt=1,4),(cfhsbf(nt),nt=1,3)
c
c	escribe en out_phot
	write(34,'(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),8(1x,F8.3),7(F7.3))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &	(ZZZ(nt),nt=1,8),(vsin(nt),nt=1,7),(ZISOL(nt),nt=1,8),
     &  (fhst(nt),nt=1,4),(cfhst(nt),nt=1,3)
	write(*,'(8(1x,F8.3))')(ZISOL(nt),nt=1,8)
c
c	escribe en out_mass
	write(40,'(A2,1x,(F5.2,1x),2(F7.4,1x),5(F8.4),3(1x,F8.3),1x,F15.8)')
     &  shapa,ZMU-1.,zetao,realeo,
     &  xmast,fmO,fmA,fmR,(1.d0-fmO),
     &  ZISOL(3),(ZISOL(3)*(fmA/fmO)),
     &  ZZZ(3),fZISOL(3)
c
c	escribe en out_contr
        write(35,'(A2,1x,(F5.2,1x),2(F7.4,1x),5(9(f6.2,1x),2x))')
c     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  shapa,ZMU-1.,zetao,realeo,
     &  (hh(1,nt),nt=1,9),(hh(2,nt),nt=1,9),(hh(3,nt),nt=1,9)
     &  ,(hh(4,nt),nt=1,9),(hh(5,nt),nt=1,9)
	do kol=1,9
	tilin(kol)=0.d0
	enddo
        do jol=1,9
        do kol=1,5
         tilin(jol)=tilin(jol)+hh(kol,jol)
        enddo
        enddo
c        write(*,'(9(f8.3,1x))')(tilin(kol),kol=1,9)
c	escribe en fort.79 el valor del CAT medido sobre espectro
        write(79,'(A2,1x,6(F7.4,1x),(1(f7.3,1x)))')
     &  shapa,ZMU-1.,zetao,realeo,zetay,reay,alpha,
     &  CATFIN
c	ENDDO
	return
	end
