c*******SUBRUTINA QUE ESCRIBE LOS RESULTADOS******************************
	SUBROUTINE e_alphaCO(alpha,izo,realeo,izy,realey)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	CHARACTER*1 shape
	CHARACTER*2 shapa
	DIMENSION vsin(50),ZZZ(8),ZISOL(8),fZISOL(8),SBF(50),vSBF(50)
     &  ,DWGI(2,8),a90(3),a90c(3),fhst(4),fhsbf(4),cfhst(4),cfhsbf(4)
      dimension tilin(9)
	common/zalpCO/ZZZ,ZISOL,fZISOL,SBF,fhst,fhsbf !common a.f, e.f
      COMMON/esc/vsin,DWGI,a90,a90c,ncols,nobs
	COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
	COMMON/lmu1/ZMU
	COMMON/lshape/shape
      COMMON/TOAGB2/hh(5,9)
	COMMON/fezsol/fez(15)
	COMMON/qua/quaqn,quaq10,quaq05,quaq15
      COMMON/masas/xmast,xmaso,xmasa,fmO,fmA,fmR !common STU.f
	COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !common a.f
cccccccccccccccccccccccccccccccccc
c	write(*,*)'realeo_e',realeo
c        write(*,*)'e_alpha.f,ZZZ(3),ZISOL(3),fZISOL(3)',
c     &  ZZZ(3),ZISOL(3),fZISOL(3)
cccccccccccccccccccccccccccccccccc
	if(shape.eq.'u')then
	 shapa='UN'
	elseif(shape.eq.'b')then
	 shapa='BI'
	elseif(shape.eq.'k')then
	 shapa='Ku'
	elseif(shape.eq.'r')then
	 shapa='Kb'
	elseif(shape.eq.'c')then
	 shapa='Ch'
	endif
	if(alpha.eq.1.0)then
	  zetay=9.99d0
	  reay=9.99d0
	else
	  zetay=fez(izy)	
	  reay=realey
	endif	
	zetao=fez(izo)
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
c	escribe en out_u
        write(29,'(A2,1x,(F5.2,1x),2(F7.4,1x),4(F14.8,1x))')
     &  shapa,ZMU-1.,zetao,realeo,
     &  fxUtm/fxUt,fxUt,
     &  fxUtms/fxUts,fxUts
c	escribe en out_sbf
	write(33,'(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),1x,4(F5.2,1x),3(F7.3))')
     &  shapa,ZMU-1.,zetao,realeo,
     &	(SBF(nsb),nsb=1,8),(vSBF(nsb),nsb=1,7),
     &  (fhsbf(nt),nt=1,4),(cfhsbf(nt),nt=1,3)
c	escribe en out_phot
	write(34,'(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),8(1x,F8.3),7(F7.3))')
     &  shapa,ZMU-1.,zetao,realeo,
     &	(ZZZ(nt),nt=1,8),(vsin(nt),nt=1,7),(ZISOL(nt),nt=1,8),
     &  (fhst(nt),nt=1,4),(cfhst(nt),nt=1,3)
c	escribe en out_mass
	write(40,'(A2,1x,(F5.2,1x),2(F7.4,1x),5(F8.4),3(1x,F8.3),1x,F15.8)')
     &  shapa,ZMU-1.,zetao,realeo,
     &  xmast,fmO,fmA,fmR,(1.d0-fmO),
     &  ZISOL(3),(ZISOL(3)*(fmA/fmO)),
     &  ZZZ(3),fZISOL(3)
c	escribe en out_contr
        write(35,'(A2,1x,(F5.2,1x),2(F7.4,1x),5(9(f6.2,1x),2x))')
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
	return
	end
