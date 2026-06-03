c*****SUBRUTINA QUE ESCRIBE LOS RESULTADOS******************************
      SUBROUTINE e_alpha(alpha,izo,realeo,izy,realey,a90m)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 shape
      CHARACTER*2 shapa
      DIMENSION vsin(50),ZZZ(8),ZISOL(8),fZISOL(8),SBF(50),vSBF(50)
      DIMENSION DWGI(2,8),fhst(4),fhsbf(4),cfhst(4),cfhsbf(4)
      DIMENSION tilin(9)
      DIMENSION a90(3),a90c(3),a90m(5)
      COMMON/esc/vsin,ZZZ,SBF,DWGI,a90,a90c,ZISOL,fZISOL,ncols,nobs,
     &fhst,fhsbf
      COMMON/tgf/t222,g222,f222
      COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
      COMMON/lshape/shape
      COMMON/TOAGB2/hh(5,9)
      COMMON/CTFIN/CATFIN
c      COMMON/NC/ncha
      COMMON/fezsol/fez(15)
      COMMON/qua/quaqn,quaq10,quaq05,quaq15
      COMMON/masas/xmast,xmaso,xmasa,fmO,fmA,fmR !common subrut STU.f
cc      COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !common a.f
      COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3
      COMMON/ab/iaba,iabaj
      COMMON/bol/abola(3)
      if(shape.eq.'u')then
	 shapa='UN'
      elseif(shape.eq.'b')then
	 shapa='BI'
      elseif(shape.eq.'t')then
	 shapa='Ba'
      elseif(shape.eq.'k')then
	 shapa='Ku'
      elseif(shape.eq.'r')then
	 shapa='Kb'
      elseif(shape.eq.'c')then
	 shapa='Ch'
      elseif(shape.eq.'f')then
	 shapa='Fe'
      elseif(shape.eq.'x')then
	 shapa='Fx'
      elseif(shape.eq.'l')then
	 shapa='LB'
      endif
      if(alpha.eq.1.0)then
       zetay=9.99d0
       reay=9.99d0
      else
       zetay=fez(izy)	
       reay=realey
      endif	
      zetao=fez(izo)
c      DO N3=1,npol !polvo
c       do J=1,8
c        AZZZ(J)=ZZZ(J)+AD(J,N3)
c       enddo
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
c     	write(*,'(f6.0,2(1x,f6.3))')(a90(nt),nt=1,3)
c     	write(*,'(f6.0,2(1x,f6.3))')(a90c(nt),nt=1,3)
c     	write(*,'(f6.0,2(1x,f6.3))')t222,g222,f222
c
c    escribe en out_li
      IF(iaba.eq.0.and.iabaj.eq.1)THEN
518    FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),32(F7.3))
       WRITE(30,518)shapa,ZMU-1.d0,zetao,realeo,(vsin(n31),n31=8,nobs)
      ENDIF
c    escribe en out_qua
517   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),7(F12.8))
      WRITE(39,517)shapa,ZMU-1.d0,zetao,realeo,quaqn/quaq10,
     &quaqn/quaq05,quaqn/quaq15,quaqn,quaq10,quaq05,quaq15
c    escribe en out_ef
516   FORMAT(A2,x,(F5.3,x),2(F7.4,x),F7.0,x,4(F7.3,x),F7.0,x,4(F7.3,x))
      WRITE(31,516)shapa,ZMU-1.d0,zetao,realeo,(a90m(nt),nt=1,5),
     &(a90c(nt),nt=1,3)
c    escribe en out_dg
515   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),2(8(F5.2,1x),2x))
      WRITE(32,515)shapa,ZMU-1.d0,zetao,realeo,(DWGI(1,nt),nt=1,8),
     &(DWGI(2,nt),nt=1,8)
c    escribe en out_u
cc514   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),7(F14.8,1x))
cc      WRITE(29,514)shapa,ZMU-1.d0,zetao,realeo,fxUtm/fxUt,fxUti/fxUt,
cc     &fxUt,fxUtms/fxUts,fxUtis/fxUts,fxUtis004/fxUts,fxUts
c    escribe en out_sbf
513   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),1x,4(F5.2,1x),3(F7.3))
      WRITE(33,513)shapa,ZMU-1.d0,zetao,realeo,(SBF(nsb),nsb=1,8),
     &(vSBF(nsb),nsb=1,7),(fhsbf(nt),nt=1,4),(cfhsbf(nt),nt=1,3)
c    escribe en out_phot
512   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),15(F7.3),8(1x,F8.3),7(F7.3))
      WRITE(34,512)shapa,ZMU-1.d0,zetao,realeo,(ZZZ(nt),nt=1,8),
     &(vsin(nt),nt=1,7),(ZISOL(nt),nt=1,8),(fhst(nt),nt=1,4),
     &(cfhst(nt),nt=1,3)
c    escribe en out_mass
511   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),6(F8.4),5(1x,F8.3),1x,
     &2(F15.8,1x))
      WRITE(40,511)shapa,ZMU-1.d0,zetao,realeo,xmast,fmO,fmA,fmR,
     &(1.d0-fmO),ZISOL(3),abola(3),(ZISOL(3)*(fmA/fmO)),
     &(abola(3)*(fmA/fmO)),ZZZ(3),abola(2),fZISOL(3),abola(1)
c    escribe en out_contr
510   FORMAT(A2,1x,(F5.2,1x),2(F7.4,1x),5(9(f6.2,1x),2x))
      WRITE(35,510)shapa,ZMU-1.d0,zetao,realeo,(hh(1,nt),nt=1,9),
     &(hh(2,nt),nt=1,9),(hh(3,nt),nt=1,9),(hh(4,nt),nt=1,9),
     &(hh(5,nt),nt=1,9)
      do kol=1,9
       tilin(kol)=0.d0
      enddo
      do jol=1,9
       do kol=1,5
        tilin(jol)=tilin(jol)+hh(kol,jol)
       enddo
      enddo
      RETURN
      END
