c Subrutina que calcula el espectros de alta resolucion
c integrados. Hace uso de las subrutinas "lbusc" y "sigmac",
c las cuales calculan los espectros estelares rangos Jones, 
c CaT y MILES para una terna de parametros atmosfericos. 
      SUBROUTINE hrsl_IRTF(k,wflujk,wflujks,al99k,fja99k,qum,qu10,qu05
     &,qu15)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c IRTF
      PARAMETER (npxk=15387) !IRTF
cccccccccccccccccccccccccccccc
      CHARACTER*41 cirtf
cccccccccccccccccccccccccccccc
      DIMENSION ssppk(npxk,2),wflujk(npxk,2),wflujks(npxk,2) !IRTF
      DIMENSION alk(3),al99k(3)
      COMMON/qsighk/xnbxsk,volpak,volk10,volk05,volk15
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
c     limites lambdas
      XK1=18500.0d0 !IRTF XK1=19400.0d0 !IRTF
	XK2=25500.0d0 !IRTF XK2=24800.0d0 !IRTF
ccccccccccccccccccccccccc
c [Fe/H] (si [Mg/Fe].ne.0) a partir de metalicidad isocrona [Z/H]=fez(k)
c solo para sigmam.f
      fiss=fez(k)
      fiss=fez(k)-DeltFeH
ccccccccccccccccccccccccc
      qum=0.0d0 !parametro calidad (IRTF) normalizado sigma_min
      qu10=0.0d0 !parametro calidad (IRTF) aceptable 1.0 sigma
      qu05=0.0d0 !parametro calidad (IRTF) aceptable 0.5 sigma
      qu15=0.0d0 !parametro calidad (IRTF) aceptable 1.5 sigma
      starmin=3.0d0 !minimo aceptable estrellas: 3 (1 x subcaja)
      fja99k=0.0d0
      do l=1,3
	 al99k(l)=0.d0 !IRTF
      enddo
      do l=1,npxk !IRTF
	 wflujk(l,2)=0.d0
	 wflujks(l,2)=0.d0
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
	 call sigmam_IRTF(tiss,giss,fiss,ssppk,alk) !IRTF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
       if(tiss.lt.10000.0d0)then
        write(cirtf(1:8),'(A2,A1,F5.0)')'iT','0',tiss
       else
        write(cirtf(1:8),'(A2,F6.0)')'iT',tiss
       endif
       if(giss.lt.0.0d0)then
        write(cirtf(9:14),'(A1,A1,F4.2)')'G','m',abs(giss)
       else
        write(cirtf(9:14),'(A1,A1,F4.2)')'G','p',abs(giss)
       endif
       if(fiss.lt.0.0d0)then
        write(cirtf(15:20),'(A1,A1,F4.2)')'Z','m',abs(fiss)
       else
        write(cirtf(15:20),'(A1,A1,F4.2)')'Z','p',abs(fiss)
       endif
c       
       if(alk(1).lt.10000.0d0)then
        write(cirtf(21:29),'(A3,A1,F5.0)')'_oT','0',alk(1)
       else
        write(cirtf(21:29),'(A3,F6.0)')'_oT',alk(1)
       endif
       if(alk(2).lt.0.0d0)then
        write(cirtf(30:35),'(A1,A1,F4.2)')'G','m',abs(alk(2))
       else
        write(cirtf(30:35),'(A1,A1,F4.2)')'G','p',abs(alk(2))
       endif
       if(alk(3).lt.0.0d0)then
        write(cirtf(36:41),'(A1,A1,F4.2)')'Z','m',abs(alk(3))
       else
        write(cirtf(36:41),'(A1,A1,F4.2)')'Z','p',abs(alk(3))
       endif
	OPEN(43,FILE=cirtf,STATUS='OLD',ERR=6976)
	CLOSE(43,STATUS='DELETE')
6976	OPEN(43,IOSTAT=IOS,FILE=cirtf,STATUS='NEW')
c	do nonoo=1,61
c	  write(43,'(A80)')cab(nonoo,4)
c	enddo
	do lupir=1,npxk
		write(43,*)(ssppk(lupir,lllll),lllll=1,2)
	enddo
	close(43)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
	 fkmil=0.0d0
       do l=1,npxk
	  rfK=0.0d0
	  if(l.eq.npxk)then
	   dl=abs(ssppk(l-1,1)-ssppk(l,1))
	  else
	   dl=abs(ssppk(l+1,1)-ssppk(l,1))
	  endif
	  if(ssppk(l,1).ge.XK1.and.ssppk(l,1).le.XK2)then
	   call respk(ssppk(l,1),rfK) !IRTF
	   fkmil=fkmil+rfK*ssppk(l,2)*dl
	  endif
       enddo
	 DO l=1,npxk !IRTF
	  ssppk(l,2)=ssppk(l,2)/fkmil
	  spfl=ssppk(l,2)*flujk
          if(m.eq.1)then
	   if(imam(k).gt.0)then
	    do mma=1,imam(k)
	     wflujk(l,2)=wflujk(l,2)+R140(k,mma)*spfl
	     wflujks(l,2)=wflujks(l,2)+R140(k,mma)*spfl*spfl
	    enddo
	   endif
	   wflujk(l,2)=wflujk(l,2)+R(k,1,14)*spfl
	   wflujks(l,2)=wflujks(l,2)+R(k,1,14)*spfl*spfl
	  else
	   wflujk(l,2)=wflujk(l,2)+R(k,m,14)*spfl
	   wflujks(l,2)=wflujks(l,2)+R(k,m,14)*spfl*spfl
	  endif
	 ENDDO
       if(m.eq.1)then
	   if(imam(k).gt.0)then
	    do mma=1,imam(k)
	     fja99k=fja99k+R140(k,mma)*flujk
	    enddo
	   endif
	   fja99k=fja99k+R(k,1,14)*flujk
	 else
	   fja99k=fja99k+R(k,m,14)*flujk
	 endif
	 do l=1,3 !IRTF
         if(m.eq.1)then
	    if(imam(k).gt.0)then
	     do mma=1,imam(k)
		al99k(l)=al99k(l)+alk(l)*R140(k,mma)*flujk
	     enddo
	    endif
	    al99k(l)=al99k(l)+alk(l)*R(k,1,14)*flujk
	   else
	    al99k(l)=al99k(l)+alk(l)*R(k,m,14)*flujk
	   endif
	 enddo
	 qum=qum+R(k,m,14)*flujk*(xnbxsk/volpak)
	 qu10=qu10+R(k,m,14)*flujk*(starmin/volk10)
	 qu05=qu05+R(k,m,14)*flujk*(starmin/volk05)
	 qu15=qu15+R(k,m,14)*flujk*(starmin/volk15)
      ENDDO
c      DO l=1,npxk !SBF
c       wflujks(l,2)=wflujks(l,2)/wflujk(l,2)
c      ENDDO
      RETURN
      END
