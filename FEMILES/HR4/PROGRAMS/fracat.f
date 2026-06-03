c Subrutina que calcula fraccion zona CaT en I Johnson
c Se pasan tiss,giss y devuelve fraccion
c Se uso Pickles library: convolucion para CaT y filtro I Johnson
c Filtro Johnson I: subrutina respij.f 
c Llamadas desde hrsl.f,STU.f
      SUBROUTINE fracat(tiss,giss,frc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tissl=dlog10(tiss)
c      af=0.157067d0
c      bf=-0.00356389d0
c      afd=-1.53103d0
c      bfd=0.472297d0
c      afg=-0.274383d0
c      bfg=+0.118478d0
c      tfnd=(-afd+af)/(bfd-bf)
c      tfng=(-afg+af)/(bfg-bf)
c      IF(giss.ge.3.5d0)THEN
c       if(tissl.lt.tfnd)then
c        frc=afd+bfd*tissl
c       else
c        frc=af+bf*tissl
c       endif
c      ELSE
c       if(tissl.lt.tfng)then
c        frc=afg+bfg*tissl
c       else
c        frc=af+bf*tissl
c       endif
c      ENDIF
c TODAS
      tisM=4.5d0 
      at0=109.916d0
      at1=-106.724d0
      at2=38.8653d0
      at3=-6.28144d0
      at4=0.380069d0
c ENANAS
      atd0=49.3947d0
      atd1=-46.5568d0
      atd2=16.4893d0
      atd3=-2.59218d0
      atd4=0.152548d0
c GIGANTES
      atg0=149.082d0
      atg1=-146.120d0
      atg2=53.6948d0
      atg3=-8.75712d0
      atg4=0.534717d0
      frcM=at0+at1*tisM+at2*(tisM*tisM)+at3*(tisM*tisM*tisM)+
     & at4*(tisM*tisM*tisM*tisM)
      frcMd=atd0+atd1*tisM+atd2*(tisM*tisM)+atd3*(tisM*tisM*tisM)+
     & atd4*(tisM*tisM*tisM*tisM)
      frcMg=atg0+atg1*tisM+atg2*(tisM*tisM)+atg3*(tisM*tisM*tisM)+
     & atg4*(tisM*tisM*tisM*tisM)
      IF(giss.gt.3.5d0)THEN
       if(tissl.gt.tisM)then
        frc=frcMd
       else
        frc=atd0+atd1*tissl+atd2*(tissl*tissl)+atd3*
     &  (tissl*tissl*tissl)+atd4*(tissl*tissl*tissl*tissl)
       endif
      ELSE
       if(tissl.gt.tisM)then
        frc=frcMg
       else
        frc=atg0+atg1*tissl+atg2*(tissl*tissl)+atg3*
     &  (tissl*tissl*tissl)+atg4*(tissl*tissl*tissl*tissl)
       endif
      ENDIF
c      if(tissl.gt.tisM)then
c       frc=frcMt
c      else
c       frc=at0+at1*tissl+at2*(tissl*tissl)+at3*(tissl*tissl*tissl)
c     & +at4*(tissl*tissl*tissl*tissl)
c      endif
c      write(61,*)tissl,giss,frc
      RETURN
      END
cENANAS
c      49.3947
c     -46.5568
c      16.4893
c     -2.59218
c     0.152548
cGIGANTES
c      149.082
c     -146.120
c      53.6948
c     -8.75712
c     0.534717
c
c	write(*,*)'fracat=',tiss,tissl,giss,frc,tfnd,tfng
c
c CALCULO FRACCION FLUJO ZONA CAT
c	IF(R(m,i,3).lt.3.582d0)THEN
c	 if(R(m,i,4).ge.3.5d0)then
c	  frx1=-0.124224437d0+0.084027545d0*R(m,i,3) 
c	 else
c	  frx1=-1.34274181d0+0.42443167d0*R(m,i,3) 
c	 endif
c	ELSE
c	  frx1=0.198935997d0-0.00599133915d0*R(m,i,3)
c	ENDIF
c !Vazdekis et al. 2003 wrong!
c	 IF(R(m,i,3).lt.3.55d0)THEN
c	  if(R(m,i,4).ge.3.5d0)then
c	   frx1=-1.62067540d0+.4969163100d0*R(m,i,3) 
c	  else
c	   frx1=-.16346262d0+.0868521340d0*R(m,i,3) 
c	  endif
c	 ELSE
c	   frx1=.15695200d0-.0033235712d0*R(m,i,3)
c	 ENDIF
c
