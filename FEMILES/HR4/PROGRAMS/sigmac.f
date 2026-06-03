cc*******SUBRUTINA SIGMAC************************************************
c This routine calculates the corresponding CaT spectrum for
c a given set of atmospheric parameters. It is called by the
c subroutine "hrsl". This routine uses "dsvdcmp" and "dsvbksb",
c "indexx" and "NICO". NICO subroutine is not included here.
c
        SUBROUTINE sigmac(tiss,giss,fiss,ala,aloec,Wcatt,fnc1,fnmgi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (nstm=1999) !miles; tambien vale para los demas arrays
        DOUBLE PRECISION ak(3,8),uk(3,8),wk(8),vk(8,8),bk(3),xk(8)
     & ,sco(8),vvk(8,8),bbk(8) !todo para metodo svd
	DIMENSION ssc(710,8),sscp(710,4,8),nsc(8),nbande(8)
	DIMENSION sscb(710,8),sscpb(710,4,8),nscb(8)
	DIMENSION ala(710,2),weigtc(8),weigc(710,8),cubstc(710,8)
	DIMENSION aloec(3),aloeci(3,8),snat(8)
	DIMENSION ilor(8),sweigt(8),silor(8),WE(8),WE0(8),
     &  silor0(1000,8),llup(1000,6)
        DIMENSION starcs(nstm,710,2) !CaT STARS ARRAY
c	dimension indxk(3)
        INTEGER ntonsc
	CHARACTER*80 ssc,sscb,starc
        COMMON/hrcat/aac(710,4),nstarc
	COMMON/clstdc/starc(710)
        COMMON/catsta/starcs !CaT STARS ARRAY: common routines: a,sigmac
c
c	tissMAX=36300.0d0
c	tissMIN=2747.0d0
c	gissMAX=5.12d0
c	gissMIN=0.00d0
c
	ddt0=0.009d0
	tlinmi=60.0d0
c	ddt0M=0.18d0 !23000,3000
c	ddt0M=0.15d0 !este es el que vale
	ddt0M=0.1696d0 !este es el que vale
c	ddt0M=0.212d0 !este es el que vale
c	ddt0M=0.365d0 !(5040/3000-5040/23000)/4
c	tlinma=4200.0d0 !33553/8=tlinma
	tlinma=3355.0d0 !20000/8=tlinma
c	tlinma=2000.0d0 !20000/8=tlinma
c	tlinma=5000.0d0 !33553/8=tlinma
	if((ddt0*tiss*tiss/5040.0d0).lt.tlinmi)then
	  ddt=abs(5040.0d0*tlinmi/(tiss*tiss))
	elseif((ddt0*tiss*tiss/5040.0d0).gt.tlinma)then
	  ddt=abs(5040.0d0*tlinma/(tiss*tiss))
	else
	  ddt=ddt0
	endif
	if(abs(ddt0M*tiss*tiss/5040.0d0).gt.tlinma)then
	  ddtmax=abs(5040.0d0*tlinma/(tiss*tiss))
	else
	  ddtmax=ddt0M
	endif
	tisst=5040.0d0/tiss
	ddg=0.18d0
	ddz=0.09d0
c
c	ddgmax=1.28d0
c	ddgmax=0.64d0
c	ddzmax=1.02d0
c	ddzmax=0.625d0
c	ddzmax=0.51d0
c	ddzmax=0.25d0
c
	ddgmax=0.512d0
	ddzmax=0.408d0
c
c	DMA=1250.0d0  
	DMA=700.0d0  !densidad limite (0.9974 frac.cum.)
	sigtim=3.0d0  !times*sigma para volumen calculo densidad
	rtimes=0.0d0
	npxstc=710
c	xs=1.0d0 !ancho minimo de la caja a usar
	xs=1.5d0 !ancho minimo de la caja a usar
        nttmax=10 !(4sigma-1sigma)/1=3
	fxs=0.5d0 !fraccion a incrementar si no hay estrellas
c	snoise=0.58d0 !limite de calidad del espectro
	snamax=367.0d0 !2.5sigma,142.,254.,367.
	snamin=22.0d0 !2sigma,13.,17.5,22. son 155*<=22.0
c	snamin=17.50d0 
c
	blc1=8474.90d0
        blc2=8807.10d0
	xnumno=0.0d0 !numero pixels en banda normalizacion
	do jij=1,npxstc 
	  ala(jij,1)=starcs(1,jij,1)
	  if(ala(jij,1).ge.blc1.and.ala(jij,1).le.blc2)then
	    xnumno=xnumno+1.00d0
	  endif	 
	  ala(jij,2)=0.0d0
	enddo
c
c CALCULO DE DENSIDAD DE ESTRELLAS Y SIGMAS 
c
       	     sigmat=0.0d0
             sigmag=0.0d0
             sigmaz=0.0d0
	DO nnnns=1,1000
	  rtimes=rtimes+1.0d0
	  starsi=0.0d0
	  dens=0.0d0
	  dddt=rtimes*sigtim*ddt
	  dddg=rtimes*sigtim*ddg
	  dddz=rtimes*sigtim*ddz
	  do l=1,nstarc
           IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	    if(fiss.gt.-0.05d0)then
	if(abs(5040.0d0/aac(l,1)-tisst).le.dddt.and.
     &  abs(aac(l,2)-giss).le.dddg.and.aac(1,3).gt.-0.15d0)then
     	      starsi=starsi+1.0d0
	endif
	    else
	if(abs(5040.0d0/aac(l,1)-tisst).le.dddt.and.
     &  abs(aac(l,2)-giss).le.dddg.and.aac(1,3).le.0.0d0)then
     	      starsi=starsi+1.0d0
	endif
	    endif
	   ELSE
	    if(abs(5040.0d0/aac(l,1)-tisst).le.dddt
     &       .and.abs(aac(l,2)-giss).le.dddg
     &       .and.abs(aac(l,3)-fiss).le.dddz)then
     	     starsi=starsi+1.0d0
	    endif
	   ENDIF
	  enddo
	  if(starsi.gt.0.0d0)then
	    dens=starsi/(2.0d0*dddt*2.0d0*dddg*2.0d0*dddz)
	    if(dens.gt.DMA)then
       	     sigmat=ddt
             sigmag=ddg
             sigmaz=ddz
            else
	     dense=((dens-DMA)/DMA)**2.      
             sigmat=ddt*exp(dense*dlog(ddtmax/ddt))
             sigmag=ddg*exp(dense*dlog(ddgmax/ddg))
             sigmaz=ddz*exp(dense*dlog(ddzmax/ddz))
c	sigmat=dens*dens-(dens/DMA)*(DMA*DMA+ddtmax-ddt)+ddtmax
c	sigmag=dens*dens-(dens/DMA)*(DMA*DMA+ddgmax-ddg)+ddgmax
c	sigmaz=dens*dens-(dens/DMA)*(DMA*DMA+ddzmax-ddz)+ddzmax
c 	     sigmat=ddtmax+(ddt-ddtmax)*dens/DMA
c	     sigmag=ddgmax+(ddg-ddgmax)*dens/DMA
c	     sigmaz=ddzmax+(ddz-ddzmax)*dens/DMA
           endif      
	    goto 222
	  endif
	ENDDO
222	continue
c
c BUSQUEDA DE ESTRELLAS DE ACUERDO A ESTAS SIGMAS PARA 8 CUBOS
c
	do lp=1,8
          nsc(lp)=0 !numero de estrellas por cubo
          nscb(lp)=0 !numero de estrellas por cubo provisional
	  nbande(lp)=0 !bandera: 0=no estrellas 1=si estrellas
	enddo
	rtt=0.0d0 !iteracion que servira para multiplicar por sigma
	rntt=-1.0d0 !contador necesario para aumentar tamanyo caja
	rbandi=0.0d0 !bandera para detener proceso si hay * en 8cubos
        DO ntt=1,1000 !de la iteracion de sigma
	 rntt=rntt+1.0d0
	 rtt=rtt+1.0d0
c	 xxs=0.7d0*(xs+rntt*fxs) !factor que multiplica a sigma
	 xxs=xs+rntt*fxs !factor que multiplica a sigma
         DO is=1,nstarc !de las estrellas de la libreria
	  IF(nbande(1).eq.0)THEN
          IF(aac(is,1).ge.tiss.and.aac(is,2).ge.giss.and.
     &    aac(is,3).ge.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &      abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(1)=nscb(1)+1
              sscb(nscb(1),1)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(1),kkkkl,1)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(1),1,1)=5040.0d0/sscpb(nscb(1),1,1)
	      sscpb(nscb(1),3,1)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &      abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(1)=nscb(1)+1
              sscb(nscb(1),1)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(1),kkkkl,1)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(1),1,1)=5040.0d0/sscpb(nscb(1),1,1)
	      sscpb(nscb(1),3,1)=fiss
            endif
	   endif
	  ELSE
           if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.
     &	   and.abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(1)=nscb(1)+1
              sscb(nscb(1),1)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(1),kkkkl,1)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(1),1,1)=5040.0d0/sscpb(nscb(1),1,1)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(2).eq.0)THEN
          IF(aac(is,1).lt.tiss.and.aac(is,2).ge.giss.and.
     &      aac(is,3).ge.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(2)=nscb(2)+1
              sscb(nscb(2),2)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(2),kkkkl,2)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(2),1,2)=5040.0d0/sscpb(nscb(2),1,2)
	      sscpb(nscb(2),3,2)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(2)=nscb(2)+1
              sscb(nscb(2),2)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(2),kkkkl,2)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(2),1,2)=5040.0d0/sscpb(nscb(2),1,2)
	      sscpb(nscb(2),3,2)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040./aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(2)=nscb(2)+1
              sscb(nscb(2),2)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(2),kkkkl,2)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(2),1,2)=5040.0d0/sscpb(nscb(2),1,2)
	   endif
	  ENDIF	  
	  ENDIF
	  ENDIF
c
	  IF(nbande(3).eq.0)THEN
          IF(aac(is,1).lt.tiss.and.aac(is,2).lt.giss.and.
     &      aac(is,3).ge.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(3)=nscb(3)+1
              sscb(nscb(3),3)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(3),kkkkl,3)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(3),1,3)=5040.0d0/sscpb(nscb(3),1,3)
	      sscpb(nscb(3),3,3)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(3)=nscb(3)+1
              sscb(nscb(3),3)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(3),kkkkl,3)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(3),1,3)=5040.0d0/sscpb(nscb(3),1,3)
	      sscpb(nscb(3),3,3)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040./aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(3)=nscb(3)+1
              sscb(nscb(3),3)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(3),kkkkl,3)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(3),1,3)=5040.0d0/sscpb(nscb(3),1,3)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(4).eq.0)THEN
          IF(aac(is,1).ge.tiss.and.aac(is,2).lt.giss.and.
     &      aac(is,3).ge.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(4)=nscb(4)+1
              sscb(nscb(4),4)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(4),kkkkl,4)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(4),1,4)=5040.0d0/sscpb(nscb(4),1,4)
	      sscpb(nscb(4),3,4)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(4)=nscb(4)+1
              sscb(nscb(4),4)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(4),kkkkl,4)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(4),1,4)=5040.0d0/sscpb(nscb(4),1,4)
	      sscpb(nscb(4),3,4)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz)then
              nscb(4)=nscb(4)+1
              sscb(nscb(4),4)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(4),kkkkl,4)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(4),1,4)=5040.0d0/sscpb(nscb(4),1,4)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(5).eq.0)THEN
          IF(aac(is,1).ge.tiss.and.aac(is,2).ge.giss.and.
     &      aac(is,3).lt.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(5)=nscb(5)+1
              sscb(nscb(5),5)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(5),kkkkl,5)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(5),1,5)=5040.0d0/sscpb(nscb(5),1,5)
	      sscpb(nscb(5),3,5)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(5)=nscb(5)+1
              sscb(nscb(5),5)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(5),kkkkl,5)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(5),1,5)=5040.0d0/sscpb(nscb(5),1,5)
	      sscpb(nscb(5),3,5)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(5)=nscb(5)+1
              sscb(nscb(5),5)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(5),kkkkl,5)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(5),1,5)=5040.0d0/sscpb(nscb(5),1,5)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(6).eq.0)THEN
          IF(aac(is,1).lt.tiss.and.aac(is,2).ge.giss.and.
     &      aac(is,3).lt.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(6)=nscb(6)+1
              sscb(nscb(6),6)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(6),kkkkl,6)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(6),1,6)=5040.0d0/sscpb(nscb(6),1,6)
	      sscpb(nscb(6),3,6)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(6)=nscb(6)+1
              sscb(nscb(6),6)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(6),kkkkl,6)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(6),1,6)=5040.0d0/sscpb(nscb(6),1,6)
	      sscpb(nscb(6),3,6)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(6)=nscb(6)+1
              sscb(nscb(6),6)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(6),kkkkl,6)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(6),1,6)=5040.0d0/sscpb(nscb(6),1,6)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(7).eq.0)THEN
          IF(aac(is,1).lt.tiss.and.aac(is,2).lt.giss.and.
     &      aac(is,3).lt.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(7)=nscb(7)+1
              sscb(nscb(7),7)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(7),kkkkl,7)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
	      sscpb(nscb(7),3,7)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(7)=nscb(7)+1
              sscb(nscb(7),7)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(7),kkkkl,7)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
	      sscpb(nscb(7),3,7)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(7)=nscb(7)+1
              sscb(nscb(7),7)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(7),kkkkl,7)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(8).eq.0)THEN
          IF(aac(is,1).ge.tiss.and.aac(is,2).lt.giss.and.
     &      aac(is,3).lt.fiss)THEN
          IF(tiss.le.4075.0d0.or.tiss.ge.9000.0d0)THEN
	   if(fiss.gt.-0.05d0)then
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).gt.-0.15d0)then
              nscb(8)=nscb(8)+1
              sscb(nscb(8),8)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(8),kkkkl,8)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(8),1,8)=5040.0d0/sscpb(nscb(8),1,8)
	      sscpb(nscb(8),3,8)=fiss
	    endif
	   else
            if(abs(5040.0d0/aac(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aac(is,2)-giss).le.sigmag*xxs
     &      .and.aac(1,3).le.0.0d0)then
              nscb(8)=nscb(8)+1
              sscb(nscb(8),8)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(8),kkkkl,8)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(8),1,8)=5040.0d0/sscpb(nscb(8),1,8)
	      sscpb(nscb(8),3,8)=fiss
	    endif
	   endif
	  ELSE
           if(abs(5040./aac(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aac(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aac(is,3)-fiss).le.sigmaz*xxs)then
              nscb(8)=nscb(8)+1
              sscb(nscb(8),8)=starc(is)
              do kkkkl=1,4
                sscpb(nscb(8),kkkkl,8)=aac(is,kkkkl)
              enddo
	      sscpb(nscb(8),1,8)=5040.0d0/sscpb(nscb(8),1,8)
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
	 ENDDO !de las estrellas de la libreria
c
c descarte de * que degeneran resolucion parametrica: cubicos extremos
c xxs=xs+rntt*fxs si ntt=1 ===> rntt=0. fxs=0.5 (frac sigma)
         if(rntt.lt.0.50d0)then
c	   xxs2=xs+rntt*fxs
	   xxs2=1.0d0
	 else
           xxs2=xs+(rntt-1.0d0)*fxs
c	   xxs2=1.0d0
	 endif
	 DO lp1=1,8
	  IF(nbande(lp1).eq.0)THEN	   
           IF(nscb(lp1).gt.1)THEN
            do nlo=1,nscb(lp1)
        if(abs(sscpb(nlo,1,lp1)-tisst).gt.abs(sigmat*xxs2).
     &  and.abs(sscpb(nlo,2,lp1)-giss).gt.abs(sigmag*xxs2).and
     &  .abs(sscpb(nlo,3,lp1)-fiss).gt.abs(sigmaz*xxs2))then
 		  inosi=2
        elseif(abs(sscpb(nlo,1,lp1)-tisst).gt.abs(sigmat*xxs2).
     &  and.abs(sscpb(nlo,2,lp1)-giss).gt.abs(sigmag*xxs2))then
		  inosi=2
        elseif(abs(sscpb(nlo,1,lp1)-tisst).gt.abs(sigmat*xxs2).
     &  and.abs(sscpb(nlo,3,lp1)-fiss).gt.abs(sigmaz*xxs2))then
		  inosi=2
        elseif(abs(sscpb(nlo,2,lp1)-giss).gt.abs(sigmag*xxs2).
     &  and.abs(sscpb(nlo,3,lp1)-fiss).gt.abs(sigmaz*xxs2))then
		  inosi=2
        else
		  inosi=0
        endif
c	if(abs(sscpb(nlo,1,lp1)-tisst).gt.abs(sigmat*xxs2))then
c	           inosi=2
c	else
c		   inosi=0
c       endif
	     if(inosi.lt.1)then
               nsc(lp1)=nsc(lp1)+1 !nuevo contador de * en subcaja lp1
               ssc(nsc(lp1),lp1)=sscb(nlo,lp1)
               do kkkkl=1,4
                sscp(nsc(lp1),kkkkl,lp1)=sscpb(nlo,kkkkl,lp1)
               enddo
	     endif
            enddo
	    IF(nsc(lp1).gt.1)THEN
               nbande(lp1)=1
	       rbandi=rbandi+1.
            ELSEIF(nsc(lp1).eq.1)THEN
               if(sscp(1,4,lp1).gt.snamin)then
	        nbande(lp1)=1
	        rbandi=rbandi+1.
	       elseif(sscp(1,4,lp1).le.snamin)then
	        nscb(lp1)=0
		nsc(lp1)=0	        
	       endif
	    ENDIF
           ELSEIF(nscb(lp1).eq.1)THEN
            if(sscpb(1,4,lp1).gt.snamin)then
                nsc(lp1)=1
                ssc(1,lp1)=sscb(1,lp1)
                do kkkkl=1,4
                 sscp(1,kkkkl,lp1)=sscpb(1,kkkkl,lp1)
                enddo
	        nbande(lp1)=1
	        rbandi=rbandi+1.
	    else
	        nscb(lp1)=0
		nsc(lp1)=0
	    endif
           ENDIF
	  ENDIF
	 ENDDO
	 if(rbandi.gt.7.5)then
	   goto 258
	 elseif(rbandi.gt.2.5.and.rbandi.lt.7.5)then
c        garantizamos que al menos tenemos 3 subcajas
	  if(ntt.ge.nttmax)then
	   goto 258
	  endif
	 endif
        ENDDO !de la iteracion aumentadora de sigma
258	continue
c a partir de ahora rbandi ya refleja el numero de cubos con *
c
c OBTENCION DE ESTRELLA POR CUBO Y ESTRELLA FINAL NO NORMALIZADA
c
	DO lp=1,8 !numero de cubos
	  do jij=1,npxstc !inicializo la estrella final cubo lp
	    cubstc(jij,lp)=0. !solo interesa flujo y NO lambda
	  enddo
	  aloeci(1,lp)=0.0d0
	  aloeci(2,lp)=0.0d0
	  aloeci(3,lp)=0.0d0
	  snat(lp)=0.0d0
	  sweigt(lp)=0.0d0
	  weigtc(lp)=0.0d0 !peso total de estrellas en cubo lp
	  sumnor=0.0d0
	  IF(nsc(lp).gt.0)THEN
	  DO lp2=1,nsc(lp) !numero de estrellas en cubo lp
	   weigc(lp2,lp)=0.0d0 !peso estrella lp2 en cubo lp
	   do is=1,nstarc !estrellas libreria
	    if(ssc(lp2,lp).eq.starc(is))then
c	     if(aac(is,4).ge.snamax)then
             if(sscp(lp2,4,lp).ge.snamax)then
	      snaton=1.0d0
             else
c	      snaton=(aac(is,4)*aac(is,4))/(snamax*snamax)
              snaton=(sscp(lp2,4,lp)*sscp(lp2,4,lp))/(snamax*snamax)
             endif
c	     weigc(lp2,lp)=(exp(-((log10(aac(is,1)/tiss))/
             weigc(lp2,lp)=exp(-(((sscp(lp2,1,lp)-tisst)/sigmat)**2.))
     &       *exp(-(((sscp(lp2,2,lp)-giss)/sigmag)**2.))     
     &       *exp(-(((sscp(lp2,3,lp)-fiss)/sigmaz)**2.))*snaton
             weigtc(lp)=weigtc(lp)+weigc(lp2,lp)    
c para evitar numeros negativos en el espectro estelar debido a baja S/R
	     atonmi=9.e19
	     do jij=1,npxstc
	        if(starcs(is,jij,2).gt.0.0d0.and.
     &             starcs(is,jij,2).lt.atonmi)atonmi=starcs(is,jij,2)	     
	     enddo
c completada evaluacion valores negativos espectro estelar
	     do jij=1,npxstc 
c para evitar valores negativos espectro estelar
	      if(starcs(is,jij,2).le.0.0d0) starcs(is,jij,2)=atonmi
c "if" para evitar dividir al final por pesos ultra0 cuando solo hay 1 * en cubo
	      if(nsc(lp).eq.1)then 
	       cubstc(jij,lp)=starcs(is,jij,2)
	      else
	  cubstc(jij,lp)=cubstc(jij,lp)+weigc(lp2,lp)*starcs(is,jij,2)
	      endif	  
	     enddo
c        write(48,'(i1,x,A6,x,2(f6.0,x),3(f6.4,x),f6.0,x,
c     &  2(f6.3,x),f6.0,x,f7.5,2(x,i3))')lp,starc(is),
c     &  5040./(tisst+sigmat),5040./(tisst-sigmat),sigmat,
c     &  sigmag,sigmaz,5040./sscp(lp2,1,lp),sscp(lp2,2,lp),
c     &  sscp(lp2,3,lp),dens,weigc(lp2,lp),nscb(lp),nsc(lp)
        aloeci(1,lp)=aloeci(1,lp)+sscp(lp2,1,lp)*weigc(lp2,lp)
	aloeci(2,lp)=aloeci(2,lp)+sscp(lp2,2,lp)*weigc(lp2,lp)
	aloeci(3,lp)=aloeci(3,lp)+sscp(lp2,3,lp)*weigc(lp2,lp)
        snat(lp)=snat(lp)+snaton*weigc(lp2,lp)
	    endif
	   enddo
	  ENDDO
          aloeci(1,lp)=aloeci(1,lp)/weigtc(lp)
          aloeci(2,lp)=aloeci(2,lp)/weigtc(lp)
          aloeci(3,lp)=aloeci(3,lp)/weigtc(lp)
	  snat(lp)=snat(lp)/weigtc(lp)
	  do jij=1,npxstc
	   cubstc(jij,lp)=cubstc(jij,lp)/weigtc(lp)	 
	   if(ala(jij,1).ge.blc1.and.ala(jij,1).le.blc2)then
	    sumnor=sumnor+cubstc(jij,lp)
	   endif	 
	  enddo
c	  estrella final del cubo lp renormalizada
	  sumnor=sumnor/xnumno	
	  do jij=1,npxstc
	   cubstc(jij,lp)=cubstc(jij,lp)/sumnor
          enddo 
	  ENDIF
c	  if(aloeci(1,lp).gt.0.)then
c         write(48,'(i1,1x,i3,1x,2(f6.0,1x,2(f6.3,1x)))')lp,nsc(lp)
c     & ,tiss,giss,fiss,5040./aloeci(1,lp),aloeci(2,lp),aloeci(3,lp)
c          endif
	ENDDO
c
c ESTRELLA FINAL NORMALIZADA
c
	na0=1 !medidor de weigtc=0, si no hubiese ninguno = 1
	DO lp=1,8 !weigtc genericos por cubo, los anteriores se anulan
	 IF(nsc(lp).gt.0)THEN
	  weigtc(lp)=exp(-(((aloeci(1,lp)-tisst)/sigmat)**2.))
     &               *exp(-(((aloeci(2,lp)-giss)/sigmag)**2.))	
     &               *exp(-(((aloeci(3,lp)-fiss)/sigmaz)**2.))
c     &    *snat(lp) !puedo eliminar snat
c	  weigtc(lp)=1.0d0
c         weigtc(lp)=snat(lp)   
	 ELSE
	  weigtc(lp)=0.0d0
	  na0=na0+1
	 ENDIF
         sweigt(lp)=weigtc(lp)*snat(lp) !error famoso al estar dentro de IF
        ENDDO 
c Ajuste weigtc mas bajos con coeficientes sco para obtener * exacta
c Definimos dimensiones fisicas y logicas matrices: m=3 n=8 mp=8 np=8
      call indexx(8,sweigt,ilor) !ilor=indica orden creciente sweigt
c      write(50,*)''
c      write(50,'(A2,x,8(f13.10,x))')'sw',(sweigt(lp),lp=1,8)
c      write(50,'(A2,x,8(i2,12x))')'il',(ilor(lp),lp=1,8)
c      write(50,'(A2,x,8(f13.10,x))')'so',(sweigt(ilor(lp)),lp=1,8)
      aloec(1)=0.0d0
      aloec(2)=0.0d0
      aloec(3)=0.0d0
      aloec1=0.0d0
      aloec2=0.0d0
      aloec3=0.0d0
      weifin=0.0d0
      weifi0=0.0d0
      taloec=9999999999.0d0 !para buscar combinacion minimiza aloec-tiss
      do iak=1,8
	  WE(iak)=0.0d0
	  WE0(iak)=0.0d0
      enddo
      ilup=1
      do iak=1,8
         if(iak.lt.na0)then
	     silor0(1,ilor(iak))=0.0d0
	 else
	     silor0(1,ilor(iak))=1.0d0
	 endif
      enddo    

      llup(ilup,1)=0
      llup(ilup,2)=0
      llup(ilup,3)=0
      llup(ilup,4)=0
      llup(ilup,5)=0
      llup(ilup,6)=0
      do lup1=na0,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=0
        llup(ilup,3)=0
        llup(ilup,4)=0
        llup(ilup,5)=0
        llup(ilup,6)=0
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1)then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
      enddo
      IF(na0.le.5)THEN
       do lup1=na0,8
       do lup2=lup1+1,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=lup2
        llup(ilup,3)=0
        llup(ilup,4)=0
        llup(ilup,5)=0
        llup(ilup,6)=0
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1.or.iak.eq.lup2)then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
       enddo
       enddo
      IF(na0.le.4)THEN
       do lup1=na0,8
       do lup2=lup1+1,8
       do lup3=lup2+1,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=lup2
        llup(ilup,3)=lup3
        llup(ilup,4)=0
        llup(ilup,5)=0
        llup(ilup,6)=0
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1.or.iak.eq.lup2.or.iak.eq.lup3
     &   )then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
       enddo
       enddo
       enddo
      IF(na0.le.3)THEN
       do lup1=na0,8
       do lup2=lup1+1,8
       do lup3=lup2+1,8
       do lup4=lup3+1,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=lup2
        llup(ilup,3)=lup3
        llup(ilup,4)=lup4
        llup(ilup,5)=0
        llup(ilup,6)=0
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1.or.iak.eq.lup2.or.iak.eq.lup3
     &   .or.iak.eq.lup4)then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
       enddo
       enddo
       enddo
       enddo
      IF(na0.le.2)THEN
       do lup1=na0,8
       do lup2=lup1+1,8
       do lup3=lup2+1,8
       do lup4=lup3+1,8
       do lup5=lup4+1,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=lup2
        llup(ilup,3)=lup3
        llup(ilup,4)=lup4
        llup(ilup,5)=lup5
        llup(ilup,6)=0
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1.or.iak.eq.lup2.or.iak.eq.lup3
     &   .or.iak.eq.lup4.or.iak.eq.lup5)then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
       enddo
       enddo
       enddo
       enddo
       enddo
      IF(na0.eq.1)THEN
       do lup1=na0,8
       do lup2=lup1+1,8
       do lup3=lup2+1,8
       do lup4=lup3+1,8
       do lup5=lup4+1,8
       do lup6=lup5+1,8
	ilup=ilup+1
        llup(ilup,1)=lup1
        llup(ilup,2)=lup2
        llup(ilup,3)=lup3
        llup(ilup,4)=lup4
        llup(ilup,5)=lup5
        llup(ilup,6)=lup6
	do iak=1,8
         if(iak.lt.na0.or.iak.eq.lup1.or.iak.eq.lup2.or.iak.eq.lup3
     &   .or.iak.eq.lup4.or.iak.eq.lup5.or.iak.eq.lup6)then
		silor0(ilup,ilor(iak))=0.0d0
	 else
		silor0(ilup,ilor(iak))=1.0d0
	 endif
	enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      ENDIF
      kolan=0
      DO iakc=1,ilup !si todos weigtc's ne 0 ==> na0=1 y el caso general es para na0-1
c Definimos la matriz ak (pesos*diferencias) y el termino independiente (0)
	do iaku=1,8
	  silor(iaku)=silor0(iakc,iaku)
	enddo
	do iaku=1,3
	 do lpu=1,8
	   ak(iaku,lpu)=0.0d0
	 enddo
	enddo
	do lpu=1,8
	 ak(1,lpu)=weigtc(lpu)*(aloeci(1,lpu)-tisst)*silor(lpu)
	 ak(2,lpu)=weigtc(lpu)*(aloeci(2,lpu)-giss)*silor(lpu)
	 ak(3,lpu)=weigtc(lpu)*(aloeci(3,lpu)-fiss)*silor(lpu)
	enddo
	bk(1)=0.0d0
	bk(2)=0.0d0
	bk(3)=0.0d0
c Salvo la matriz original ak y trabajo con uk que sera destruida
	do iksvd=1,3
	 do jksvd=1,8
	   uk(iksvd,jksvd)=ak(iksvd,jksvd)
	 enddo
	enddo
c	write(*,*)5040.d0/tisst,giss,fiss
	call dsvdcmp(uk,3,8,3,8,wk,vk,ntonsc) !resolvemos sist. ec. lin.
	if(ntonsc.eq.1) goto 5045
	wmax=0.0D0
	do jksvd=1,8
	 if(wk(jksvd).gt.wmax) wmax=wk(jksvd)
	enddo
	wmin=wmax*1.0d-12 !valor original
	do jksvd=1,8
	 if(wk(jksvd).lt.wmin) wk(jksvd)=0.0D0
	enddo
	call dsvbksb(uk,wk,vk,3,8,3,8,bk,xk)
c calculo de la solucion mas cercana a 1. 1. 1. 1. 1. 1. 1. 1.
c la solucion del sistema anterior es d_sol + lambda*vk(i,j) si wk(i)=0.
	iiak=0
	do iak1=1,8
	  if(wk(iak1).eq.0.0d0)then
	    iiak=iiak+1
	    do iak2=1,8
	      vvk(iiak,iak2)=vk(iak2,iak1)
	    enddo
	  endif
	enddo
c
c	nnp=iiak mmp=iiak nn=iiak mm=iiak
c
c Empezamos un bucle de busqueda de solucion mas proxima a 1 y si no vamos
c aumentando hasta llegar a una solucion positiva:
	do iak1=1,iiak
	  bbk(iak1)=0.0d0
	  do iak2=1,8
c	    bbk(iak1)=bbk(iak1)+vvk(iak1,iak2)*(xk(iak2)-(1.0d0+dibu))
	    bbk(iak1)=bbk(iak1)+vvk(iak1,iak2)*(xk(iak2)-(1.0d0))
	  enddo
	  bbk(iak1)=-bbk(iak1) !faltaba el signo - en las formulas
c	  xxk(iak1)=bbk(iak1)
	enddo
c calculo del vector solucion mas proximo a 11111111
5045	DO lpu=1,8
	 if(ntonsc.eq.1)then !esta opcion es por si falla dsvdcmp
	  sco(lpu)=1.0d0
	 else
	  sco(lpu)=xk(lpu)
	  do iak1u=1,iiak
	 	sco(lpu)=sco(lpu)+bbk(iak1u)*vvk(iak1u,lpu)
	  enddo
	 endif
	 WE(lpu)=sco(lpu)*silor(lpu)*weigtc(lpu)
	ENDDO
c Ahora averiguamos si es mayor que 1 o menor que 0
       abigMA=MAX(WE(1),WE(2),WE(3),WE(4),WE(5),WE(6),WE(7),WE(8))
       abigMI=MIN(WE(1),WE(2),WE(3),WE(4),WE(5),WE(6),WE(7),WE(8))
c
c       write(48,'(A2,x,8(f13.10,x))')'S0',(silor0(iakc,lupu),lupu=1,8)
c       write(48,'(A2,x,8(f13.10,x))')'SI',(silor(lupu),lupu=1,8)
cccc       write(48,'(A2,1x,8(f13.10,1x),6(I2,1x))')'WE',
cccc     & (WE(lupu),lupu=1,8),(llup(iakc,lupu),lupu=1,6)
c
c      esto es por si la solucion es unica y vale 000       
       kolawn=0
       kolan=0
       if(abs(abigMA).lt.0.10d-10.and.abs(abigMI).lt.0.10d-10)then
         abigMA=MAX(weigtc(1),weigtc(2),weigtc(3),weigtc(4),
     &              weigtc(5),weigtc(6),weigtc(7),weigtc(8))
         do lp=1,8
	     WE(lp)=weigtc(lp) !si multiplico por silor(lp) rompe programa
         enddo
	 kolawn=1
       endif
       do lp=1,8
          if(WE(lp).le.0.0d0) then
	     WE(lp)=0.0d0
	  else
	     WE(lp)=WE(lp)/abigMA
	  endif
       enddo
       aloec(1)=0.0d0
       aloec(2)=0.0d0
       aloec(3)=0.0d0
       weifin=0.0d0
       do lp=1,8
          weifin=weifin+WE(lp) 
	  aloec(1)=WE(lp)*aloeci(1,lp)+aloec(1)
	  aloec(2)=WE(lp)*aloeci(2,lp)+aloec(2)
	  aloec(3)=WE(lp)*aloeci(3,lp)+aloec(3)
       enddo
       aloec(1)=aloec(1)/weifin
       aloec(2)=aloec(2)/weifin
       aloec(3)=aloec(3)/weifin
       IF(abigMI.ge.0.0d0.and.kolawn.eq.0)THEN !no vale kolawn=1 al ser sol. 000
	  kolan=2
          goto 2537
       ELSE
c        etaloe=exp(-((aloec(1)-tisst)/sigmat)**2.)
c     &         *exp(-((aloec(2)-giss)/sigmag)**2.)    
c     &         *exp(-((aloec(3)-fiss)/sigmaz)**2.)
        etaloe=abs(aloec(1)-tisst)
	if(etaloe.lt.taloec)then
	  taloec=etaloe
	  kolan=0
	  aloec1=aloec(1)
	  aloec2=aloec(2)
	  aloec3=aloec(3)
	  weifi0=weifin
	  do lup=1,8
	     WE0(lup)=WE(lup)
	  enddo
	endif
       ENDIF
      ENDDO !bucle iakc
2537   continue
       IF(kolan.lt.2)THEN
	weifin=weifi0
	aloec(1)=aloec1
	aloec(2)=aloec2
	aloec(3)=aloec3
	do lp=1,8
	 WE(lp)=WE0(lp)
	enddo
       ENDIF
cccc       write(48,'(A3,1x,8(f13.11,1x))')'wei',(weigtc(lp),lp=1,8)
cccc       write(48,'(A3,1x,8(f13.11,1x))')'WEF',(WE(lp),lp=1,8)
c        write(51,*)kolawn,kolan
       aloec(1)=5040.0d0/aloec(1) !paso a Teff
c calculamos espectro estrella final
       do lp=1,8
	if(WE(lp).gt.0.0d0)then
	 do jij=1,npxstc
          ala(jij,2)=ala(jij,2)+WE(lp)*cubstc(jij,lp)
	 enddo
	endif
       enddo
c normalizamos espectro estrella final	
       c1li=8473.5d0
       c1lf=8484.5d0
       cmgii=8781.0d0
       cmgif=8846.5d0
       fnc1=0.0d0
       fnmgi=0.0d0
       sumnor=0.0d0
       do jij=1,npxstc
	 ala(jij,2)=ala(jij,2)/weifin
	 if(jij.eq.npxstc)then
	  dlcat=abs(ala(jij-1,1)-ala(jij,1))
	 else
	  dlcat=abs(ala(jij+1,1)-ala(jij,1))
	 endif
	 if(ala(jij,1).ge.blc1.and.ala(jij,1).le.blc2)then
	    sumnor=sumnor+ala(jij,2)*dlcat
	 endif
c calculo del flujo en el continuo c1
	 if(ala(jij,1).gt.c1li.and.ala(jij,1).lt.c1lf)then
	    fnc1=fnc1+ala(jij,2)*dlcat
	 endif
	 if(ala(jij,1).gt.cmgii.and.ala(jij,1).lt.cmgif)then
	    fnmgi=fnmgi+ala(jij,2)*dlcat
	 endif
       enddo
       fnc1=fnc1/sumnor
       fnmgi=fnmgi/sumnor
c       write(16,*)fnc1,sumnor
       sumnor=sumnor/xnumno
c estrella final renormalizada
       do jij=1,npxstc
	   ala(jij,2)=ala(jij,2)/sumnor	 	
       enddo
       call NICO(ala,Wcatt) !valor del CaT en * final
c       write(48,'(A4,i2,1x,2(f6.0,1x,2(f6.3,1x)),f10.4)')'FIN='
c     &          ,kolan,tiss,giss,fiss,(aloec(ka),ka=1,3),Wcatt
c       write(91,'(2(f6.0,1x,2(f6.3,1x)),f10.4))')
c     &   (aloec(ka),ka=1,3),tiss,giss,fiss,Wcatt !comparar con fort.88
        return
	end	


