      SUBROUTINE sigmam_MIUSC(tiss,giss,fiss,ala,aloec)
CC Version 22-January-03-2014.
C------------------------------------------------------------------------------
C 	Copyright A. Vazdekis 
C 	Instituto de Astrofisica de Canarias
C 	E-mail: vazdekis@iac.es
C------------------------------------------------------------------------------
c This program computes the SED of a star, in the spectral 
c range 3540.5-7409.6A, for a given set of atmospheric parameters
C------------------------------------------------------------------------------
C INPUT FILES:
C	STAR_HEADER_MILES: header file
C	PARAM_MILES: Atmopheric parameters + SN file
C	m#V: normalized Miles stellar spectra ascii files
C------------------------------------------------------------------------------
C Subroutines:
C	sigmam: interpolator
C	respv: normalization according to the V filter
C	indexx: Numerical Recipes a_interp.f
C	dsvdcmp: Numerical Recipes
C	dsvbksb: Numerical Recipes
C	hunt: Numerical Recipes
C	polint: Numerical Recipes
C------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c MIUSC
      PARAMETER (npxi=6672) !MIUSC
      PARAMETER (nsti=1100)
      DIMENSION ssc(nsti,8),sscp(nsti,4,8),sscb(nsti,8),sscpb(nsti,4,8)
	DIMENSION ala(npxi,2),weigc(nsti,8),cubstc(npxi,8)
c
	DIMENSION ak(3,8),uk(3,8),wk(8),vk(8,8),bk(3),xk(8)
     & ,sco(8),vvk(8,8),bbk(8) !todo para metodo svd
      DIMENSION nsc(8),nbande(8),nscb(8)
      DIMENSION weigtc(8),sigbox(8,3)
      DIMENSION aloec(3),aloeci(3,8),snat(8)
      DIMENSION ilor(8),sweigt(8),silor(8),WE(8),WE0(8),silor0(1000,8)
     &,llup(1000,6)
      INTEGER ilup,iiak
      CHARACTER*80 ssc,sscb
c MIUSC
      CHARACTER*80 stari
      COMMON/aamcoi/aai(nsti,4) !MIUSC
      COMMON/hri1/nstari !MIUSC
      COMMON/hri2/stari(nsti) !MIUSC
      COMMON/iilsta/staris(nsti,npxi,2) !MIUSC STARS ARRAY: common routines: a,sigmam
      COMMON/qsighi/xnbxsi,volpai,voli10,voli05,voli15
c
      iiak=0
      nbxstt=0 !numero total estrellas en las 8 subcajas
      re10=0.6827d0*0.5d0 !frac. subbox max. norm.Q (0.34=>1sigma)
      re05=0.3829d0*0.5d0 !frac. subbox max. norm.Q (0.19=>0.5sigma)
      re15=0.8664d0*0.5d0 !frac. subbox max. norm.Q (0.43=>1.5sigma)
	ultra0=1.0d-12
	t50=5040.0d0
c	tz0_i=3500.0d0 !temp mas baja para la que distinguimos[M/H]
	tz0_i=3500.0d0 !temp mas baja para la que distinguimos[M/H]
	tz0_f=9000.0d0 !temp mas alta para la que distinguimos [M/H]
	blc1=4900.0d0 !lambda inicial filtro V donde supera 0.1
        blc2=6300.0d0 !lambda final filtro V donde supera 0.1
	ddt0=0.009d0
	tlinmi=60.0d0
	ddt0M=0.1696d0
	tlinma=3355.0d0
	if((ddt0*tiss*tiss/t50).lt.tlinmi)then
	  ddt=abs(t50*tlinmi/(tiss*tiss))
	elseif((ddt0*tiss*tiss/t50).gt.tlinma)then
	  ddt=abs(t50*tlinma/(tiss*tiss))
	else
	  ddt=ddt0
	endif
	if(abs(ddt0M*tiss*tiss/t50).gt.tlinma)then
	  ddtmax=abs(t50*tlinma/(tiss*tiss))
	else
	  ddtmax=ddt0M
	endif
	tisst=t50/tiss
	ddg=0.18d0
	ddz=0.09d0
	ddgmax=0.512d0
	ddzmax=0.408d0
	do isbox=1,8
	      sigbox(isbox,1)=ddtmax
	      sigbox(isbox,2)=ddgmax
	      sigbox(isbox,3)=ddzmax
	enddo
	DMA=800.0d0 !densidad limite (0.9974 frac.cum.)(CAT=700)
	sigtim=3.0d0 !times*sigma para volumen calculo densidad
	rtimes=0.0d0
	xs=1.5d0 !ancho minimo de la caja a usar
        nttmax=10 !(4sigma-1sigma)/1=3
	fxs=0.5d0 !fraccion a incrementar si no hay estrellas
	snamax=390.0d0 !0.970408 CAT(367)CAT 2.5sigma,142.,254.,367.
	snamin=40.0d0 !0.03 CAT(22) !2sigma,13.,17.5,22. son 155*<=22.0
c
	xnumno=0.0d0 !numero pixels en banda normalizacion
	do jij=1,npxi 
	  ala(jij,1)=staris(1,jij,1)
	  if(ala(jij,1).ge.blc1.and.ala(jij,1).le.blc2)then
	    xnumno=xnumno+1.00d0
	  endif	 
	  ala(jij,2)=0.0d0
	enddo
	CLOSE(99)
c CALCULO DE DENSIDAD DE ESTRELLAS Y SIGMAS 
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
	  do l=1,nstari
           IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	    if(fiss.gt.-0.05d0)then
	if(abs(t50/aai(l,1)-tisst).le.dddt.and.
c     &  abs(aai(l,2)-giss).le.dddg.and.aai(1,3).gt.-0.15d0)then
     &  abs(aai(l,2)-giss).le.dddg)then
     	      starsi=starsi+1.0d0
	endif
c	    else
c	if(abs(t50/aai(l,1)-tisst).le.dddt.and.
c     &  abs(aai(l,2)-giss).le.dddg.and.aai(1,3).le.0.0d0)then
c     	      starsi=starsi+1.0d0
c	endif
c	    endif
	   ELSE
	    if(abs(t50/aai(l,1)-tisst).le.dddt
     &       .and.abs(aai(l,2)-giss).le.dddg
     &       .and.abs(aai(l,3)-fiss).le.dddz)then
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
           endif      
	    goto 222
	  endif
	ENDDO
222	continue
c	   WRITE(*,*)dddt,dddg,dddz
c	   WRITE(*,*)ddt,ddg,ddz
c	   WRITE(*,*)starsi,dens
c	   WRITE(*,*)sigmat,sigmag,sigmaz
c BUSQUEDA ESTRELLAS DE ACUERDO A ESTAS SIGMAS PARA 8 CUBOS
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
	 xxs=xs+rntt*fxs !factor que multiplica a sigma
         DO is=1,nstari !de las estrellas de la libreria
c       __________________________________________________
c	 if (is.ne.nist) then !exclusion estrella a ajustar
c       __________________________________________________
	  IF(nbande(1).eq.0)THEN
          IF(aai(is,1).ge.tiss.and.aai(is,2).ge.giss.and.
     &    aai(is,3).ge.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &      abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(1)=nscb(1)+1
              sscb(nscb(1),1)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(1),kkkkl,1)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(1),1,1)=t50/sscpb(nscb(1),1,1)
	      sscpb(nscb(1),3,1)=fiss
	      sigbox(1,1)=sigmat*xxs
	      sigbox(1,2)=sigmag*xxs
	      sigbox(1,3)=ddzmax
	    endif
c	   else
c            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c              nscb(1)=nscb(1)+1
c              sscb(nscb(1),1)=stari(is)
c              do kkkkl=1,4
c                sscpb(nscb(1),kkkkl,1)=aai(is,kkkkl)
c              enddo
c	      sscpb(nscb(1),1,1)=t50/sscpb(nscb(1),1,1)
c	      sscpb(nscb(1),3,1)=fiss
c            endif
c	   endif
	  ELSE
           if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.
     &	   and.abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(1)=nscb(1)+1
              sscb(nscb(1),1)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(1),kkkkl,1)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(1),1,1)=t50/sscpb(nscb(1),1,1)
	      sigbox(1,1)=sigmat*xxs
	      sigbox(1,2)=sigmag*xxs
	      sigbox(1,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(2).eq.0)THEN
          IF(aai(is,1).lt.tiss.and.aai(is,2).ge.giss.and.
     &      aai(is,3).ge.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(2)=nscb(2)+1
              sscb(nscb(2),2)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(2),kkkkl,2)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(2),1,2)=t50/sscpb(nscb(2),1,2)
	      sscpb(nscb(2),3,2)=fiss
	      sigbox(2,1)=sigmat*xxs
	      sigbox(2,2)=sigmag*xxs
	      sigbox(2,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(2)=nscb(2)+1
c	       sscb(nscb(2),2)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(2),kkkkl,2)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(2),1,2)=t50/sscpb(nscb(2),1,2)
c	       sscpb(nscb(2),3,2)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(5040./aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(2)=nscb(2)+1
              sscb(nscb(2),2)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(2),kkkkl,2)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(2),1,2)=t50/sscpb(nscb(2),1,2)
	      sigbox(2,1)=sigmat*xxs
	      sigbox(2,2)=sigmag*xxs
	      sigbox(2,3)=sigmaz*xxs
	   endif
	  ENDIF	  
	  ENDIF
	  ENDIF
c
	  IF(nbande(3).eq.0)THEN
          IF(aai(is,1).lt.tiss.and.aai(is,2).lt.giss.and.
     &      aai(is,3).ge.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(3)=nscb(3)+1
              sscb(nscb(3),3)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(3),kkkkl,3)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(3),1,3)=t50/sscpb(nscb(3),1,3)
	      sscpb(nscb(3),3,3)=fiss
	      sigbox(3,1)=sigmat*xxs
	      sigbox(3,2)=sigmag*xxs
	      sigbox(3,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(3)=nscb(3)+1
c	       sscb(nscb(3),3)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(3),kkkkl,3)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(3),1,3)=t50/sscpb(nscb(3),1,3)
c	       sscpb(nscb(3),3,3)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(5040./aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(3)=nscb(3)+1
              sscb(nscb(3),3)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(3),kkkkl,3)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(3),1,3)=t50/sscpb(nscb(3),1,3)
	      sigbox(3,1)=sigmat*xxs
	      sigbox(3,2)=sigmag*xxs
	      sigbox(3,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(4).eq.0)THEN
          IF(aai(is,1).ge.tiss.and.aai(is,2).lt.giss.and.
     &      aai(is,3).ge.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(4)=nscb(4)+1
              sscb(nscb(4),4)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(4),kkkkl,4)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(4),1,4)=t50/sscpb(nscb(4),1,4)
	      sscpb(nscb(4),3,4)=fiss
	      sigbox(4,1)=sigmat*xxs
	      sigbox(4,2)=sigmag*xxs
	      sigbox(4,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(4)=nscb(4)+1
c	       sscb(nscb(4),4)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(4),kkkkl,4)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(4),1,4)=t50/sscpb(nscb(4),1,4)
c	       sscpb(nscb(4),3,4)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz)then
              nscb(4)=nscb(4)+1
              sscb(nscb(4),4)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(4),kkkkl,4)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(4),1,4)=t50/sscpb(nscb(4),1,4)
	      sigbox(4,1)=sigmat*xxs
	      sigbox(4,2)=sigmag*xxs
	      sigbox(4,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(5).eq.0)THEN
          IF(aai(is,1).ge.tiss.and.aai(is,2).ge.giss.and.
     &      aai(is,3).lt.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(5)=nscb(5)+1
              sscb(nscb(5),5)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(5),kkkkl,5)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(5),1,5)=t50/sscpb(nscb(5),1,5)
	      sscpb(nscb(5),3,5)=fiss
	      sigbox(5,1)=sigmat*xxs
	      sigbox(5,2)=sigmag*xxs
	      sigbox(5,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(5)=nscb(5)+1
c	       sscb(nscb(5),5)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(5),kkkkl,5)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(5),1,5)=t50/sscpb(nscb(5),1,5)
c	       sscpb(nscb(5),3,5)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(5)=nscb(5)+1
              sscb(nscb(5),5)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(5),kkkkl,5)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(5),1,5)=t50/sscpb(nscb(5),1,5)
	      sigbox(5,1)=sigmat*xxs
	      sigbox(5,2)=sigmag*xxs
	      sigbox(5,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(6).eq.0)THEN
          IF(aai(is,1).lt.tiss.and.aai(is,2).ge.giss.and.
     &      aai(is,3).lt.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(6)=nscb(6)+1
              sscb(nscb(6),6)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(6),kkkkl,6)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(6),1,6)=t50/sscpb(nscb(6),1,6)
	      sscpb(nscb(6),3,6)=fiss
	      sigbox(6,1)=sigmat*xxs
	      sigbox(6,2)=sigmag*xxs
	      sigbox(6,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(6)=nscb(6)+1
c	       sscb(nscb(6),6)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(6),kkkkl,6)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(6),1,6)=t50/sscpb(nscb(6),1,6)
c	       sscpb(nscb(6),3,6)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(6)=nscb(6)+1
              sscb(nscb(6),6)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(6),kkkkl,6)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(6),1,6)=t50/sscpb(nscb(6),1,6)
	      sigbox(6,1)=sigmat*xxs
	      sigbox(6,2)=sigmag*xxs
	      sigbox(6,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(7).eq.0)THEN
          IF(aai(is,1).lt.tiss.and.aai(is,2).lt.giss.and.
     &      aai(is,3).lt.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(7)=nscb(7)+1
              sscb(nscb(7),7)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(7),kkkkl,7)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
	      sscpb(nscb(7),3,7)=fiss
	      sigbox(7,1)=sigmat*xxs
	      sigbox(7,2)=sigmag*xxs
	      sigbox(7,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(7)=nscb(7)+1
c	       sscb(nscb(7),7)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(7),kkkkl,7)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
c	       sscpb(nscb(7),3,7)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(7)=nscb(7)+1
              sscb(nscb(7),7)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(7),kkkkl,7)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(7),1,7)=5040./sscpb(nscb(7),1,7)
	      sigbox(7,1)=sigmat*xxs
	      sigbox(7,2)=sigmag*xxs
	      sigbox(7,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c
	  IF(nbande(8).eq.0)THEN
          IF(aai(is,1).ge.tiss.and.aai(is,2).lt.giss.and.
     &      aai(is,3).lt.fiss)THEN
          IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
c	   if(fiss.gt.-0.05d0)then
            if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
     &	    abs(aai(is,2)-giss).le.sigmag*xxs)then
c     &      .and.aai(1,3).gt.-0.15d0)then
              nscb(8)=nscb(8)+1
              sscb(nscb(8),8)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(8),kkkkl,8)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(8),1,8)=t50/sscpb(nscb(8),1,8)
	      sscpb(nscb(8),3,8)=fiss
	      sigbox(8,1)=sigmat*xxs
	      sigbox(8,2)=sigmag*xxs
	      sigbox(8,3)=ddzmax
	    endif
c	    else
c	     if(abs(t50/aai(is,1)-tisst).le.sigmat*xxs.and.
c     &      abs(aai(is,2)-giss).le.sigmag*xxs
c     &      .and.aai(1,3).le.0.0d0)then
c	       nscb(8)=nscb(8)+1
c	       sscb(nscb(8),8)=stari(is)
c	       do kkkkl=1,4
c		 sscpb(nscb(8),kkkkl,8)=aai(is,kkkkl)
c	       enddo
c	       sscpb(nscb(8),1,8)=t50/sscpb(nscb(8),1,8)
c	       sscpb(nscb(8),3,8)=fiss
c	     endif
c	    endif
	  ELSE
           if(abs(5040./aai(is,1)-tisst).le.sigmat*xxs.and.
     &	   abs(aai(is,2)-giss).le.sigmag*xxs.and.
     &	   abs(aai(is,3)-fiss).le.sigmaz*xxs)then
              nscb(8)=nscb(8)+1
              sscb(nscb(8),8)=stari(is)
              do kkkkl=1,4
                sscpb(nscb(8),kkkkl,8)=aai(is,kkkkl)
              enddo
	      sscpb(nscb(8),1,8)=t50/sscpb(nscb(8),1,8)
	      sigbox(8,1)=sigmat*xxs
	      sigbox(8,2)=sigmag*xxs
	      sigbox(8,3)=sigmaz*xxs
	   endif
	  ENDIF
	  ENDIF
	  ENDIF
c       __________________________________________________
c	 endif !exclusion estrella a ajustar
c       __________________________________________________
	 ENDDO !de las estrellas de la libreria
c
c descarte de * que degeneran resolucion parametrica: cubicos extremos
c xxs=xs+rntt*fxs si ntt=1 ===> rntt=0. fxs=0.5 (frac sigma)
         if(rntt.lt.0.50d0)then
	   xxs2=1.0d0
	 else
           xxs2=xs+(rntt-1.0d0)*fxs
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
        ENDDO !iteracion incrementadora de sigma
258	continue
c a partir de ahora rbandi ya refleja el numero de cubos con *
c
c OBTENCION DE ESTRELLA POR CUBO Y ESTRELLA FINAL NO NORMALIZADA
c
	DO lp=1,8 !numero de cubos
	  do jij=1,npxi !inicializo la estrella final cubo lp
	    cubstc(jij,lp)=0.0d0 !solo interesa flujo y NO lambda
	  enddo
	  aloeci(1,lp)=0.0d0
	  aloeci(2,lp)=0.0d0
	  aloeci(3,lp)=0.0d0
	  snat(lp)=0.0d0
	  sweigt(lp)=0.0d0
	  weigtc(lp)=0.0d0 !peso total de estrellas en cubo lp
	  sumnor=0.0d0
229	FORMAT(2(I2,1x),A21,1x,f7.0,1x,f4.2,1(1x,f5.2),1x,f5.3,1x,G20.10)
	  IF(nsc(lp).gt.0)THEN
	  DO lp2=1,nsc(lp) !numero de estrellas en cubo lp
	   weigc(lp2,lp)=0.0d0 !peso estrella lp2 en cubo lp
	   do is=1,nstari !estrellas libreria
	    if(ssc(lp2,lp).eq.stari(is))then
            if(sscp(lp2,4,lp).ge.snamax)then
	       snaton=1.0d0
            else
             snaton=(sscp(lp2,4,lp)*sscp(lp2,4,lp))/(snamax*snamax)
            endif
            weigc(lp2,lp)=exp(-(((sscp(lp2,1,lp)-tisst)/sigmat)**2.))
     &      *exp(-(((sscp(lp2,2,lp)-giss)/sigmag)**2.))     
     &      *exp(-(((sscp(lp2,3,lp)-fiss)/sigmaz)**2.))*snaton
            weigtc(lp)=weigtc(lp)+weigc(lp2,lp)
c para evitar numeros negativos en el espectro estelar debido a baja S/R
	     atonmi=9.e19
	     do jij=1,npxi
	      if(staris(is,jij,2).gt.0.0d0.and.
     &      staris(is,jij,2).lt.atonmi)atonmi=staris(is,jij,2)	     
	     enddo
c completada evaluacion valores negativos espectro estelar
	     do jij=1,npxi 
c para evitar valores negativos espectro estelar
	      if(staris(is,jij,2).le.0.0d0) staris(is,jij,2)=atonmi
c "if" para evitar dividir al final por pesos ultra0 cuando solo hay 1 * en cubo
	      if(nsc(lp).eq.1)then 
	       cubstc(jij,lp)=staris(is,jij,2)
	      else
	       cubstc(jij,lp)=cubstc(jij,lp)+weigc(lp2,lp)*
     &       staris(is,jij,2)
	      endif	  
	     enddo
	     if(nsc(lp).eq.1)then
            aloeci(1,lp)=sscp(lp2,1,lp)
	      aloeci(2,lp)=sscp(lp2,2,lp)
	      aloeci(3,lp)=sscp(lp2,3,lp)
            snat(lp)=snaton
	     else
            aloeci(1,lp)=aloeci(1,lp)+sscp(lp2,1,lp)*weigc(lp2,lp)
	      aloeci(2,lp)=aloeci(2,lp)+sscp(lp2,2,lp)*weigc(lp2,lp)
	      aloeci(3,lp)=aloeci(3,lp)+sscp(lp2,3,lp)*weigc(lp2,lp)
            snat(lp)=snat(lp)+snaton*weigc(lp2,lp)
	     endif	  
	    endif
	   enddo
c DESCOMENTAR PARA ESCUPIR LO QUE HACE PARA CADA ESTRELLA
c      write(555,229)lp,lp2,ssc(lp2,lp),t50/sscp(lp2,1,lp),
c     &(sscp(lp2,la,lp),la=2,3),snaton,weigc(lp2,lp)
	  ENDDO
	  if(weigtc(lp).lt.ultra0) weigtc(lp)=ultra0
	  if(nsc(lp).gt.1)then
          aloeci(1,lp)=aloeci(1,lp)/weigtc(lp)
          aloeci(2,lp)=aloeci(2,lp)/weigtc(lp)
          aloeci(3,lp)=aloeci(3,lp)/weigtc(lp)
	   snat(lp)=snat(lp)/weigtc(lp)
	  endif
c DESCOMENTAR PARA ESCUPIR LO QUE HACE PARA CADA ESTRELLA
c      write(555,*)'final=',lp,5040.d0/aloeci(1,lp),
c     &aloeci(2,lp),aloeci(3,lp),snat(lp),nsc(lp),weigtc(lp)
	  do jij=1,npxi
	   cubstc(jij,lp)=cubstc(jij,lp)/weigtc(lp)	 
	   if(ala(jij,1).ge.blc1.and.ala(jij,1).le.blc2)then
	    sumnor=sumnor+cubstc(jij,lp)
c esto habra que aplicarlo cuando se trate de normalizar
c la estrella final con el objetivo de integrar la SSP.
c Sin embargo para calcular la estrella resultante en base
c a las estrellas de los cubos uno no tiene porque proceder
c de esta forma.
c
c call respv(ala(jij,1),respon)
c sumnor=sumnor+respon*cubstc(jij,lp)*abs(ala(jij+1,1)-ala(jij,1))
	   endif	 
	  enddo
c	  estrella final del cubo lp renormalizada
	  sumnor=sumnor/xnumno !si normalizaramos a V se eliminaria 
	  do jij=1,npxi
	   cubstc(jij,lp)=cubstc(jij,lp)/sumnor
         enddo 
	  ENDIF
	ENDDO
c
c ESTRELLA FINAL NORMALIZADA
c
        na0=1 !medidor de weigtc=0, si no hubiese ninguno = 1
	DO lp=1,8 !weigtc genericos por cubo, los anteriores se anulan
	 IF(nsc(lp).gt.0)THEN
	  weigtc(lp)=exp(-(((aloeci(1,lp)-tisst)/sigmat)**2.))
     &  *exp(-(((aloeci(2,lp)-giss)/sigmag)**2.))
     &  *exp(-(((aloeci(3,lp)-fiss)/sigmaz)**2.))
	 ELSE
	  weigtc(lp)=0.0d0
	  na0=na0+1
	 ENDIF
         sweigt(lp)=weigtc(lp)*snat(lp) !error famoso al estar dentro de IF
        ENDDO 
c       WRITE(48,'(A2,1x,8(f13.10,1x))')'sw',(sweigt(lupu),lupu=1,8)
c       WRITE(48,'(A2,1x,8(I3,1x))')'ns',(nsc(lupu),lupu=1,8)
c Ajuste weigtc mas bajos con coeficientes sco para obtener * exacta
c Definimos dimensiones fisicas y logicas matrices: m=3 n=8 mp=8 np=8
      call indexx(8,sweigt,ilor) !ilor=indica orden creciente sweigt
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
      DO iakc=1,ilup 
C si todos weigtc's ne 0 ==> na0=1 y el caso general es para na0-1
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
	call dsvdcmp(uk,3,8,3,8,wk,vk,ntons) !resolvemos sist. ec. lin.
	if(ntons.eq.1) goto 5015
	wmax=0.0D0
	do jksvd=1,8
	 if(wk(jksvd).gt.wmax) wmax=wk(jksvd)
	enddo
	wmin=wmax*1.0d-12 !valor original
	do jksvd=1,8
	 if(wk(jksvd).lt.wmin) wk(jksvd)=0.0D0
	enddo
	call dsvbksb(uk,wk,vk,3,8,3,8,bk,xk)
c calculo solucion mas cercana a 1. 1. 1. 1. 1. 1. 1. 1.
c solucion sistema anterior: d_sol+lambda*vk(i,j) si wk(i)=0
	iiak=0
	do iak1=1,8
	  if(wk(iak1).lt.ultra0)then
	    iiak=iiak+1
	    do iak2=1,8
	      vvk(iiak,iak2)=vk(iak2,iak1)
	    enddo
	  endif
	enddo
c Bucle busqueda solucion mas proxima a 1 y si no vamos
c aumentando hasta llegar a una solucion positiva:
	do iak1=1,iiak
	  bbk(iak1)=0.0d0
	  do iak2=1,8
	    bbk(iak1)=bbk(iak1)+vvk(iak1,iak2)*(xk(iak2)-(1.0d0))
	  enddo
	  bbk(iak1)=-bbk(iak1) !faltaba el signo - en las formulas
	enddo
c calculo del vector solucion mas proximo a 11111111
5015	DO lpu=1,8
	  if(ntons.eq.1)then !esta opcion es por si falla dsvdcmp
	    sco(lpu)=1.0d0
	  else
	    sco(lpu)=xk(lpu)
	    do iak1u=1,iiak
	 	sco(lpu)=sco(lpu)+bbk(iak1u)*vvk(iak1u,lpu)
	    enddo
	    WE(lpu)=sco(lpu)*silor(lpu)*weigtc(lpu)
	  endif
	ENDDO
       abigMA=MAX(WE(1),WE(2),WE(3),WE(4),WE(5),WE(6),WE(7),WE(8))
       abigMI=MIN(WE(1),WE(2),WE(3),WE(4),WE(5),WE(6),WE(7),WE(8))
c    esto es por si la solucion es unica y vale 000       
       kolawn=0
       if(abigMA.lt.ultra0)then
	  kolawn=1
         abigMA=MAX(weigtc(1),weigtc(2),weigtc(3),weigtc(4),
     &              weigtc(5),weigtc(6),weigtc(7),weigtc(8))
         do lp=1,8
	     WE(lp)=weigtc(lp) !si multiplico por silor(lp) rompe programa
         enddo
       endif
c      Normalizacion de los pesos
       nwe=0
       do lp=1,8
         if(WE(lp).lt.ultra0)then
	     WE(lp)=0.0d0
	  else
	     WE(lp)=WE(lp)/abigMA
	     if(WE(lp).gt.0.001d0) nwe=nwe+1
	  endif
       enddo
c	if((kolawn.eq.1.and.nwe.gt.1).or.(kolawn.eq.0.and.nwe.gt.2))then
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
c	endif
       IF(abigMI.ge.0.0d0.and.kolawn.eq.0)THEN !no vale kolawn=1 al ser sol. 000
	  kolan=2
	  if(nwe.gt.2) goto 2537
       ELSE
         etaloe=abs(aloec(1)-tisst)
c  Para evitar perder solucion original en caso que abigMA>ultra0 y nwe<2	  
         if(etaloe.lt.taloec.and.nwe.gt.1)then 
	    kolan=0
	    taloec=etaloe
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
2537  continue
      IF(kolan.lt.2)THEN
	 weifin=weifi0
	 aloec(1)=aloec1
	 aloec(2)=aloec2
	 aloec(3)=aloec3
	 do lp=1,8
	   WE(lp)=WE0(lp)
	 enddo
      ENDIF
      IF((kolawn.eq.1.and.nwe.lt.2).or.(kolawn.eq.0.and.nwe.lt.3))THEN
	 kolawn=9
        aloec(1)=0.0d0
        aloec(2)=0.0d0
        aloec(3)=0.0d0
        weifin=0.0d0
        abigMA=MAX(weigtc(1),weigtc(2),weigtc(3),weigtc(4),
     &             weigtc(5),weigtc(6),weigtc(7),weigtc(8))
        do lp=1,8
	  WE(lp)=weigtc(lp)/abigMA
         weifin=weifin+WE(lp) 
	  aloec(1)=WE(lp)*aloeci(1,lp)+aloec(1)
	  aloec(2)=WE(lp)*aloeci(2,lp)+aloec(2)
	  aloec(3)=WE(lp)*aloeci(3,lp)+aloec(3)
        enddo
        aloec(1)=aloec(1)/weifin
        aloec(2)=aloec(2)/weifin
        aloec(3)=aloec(3)/weifin
      ENDIF 
      aloec(1)=t50/aloec(1) !paso a Teff
c calculamos espectro estrella final
      nwef=0
      do lp=1,8
	 if(WE(lp).gt.0.001d0)then
	  nwef=nwef+1
	  do jij=1,npxi
          ala(jij,2)=ala(jij,2)+WE(lp)*cubstc(jij,lp)
	  enddo
	 endif
      enddo
c fichero para testear lo que hace:
c      write(555,'(A3,1x,8(f11.9,1x))')'wei',(weigtc(lp),lp=1,8)
c      write(555,'(A3,1x,8(f11.9,1x))')'WEF',(WE(lp),lp=1,8)
c      write(555,*)'kolan=',kolan,' kolawn=',kolawn
c      write(555,*)'nwe=',nwe,' nwef=',nwef
c      write(555,*)'weifin=',weifin
c      if(weifin.le.1.0d0)WRITE(*,*)'WEIFIN MIUSC LT ZERO'
c      if(weifin.le.1.0d0)write(555,*)'WEIFIN LT ZERO'
c      if(nwef.eq.1) then
c      write(555,'(A)')'ALMOST 1 BOX WAS USED FOR THE INTERPOLATION'
c      elseif(nwef.eq.2)then
c      write(555,'(A)')'ALMOST 2 BOXES WERE USED FOR THE INTERPOLATION'
c      endif
C  Calcular volumenes calidad modelo  
      volpai=0.d0 !volumen escalado a sigma minimo
      do lp=1,8
       nbxstt=nbxstt+nsc(lp)
       volpai=volpai+sqrt((sigbox(lp,1)/(xs*ddt))**2.+
     & (sigbox(lp,2)/(xs*ddg))**2.+(sigbox(lp,3)/(xs*ddz))**2.)
      enddo
      xnbxsi=dble(nbxstt)
      voli10=0.d0 !volumen maximo aceptable (sigma=1.0)
      voli10=8.0d0*sqrt(((re10*ddtmax)/(xs*ddt))**2.+
     & ((re10*ddgmax)/(xs*ddg))**2.+((re10*ddzmax)/(xs*ddz))**2.) 
      voli05=0.d0 !volumen maximo aceptable (sigma=0.5)
      voli05=8.0d0*sqrt(((re05*ddtmax)/(xs*ddt))**2.+
     & ((re05*ddgmax)/(xs*ddg))**2.+((re05*ddzmax)/(xs*ddz))**2.) 
      voli15=0.d0 !volumen maximo aceptable (sigma=1.5)
      voli15=8.0d0*sqrt(((re15*ddtmax)/(xs*ddt))**2.+
     & ((re15*ddgmax)/(xs*ddg))**2.+((re15*ddzmax)/(xs*ddz))**2.) 
      RETURN
      END	
