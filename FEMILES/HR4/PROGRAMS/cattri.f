c*******SUBRUTINA CaIIT************************************************
        SUBROUTINE cattri(i,WCATS,WPAT,WCAT,WSTIO,WMGI
     &  ,stars1,stars2,stars3,stars4,stars5,t2,g2,f2)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
	common/cfnc/fnc11(2999),fnmg(2999)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c     feh(i)=fez(i)+dlog10(ra(nra,15)/ra(1,15))
	WCATS=0.0d0
	WPAT=0.0d0
	WCAT=0.0d0
	WSTIO=0.0d0
	WMGI=0.0d0
	stars1=0.0d0
	stars2=0.0d0
	stars3=0.0d0
	stars4=0.0d0
	stars5=0.0d0
	t2=0.0d0
	f2=0.0d0
	g2=0.0d0
	if(jotai(i).eq.0)goto 40491
	DO j=1,jotai(i)
	 at=10**R(i,j,3)
	 tetas=5040.d0/at
	 ggrav=R(i,j,4)
	 zfeff=feh(i)
         call nIR(tetas,ggrav,zfeff,ZCATS,ecats,ZPAT,epat,ZCAT,
     +   ecat,ZSTIO,estio,ZMGI,emgi,iflag,iflag2,iflag3,iflag4)
	 vikc=(-2.5d0)*dlog10(R(i,j,7)/R(i,j,9))
	 if(vikc.le.1.5d0)then
	    vij=-0.005d0+1.273d0*vikc
	 else
	    vij=0.723d0+0.486d0*vikc+0.215d0*vikc*vikc
	 endif
	 flujv=10**((-0.4d0)*((-2.5d0)*dlog10(R(i,j,7))-3.762d0))
	 flujcj=(flujv*(10**(+0.4d0*(vij))))
c MgI todo   frxmgi=.044519692 -.0041583698 *R(m,i,3)
c MgI d      frxmgi=.43927312  -.218387380  *R(m,i,3)+.028949539*R(m,i,3)*R(m,i,3)
c MgI g      frxmgi=.96993812  -.497936560  *R(m,i,3)+.065761539*R(m,i,3)*R(m,i,3)
c
c sTiO todo  frxsTiO= .0051969893-.00016500045*R(m,i,3)
c sTiO d     frxsTiO=-.38456164  +.213282890  *R(m,i,3)-.029221114*R(m,i,3)*R(m,i,3)
c sTiO g     frxsTiO=-.097070614 +.028402666  *R(m,i,3)
c
c CaT todo   frxCaT=.15695200  -.0033235712 *R(m,i,3)
c CaT d      frxCaT=-.16346262  +.0868521340 *R(m,i,3)
c CaT g      frxCaT=-1.62067540  +.4969163100 *R(m,i,3)
	 IF(R(i,j,3).lt.3.75d0)THEN
	  if(ggrav.ge.3.5d0)then
       frxmgi=.43927312d0-.218387380d0*R(i,j,3)+
     & .028949539d0*R(i,j,3)*R(i,j,3)
	  else
       frxmgi=.96993812d0-.497936560d0*R(i,j,3)+
     & .065761539d0*R(i,j,3)*R(i,j,3)
	  endif
	 ELSEIF(R(i,j,3).ge.3.75d0.and.R(i,j,3).lt.3.90d0)THEN	   
       frxmgi=.96993812d0-.497936560d0*R(i,j,3)+
     & .065761539d0*R(i,j,3)*R(i,j,3)
	 ELSE
       frxmgi=.044519692d0-.0041583698d0 *R(i,j,3)
	 ENDIF
	 IF(R(i,j,3).lt.3.57d0)THEN
	  if(ggrav.ge.3.5d0)then
       frxstio=-.38456164d0+.213282890d0*R(i,j,3)-
     & .029221114d0*R(i,j,3)*R(i,j,3)
	  else
       frxstio=-.097070614d0+.028402666d0*R(i,j,3) 
	  endif
	 ELSEIF(R(i,j,3).ge.3.57d0.and.R(i,j,3).lt.3.62d0)THEN 	   
       frxstio=-.38456164d0+.213282890d0*R(i,j,3)-
     & .029221114d0*R(i,j,3)*R(i,j,3) 
	 ELSE
       frxstio=.0051969893d0-.00016500045d0*R(i,j,3) 
	 ENDIF
c Vazdekis et al. 2003!
	IF(R(i,j,3).lt.3.55d0)THEN
	 if(ggrav.ge.3.5d0)then
	  frxcat=-1.62067540d0+.4969163100d0*R(i,j,3) 
	 else
	  frxcat=-.16346262d0+.0868521340d0*R(i,j,3) 
	 endif
	ELSE
	  frxcat=.15695200d0-.0033235712d0*R(i,j,3)
	ENDIF
	 flujc=(frxcat*flujcj)/(8807.d0-8475.d0)
	 flujst=fnc11(j)*frxcat*flujcj
	 flujmg=fnmg(j)*frxcat*flujcj
        if(j.eq.1)then
	  if(imam(i).gt.0)then
	   do mma=1,imam(i)
	 WCATS=WCATS+ZCATS*R140(i,mma)*flujc
	 WPAT=WPAT+ZPAT*R140(i,mma)*flujc
	 WCAT=WCAT+ZCAT*R140(i,mma)*flujc
	 WSTIO=WSTIO+ZSTIO*R140(i,mma)*flujst
	 WMGI=WMGI+ZMGI*R140(i,mma)*flujmg
	 stars1=stars1+R140(i,mma)*flujc
	 stars2=stars2+R140(i,mma)*flujc
	 stars3=stars3+R140(i,mma)*flujc
	 stars4=stars4+R140(i,mma)*flujst
	 stars5=stars5+R140(i,mma)*flujmg
	 t2=t2+at*R140(i,mma)*flujc
	 g2=g2+ggrav*R140(i,mma)*flujc
	 f2=f2+zfeff*R140(i,mma)*flujc
	   enddo
	  endif
	 WCATS=WCATS+ZCATS*R(i,j,14)*flujc
	 WPAT=WPAT+ZPAT*R(i,j,14)*flujc
	 WCAT=WCAT+ZCAT*R(i,j,14)*flujc
	 WSTIO=WSTIO+ZSTIO*R(i,j,14)*flujst
	 WMGI=WMGI+ZMGI*R(i,j,14)*flujmg
	 stars1=stars1+R(i,j,14)*flujc
	 stars2=stars2+R(i,j,14)*flujc
	 stars3=stars3+R(i,j,14)*flujc
	 stars4=stars4+R(i,j,14)*flujst
	 stars5=stars5+R(i,j,14)*flujmg
	 t2=t2+at*R(i,j,14)*flujc
	 g2=g2+ggrav*R(i,j,14)*flujc
	 f2=f2+zfeff*R(i,j,14)*flujc
	else
	 WCATS=WCATS+ZCATS*R(i,j,14)*flujc
	 WPAT=WPAT+ZPAT*R(i,j,14)*flujc
	 WCAT=WCAT+ZCAT*R(i,j,14)*flujc
	 WSTIO=WSTIO+ZSTIO*R(i,j,14)*flujst
	 WMGI=WMGI+ZMGI*R(i,j,14)*flujmg
	 stars1=stars1+R(i,j,14)*flujc
	 stars2=stars2+R(i,j,14)*flujc
	 stars3=stars3+R(i,j,14)*flujc
	 stars4=stars4+R(i,j,14)*flujst
	 stars5=stars5+R(i,j,14)*flujmg
	 t2=t2+at*R(i,j,14)*flujc
	 g2=g2+ggrav*R(i,j,14)*flujc
	 f2=f2+zfeff*R(i,j,14)*flujc
	endif
	ENDDO
40491   continue
	return
	END
