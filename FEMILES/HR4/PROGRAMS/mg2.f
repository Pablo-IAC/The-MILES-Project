c*******SUBRUTINA DEL Mg2**************************************************
c        SUBROUTINE mg2(nra,i,W,stars)
        SUBROUTINE mg2(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,8)/ra(1,8))
	spend=(5177.-4400.)/(5500.-4400.)
	ancho=42.5
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40012
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=1.0425-4.4418e-4*at+7.7682e-8*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=0.9232-0.2130*R(i,j,1)+0.02454*R(i,j,1)**2
			else
				goto 00012
			endif
	elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=.0762+47.2005*t*t-146.1148*t**3+.0836*feh(i)-
     &  .0077*feh(i)**3+.0102*R(i,j,4)**2+.8210*t*feh(i)+
     &  1.0168*t*R(i,j,4)-12.0608*R(i,j,4)*t*t
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.3.8131))then
	t=10**t
	ZMAG=1.2804-1.1895*t+.0468*feh(i)-.2596*R(i,j,4)+
     &  .0017*R(i,j,4)**3+.2622*R(i,j,4)*t*t
	elseif((R(i,j,3).gt.3.8131).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=0.1212*t**3+.0105*feh(i)
	else
			goto 00012
	endif
		ZMAG=ancho*(1.0d0-10**((-0.4d0)*ZMAG))
        if(j.eq.1)then
	  if(imam(i).gt.0)then
	   do mma=1,imam(i)
	 W=W+ZMAG*R140(i,mma)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	 stars=stars+R140(i,mma)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	   enddo
	  endif
	 W=W+ZMAG*R(i,j,14)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	else
	 W=W+ZMAG*R(i,j,14)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,6)+(R(i,j,7)-R(i,j,6))*spend)
	endif
	ENDIF
00012	continue
	ENDDO	
40012   continue
	return
	END
