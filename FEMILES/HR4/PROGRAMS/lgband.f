c*******SUBRUTINA DE G-band************************************************
c        SUBROUTINE lgband(nra,i,W,stars)
        SUBROUTINE lgband(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,3)/ra(1,3))
	spend=(4300.-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40004
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-11.030+2.9856e-3*at+4.1073e-7*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=12.202-3.3556*R(i,j,1)+0.27811*R(i,j,1)**2
			else
				goto 00004
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=6.1924+33.0914*t-220.4816*t*t-0.6158*feh(i)**2
     &  -5.5839*R(i,j,4)*t
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=40.874-210.2997*t+311.2509*t*t-135.3966*t**3+1.74*t*feh(i)
		else
			goto 00004
		endif
        if(j.eq.1)then
	  if(imam(i).gt.0)then
	   do mma=1,imam(i)
	 W=W+ZMAG*R140(i,mma)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	 stars=stars+R140(i,mma)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	   enddo
	  endif
	 W=W+ZMAG*R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	else
	 W=W+ZMAG*R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
	endif
	ENDIF
00004	continue
	ENDDO	
40004   continue
	return
	END
