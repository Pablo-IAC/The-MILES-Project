c*******SUBRUTINA DE Fe4383************************************************
c        SUBROUTINE fe4383(nra,i,W,stars)
        SUBROUTINE fe4383(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,3)/ra(1,3))
	spend=(4383.-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40005
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=216.13-1.4386e-1*at+2.3497e-5*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=27.190-9.0343*R(i,j,1)+0.82927*R(i,j,1)**2
			else
				goto 00005
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=2.9642+95.347*t-377.3676*t*t+3.2305*feh(i)-0.1731*feh(i)**3
     &  +0.9175*R(i,j,4)-12.4649*t*R(i,j,4)
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.3.8131))then
	t=10**t
	ZMAG=-15.2137+21.2603*t+0.7105*feh(i)**2+2.9105*t*feh(i)
		elseif((R(i,j,3).gt.3.8131).and.(R(i,j,3).le.4.1225))then
	t=10**t	
	ZMAG=11.4367-37.4777*t+36.3357*t*t+0.9215*feh(i)-0.8214*R(i,j,4)
		else
			goto 00005
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
00005	continue
	ENDDO	
40005   continue
	return
	END
