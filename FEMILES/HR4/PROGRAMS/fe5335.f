c*******SUBRUTINA DEL Fe5335***********************************************
c        SUBROUTINE fe5335(nra,i,W,stars)
        SUBROUTINE fe5335(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,19)/ra(1,19))
	spend=(5335.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40015
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5824)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-6.9412+5.2087e-3*at-6.1016e-7*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=13.303-3.8264*R(i,j,1)+0.32262*R(i,j,1)**2
			else
				goto 00015
			endif
		elseif((R(i,j,3).ge.3.5824).and.(R(i,j,3).le.3.7076))then
	ZMAG=1.9450+19.2065*t+1.6282*feh(i)+.2817*feh(i)**2
     &  +.0476*R(i,j,4)**2
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=-1.0386+3.4925*t*t+.3224*feh(i)**2+1.3286*t*feh(i)
		else
			goto 00015
		endif
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
00015	continue
	ENDDO	
40015   continue
	return
	end
