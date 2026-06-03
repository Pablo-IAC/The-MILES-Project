c*******SUBRUTINA DE Hbeta************************************************
c        SUBROUTINE hbeta(nra,i,W,stars)
        SUBROUTINE hbeta(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,1)/ra(1,1))
	spend=(4861.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40009
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-36.196+2.9880e-2*at-5.2742e-6*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=15.411-7.7717*R(i,j,1)+0.88109*R(i,j,1)**2
			else
				goto 00009
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=2.1471-307.5695*t*t+2075.374*t**3+.2079*feh(i)-
     &  .1680*R(i,j,4)-.5883*t*R(i,j,4)**2
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.3.8131))then
	t=10**t
	ZMAG=29.0304-47.1204*t+19.8378*t*t-.0038*R(i,j,4)**3
		elseif((R(i,j,3).gt.3.8131).and.(R(i,j,3).le.4.1225))then
	t=10**t	
	ZMAG=-48.3032+260.1385*t-373.5425*t**2+161.9056*t**3
		else
			goto 00009
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
00009	continue
	ENDDO	
40009   continue
	return
	END
