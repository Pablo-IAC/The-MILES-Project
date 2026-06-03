c*******SUBRUTINA DEL Mgb**************************************************
c        SUBROUTINE mgb(nra,i,W,stars)
        SUBROUTINE mgb(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,8)/ra(1,8))
	spend=(5177.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40013
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-60.119+5.9098e-2*at-1.1004e-5*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=14.228-5.0136*R(i,j,1)+0.62793*R(i,j,1)**2
			else
				goto 00013
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=-1.1228+46.2019*t+681.4774*t**3+.7789*feh(i)-.0607*feh(i)**3
     &  +1.3514*R(i,j,4)+10.1217*t*feh(i)-115.2988*R(i,j,4)*t*t
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.3.8131))then
	t=10**t
	ZMAG=8.7334+.8537*feh(i)-4.3748*R(i,j,4)+.7815*t*R(i,j,4)**2
		elseif((R(i,j,3).gt.3.8131).and.(R(i,j,3).le.4.1225))then
	t=10**t	
	ZMAG=-0.3996+3.3480*t**3
		else
			goto 00013
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
00013	continue
	ENDDO	
40013   continue
	return
	END
