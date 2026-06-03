c*******SUBRUTINA DE Fe4668************************************************
c        SUBROUTINE fe4668(nra,i,W,stars)
        SUBROUTINE fe4668(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,3)/ra(1,3))
	spend=(4668.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 50006
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-84.851+8.1794e-2*at-1.4932e-5*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=2.8569-1.9164*R(i,j,1)+0.47471*R(i,j,1)**2
			else
				goto 10006
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=154.757*t-1508.7241*t*t+7517.5278*t**3+7.281*feh(i)+
     &  2.9184*feh(i)**2+0.4048*feh(i)**3+3.3396*R(i,j,4)
     &  -0.4802*R(i,j,4)**2-25.0169*t*R(i,j,4)
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=3.243-19.6869*t+21.6173*t*t+1.1699*feh(i)**2+
     &  4.9505*feh(i)*t*t
		else
			goto 10006
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
10006	continue
	ENDDO	
50006   continue
	return
	END
