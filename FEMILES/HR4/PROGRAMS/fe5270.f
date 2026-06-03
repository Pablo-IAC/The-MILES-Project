c*******SUBRUTINA DEL Fe5270***********************************************
c        SUBROUTINE fe5270(nra,i,W,stars)
        SUBROUTINE fe5270(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,19)/ra(1,19))
	spend=(5270.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 00014
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-24.153+1.6448e-2*at-2.3709e-6*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=13.998-4.1814*R(i,j,1)+0.38014*R(i,j,1)**2
			else
				goto 40014
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=3.4111+17.2923*t-596.5421*t**3+1.2899*feh(i)-.5078*R(i,j,4)+
     &  .1066*R(i,j,4)**2
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=-.8921+3.6501*t*t+.8108*feh(i)
		else
			goto 40014
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
40014	continue
	ENDDO	
00014   continue
	return
	end
