c*******SUBRUTINA DE HgammaA************************************************
c        SUBROUTINE hgamaa(nra,i,W,stars)
        SUBROUTINE hgamaa(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,1)/ra(1,1))
	spend=(4341.625-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 4042
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	tetas=5040./at
	IF((lc.eq.5).or.(lc.eq.3))THEN
	  if(R(i,j,3).gt.3.826075.and.R(i,j,3).le.4.13)then
	ZMAG=-46.529+213.84*tetas-193.0285*tetas*tetas
	  elseif(R(i,j,3).ge.3.7076.and.R(i,j,3).le.3.826075)then
	ZMAG=99.846-180.61*tetas+74.564*tetas*tetas
     &	-0.02066*R(i,j,4)**3-2.56*tetas*tetas*feh(i)
	  elseif(R(i,j,3).lt.3.7076)then
		if(R(i,j,3).lt.3.5798.and.lc.eq.3)then
	ZMAG=-35.04+0.02845*at-5.728e-6*at*at
	  	elseif(R(i,j,3).lt.3.6335.and.lc.eq.5)then
	ZMAG=-8.215+8.556e-3*at-2.118e-6*at*at
	  	else
	ZMAG=-7.3344-47.269*t+1565.6*t*t*t-3.7692*feh(i)
     &	+0.5623*feh(i)*R(i,j,4)
	  	endif
	  else
	  	goto 5042
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
5042	continue
	ENDDO
4042    continue
	return
	END
