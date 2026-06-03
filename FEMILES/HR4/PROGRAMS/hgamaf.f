c*******SUBRUTINA DE Hgammaf************************************************
c        SUBROUTINE hgamaf(nra,i,W,stars)
        SUBROUTINE hgamaf(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,1)/ra(1,1))
	spend=(4341.75-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 4040
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	tetas=5040./at
	IF((lc.eq.5).or.(lc.eq.3))THEN
	  if(R(i,j,3).gt.3.826075.and.R(i,j,3).le.4.13)then
	ZMAG=-20.397+106.87*tetas-98.0065*tetas*tetas-0.2997*feh(i)
	  elseif(R(i,j,3).ge.3.7076.and.R(i,j,3).le.3.826075)then
	ZMAG=44.554-59.657*tetas+13.508*tetas*tetas*tetas
     &	-0.01438*R(i,j,4)**3-0.954*tetas*tetas*feh(i)
	  elseif(R(i,j,3).lt.3.7076)then
		if(R(i,j,3).lt.3.5798.and.lc.eq.3)then
	ZMAG=-9.007+5.729e-3*at-1.049e-6*at*at
	  	elseif(R(i,j,3).lt.3.6335.and.lc.eq.5)then
	ZMAG=-14.155+7.846e-3*at-1.242e-6*at*at
	  	else
	ZMAG=-1.3462-24.41*t+1100.3*t*t*t-1.3551*feh(i)
     &	-0.2741*R(i,j,4)+0.1976*feh(i)*R(i,j,4)
	  	endif
	  else
	  	goto 5040
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
5040	continue
	ENDDO
4040    continue
	return
	END
