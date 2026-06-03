c*******SUBRUTINA DE HdeltaA************************************************
c        SUBROUTINE hdelta(nra,i,W,stars)
        SUBROUTINE hdelta(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,1)/ra(1,1))
	spend=(4102.875-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 4043
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	tetas=5040./at
	IF((lc.eq.5).or.(lc.eq.3))THEN
	  if(R(i,j,3).gt.3.826075.and.R(i,j,3).le.4.13)then
	ZMAG=-21.3+87.328*tetas-89.3136*tetas*tetas*tetas
	  elseif(R(i,j,3).ge.3.7076.and.R(i,j,3).le.3.826075)then
	ZMAG=35.982-39.599*tetas-0.4963*feh(i)*feh(i)
     &	-0.01241*R(i,j,4)**3-2.8349*tetas*tetas*feh(i)
	  elseif(R(i,j,3).lt.3.7076)then
		if(R(i,j,3).lt.3.5798.and.lc.eq.3)then
	ZMAG=-44.638+0.02945*at-4.949e-6*at*at
	  	elseif(R(i,j,3).lt.3.6335.and.lc.eq.5)then
	ZMAG=-9.236+0.01074*at-2.312e-6*at*at
	  	else
	ZMAG=-4.2019-70.59*t+382.65*t*t-4.8387*feh(i)
     &	+0.1443*feh(i)*R(i,j,4)*R(i,j,4)+67.602*t*t*R(i,j,4)
	  	endif
	  else
	  	goto 5043
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
5043	continue
	ENDDO
4043    continue
	return
	END

