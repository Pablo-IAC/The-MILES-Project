c*******SUBRUTINA DEL TiO2*************************************************
c        SUBROUTINE ltio2(nra,i,W,stars)
        SUBROUTINE ltio2(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,5)/ra(1,5))
	spend=(6232.6-5500.)/(6400.-5500.)
	ancho=82.5
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40021
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=-3.3671+3.2580e-3*at-5.9700e-7*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=-1.0903+.4335*R(i,j,1)-0.02141*R(i,j,1)**2
			else
				goto 00021
			endif
	elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=-.5837+exp(-.5296+286.8986*t**3+.3643*t*R(i,j,4)+
     &  10.6302*feh(i)*t*t-7.5956*R(i,j,4)*t*t) 
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=-.0115+.0144*t*t
	else
			goto 00021
	endif
		ZMAG=ancho*(1.0d0-10**((-0.4d0)*ZMAG))
        if(j.eq.1)then
	  if(imam(i).gt.0)then
	   do mma=1,imam(i)
	 W=W+ZMAG*R140(i,mma)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	 stars=stars+R140(i,mma)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	   enddo
	  endif
	 W=W+ZMAG*R(i,j,14)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	else
	 W=W+ZMAG*R(i,j,14)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	 stars=stars+R(i,j,14)*(R(i,j,7)+(R(i,j,8)-R(i,j,7))*spend)
	endif
	ENDIF
00021	continue
	ENDDO
40021   continue
	return
	END
