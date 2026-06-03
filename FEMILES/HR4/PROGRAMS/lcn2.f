c*******SUBRUTINA DEL CN2**************************************************
c        SUBROUTINE lcn2(nra,i,W,stars)
        SUBROUTINE lcn2(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,3)/ra(1,3))
	ancho=35.0
	spend=(4161.-3600.)/(4400.-3600.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40002
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		IF(feh(i).lt.-1.0)THEN
		if((R(i,j,3).le.3.7076).and.(R(i,j,3).ge.3.60))then
	ZMAG=0.170592*R(i,j,3)**2-1.96097e-2*R(i,j,4)**2-
     &  0.186713*feh(i)**2-.192425*R(i,j,3)-5.203332e-2*R(i,j,4)
     &  -10.52774*feh(i)+7.33773e-2*R(i,j,3)*R(i,j,4)+2.65265*R(i,j,3)
     &  *feh(i)+5.41124e-2*R(i,j,4)*feh(i)-2.5165
		elseif((R(i,j,3).le.4.045).and.(R(i,j,3).gt.3.7076))then
	ZMAG=0.187522*R(i,j,3)**2-1.48142e-2*R(i,j,4)**2
     &  +1.69581e-2*feh(i)**2-4.24338*R(i,j,3)-2.62516*R(i,j,4)-0.478391
     &  *feh(i)+0.702715*R(i,j,3)*R(i,j,4)+
     &  0.201885*R(i,j,3)*feh(i)-4.85303e-2*R(i,j,4)*feh(i)+13.5325
		else
			goto 00002
		endif
		ELSE
		if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=1.3614-1.2333e-3*at+2.4043e-7*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=0.5443-0.2566*R(i,j,1)+0.02273*R(i,j,1)**2
			else
				goto 00002
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=.1901+5.2625*t-39.4219*t*t+.2873*feh(i)+.206*feh(i)**2+
     &  .0535*feh(i)**3-.0062*R(i,j,4)**2-.4206*t*R(i,j,4)
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=0.8763-4.6267*t+6.5272*t*t-2.7719*t**3+0.0222*feh(i)*t*t
		else
			goto 00002
		endif
		ENDIF
		ZMAG=ancho*(1.0d0-10**((-0.4d0)*ZMAG))
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
00002	continue
	ENDDO
40002   continue
	return
	END
