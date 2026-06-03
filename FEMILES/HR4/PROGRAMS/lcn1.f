c*******SUBRUTINA DEL CN1**************************************************
c        SUBROUTINE lcn1(nra,i,W,stars)
        SUBROUTINE lcn1(i,W,stars)
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
	if(jotai(i).eq.0)goto 40001
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
		IF(feh(i).lt.-1.0)THEN
		if((R(i,j,3).le.3.7076).and.(R(i,j,3).ge.3.60))then
	ZMAG=(-2.62644)*R(i,j,3)**2-4.27985e-2*R(i,j,4)**2+0.42298*
     &  feh(i)**2+4.52476*R(i,j,3)-7.66026*R(i,j,4)+16.2199*feh(i)
     &  +2.11758*R(i,j,3)*R(i,j,4)-3.9428*R(i,j,3)*feh(i)-
     &  0.118944*R(i,j,4)*feh(i)+19.9478
		elseif((R(i,j,3).le.4.045).and.(R(i,j,3).gt.3.7076))then
	ZMAG=(-5.58056e-2)*R(i,j,3)**2-8.90992e-3*R(i,j,4)**2+
     &  8.17567e-2*feh(i)**2-2.95409*R(i,j,3)-2.85048*R(i,j,4)+
     &  7.01824e-2*feh(i)+.766979*R(i,j,3)*R(i,j,4)+
     &  9.45437e-2*R(i,j,3)*feh(i)-2.0301e-2*R(i,j,4)*feh(i)+12.1739
		else
			goto 00001
		endif
		ELSE
		if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=2.0379-1.6262e-3*at+2.931e-7*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=0.1234-0.0899*R(i,j,1)+0.00689*R(i,j,1)**2
			else
				goto 00001
			endif
		elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=.1723+3.766*t-33.1882*t*t+.2709*feh(i)+.1888*feh(i)**2+
     &  .0467*feh(i)**3-.0066*R(i,j,4)**2-.0622*t*R(i,j,4)**2
		elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=1.3308-7.2601*t+10.5143*t*t-4.6009*t**3+.0204*feh(i)
		else
			goto 00001
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
00001	continue
	ENDDO
40001   continue
	return
	END
