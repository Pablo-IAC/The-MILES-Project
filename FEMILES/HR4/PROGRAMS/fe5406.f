c*******SUBRUTINA DEL Fe5406***********************************************
c        SUBROUTINE fe5406(nra,i,W,stars)
        SUBROUTINE fe5406(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,19)/ra(1,19))
	spend=(5406.-4400.)/(5500.-4400.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40016
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5824)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=16.020-1.3358e-2*at+2.6283e-6*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=7.8350-1.9903*R(i,j,1)+0.13898*R(i,j,1)**2
			else
				goto 00016
			endif
	elseif((R(i,j,3).ge.3.5824).and.(R(i,j,3).le.3.7076))then
	ZMAG=1.3752+14.5727*t+1.0494*feh(i)-0.0584*feh(i)**3
     &  +.0211*R(i,j,4)**2
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=-0.4650+1.9841*t**3+0.2312*t*feh(i)**2+0.9676*t*t*feh(i)
	else
			goto 00016
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
00016	continue
	ENDDO	
40016   continue
	return
	end
