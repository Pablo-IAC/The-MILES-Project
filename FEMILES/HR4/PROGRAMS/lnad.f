c*******SUBRUTINA DE NaD***************************************************
c        SUBROUTINE lnad(nra,i,W,stars)
        SUBROUTINE lnad(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,7)/ra(1,7))
c	El filtro en R es Cousins
	spend=(5895.-5500.)/(6400.-5500.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40019
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5527)then
			if((lc.eq.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=57.265-3.0626e-2*at+4.4396e-6*at**2
			elseif((lc.eq.5).and.(R(i,j,1).le.7.0))then
	ZMAG=3.6794+1.9365*R(i,j,1)-0.07710*R(i,j,1)**2
			else
				goto 00019
			endif
	elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=1.6582+60.3902*t-752.5432*t**3+2.2059*feh(i)+
     &  .7747*feh(i)**2+.0171*R(i,j,4)**3+15.8952*t*feh(i)-
     &  28.7499*t*R(i,j,4)+6.2953*t*R(i,j,4)**2
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=2.3060*t**3-0.670*feh(i)+1.5383*feh(i)*t*t
	else
		goto 00019
	endif
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
00019	continue
	ENDDO
40019   continue	
	return
	END
