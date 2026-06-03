c*******SUBRUTINA DEL Fe5782***********************************************
c        SUBROUTINE fe5782(nra,i,W,stars)
        SUBROUTINE fe5782(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,19)/ra(1,19))
	spend=(5782.-5500.)/(6400.-5500.)
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40018
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=17.717-1.3024e-2*at+2.2894e-6*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=4.1657-1.2528*R(i,j,1)+0.09038*R(i,j,1)**2
			else
				goto 00018
			endif
	elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
      ZMAG=.7642+10.1171*t-345.4548*t**3+.5284*feh(i)
     &-.1744*t*R(i,j,4)**2
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=0.6296*t**3+0.2159*feh(i)
	else
			goto 00018
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
00018	continue
	ENDDO
40018   continue
	return
	end
