c*******SUBRUTINA DEL Mg1**************************************************
c        SUBROUTINE mg1(nra,i,W,stars)
        SUBROUTINE mg1(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
      COMMON/MAMA/R140(12,99000),imam(12)
	feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,8)/ra(1,8))
	spend=(5101.5-4400.)/(5500.-4400.)
	ancho=65.0
	W=0.0
	stars=0.0
	if(jotai(i).eq.0)goto 40011
	DO j=1,jotai(i)
	lc=int(R(i,j,13)+0.15)
	t=3.70243-R(i,j,3)
	at=10**R(i,j,3)
	IF((lc.eq.5).or.(lc.eq.3))THEN
	if(R(i,j,3).lt.3.5527)then
			if((lc.le.3).and.(R(i,j,3).ge.3.3617))then
	ZMAG=4.3350-3.0777e-3*at+5.2561e-7*at**2
			elseif((lc.ge.4).and.(R(i,j,1).le.7.0))then
	ZMAG=0.5806-0.0863*R(i,j,1)+0.00524*R(i,j,1)**2
			else
				goto 00011
			endif
	elseif((R(i,j,3).ge.3.5527).and.(R(i,j,3).le.3.7076))then
	ZMAG=0.0345-2.8343*t+75.3661*t*t-327.9276*t**3+0.0415*feh(i)
     &  -0.0050*feh(i)**3+0.0027*R(i,j,4)**2+
     &  0.5013*t*feh(i)+1.1371*t*R(i,j,4)-8.1113*R(i,j,4)*t*t
	elseif((R(i,j,3).gt.3.7076).and.(R(i,j,3).le.4.1225))then
	t=10**t
	ZMAG=.0536-.4697*t*t+.4868*t**3-.0815*t*feh(i)+.1122*t*t*feh(i)
	else
			goto 00011
	endif
		ZMAG=ancho*(1.0d0-10**((-0.4d0)*ZMAG))
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
00011	continue
	ENDDO	
40011   continue
	return
	END
