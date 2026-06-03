c      character*42 tonto
      character*6 miles !toda la linea de la estrella MILES
      character*14 hd
      character*1 s
      REAL*4 a
      PARAMETER (nmil=1154)
      dimension a(nmil,21),miles(nmil),hd(nmil)
      cero=0.0
      f=2.0
      dmint=60.*f
      dming=0.18*f
      dminf=0.09*f
44    FORMAT(A6,9(1X,F8.2),2X,A1,1X,A14)
      open(99,file='PARAM_MILES',status='old')
      DO i=1,999999
        read(99,*,end=322)miles(i),(a(i,j),j=1,21),hd(i)
c	write(*,*)miles(i),a(i,1),a(i,2),a(i,3),a(i,18),a(i,19),
c     &a(i,20),a(i,21)
        adt=0.0
	adg=0.0
	adf=0.0
	s='-'
        IF(a(i,21).gt.0.9)THEN
	 s='s'
	 dt=ABS(a(i,1)-a(i,18))
	 dg=ABS(a(i,2)-a(i,19))
	 df=ABS(a(i,3)-a(i,20))
	 if(dt.gt.dmint)adt=dt
	 if(dg.gt.dming)adg=dg
	 if(df.gt.dminf)adf=df
	ENDIF
	write(93,44)miles(i),a(i,1),a(i,2),a(i,3),a(i,18),a(i,19),
     &a(i,20),adt,adg,adf,s,hd(i)
      ENDDO
322   close(99)
      END
