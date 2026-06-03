      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 miuk
5478	write(*,*)'MILES(m),MIUSC(i),NGSL(u),IRTF(k)?'
	read(*,*)miuk
	IF(miuk.eq.'m'.or.miuk.eq.'M')THEN
	  miukn=1
	  npxm=4300
	  nstm=1999
c	  write(*,*)'IN ADDITION IT WILL BE COMPUTED: B,R & C'
c	  dirout='ls ./OUT'
c	  mkout='mkdir ./OUT ./OUT/B ./OUT/C ./OUT/R ./OUT/M'
	ELSEIF(miuk.eq.'i'.or.miuk.eq.'I')THEN
	  miukn=2
	  npxm=6672
	  nstm=1999
	ELSEIF(miuk.eq.'u'.or.miuk.eq.'U')THEN
	  miukn=3
	  npxm=6192
	  nstm=1100
	ELSEIF(miuk.eq.'k'.or.miuk.eq.'K')THEN
	  miukn=4
	  npxm=15387
	  nstm=300
	ELSE
	  goto 5478
	ENDIF
	write(*,*)npxm,nstm
	call ano(npxm,nstm)
	end
	
	subroutine ano(npxm,nstm)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	INTEGER a
	DIMENSION a(npxm,nstm)
	a(1,1)=1.
	a(npxm,nstm)=99.
	write(*,*)a(1,1),a(npxm,nstm)
c	nnn=npxm*nstm
	return
	end
