	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	fimm=0.d0
	ZMU=2.3d0
	xxmasa=0.0d0
	dxm=0.1d0
      xnorm=4.d0/( 2.653238076d0+(1.d0-(.2d0)**(-ZMU))/(1.d0-ZMU))
      write(*,*)xnorm
	xnorm=1.0d0
      DO i=1,20
	 fimm=0.0d0
	 fimd=0.0d0
	 xxmasa=xxmasa+dxm
       x=xxmasa/0.5d0
	 xx=0.5**(-2.3)
 	 fimm=xxmasa**(-ZMU)
       IF(x .LT. 1.0d0) THEN
         fimd = xnorm*(x**(-ZMU))*xx
c          fimd2 = xxmasa**(-ZMU)
      ELSE
         fimd = xnorm*(x**(-2.3d0))*xx
       ENDIF
       WRITE(*,*)xxmasa,fimm,fimd,fimm/fimd
      ENDDO
	END



