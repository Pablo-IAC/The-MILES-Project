c************************************************************************
	SUBROUTINE NICO(xccat,FINDEX)
	IMPLICIT double precision (A-H,O-Z)
	PARAMETER (NXMAX=2042)
	DIMENSION S(NXMAX),xccat(710,2)
c	DIMENSION S(710),xccat(710,2)
	do iz=1,710
	  S(iz)=xccat(iz,2)
	enddo
	STWV=xccat(1,1)
c	write(*,*)STWV,xccat(1,1)
	DISP=0.85d0
	RVEL=0.0d0
	NPIXELS=710
	call INDICES(40,NPIXELS,S,STWV,DISP,RVEL,FINDEX)
	return
	end
C******************************************************************************
	SUBROUTINE INDICES(NINDEX,NPIXELS,S,STWV,DISP,RVEL,FINDEX)
	IMPLICIT double precision (A-H,O-Z)
C Version 19-January-2000
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C******************************************************************************
C This subroutine computes a line-strength index in a single spectrum.
C INPUT:
C        NINDEX: number of index to be computed (see list below)
C        S(NPIXELS): spectrum to be measured
C        STWV: central wavelength of first pixel
C        DISP: dispersion (Angs./pixel)   Note: the scale must be linear
C        RVEL: spectrum radial velocity (km/sec)
C OUTPUT:
C        FINDEX: measured index
C------------------------------------------------------------------------------
C List of defined line-strength indices (see also subroutine INDEXDEF)
C (01) D4000...: 4000-Angs. break (Bruzual 1983)
C (02) B4000...: modified 4000-Angs. break (Gorgas et al. 1998)
C (03) HdA.....: Hdelta A (Worthey & Ottaviani 1997)
C (04) HdF.....: Hdelta F (Worthey & Ottaviani 1997)
C (05) CN1.....: Lick
C (06) CN2.....: Lick
C (07) Ca4227..: Lick
C (08) G4300...: Lick
C (09) HgA.....: Hgamma A (Worthey & Ottaviani 1997)
C (10) HgF.....: Hgamma F (Worthey & Ottaviani 1997)
C (11) Fe4383..: Lick
C (12) Ca4455..: Lick
C (13) Fe4531..: Lick
C (14) Fe4668..: Lick
C (15) Hbeta...: Lick
C (16) Hbeta_p.: Hbeta plus from Gonzalez thesis (p116) 
C (17) OIII_1..: OIII_1 from Gonzales thesis (p116)
C (18) OIII_2..: OIII_2 from Gonzales thesis (p116)
C (19) Fe5015..: Lick
C (20) Mg1.....: Lick
C (21) Mg2.....: Lick
C (22) Mgb5177.: Lick
C (23) Fe5270..: Lick
C (24) Fe5335..: Lick
C (25) Fe5406..: Lick
C (26) Fe5709..: Lick
C (27) Fe5782..: Lick
C (28) Na5895..: Lick
C (29) TiO1....: Lick
C (30) TiO2....: Lick
C (31) Ca1(DTT): Diaz, Terlevich & Terlevich (1989)
C (32) Ca2(DTT): Diaz, Terlevich & Terlevich (1989)
C (33) Ca3(DTT): Diaz, Terlevich & Terlevich (1989)
C (34) MgI(DTT): Diaz, Terlevich & Terlevich (1989)
C (35) CaT.....: CGCVP et al. (2000)
C (36) CaT1....: CGCVP et al. (2000)
C (37) CaT2....: CGCVP et al. (2000)
C (38) CaT3....: CGCVP et al. (2000)
C (39) PaT.....: CGCVP et al. (2000)
C (40) CaT*....: CGCVP et al. (2000)
C------------------------------------------------------------------------------
C parametros dimensionales
	INTEGER NXMAX              !maximum dimension in the spectral direction
	PARAMETER (NXMAX=2042)
C
        INTEGER NBDMAX        !maximum no. of bands allowed for a generic index
        PARAMETER (NBDMAX=198)                       !note that NBDMAX=NWVMAX/2
C
        INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
        PARAMETER (NWVMAX=396)
C constantes
c        REAL C                                           !speed of light (km/s)
        PARAMETER (C=2.9979246E+5)
C------------------------------------------------------------------------------
C parametros de la subrutina
	INTEGER NINDEX
	INTEGER NPIXELS
	DIMENSION S(NXMAX)
c	REAL STWV,DISP
c	REAL RVEL
c	REAL FINDEX
C------------------------------------------------------------------------------
C variables locales:
	INTEGER ITI                                             !tipo de indice
	INTEGER J
	INTEGER NBAND,NB
	INTEGER J1(NBDMAX),J2(NBDMAX)
	INTEGER J1MIN,J2MAX      !calcula los limites reales del indice a medir
	INTEGER NCEFF
	INTEGER NCONTI,NABSOR
	DIMENSION WV(NWVMAX)                                   !l.d.o. de las bandas
	DIMENSION  FWV(NWVMAX/4)   !constantes para multiplicar se\~{n}al de asbsorc.
	DIMENSION  SS(NXMAX)
c	REAL RCVEL,RCVEL1                                            !v/c,(1+z)
c	REAL CA,CB,C3,C4
	DIMENSION D1(NBDMAX),D2(NBDMAX),RL(NBDMAX),RG(NBDMAX)
c       REAL SUMRL
c	REAL WLA
c	REAL TC                                                          !flujo
	DIMENSION FX(NBDMAX)                                    !flujo en las bandas
c	REAL SB,SR
	DIMENSION  SC(NXMAX)
	DOUBLE PRECISION MWB,MWR
c	REAL FSMEAN
	DIMENSION  WL(NXMAX)                                     !pesos para el D4000
c	REAL WLMIN              !limite inferior en l.d.o. (espectro observado)
c	REAL F
c	REAL AMC,BMC       !parametros de recta en ajuste por minimos cuadrados
c	REAL SUMX,SUMY,SUM0,SUMXY,SUMXX  !sumatorios para ajuste por min. cuad.
c	REAL DETER            !variable local para ajuste por minimos cuadrados
	CHARACTER*8 INDEXNAME
	DOUBLE PRECISION SMEAN
	LOGICAL IFCHAN(NXMAX)   !indica los canales usados para medir el indice
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C proteccion
	IF(NPIXELS.GT.NXMAX)THEN
	  WRITE(*,101) 'FATAL ERROR: NPIXELS.GT.NXMAX'
	  WRITE(*,101) '=> Enlarge NXMAX and recompile the program.'
	  STOP
	END IF
C------------------------------------------------------------------------------
	FINDEX=0.                          !hasta que se demuestre lo contrario
C Definimos el indice que queremos medir
	CALL INDEXDEF(NINDEX,ITI,INDEXNAME,WV,FWV)
C Numero de bandas, limites y anchura teniendo en cuenta la velocidad radial
	NBAND=3               !indices atomicos y moleculares (lo mas probable)
	IF(ITI.EQ.3) NBAND=2                    !el D4000 solo tiene dos bandas
	IF(ITI.EQ.4) NBAND=2                    !el B4000 solo tiene dos bandas
	IF((ITI.GE.101).AND.(ITI.LE.9999))THEN     !para los indices genericos:
	  NCONTI=(ITI/100)                        !numero de bandas de continuo
	  NABSOR=ITI-NCONTI*100                !numero de bandas de absorciones
	  NBAND=NCONTI+NABSOR
	END IF
C..............................................................................
	WLMIN=STWV-DISP/2.
	RCVEL=RVEL/C
	RCVEL1=1.+RCVEL
	RCVEL1=RCVEL1/SQRT(1.-RCVEL*RCVEL)              !correccion relativista
	DO NB=1,NBAND                                          !para cada banda
	  CA=WV(2*NB-1)*RCVEL1                           !redshifted wavelength
	  CB=WV(2*NB)*RCVEL1                             !redshifted wavelength
	  C3=(CA-WLMIN)/DISP+1.d0                           !band limit (channel)
	  C4=(CB-WLMIN)/DISP  
c	  write(*,*)'C3,C4,CA,CB',C3,C4,CA,CB,STWV  !band limit (channel)
	  IF((C3.LT.1.d0).OR.(C4.GT.dble(NPIXELS)))THEN       !index out of range
	    WRITE(*,101) 'ERROR: index out of range in subroutine '//
     >       'INDICES.'
	    WRITE(*,100) 'Press Control Return to continue...'
	    READ(*,*)
	    RETURN
	  END IF
	  J1(NB)=INT(C3)                          !band limit: integer(channel)
	  J2(NB)=INT(C4)                          !band limit: integer(channel)
	  D1(NB)=C3-dble(J1(NB))            !fraction (excess) of first channel
	  D2(NB)=C4-dble(J2(NB))             !fraction (defect) of last channel
	  RL(NB)=CB-CA                       !redshifted band width (angstroms)
	  RG(NB)=C4-C3                        !redshifted band width (channels)
	END DO
C..............................................................................
	J1MIN=J1(1)     !calculamos limites por si las bandas no estan en orden
	J2MAX=J2(1)
	DO NB=2,NBAND
	  IF(J1MIN.GT.J1(NB)) J1MIN=J1(NB)
	  IF(J2MAX.LT.J2(NB)) J2MAX=J2(NB)
	END DO
C------------------------------------------------------------------------------
C Fijamos los canales a usar para medir el indice (utilizando la variable
C logica evitamos el problema de la posible superposicion de las bandas)
	DO J=1,NPIXELS
	  IFCHAN(J)=.FALSE.              !inicializamos: ningun canal utilizado
	END DO
C
	DO NB=1,NBAND                                    !recorremos las bandas
	  DO J=J1(NB),J2(NB)+1
	    IFCHAN(J)=.TRUE.
	  END DO
	END DO
C
	NCEFF=0             !contamos el numero de canales effectivo a utilizar
	DO J=1,NPIXELS
	  IF(IFCHAN(J)) NCEFF=NCEFF+1
	END DO
C------------------------------------------------------------------------------
C Normalizamos espectro a medir, usando la se\~{n}al solo en la region
C del indice
	SMEAN=0.D0
	DO J=1,NPIXELS
	  IF(IFCHAN(J)) SMEAN=SMEAN+DBLE(S(J))      !solo canales en las bandas
	END DO
	SMEAN=SMEAN/DBLE(NCEFF)              !valor promedio (DOUBLE PRECISION)
	FSMEAN=dble(SMEAN)                               !valor promedio (REAL)
C
	DO J=1,NPIXELS
	  SS(J)=S(J)/FSMEAN                              !normalizamos espectro
	END DO
C------------------------------------------------------------------------------
C Calculamos pseudo continuo en indices moleculares y atomicos (ITI=1,2)
C (formulas en Tesis de JJGG, pag. 35)
C
C NOTA: al transformar las integrales en sumatorios, habria que multiplicar
C cada valor de la funcion a integrar por el incremento en longitud de onda
C (que en el sumatorio coincide con DISP). Sin embargo, este valor es factor
C comun y puede salir fuera del sumatorio, por lo que solo hace falta
C introducirlo al final. Asimismo, como al calcular los limites de las bandas,
C J1() y J2(), hemos tenido presente la vel. radial, la anchura de las bandas,
C RL(), es mayor en un factor (1+z) que la anchura cuando el objeto no presenta
C velocidad radial. Es decir, el sumatorio se extiende sobre una region (en 
C longitud de onda) algo mayor. Sin embargo el efecto queda anulado al dividir
C por RL() que se encuentra ensanchado en el mismo factor.
C
	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN
C..............................................................................
	  SB=0.d0                              !cuentas promedio en la banda azul
	  DO J=J1(1),J2(1)+1                          !recorremos la banda azul
	    IF(J.EQ.J1(1))THEN
	      F=1.-D1(1)
	    ELSEIF(J.EQ.J2(1)+1)THEN
	      F=D2(1)
	    ELSE
	      F=1.d0
	    END IF
	    SB=SB+F*SS(J)
	  END DO
	  SB=SB*DISP                                     !completamos sumatorio
	  SB=SB/RL(1)                   !dividimos por anchura de la banda azul
C..............................................................................
	  SR=0.d0                              !cuentas promedio en la banda roja
	  DO J=J1(3),J2(3)+1                          !recorremos la banda roja
	    IF(J.EQ.J1(3))THEN
	      F=1.-D1(3)
	    ELSEIF(J.EQ.J2(3)+1)THEN
	      F=D2(3)
	    ELSE
	      F=1.d0
	    END IF
	    SR=SR+F*SS(J)
	  END DO
	  SR=SR*DISP                                     !completamos sumatorio
	  SR=SR/RL(3)                   !dividimos por anchura de la banda roja
C..............................................................................
C Trabajamos en la escala en l.d.o. sin corregir de Vrad (es decir, con las
C bandas de los indices desplazadas a Vrad correspondiente). Se obtiene lo
C mismo si mantenemos las bandas de los indices a Vrad=0 pero al calcular WLA
C para cada canal J dividimos por RCVEL1 para obtener la escala en l.d.o.
C corregida de Vrad.
	  MWB=(WV(1)+WV(2))/2.d0             !mean wavelength blue band at Vrad=0
	  MWB=MWB*RCVEL1                               !idem a Vrad considerado
	  MWR=(WV(5)+WV(6))/2.d0              !mean wavelength red band at Vrad=0
	  MWR=MWR*RCVEL1                               !idem a Vrad considerado
C Calculamos el valor del pseudo continuo desde el borde mas azul de todas las
C bandas al borde mas rojo de todas las bandas (no importa que las banda
C "azul" no sea la mas azul, etc.)
	  DO J=J1MIN,J2MAX+1
	    WLA=dble(J-1)*DISP+STWV                         !l.d.o. del canal J
	    SC(J)=(SB*(MWR-WLA)+SR*(WLA-MWB))/(MWR-MWB)
	  END DO
C..............................................................................
	END IF
C------------------------------------------------------------------------------
C Pesos para el D4000: debido a que la integral es  F(ldo)*d(nu), hay que
C multiplicar el flujo por el cuadrado de la longitud de onda, y de este
C modo la integral se transforma en F(ldo)*d(ldo).
C
C NOTA: estamos corrigiendo WLA de Vrad, aunque luego medimos el indice con
C los limites de las bandas multiplicados por (1+z). Esto no es importante
C porque solo hay un factor Cte=(1+z) en los pesos que no influye en el
C indice. Sin embargo, vamos a trabajar con los pesos asi porque de esta
C manera la normalizacion a 4000 siempre nos proporciona pesos cercanos a uno
C independientemente del valor de Vrad.
	IF(ITI.EQ.3)THEN
	  DO J=1,NPIXELS
	    WLA=dble(J-1)*DISP+STWV                         !l.d.o. del canal J 
	    WLA=WLA/RCVEL1               !l.d.o. del canal J corrigiendo a Vrad
	    WLA=WLA/4000.d0      !normalizamos a 4000 para tener todo proximo a 1
	    WL(J)=WLA*WLA
	  END DO
	END IF
C------------------------------------------------------------------------------
C Para el B4000 usaremos codigo comun con el D4000, donde los pesos son iguales
C a uno.
	IF(ITI.EQ.4)THEN
	  DO J=1,NPIXELS
	    WL(J)=1.0d0
	  END DO
	END IF
C------------------------------------------------------------------------------
C Para los indices genericos calculamos el pseudocontinuo ajustando por minimos
C cuadrados (pesando con errores si procede) a todos los pixels de las bandas 
C de continuo.
	IF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C..............................................................................
C calculamos la recta del continuo mediante minimos cuadrados (la recta es
C de la forma y= amc * x + bmc
C (NOTA: para la variable "x" utilizamos el valor del numero de pixel en lugar
C de la longitud de onda porque, en principio, son numeros mas pequenhos)
	  SUM0=0.d0
	  SUMX=0.d0
	  SUMY=0.d0
	  SUMXY=0.d0
	  SUMXX=0.d0
	  DO NB=1,NCONTI               !recorremos todas las bandas de continuo
	    DO J=J1(NB),J2(NB)+1  !recorremos todos los pixels de dichas bandas
	      IF(J.EQ.J1(NB))THEN        !comprobamos efecto de borde izquierdo
	        F=1.-D1(NB)
	      ELSEIF(J.EQ.J2(NB)+1)THEN                !efecto de borde derecho
	        F=D2(NB)
	      ELSE                                  !pixels sin efecto de borde
	        F=1.d0
	      END IF
	      SUM0=SUM0+F
	      SUMX=SUMX+F*dble(J)
	      SUMY=SUMY+F*SS(J)
	      SUMXY=SUMXY+F*dble(J)*SS(J)
	      SUMXX=SUMXX+F*dble(J)*dble(J)
	    END DO
	  END DO
	  DETER=SUM0*SUMXX-SUMX*SUMX
	  AMC=(SUM0*SUMXY-SUMX*SUMY)/DETER
	  BMC=(SUMXX*SUMY-SUMX*SUMXY)/DETER
C..............................................................................
C calculamos el pseudocontinuo desde el borde mas azul de todas las bandas
C hasta el borde mas rojo
	  DO J=J1MIN,J2MAX+1
	    SC(J)=AMC*dble(J)+BMC
	  END DO
C..............................................................................
	END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C medimos indices
	IF((ITI.EQ.1).OR.(ITI.EQ.2))THEN        !indices moleculares y atomicos
	  TC=0.d0
	  DO J=J1(2),J2(2)+1                       !recorremos la banda central
	    IF(J.EQ.J1(2))THEN
	      F=1.-D1(2)
	    ELSEIF(J.EQ.J2(2)+1)THEN
	      F=D2(2)
	    ELSE
	      F=1.d0
	    END IF
	    TC=TC+F*SS(J)/SC(J)
	  END DO
	  TC=TC*DISP                                     !completamos sumatorio
	  IF(ITI.EQ.1)THEN                                    !indice molecular
c	    FINDEX=-2.5d0*ALOG10(TC/RL(2))
	    FINDEX=(-2.5d0)*DLOG10(TC/RL(2))
	  ELSEIF(ITI.EQ.2)THEN                                  !indice atomico
	    FINDEX=RL(2)-TC
	    FINDEX=FINDEX/RCVEL1                       !correccion por redshift
	  END IF
C------------------------------------------------------------------------------
	ELSEIF((ITI.EQ.3).OR.(ITI.EQ.4))THEN                             !D4000
C NOTA: el sumatorio TC no incluye el factor DISP, que corresponderia
C con el incremento (diferencial en la integral), dado que al ser un factor
C constante tampoco altera el resultado a la hora de computar el D4000 (o el
C B4000).
	  DO NB=1,NBAND                                        !bucle en bandas
	    TC=0.d0
	    DO J=J1(NB),J2(NB)+1                        !recorremos la banda NB
	      IF(J.EQ.J1(NB))THEN
		F=1.-D1(NB)
	      ELSEIF(J.EQ.J2(NB)+1)THEN
		F=D2(NB)
	      ELSE
		F=1
	      END IF
	      TC=TC+F*SS(J)*WL(J)
	    END DO
	    FX(NB)=TC
	  END DO
	  FINDEX=FX(2)/FX(1)
C------------------------------------------------------------------------------
	ELSEIF((ITI.GE.101).AND.(ITI.LE.9999))THEN
C recorremos las bandas con absorciones
	  TC=0.d0
	  SUMRL=0.d0
	  DO NB=1,NABSOR
	    DO J=J1(NCONTI+NB),J2(NCONTI+NB)+1     !recorremos todos los pixels
	      IF(J.EQ.J1(NCONTI+NB))THEN !comprobamos efecto de borde izquierdo
	        F=1.-D1(NCONTI+NB)
	      ELSEIF(J.EQ.J2(NCONTI+NB)+1)THEN         !efecto de borde derecho
	        F=D2(NCONTI+NB)
	      ELSE                                  !pixels sin efecto de borde
	        F=1.d0
	      END IF
	      TC=TC+F*FWV(NB)*SS(J)/SC(J)   !multiplicamos absorcion por factor
	    END DO
	    SUMRL=SUMRL+FWV(NB)*RL(NCONTI+NB)
	  END DO
	  TC=TC*DISP                                  !completamos el sumatorio
	  FINDEX=SUMRL-TC         !medimos el indice generico como los atomicos
	  FINDEX=FINDEX/RCVEL1                         !correccion por redshift
C------------------------------------------------------------------------------
	END IF
C------------------------------------------------------------------------------
C Fin de subrutina
c100	FORMAT(A,$)
101	FORMAT(A)
100	FORMAT(A)
	END
C
C******************************************************************************
C Retorna las "bandpasses" del indice numero NINDEX, a la vez que el tipo de
C indice:
C ITI=3 o 4: D4000 o B4000
C ITI=2: indice atomico
C ITI=1: indice molecular
C------------------------------------------------------------------------------
	SUBROUTINE INDEXDEF(NINDEX,ITI,NAME,WV,FWV)
	IMPLICIT double precision (A-H,O-Z)
c	IMPLICIT NONE
C parametros generales
        INTEGER NWVMAX      !2*maximum no. of bands allowed for a generic index
        PARAMETER (NWVMAX=396)                       !see subroutine selindex.f
C parametros de la subrutina
	INTEGER NINDEX
	INTEGER ITI
	CHARACTER*8 NAME
	DIMENSION WV(NWVMAX)                                   !l.d.o. de las bandas
	DIMENSION FWV(NWVMAX/4)   !constantes para multiplicar se\~{n}al de asbsorc.
C------------------------------------------------------------------------------
	IF(NINDEX.EQ.1)THEN
	  NAME='D4000'
	  ITI=3
	  WV(1)=3750.000d0
	  WV(2)=3950.000d0
	  WV(3)=4050.000d0
	  WV(4)=4250.000d0
	ELSEIF(NINDEX.EQ.2)THEN
	  NAME='B4000'
	  ITI=4
	  WV(1)=3750.000d0
	  WV(2)=3950.000d0
	  WV(3)=4050.000d0
	  WV(4)=4250.000d0
	ELSEIF(NINDEX.EQ.3)THEN
	  NAME='HdA'
	  ITI=2
	  WV(1)=4041.600d0
	  WV(2)=4079.750d0
	  WV(3)=4083.500d0
	  WV(4)=4122.250d0
	  WV(5)=4128.500d0
	  WV(6)=4161.000d0
	ELSEIF(NINDEX.EQ.4)THEN
	  NAME='HdF'
	  ITI=2
	  WV(1)=4057.250d0
	  WV(2)=4088.500d0
	  WV(3)=4091.000d0
	  WV(4)=4112.250d0
	  WV(5)=4114.750d0
	  WV(6)=4137.250d0
	ELSEIF(NINDEX.EQ.5)THEN
	  NAME='CN1'
	  ITI=1
	  WV(1)=4080.125d0
	  WV(2)=4117.625d0
	  WV(3)=4142.125d0
	  WV(4)=4177.125d0
	  WV(5)=4244.125d0
	  WV(6)=4284.125d0
	ELSEIF(NINDEX.EQ.6)THEN
	  NAME='CN2'
	  ITI=1
	  WV(1)=4083.875d0
	  WV(2)=4096.375d0
	  WV(3)=4142.125d0
	  WV(4)=4177.125d0
	  WV(5)=4244.125d0
	  WV(6)=4284.125d0
	ELSEIF(NINDEX.EQ.7)THEN
	  NAME='Ca4227'
	  ITI=2
	  WV(1)=4211.000d0
	  WV(2)=4219.750d0
	  WV(3)=4222.250d0
	  WV(4)=4234.750d0
	  WV(5)=4241.000d0
	  WV(6)=4251.000d0
	ELSEIF(NINDEX.EQ.8)THEN
	  NAME='G4300'
	  ITI=2
	  WV(1)=4266.375d0
	  WV(2)=4282.625d0
	  WV(3)=4281.375d0
	  WV(4)=4316.375d0
	  WV(5)=4318.875d0
	  WV(6)=4335.125d0
	ELSEIF(NINDEX.EQ.9)THEN
	  NAME='HgA'
	  ITI=2
	  WV(1)=4283.500d0
	  WV(2)=4319.750d0
	  WV(3)=4319.750d0
	  WV(4)=4363.500d0
	  WV(5)=4367.250d0
	  WV(6)=4419.750d0
	ELSEIF(NINDEX.EQ.10)THEN
	  NAME='HgF'
	  ITI=2
	  WV(1)=4283.500d0
	  WV(2)=4319.750d0
	  WV(3)=4331.250d0
	  WV(4)=4352.250d0
	  WV(5)=4354.750d0
	  WV(6)=4384.750d0
	ELSEIF(NINDEX.EQ.11)THEN
	  NAME='Fe4383'
	  ITI=2
	  WV(1)=4359.125d0
	  WV(2)=4370.375d0
	  WV(3)=4369.125d0
	  WV(4)=4420.375d0
	  WV(5)=4442.875d0
	  WV(6)=4455.375d0
	ELSEIF(NINDEX.EQ.12)THEN
	  NAME='Ca4455'
	  ITI=2
	  WV(1)=4445.875d0
	  WV(2)=4454.625d0
	  WV(3)=4452.125d0
	  WV(4)=4474.625d0
	  WV(5)=4477.125d0
	  WV(6)=4492.125d0
	ELSEIF(NINDEX.EQ.13)THEN
	  NAME='Fe4531'
	  ITI=2
	  WV(1)=4504.250d0
	  WV(2)=4514.250d0
	  WV(3)=4514.250d0
	  WV(4)=4559.250d0
	  WV(5)=4560.500d0
	  WV(6)=4579.250d0
	ELSEIF(NINDEX.EQ.14)THEN
	  NAME='Fe4668'
	  ITI=2
	  WV(1)=4611.500d0
	  WV(2)=4630.250d0
	  WV(3)=4634.000d0
	  WV(4)=4720.250d0
	  WV(5)=4742.750d0
	  WV(6)=4756.500d0
	ELSEIF(NINDEX.EQ.15)THEN
	  NAME='Hbeta'
	  ITI=2
	  WV(1)=4827.875d0
	  WV(2)=4847.875d0
	  WV(3)=4847.875d0
	  WV(4)=4876.625d0
	  WV(5)=4876.625d0
	  WV(6)=4891.625d0
	ELSEIF(NINDEX.EQ.16)THEN
	  NAME='Hbeta_p'
	  ITI=2
	  WV(1)=4815.000d0
	  WV(2)=4845.000d0
	  WV(3)=4851.320d0
	  WV(4)=4871.320d0
	  WV(5)=4880.000d0
	  WV(6)=4930.000d0
	ELSEIF(NINDEX.EQ.17)THEN
	  NAME='OIII_1'
	  ITI=2
	  WV(1)=4885.000d0
	  WV(2)=4935.000d0
	  WV(3)=4948.920d0
	  WV(4)=4978.920d0
	  WV(5)=5030.000d0
	  WV(6)=5070.000d0
	ELSEIF(NINDEX.EQ.18)THEN
	  NAME='OIII_2'
	  ITI=2
	  WV(1)=4885.000d0
	  WV(2)=4935.000d0
	  WV(3)=4996.850d0
	  WV(4)=5016.850d0
	  WV(5)=5030.000d0
	  WV(6)=5070.000d0
	ELSEIF(NINDEX.EQ.19)THEN
	  NAME='Fe5015'
	  ITI=2
	  WV(1)=4946.500d0
	  WV(2)=4977.750d0
	  WV(3)=4977.750d0
	  WV(4)=5054.000d0
	  WV(5)=5054.000d0
	  WV(6)=5065.250d0
	ELSEIF(NINDEX.EQ.20)THEN
	  NAME='Mg1'
	  ITI=1
	  WV(1)=4895.125d0
	  WV(2)=4957.625d0
	  WV(3)=5069.125d0
	  WV(4)=5134.125d0
	  WV(5)=5301.125d0
	  WV(6)=5366.125d0
	ELSEIF(NINDEX.EQ.21)THEN
	  NAME='Mg2'
	  ITI=1
	  WV(1)=4895.125d0
	  WV(2)=4957.625d0
	  WV(3)=5154.125d0
	  WV(4)=5196.625d0
	  WV(5)=5301.125d0
	  WV(6)=5366.125d0
	ELSEIF(NINDEX.EQ.22)THEN
	  NAME='Mgb5177'
	  ITI=2
	  WV(1)=5142.625d0
	  WV(2)=5161.375d0
	  WV(3)=5160.125d0
	  WV(4)=5192.625d0
	  WV(5)=5191.375d0
	  WV(6)=5206.375d0
	ELSEIF(NINDEX.EQ.23)THEN
	  NAME='Fe5270'
	  ITI=2
	  WV(1)=5233.150d0
	  WV(2)=5248.150d0
	  WV(3)=5245.650d0
	  WV(4)=5285.650d0
	  WV(5)=5285.650d0
	  WV(6)=5318.150d0
	ELSEIF(NINDEX.EQ.24)THEN
	  NAME='Fe5335'
	  ITI=2
	  WV(1)=5304.625d0
	  WV(2)=5315.875d0
	  WV(3)=5312.125d0
	  WV(4)=5352.125d0
	  WV(5)=5353.375d0
	  WV(6)=5363.375d0
	ELSEIF(NINDEX.EQ.25)THEN
	  NAME='Fe5406'
	  ITI=2
	  WV(1)=5376.250d0
	  WV(2)=5387.500d0
	  WV(3)=5387.500d0
	  WV(4)=5415.000d0
	  WV(5)=5415.000d0
	  WV(6)=5425.000d0
	ELSEIF(NINDEX.EQ.26)THEN
	  NAME='Fe5709'
	  ITI=2
	  WV(1)=5672.875d0
	  WV(2)=5696.625d0
	  WV(3)=5696.625d0
	  WV(4)=5720.375d0
	  WV(5)=5722.875d0
	  WV(6)=5736.625d0
	ELSEIF(NINDEX.EQ.27)THEN
	  NAME='Fe5782' 
	  ITI=2
	  WV(1)=5765.375d0
	  WV(2)=5775.375d0
	  WV(3)=5776.625d0
	  WV(4)=5796.625d0
	  WV(5)=5797.875d0
	  WV(6)=5811.625d0
	ELSEIF(NINDEX.EQ.28)THEN
	  NAME='Na5895'
	  ITI=2
	  WV(1)=5860.625d0
	  WV(2)=5875.625d0
	  WV(3)=5876.875d0
	  WV(4)=5909.375d0
	  WV(5)=5922.125d0
	  WV(6)=5948.125d0
	ELSEIF(NINDEX.EQ.29)THEN
	  NAME='TiO1'
	  ITI=1
	  WV(1)=5816.625d0
	  WV(2)=5849.125d0
	  WV(3)=5936.625d0
	  WV(4)=5994.125d0
	  WV(5)=6038.625d0
	  WV(6)=6103.625d0
	ELSEIF(NINDEX.EQ.30)THEN
	  NAME='TiO2'
	  ITI=1
	  WV(1)=6066.625d0
	  WV(2)=6141.625d0
	  WV(3)=6189.625d0
	  WV(4)=6272.125d0
	  WV(5)=6372.625d0
	  WV(6)=6415.125d0
	ELSEIF(NINDEX.EQ.31)THEN
	  NAME='Ca1(DTT)'
	  ITI=2
	  WV(1)=8447.500d0
	  WV(2)=8462.500d0
	  WV(3)=8483.000d0
	  WV(4)=8513.000d0
	  WV(5)=8842.500d0
	  WV(6)=8857.500d0
	ELSEIF(NINDEX.EQ.32)THEN
	  NAME='Ca2(DTT)'
	  ITI=2
	  WV(1)=8447.500d0
	  WV(2)=8462.500d0
	  WV(3)=8527.000d0
	  WV(4)=8557.000d0
	  WV(5)=8842.500d0
	  WV(6)=8857.500d0
	ELSEIF(NINDEX.EQ.33)THEN
	  NAME='Ca3(DTT)'
	  ITI=2
	  WV(1)=8447.500d0
	  WV(2)=8462.500d0
	  WV(3)=8647.000d0
	  WV(4)=8677.000d0
	  WV(5)=8842.500d0
	  WV(6)=8857.500d0
	ELSEIF(NINDEX.EQ.34)THEN
	  NAME='MgI(DTT)'
	  ITI=2
	  WV(1)=8775.000d0
	  WV(2)=8787.000d0
	  WV(3)=8799.500d0
	  WV(4)=8814.500d0
	  WV(5)=8845.000d0
	  WV(6)=8855.000d0
	ELSEIF(NINDEX.EQ.35)THEN
	  NAME='CaT'
	  ITI=503
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8484.00d0
	  WV(12)=8513.00d0
	  FWV(1)=1.0d0
          WV(13)=8522.00d0
	  WV(14)=8562.00d0
	  FWV(2)=1.0d0
          WV(15)=8642.00d0
	  WV(16)=8682.00d0
	  FWV(3)=1.0d0                   
	ELSEIF(NINDEX.EQ.36)THEN
	  NAME='CaT1'
	  ITI=501
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8484.00d0
	  WV(12)=8513.00d0
	  FWV(1)=1.0d0
	ELSEIF(NINDEX.EQ.37)THEN
	  NAME='CaT2'
	  ITI=501
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8522.00d0
	  WV(12)=8562.00d0
	  FWV(1)=1.0d0
	ELSEIF(NINDEX.EQ.38)THEN
	  NAME='CaT3'
	  ITI=501
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8642.00d0
	  WV(12)=8682.00d0
	  FWV(1)=1.0d0
	ELSEIF(NINDEX.EQ.39)THEN
	  NAME='PaT'
	  ITI=503
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8461.00d0
	  WV(12)=8474.00d0
	  FWV(1)=1.0d0
          WV(13)=8577.00d0
	  WV(14)=8619.00d0
	  FWV(2)=1.0d0
          WV(15)=8730.00d0
	  WV(16)=8772.00d0
	  FWV(3)=1.0d0
	ELSEIF(NINDEX.EQ.40)THEN
	  NAME='CaT*'
	  ITI=506
          WV(1)=8474.000d0
	  WV(2)=8484.000d0
          WV(3)=8563.000d0
	  WV(4)=8577.000d0
          WV(5)=8619.000d0
	  WV(6)=8642.000d0
          WV(7)=8700.000d0
	  WV(8)=8725.000d0
          WV(9)=8776.000d0
	  WV(10)=8792.00d0
          WV(11)=8484.00d0
	  WV(12)=8513.00d0
	  FWV(1)=1.0d0
          WV(13)=8522.00d0
	  WV(14)=8562.00d0
	  FWV(2)=1.0d0
          WV(15)=8642.00d0
	  WV(16)=8682.00d0
	  FWV(3)=1.0d0
          WV(17)=8461.00d0
	  WV(18)=8474.00d0
	  FWV(4)=-0.93d0
          WV(19)=8577.00d0
	  WV(20)=8619.00d0
	  FWV(5)=-0.93d0
          WV(21)=8730.00d0
	  WV(22)=8772.00d0 
	  FWV(6)=-0.93d0
	ELSE
	  WRITE(*,100) 'NINDEX='
	  WRITE(*,*) NINDEX
	  WRITE(*,101) 'FATAL ERROR: index number out of range.'
	  WRITE(*,101) '(subroutine INDEXDEF)'
	  STOP
	END IF
C
c100	FORMAT(A,$)
100	FORMAT(A)
101	FORMAT(A)
	END
