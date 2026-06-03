      INTEGER nmil,nxsl
      PARAMETER (nmil=1154)
      PARAMETER (nxsl=192)
c      CHARACTER*13 miles(nmil),m(nmil),nadam
      CHARACTER*13 miles(nmil),m(nmil),nadam
      CHARACTER*500 tonton,milon
      CHARACTER*6 nada
      CHARACTER*26 xsl(nxsl),xslm
      CHARACTER*279 tonto
      CHARACTER*332 a(nmil)
      nada  = '------'
      nadam = '-------------'
      OPEN(99, FILE='PARAM_MILES_star_names', STATUS='OLD')
      DO i = 1, nmil
         READ(99,'(A)')milon
         miles(i)=milon(1:13)
         a(i)= milon(14:)
c         mm=LNBLNK(miles(i))
         m(i)=miles(i)(1:LNBLNK(miles(i)))
         WRITE(90,*)m(i)
      ENDDO
      CLOSE(99)
c
      OPEN(99, FILE='PARAM_XSL_arentsen_parameters', STATUS='OLD')
      DO i = 1, 999999
         READ(99,*,END=999)xsl(i)
      ENDDO
999   CLOSE(99)
      OPEN(99, FILE='PARAM_XSL_arentsen_parameters', STATUS='OLD')
      DO i = 1, 999999
         READ(99,'(A)',END=888)tonton
c         xsl=tonton(1:26)
c         tonto=tonton(27:)
         tonto=tonton
c         mx=LNBLNK(xsl)
c	 xslm=xsl(1:LNBLNK(xsl))
c	 xslm=xsl(i)
c         WRITE(xslm(1:mx),'(A)')xsl
c         WRITE(89,*)xslm,LNBLNK(xsl),xsl
         DO j = 1, nmil
            IF (xsl(i) .EQ. m(j)) THEN
               WRITE(91,'(A26,A304,A332)')xsl(i),tonto,a(j)
c               WRITE(*,'(A26,A13)') xsl, miles(j)
              WRITE(*,*)xsl(i),m(j)
               GOTO 233
            ENDIF
         ENDDO
         WRITE(91,'(A26,A279)')xsl(i),tonto
233      CONTINUE
      ENDDO
888   CONTINUE
      CLOSE(99)
      END
