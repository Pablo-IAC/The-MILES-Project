      INTEGER nmil,nxsl
      PARAMETER (nmil=1154)
      PARAMETER (nxsl=192)
      DIMENSION pm(nmil,5),px(nxsl,5)
      CHARACTER*6 nm(nmil),nxm(nxsl),nx(nxsl)
c ORDENAR FICHERO PRIMERO      
      OPEN(99, FILE='PARAM_MILES_FINAL_VERSION_OK', STATUS='OLD')
      DO i = 1, nmil
         READ(99,*)nm(i),(pm(i,j),j=1,5)
      ENDDO
      CLOSE(99)
      OPEN(99, FILE='PARAM_XSL_UPDATED_MILES_bis', STATUS='OLD')
      DO i = 1, nxsl
         READ(99,*)nxm(i),nx(i),(px(i,j),j=1,5)
         DO j=1,nmil
	  IF (nxm(i).eq.nm(j)) THEN
	    IF (px(i,1).ne.pm(j,1).OR.px(i,2).ne.pm(j,2).OR. 
     &px(i,3).ne.pm(j,3).OR.px(i,5).ne.pm(j,5)) 
     &WRITE(77,*)nxm(i),nm(j),nx(i),'REVISE'
           ELSE
	    GOTO 22
           ENDIF
22	   CONTINUE
	 ENDDO
      ENDDO
      CLOSE(99)
      END
