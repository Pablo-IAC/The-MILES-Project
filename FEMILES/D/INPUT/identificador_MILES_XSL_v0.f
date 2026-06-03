      character*13 miles,m,nadam
      character*6 nada
      character*26 xsl,xslm
      character*279 tonto
      character*400 tonton,milon
      PARAMETER (nmil=1154)
      PARAMETER (nxsl=536)
c      character*346 a !toda la linea de la estrella MILES
      character*332 a !toda la linea de la estrella MILES
c      dimension a(1154),miles(1154)
      dimension a(nmil),miles(nmil),m(nmil)
      nada='------'
      nadam='-------------'
      open(99,file='PARAM_MILES_star_names',status='old')
      do i=1,nmil
        read(99,'(A)')milon
        miles(i)=milon(1:13)
        a=milon(14:)
	mm=lnblnk(miles(i))
        write(m(i)(1:mm),'(A)')miles(i)        
      enddo
      close(99)
      open(99,file='PARAM_XSL_arentsen_parameters',status='old')
      DO i=1,nxsl
c       read(99,*)xsl,tonto
       read(99,'(A)')tonton
       xsl=tonton(1:26)
       tonto=tonton(27:)
	mx=lnblnk(xsl)
        write(xslm(1:mx),'(A)')xsl        
       do j=1,nmil
        if(xslm.eq.m(j))then
         write(91,'(A26,A13,A279,A332)')xsl,miles(j),tonto,a(j)
c         write(*,'(A26,A279,A332)')xsl,tonto,a(j)
         write(*,'(A26,A13)')xsl,miles(j)
         goto 233
        endif
       enddo
       write(91,'(A26,A13,A279)')xsl,nadam,tonto
c       write(*,'(A26,A279,A6)')xsl,tonto,nada
233	 continue
      ENDDO
      close(99)
      END
