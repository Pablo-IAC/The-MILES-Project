      character*18 miles,ngsl
      character*260 tonto
      character*5 xx,YY
      character*346 a !toda la linea de la estrella MILES
c      dimension a(1130),miles(1130)
      dimension a(1293),miles(1293)
      xx='MILES'
      YY='NGSL '
c      open(99,file='PARAM_MILES_extended_work7',status='old')
      open(99,file='MILESIRTFGTC.names',status='old')
      do i=1,1293
c        read(99,'(A13,A346)')miles(i),a(i)
        read(99,*)miles(i)
c	write(*,*)miles(i),a(i)
      enddo
      close(99)
      open(99,file='PARAM_NGSL.names',status='old')
      DO i=1,424 !no incluyen estrellas con nombres no identificados por simbad
c       read(99,'(A18,A260)')ngsl,tonto
       read(99,*)ngsl
       do j=1,1293
        if(ngsl.eq.miles(j))then
c         write(87,'(A18,A18,A260,A346)')ngsl,miles(j),tonto,a(j)
         write(87,'(A5,X,A18,A18)')xx,ngsl,miles(j)
         goto 233
        endif
       enddo
c       write(87,'(A18,A260)')ngsl,tonto
c       write(87,'(A4,X,A18,A260)')YY,ngsl,tonto
       write(87,'(A5,X,A18,A18)')YY,ngsl,ngsl
233	 continue
      ENDDO
      close(99)
      END
