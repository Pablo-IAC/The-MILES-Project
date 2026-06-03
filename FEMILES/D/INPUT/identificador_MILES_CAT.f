      character*18 miles
      character*13 cat
      character*6 ws
      character*260 tonto
      character*5 xx,YY
      character*346 a !toda la linea de la estrella MILES
c      dimension a(1130),miles(1130)
      dimension a(1763),miles(1763)
      xx='MILES'
      YY='CAT  '
c      open(99,file='PARAM_MILES_extended_work7',status='old')
      open(99,file='MILESIRTFNGSLMIUSCGTC.names',status='old')
      do i=1,1763
c        read(99,'(A13,A346)')miles(i),a(i)
        read(99,*)miles(i)
c	write(*,*)miles(i),a(i)
      enddo
      close(99)
      open(99,file='PARAM_CAT.names',status='old')
      DO i=1,706 !no incluyen estrellas con nombres no identificados por simbad
c       read(99,'(A18,A260)')ngsl,tonto
       read(99,*)cat,ws
       do j=1,1763
        if(cat.eq.miles(j))then
c         write(87,'(A18,A18,A260,A346)')ngsl,miles(j),tonto,a(j)
         write(84,'(A5,X,A18,A7,A18)')xx,cat,ws,miles(j)
         goto 233
        endif
       enddo
c       write(87,'(A18,A260)')ngsl,tonto
c       write(87,'(A4,X,A18,A260)')YY,ngsl,tonto
       write(84,'(A5,X,A18,A7,A18)')YY,cat,ws,cat
233	 continue
      ENDDO
      close(99)
      END
