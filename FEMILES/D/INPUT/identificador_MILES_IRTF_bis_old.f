      character*13 miles,irtf
      character*42 tonto
      character*346 a !toda la linea de la estrella MILES
      dimension a(1154),miles(1154)
      open(99,file='PARAM_MILES',status='old')
      do i=1,1130
        read(99,'(A13,A346)')miles(i),a(i)
	write(*,*)miles(i),a(i)
      enddo
      close(99)
      open(99,file='PARAM_IRTF',status='old')
      DO i=1,380
       read(99,'(A13,A42)')irtf,tonto
       do j=1,1130
        if(irtf.eq.miles(j))then
         write(94,'(A13,A42,A346)')irtf,tonto,a(j)
         goto 233
        endif
       enddo
       write(94,'(A13,A42)')irtf,tonto
233	 continue
      ENDDO
      close(99)
      END
