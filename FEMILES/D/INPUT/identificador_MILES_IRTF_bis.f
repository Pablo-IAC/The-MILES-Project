      PARAMETER (nm=1154)
      PARAMETER (ny=238)
c      character*42 tonto
      character*325 am(nm) !toda la linea de la estrella MILES
      DIMENSION pm(nm,5)
      CHARACTER*6 m(nm),ym,y !,ym(ny),y(ny)
      open(99,file='PARAM_MILES',status='old')
      do i=1,nm
c        read(99,*)m(i),(pm(i,j),j=1,5),am(i)
        read(99,'(A6,X,A325)')m(i),am(i)
      enddo
      close(99)
      open(99,file='PARAM_IRTF',status='old')
c      DO i=1,380
      DO i=1,ny
       read(99,*)ym,y
       do j=1,nm
        if(ym.eq.m(j))then
         write(76,'(3(A6,X),A325)')ym,y,m(j),am(j)
c         write(76,'(3(A6,X),F8.2,X,4(F6.2,X),A291)')
c     &ym,y,m(j),pm(j,1),(pm(j,l),l=2,5),am(j)
         goto 233
        endif
       enddo
c       write(94,'(A13,A42)')irtf,tonto
233	 continue
      ENDDO
      close(99)
      END
