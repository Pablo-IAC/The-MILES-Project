	character*3 iii
	character*6 a,aa
	dimension a(852)
	iii='non'
      open(99,file='../MILES_ALPHA/PARAM_MILES_mgfe_MILONE',
     & status='old')
	do i=1,852
	 read(99,*)a(i)
	enddo
      close(99)
	open(99,file='PARAM_MILES_mgfe',status='old')
	DO i=1,925
	 n=0
	 read(99,'(A6)')aa
	 do j=1,852
	  if(aa.eq.a(j))then
	   n=1
	   goto 233
	  endif
	 enddo
233	 continue
	 if(n.eq.0)then
	   write(66,'(A3,1(1X,A6))')iii,aa
	 else
	   write(66,'(A3,2(1X,A6))')'___',aa,a(j)
	 endif
	ENDDO
	close(99)
	end
