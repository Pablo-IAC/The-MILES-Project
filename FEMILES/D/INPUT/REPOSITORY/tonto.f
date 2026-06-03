	dimension x(75)
	l=9999
	open(99,file='tiemposP',status='old')
	do i=1,75
	 read(99,*)t
	 x(i)=(10**t)/1.e9
	enddo
	do i=1,75
	 write(98,'(3(I4,1X),F8.5,1X,2(I4,1X))')l,l,l,x(76-i),l,l
	enddo
	end
