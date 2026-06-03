	dimension f1(311,3),f2(311,7),ng(311),hd(311),rr(311,3),r(311,3)
	dimension M(311)
 	character*6 ng
	character*20 hd
	character*30 M
      open(99,file='PARAM_NGSL_13_04_2014',status='old')
      DO i=1,311
	 read(99,*)(f1(i,j),j=1,3),ng(i),(f2(i,jj),jj=1,7),hd(i)
	 if(i.le.119)then
	   M(i)='MIL                         '
	 elseif(i.gt.119.and.i.le.135)then
	   M(i)='CAT                         '
	 else
	   if(f1(i,2).ge.3.0)then
	    M(i)='Teff-(13.1784+(z*(-41.9432))'
	   else
	    M(i)='Teff-(54.7807+(z*(-115.537))'
	   endif
	 endif
c trasnformamos todas las * con el polinomio de NGSL_params.pro
	 if(f1(i,2).ge.3.0)then
	    rr(i,1)=f1(i,1)-(13.1784+(f1(i,3))*(-41.9432))
	    rr(i,2)=f1(i,2)
	    rr(i,3)=f1(i,3)
	 else
	    rr(i,1)=f1(i,1)-(54.7807+(f1(i,3))*(-115.537))
	    rr(i,2)=f1(i,2)
	    rr(i,3)=f1(i,3)
	 endif
      END DO
      close(99)
c Si es MILES/CAT lo dejo como estaba
      DO i=1,311
	 if(i.le.135)then
	  r(i,1)=f2(i,1)
	  r(i,2)=f2(i,2)
	  r(i,3)=f2(i,3)
	 else
	  r(i,1)=rr(i,1)
	  r(i,2)=rr(i,2)
	  r(i,3)=rr(i,3)
	 endif
      END DO
44    FORMAT(A6,X,F6.0,2(X,F5.2),4(X,F5.0),X,A20,X,A30,
     &2(X,F6.0,2(X,F5.2)))
      c=99.
	open(99,file='PARAM_NGSL_13_04_2014_OK',status='new')
	DO i=1,311
       write(99,44)ng(i),(r(i,j),j=1,3),c,(f2(i,jj),jj=5,7),hd(i),
     &M(i),(f1(i,jjj),jjj=1,3),(rr(i,jjjj),jjjj=1,3)
	END DO
	close(99)
	end
