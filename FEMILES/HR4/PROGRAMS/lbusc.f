c*******SUBRUTINA LBUSC************************************************
        SUBROUTINE lbusc(tiss,giss,fiss,viss,xxb,xxr,aloe,fra4000)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION xxb(1120,2),xxr(1120,2)
        DIMENSION mf(1000),mfv(1000)
        DIMENSION sta(650,6),sel(650,2),wn(650,2),fstar(650,3),aloe(3)
        CHARACTER*20 sel,ast,wn
        COMMON/hrbus/aa(650,6),nstar
        COMMON/chrbus/ast(650,2)
        dtiss=75.d0
        dfiss=0.15d0
        dgiss=0.5d0
        dviss=0.5d0
        call cajas(tiss,giss,fiss,cajat)
        DO nonoo=1,1000
	 IF(fiss.ge.0.15d0)THEN
	  if(giss.lt.3.5d0)then                               !giants
	   if(tiss.le.3150.d0)then
             tup=3125.d0+dble(nonoo)*25.d0
             tdo=2500.d0
           elseif(tiss.le.3300.d0.and.tiss.gt.3150.d0)then
             tup=tiss+160.d0+dble(nonoo)*18.d0
             tdo=2900.d0
           elseif(tiss.le.3600.d0.and.tiss.gt.3300.d0)then
             tup=tiss+145.d0+dble(nonoo)*18.d0
             tdo=tiss-175.d0-dble(nonoo)*18.d0 !en 3600 cajas daria 178
	   else
             tup=tiss+cajat+15.d0+dble(nonoo)*15.d0
             tdo=tiss-cajat-15.d0-dble(nonoo)*15.d0
	   endif
	  else                                              !dwarfs
	   if(tiss.lt.4800.d0.and.tiss.ge.4200.d0)then
             tup=tiss+140.d0+dble(nonoo)*10.d0
             tdo=tiss-140.d0-dble(nonoo)*10.d0 !en 4800 cajas daria 150
           elseif(tiss.lt.3900.d0)then
             tup=tiss+270.d0+dble(nonoo)*10.d0
             tdo=tiss-270.d0-dble(nonoo)*10.d0 !en 3900 cajas daria 270
	   else
             tup=tiss+cajat+dble(nonoo)*15.d0
             tdo=tiss-cajat-dble(nonoo)*15.d0
	   endif
	  endif
	 ELSEIF(fiss.lt.0.15d0.and.fiss.gt.-0.15d0)THEN
	  if(giss.lt.3.5d0)then                               !giants
	   if(tiss.le.3150.d0)then
             tup=3125.d0+dble(nonoo)*25.d0
             tdo=2800.d0
           elseif(tiss.le.3300.d0.and.tiss.gt.3150.d0)then
             tup=tiss+160.d0+dble(nonoo)*18.d0
             tdo=2900.d0
           else
             tup=tiss+cajat+dble(nonoo)*15.d0
             tdo=tiss-cajat-dble(nonoo)*15.d0
           endif
	  else                                              !dwarfs
           if(tiss.lt.3900.d0)then
             tup=tiss+270.d0+dble(nonoo)*10.d0
             tdo=tiss-270.d0-dble(nonoo)*10.d0 !en 3900 cajas daria 270
	   else
             tup=tiss+cajat+dble(nonoo)*15.d0
             tdo=tiss-cajat-dble(nonoo)*15.d0
	   endif
	  endif
	 ELSEIF(fiss.le.-0.15d0.and.fiss.gt.-0.9d0)THEN
	  if(giss.lt.3.5d0)then                               !giants
	   if(tiss.gt.5000.d0.and.fiss.lt.-0.55d0)then  
             tup=tiss+cajat+80.d0+dble(nonoo)*15.d0 !es un remiendo
             tdo=tiss-cajat-70.d0-dble(nonoo)*15.d0 !para las hot HB
           else
             tup=tiss+cajat+35.d0+dble(nonoo)*15.d0
             tdo=tiss-cajat-35.d0-dble(nonoo)*15.d0
           endif
	  else                                              !dwarfs
	   if(tiss.gt.5950.d0.and.tiss.lt.6150.d0.and.fiss.lt.-0.55d0)then  
             tup=tiss+cajat+140.d0+dble(nonoo)*15.d0 !es un remiendo
             tdo=tiss-cajat-25.d0-dble(nonoo)*15.d0  !para el Toff
           else
             tup=tiss+cajat+dble(nonoo)*15.d0
             tdo=tiss-cajat-dble(nonoo)*15.d0
           endif
	  endif
	 ELSE
             tup=tiss+cajat+dble(nonoo)*15.d0
             tdo=tiss-cajat-dble(nonoo)*15.d0
	 ENDIF
	 if(tiss.gt.28000.d0)then
	  tdo=28000.d0
	  tup=200000.d0
	 endif
         nok=0
         do iijo=1,nstar
          if(aa(iijo,1).lt.tup.and.aa(iijo,1).gt.tdo)then
           IF(tiss.gt.9000.d0)THEN
                  nok=nok+1
                  sel(nok,1)=ast(iijo,1)
                  sel(nok,2)=ast(iijo,2)
                  do kkkkl=1,6
                   sta(nok,kkkkl)=aa(iijo,kkkkl)
                  enddo
           ELSEIF(tiss.le.9000.)THEN
            IF(giss.ge.3.6d0)THEN
              if(aa(iijo,2).ge.3.3)then
                  nok=nok+1
                  sel(nok,1)=ast(iijo,1)
                  sel(nok,2)=ast(iijo,2)
                  do kkkkl=1,6
                   sta(nok,kkkkl)=aa(iijo,kkkkl)
                  enddo
              endif
            ELSEIF(giss.gt.3.1d0.and.giss.lt.3.6d0)THEN
              if(aa(iijo,2).gt.2.45d0.and.aa(iijo,2).lt.4.1d0)then
                 nok=nok+1
                 sel(nok,1)=ast(iijo,1)
                 sel(nok,2)=ast(iijo,2)
                 do kkkkl=1,6
                   sta(nok,kkkkl)=aa(iijo,kkkkl)
                 enddo
              endif
            ELSEIF(giss.le.3.1d0.and.giss.gt.1.5d0)THEN
              if(aa(iijo,2).lt.3.7d0.and.aa(iijo,2).gt.0.9d0)then
                 nok=nok+1
                 sel(nok,1)=ast(iijo,1)
                 sel(nok,2)=ast(iijo,2)
                 do kkkkl=1,6
                   sta(nok,kkkkl)=aa(iijo,kkkkl)
                 enddo
              endif
            ELSEIF(giss.le.1.5d0)THEN
              if(aa(iijo,2).lt.2.75d0)then
                 nok=nok+1
                 sel(nok,1)=ast(iijo,1)
                 sel(nok,2)=ast(iijo,2)
                 do kkkkl=1,6
                   sta(nok,kkkkl)=aa(iijo,kkkkl)
                 enddo
              endif
            ENDIF
           ENDIF
          endif
         enddo
         if(nok.ge.1) goto 2501
        ENDDO
2501    continue
        DO nonoo=1,300
         IF(fiss.ge.0.15d0)THEN 
          fup=5.0d0
          if(giss.lt.3.5d0)then                          !giants
            if(tiss.ge.9000.d0)then
              fdo=-5.0
            elseif(tiss.lt.9000.d0.and.tiss.ge.5100.d0)then
              fdo=fiss-0.2d0-dble(nonoo)*0.01d0
            elseif(tiss.lt.3300.d0)then
              fdo=fiss-0.2d0-dble(nonoo)*0.01d0
            else
             fdo=fiss-0.09d0-dble(nonoo)*0.01d0  !aqui esta el grueso
	    endif
          else                                         !dwarfs
            if(tiss.ge.12000.d0)then
              fdo=-5.0d0
            elseif(tiss.lt.12000.d0.and.tiss.ge.6800.d0)then
              fdo=fiss-0.17d0-dble(nonoo)*0.01d0
            elseif(tiss.lt.4750.d0)then
              fdo=fiss-0.19d0-dble(nonoo)*0.01d0
            else
             fdo=fiss-0.08d0-dble(nonoo)*0.01d0 !aqui esta el grueso
	    endif
          endif
	 ELSEIF(fiss.lt.0.15d0.and.fiss.gt.-0.15d0)THEN
          if(giss.lt.3.5d0)then                          !giants
            if(tiss.ge.8500.d0)then
              fdo=-5.0d0
              fup=5.0d0
            elseif(tiss.lt.8500.d0.and.tiss.ge.5400.d0)then
              fdo=fiss-0.18d0-dble(nonoo)*0.01d0
              fup=fiss+0.18d0+dble(nonoo)*0.01d0
            else
             fdo=fiss-0.13d0-dble(nonoo)*0.01d0
             fup=fiss+0.13d0+dble(nonoo)*0.01d0
            endif
          else                                         !dwarfs
            if(tiss.ge.12000.d0)then
              fdo=-5.0d0
              fup=5.0d0
            elseif(tiss.lt.12000.d0.and.tiss.ge.7200.d0)then
             fup=fiss+0.18d0+dble(nonoo)*0.02d0
             fdo=fiss-0.18d0-dble(nonoo)*0.02d0
            else             
             fdo=fiss-0.12d0-dble(nonoo)*0.01d0
             fup=fiss+0.11d0+dble(nonoo)*0.01d0
            endif
          endif
	 ELSEIF(fiss.le.-0.15d0.and.fiss.gt.-0.55d0)THEN
          if(giss.lt.3.5d0)then                          !giants
            if(tiss.ge.9500.d0)then
             fdo=-5.0d0
             fup=+5.0d0
            else
             fup=fiss+0.13d0+dble(nonoo)*0.01d0
             fdo=fiss-0.18d0-dble(nonoo)*0.02d0
            endif
          else                                         !dwarfs
            if(tiss.ge.16000.d0)then
             fdo=-5.0d0
             fup=+5.0d0
            else
             fup=fiss+0.13d0+dble(nonoo)*0.01d0
             fdo=fiss-0.18d0-dble(nonoo)*0.02d0
            endif
          endif
	 ELSEIF(fiss.le.-0.55d0.and.fiss.gt.-0.9d0)THEN
          if(giss.lt.3.5d0)then                          !giants
            if(tiss.ge.9500.d0)then
             fdo=-5.0d0
             fup=+5.0d0
	    elseif(tiss.gt.5000.d0.and.tiss.lt.9000.d0)then  
             fup=fiss+0.12d0+dble(nonoo)*0.01d0 !es un remiendo
             fdo=fiss-0.37d0-dble(nonoo)*0.02d0 !para las hot HB
            else
             fup=fiss+0.12d0+dble(nonoo)*0.01d0
             fdo=fiss-0.365d0-dble(nonoo)*0.02d0
            endif
          else                                         !dwarfs
            if(tiss.ge.16000.d0)then
             fdo=-5.0d0
             fup=+5.0d0
            elseif(tiss.lt.6150.d0.and.tiss.gt.5950.d0)then  
             fup=fiss+0.14d0+dble(nonoo)*0.01d0 !es un remiendo
             fdo=fiss-0.365d0-dble(nonoo)*0.02d0 !para el Toff
            else
             fup=fiss+0.12d0+dble(nonoo)*0.01d0
             fdo=fiss-0.365d0-dble(nonoo)*0.02d0
            endif
          endif
	 ELSEIF(fiss.le.-0.9d0)THEN
             fup=fiss+0.20d0+dble(nonoo)*0.02d0
             fdo=fiss-0.30d0-dble(nonoo)*0.02d0
      	 ELSE
             fup=fiss+0.09d0+dble(nonoo)*0.02d0
             fdo=fiss-0.09d0-dble(nonoo)*0.02d0
         ENDIF
         nokf=0
         do iijo=1,nok
          if(sta(iijo,3).lt.fup.and.sta(iijo,3).gt.fdo)then
               mf(iijo)=iijo
               nokf=nokf+1
          else
               mf(iijo)=0
          endif
         enddo
         if(nokf.ge.1)goto 2601
        ENDDO 
2601    continue
	DO nonoo=1,300
         if(giss.ge.3.5d0.and.tiss.lt.6150.d0.and.tiss.gt.5950.d0.and.
     &   fiss.lt.-0.55d0)then
	  vup=viss+1.25d0+dble(nonoo)*0.25d0
	  vdo=viss-1.25d0-dble(nonoo)*0.25d0
          gup=5.3d0
          gdo=3.6d0-dble(nonoo)*0.10d0
	 elseif(giss.lt.3.5d0.and.tiss.gt.5000.d0.and.fiss.lt.-0.55d0)then 
	  vup=viss+1.85d0+dble(nonoo)*0.25d0
	  vdo=viss-1.85d0-dble(nonoo)*0.25d0
          gup=giss+0.8d0+dble(nonoo)*0.10d0
          gdo=giss-0.8d0-dble(nonoo)*0.10d0
	 else
	  vup=viss+0.5d0+dble(nonoo)*0.25d0
	  vdo=viss-0.5d0-dble(nonoo)*0.25d0
          gup=giss+0.5d0+dble(nonoo)*0.10d0
          gdo=giss-0.5d0-dble(nonoo)*0.10d0
         endif
	 nokfv=0
       do kin=1,nok
       if(sta(kin,4).lt.vup.and.sta(kin,4).gt.vdo.and.mf(kin).ne.0)then
           if(sta(kin,2).lt.gup.and.sta(kin,2).gt.gdo)then
	       mfv(kin)=kin
	       nokfv=nokfv+1
	   endif
       else
	   mfv(kin)=0
       endif
       enddo
	 if(nokfv.ge.1)goto 27
	ENDDO	
27	continue
        mmfvg=0        
        do nonoo=1,nok
         if(mfv(nonoo).ne.0)then
           mmfvg=mmfvg+1
           wn(mmfvg,1)=sel(nonoo,1)                                               
           wn(mmfvg,2)=sel(nonoo,2)
           fstar(mmfvg,3)=sta(nonoo,3)
           fstar(mmfvg,2)=sta(nonoo,2)
           fstar(mmfvg,1)=sta(nonoo,1)
         endif
        enddo
	close(90)
c1971	format(2(f6.0,1x),2(f5.2,1x),I3,1x,f6.0,1x,2(f6.2,1x))
c        write(59,1971)tup,tdo,fup,fdo,mmfvg,tiss,giss,viss
	call lstdor(mmfvg,wn,xxb,xxr)
	aloe(1)=0.d0
	aloe(2)=0.d0
	aloe(3)=0.d0
	DO lina=1,mmfvg
		aloe(1)=aloe(1)+fstar(lina,1)
		aloe(2)=aloe(2)+fstar(lina,2)
		aloe(3)=aloe(3)+fstar(lina,3)
c	write(43,'(f6.0,1x,f6.3,1x,f6.3)')(fstar(lina,kakak),kakak=1,3)
	ENDDO
	do jijo=1,3
		aloe(jijo)=aloe(jijo)/dble(mmfvg)
	enddo
	fra4000=0.0d0
	fra4436=0.0d0
	do j4=1,1107
	 if(xxb(j4,1).gt.3855.4d0.and.xxb(j4,1).lt.3950.2d0)then
	  fra4000=fra4000+xxb(j4,2)
	 endif
	 if(xxb(j4,1).gt.4421.0d0.and.xxb(j4,1).lt.4451.5d0)then
	  fra4436=fra4436+xxb(j4,2)
	 endif
	enddo
	fra4000=fra4000/fra4436	 
c	write(43,'(A14,f6.0,1x,f6.3,1x,f6.3)')'              ',
c     &  (aloe(kakak),kakak=1,3)
	return
	end
c************************************************************************
        SUBROUTINE lstdor(mmm,wiii,xb,xr)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (nstm=1999)
	DIMENSION sb(650,1120,2),xb(1120,2)
	DIMENSION sr(650,1120,2),xr(1120,2)
c	double precision sb(650,1120,2),xb(1120,2)
c	double precision sr(650,1120,2),xr(1120,2)
	DIMENSION wiii(650,2)
        DIMENSION stars(nstm,2,1107,2) !Jones b&v STARS ARRAY
	CHARACTER*20 star,wiii
        COMMON/jonsta/stars !Jones b&v STARS ARRAY:commons:a,lstdor(en lbusc.f)
	COMMON/lstd/nsss1,nsss2
	COMMON/clstd/star(650,2)
c	No variar la longitud A10
	DO lina=1,mmm
	    do kin=1,nsss1
	      if(wiii(lina,1).eq.star(kin,1))then
c	      open(99,file=star(kin,1),status='old')
	      do jij=1,1107
c	         read(99,*)(sb(lina,jij,kikii),kikii=1,2)
	        sb(lina,jij,1)=stars(kin,1,jij,1)
	        sb(lina,jij,2)=stars(kin,1,jij,2)
	      enddo
c	      close(99)
	      endif
	    enddo
	ENDDO
	DO lina=1,mmm
	    do kin=1,nsss2
	      if(wiii(lina,2).eq.star(kin,2))then
c	      open(99,file=star(kin,2),status='old')
	      do jij=1,1107
c	         read(99,*)(sr(lina,jij,kikii),kikii=1,2)
	        sr(lina,jij,1)=stars(kin,2,jij,1)
	        sr(lina,jij,2)=stars(kin,2,jij,2)
	      enddo
c	      close(99)
	      endif
	    enddo
	ENDDO
	DO jijo=1,1107
	xb(jijo,2)=0.d0
	xb(jijo,1)=0.d0
	do lijo=1,mmm
c if(sb(l,j,1).eq.sb(lijo,jijo,1))then
	if(sb(lijo,jijo,1).gt.3855.4d0.and.sb(lijo,jijo,1).lt.4476.5d0)then
	 xb(jijo,2)=xb(jijo,2)+sb(lijo,jijo,2)
	else
	 xb(jijo,2)=1.d0
	endif
c else
c write(*,*)'AZUL_ERROR',sb(1,j,1),sb(2,j,1)
c endif
	enddo
	ENDDO	
	DO jijo=1,1107
	  xr(jijo,2)=0.d0
	  do lijo=1,mmm
	if(sr(lijo,jijo,1).gt.4794.9d0.and.sr(lijo,jijo,1).lt.5465.1d0)then
		xr(jijo,2)=xr(jijo,2)+sr(lijo,jijo,2)
	else
		xr(jijo,2)=1.d0
	endif
	  enddo
	ENDDO
	do jijo=1,1107
		xb(jijo,2)=xb(jijo,2)/dble(mmm)	
		xr(jijo,2)=xr(jijo,2)/dble(mmm)
		xb(jijo,1)=sb(1,jijo,1)
		xr(jijo,1)=sr(1,jijo,1)
	enddo	
	return
	END

cc*******SUBRUTINA CAJAS_TEFF********************************************
        SUBROUTINE cajas(X0,grav,f,Y0)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c        DOUBLE PRECISION X0,Y0
        PARAMETER (NDEG=3)
        DIMENSION A(NDEG+1)
c        DOUBLE PRECISION A(NDEG+1)
c	if(grav.lt.4.0d0.and.grav.gt.3.0d0)then
c	tiene que haber 2 grav, el pedido y el de busqueda	
	IF(grav.lt.3.5d0)THEN !giants
         if(X0.ge.5700.d0)then
          t1=5700.d0
          d1=334.d0
          t2=30000.d0
          d2=10000.d0
	  Y0=d1+(X0-t1)*((d2-d1)/(t2-t1))
	 else
          A(01)=  5054.31198d0
          A(02)= -2.24320061d0
          A(03)=  0.000281794717d0
          A(04)= -5.88359573d-9
          call orden3(A,X0,Y0)	 
         endif
c	  call recta(5700.d0,334.d0,30000.d0,10000.d0,X0,Y0)
c         endif
        ELSE !dwarfs
         if(f.lt.0.15d0.and.f.gt.-0.15d0)then
          if(X0.gt.6750.d0)then
           t1=6750.d0
           d1=190.d0
           t2=30000.d0
           d2=10000.d0
	   Y0=d1+(X0-t1)*((d2-d1)/(t2-t1))
          else
           A(01)=  2146.65379d0
           A(02)= -0.744303943d0
           A(03)=  7.58224914d-5
           A(04)= -1.2563799d-9
           call orden3(A,X0,Y0)
          endif
         elseif(f.lt.-1.0d0)then
          if(X0.gt.6050.d0)then
           t1=6050.d0
           d1=250.d0
           t2=30000.d0
           d2=10000.d0
	   Y0=d1+(X0-t1)*((d2-d1)/(t2-t1))
          else
           A(01)=  3547.11311d0
           A(02)= -1.23146692d0
           A(03)=  0.000127678143d0
           A(04)= -2.35093402d-9
           call orden3(A,X0,Y0)
          endif
         else
          if(X0.gt.6300.d0)then
           t1=6300.d0
           d1=247.d0
           t2=30000.d0
           d2=10000.d0
	   Y0=d1+(X0-t1)*((d2-d1)/(t2-t1))
          else
           A(01)=  3697.24345d0
           A(02)= -1.2970467d0
           A(03)=  0.000134677364d0
           A(04)= -2.49797294d-9
           call orden3(A,X0,Y0)
          endif
	 endif
	ENDIF
	Y0=Y0*0.5d0
	return
	end

	SUBROUTINE orden3(A,X0,Y0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c        DOUBLE PRECISION X0,Y0
        PARAMETER (NDEG=3)
        DIMENSION A(NDEG+1)
c        DOUBLE PRECISION A(NDEG+1)
        Y0=A(NDEG+1)
c       IF(NDEG.GT.0)THEN
          DO K=NDEG,1,-1
            Y0=Y0*X0+A(K)
          END DO
c       END IF
	return
        END
