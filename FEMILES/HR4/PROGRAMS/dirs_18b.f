	SUBROUTINE dirs(aaNa,aaFe)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 ctonti,ctun,ctunn
      CHARACTER*80 jmiku
      COMMON/ab/iaba
	COMMON/jmik/jmiku(20),lenj(20) !nombre dirs estrellas
	COMMON/jmik2/ctunn,lctunn
	jmiku(1)='../D/SL_BASE/STARS_4/'
	jmiku(2)='../D/SL_BASE/STARS_5/'
	IF(iaba.eq.1.or.iaba.eq.2)THEN
	 if(aaNa.ge.0.0d0)then
	  ctun='p'
	 else
	  ctun='m'
	 endif
	 write(ctonti,'(F5.2)')aaNa
	 write(ctun(2:2),'(A1)')ctonti(2:2)
	 write(ctun(3:3),'(A1)')ctonti(4:4)
	 write(ctonti,'(F5.2)')aaFe
	 if(aaFe.ge.0.0d0)then
	  write(ctun(4:8),'(A5)')'_aFep'
	 else
	  write(ctun(4:8),'(A5)')'_aFem'
	 endif
	 write(ctonti,'(F5.2)')aaFe
	 write(ctun(9:9),'(A1)')ctonti(2:2)
	 write(ctun(10:10),'(A1)')ctonti(4:4)
	 if(iaba.eq.2)then
	  write(ctun(11:21),'(A10)')'_MgFeCORR/'
	  k8=30
	  k9=60
60	  format(A11)
	 elseif(iaba.eq.1)then
	  write(ctun(11:11),'(A1)')'/'
	  k8=30
	  k9=40
61	  format(A30)
	 endif
	 jmiku(3)='../D/SL_NaFe_aFe/STARS_M_NaFe'
	 jmiku(4)='../D/SL_NaFe_aFe/STARS_C_NaFe'
	 jmiku(5)='../D/SL_NaFe_aFe/STARS_I_NaFe'
	 jmiku(6)='../D/SL_NaFe_aFe/STARS_K_NaFe'
	 jmiku(7)='../D/SL_NaFe_aFe/STARS_U_NaFe'
	 jmiku(7)='../D/SL_NaFe_aFe/STARS_U_NaFe'
	 IF(iaba.eq.1)then
	  do iku=3,7
	   write(jmiku(iku)(k8:k9),60)ctun
	  enddo
	 ELSEIF(iaba.eq.2)THEN
	  do iku=3,7
	   write(jmiku(iku)(k8:k9),61)ctun
	  enddo
	 ENDIF
	 lctun=lnblnk(ctun)
	 ctunn='_NaFe'//ctun(1:lctun-1)
c	 write(ctunn(1:lctun-1)
	 lctunn=lnblnk(ctunn)
	ELSE !aplica a iaba=0,9
	 jmiku(3)='../D/SL_BASE/STARS_M/'
	 jmiku(4)='../D/SL_BASE/STARS_C/'
	 jmiku(5)='../D/SL_BASE/STARS_I/'
	 jmiku(6)='../D/SL_BASE/STARS_K/'
	 jmiku(7)='../D/SL_BASE/STARS_U/'
	 lctun=0
	 lctunn=0
	 ctunn=''
	ENDIF
	do i6=1,20
	 lenj(i6)=lnblnk(jmiku(i6))
	enddo
      RETURN
      END
