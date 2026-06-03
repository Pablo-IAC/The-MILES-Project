	SUBROUTINE a_IRTF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c IRTF
	PARAMETER (npxk=15387) !IRTF #pixels in the spectra
	PARAMETER (nstk=300) !IRTF #stars in the library
	DIMENSION wfk(npxk,2),wwffko(npxk,2),wwffky(npxk,2) !IRTF
	DIMENSION a9ko(3),a9ky(3),a90k(3) !IRTF
	DIMENSION cabk(150,4) !cabecera IRTF
c
	DIMENSION bb(4),sson(4)
	DIMENSION z(12,150000,16) !para TERAMO, Z00(15) common/eststu
	CHARACTER*7 cmgfe
	CHARACTER*8 cisoc
c IRTF
	CHARACTER*80 stark
	CHARACTER*80 cabk
c
	CHARACTER*80 chatom,chato
	CHARACTER*1 shape
	CHARACTER*9 anom
	CHARACTER*80 jmiku,ctunn
	COMMON/aaaa/izold(15),izyng(15),ttold(99),ttyng(99),pizold(15)
     & ,pizyng(15),zi(150),zalpha(22),nnIMF,nalpha,nzold,ntold
     & ,nzyng,ntyng
c IRTF
	COMMON/aamcok/aak(nstk,4) !IRTF
      COMMON/hrk1/nstark !IRTF
      COMMON/hrk2/stark(nstk) !IRTF
      COMMON/kilsta/starks(nstk,npxk,2) !IRTF STARS ARRAY: common routines:a,sigmam
      COMMON/escK/a90k !a,e
c
	COMMON/labeli/cisoc,cmgfe !aaa,a
	COMMON/matina/a(4,4)
	COMMON/fezsol/fez(15)
      COMMON/tefsca/nghb
	COMMON/lmu0/bicl,bicp,bich
	COMMON/lmu1/ZMU
	COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3
	COMMON/lshape/shape
	COMMON/lson/sson
	COMMON/stus/ZMASA,tiso(75),num(15),numtis    
	COMMON/zsolo/z
	COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
	COMMON/KEYS/keypT(18)
      COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
	COMMON/REMNAN/re(15,2,99),nr(15)
      COMMON/NC/ncha !indice tipo isocronas
	COMMON/qua/quaqn,quaq10,quaq05,quaq15
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk
      COMMON/ab/iaba
	COMMON/jmik/jmiku(20),lenj(20) !nombre dirs estrellas
	COMMON/jmik2/ctunn,lctunn
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
246   format(A6)
	OPEN(98,FILE='../D/INPUT/PARAM_IRTF',STATUS='OLD') !IRTF
      nstark=0
	do k=1,99999
       read(98,*,end=25)anom,(aak(k,l),l=1,4)
c      Aplicamos escala GonzalezHernandezBonifacio09?
       if(nghb.eq.1)then
        if(aak(k,1).gt.3750.0d0.and.aak(k,1).lt.7500.0d0)then
         aak(k,1)=aak(k,1)+(-116.d0+0.0312d0*aak(k,1))
        elseif(aak(k,1).le.3750.0d0)then
         aak(k,1)=aak(k,1)+(-116.d0+0.0312d0*3750.0d0)
        elseif(aak(k,1).ge.7500.0d0)then
         aak(k,1)=aak(k,1)+(-116.d0+0.0312d0*7500.0d0)
        endif
       endif
       stark(k)='./STARS_IRTF/'//anom !IRTF
       nstark=nstark+1
       stark(k)=jmiku(6)
	 write(stark(k)(lenj(6)+1:lenj(6)+1+6),246)anom
      enddo
25    CLOSE(98)
c LECTURA DE FICHEROS ASCII DE LIBRERIA ESPECTRAL
c SE CARGA EN MATRIZ:starks(nstark,npxk,2)
	DO ios=1,nstark !IRTF
	 open(99,file=stark(ios),status='old')
	 do iop=1,npxk
	   read(99,*)starks(ios,iop,1),starks(ios,iop,2)
	   if(ios.eq.1)wfk(iop,1)=starks(ios,iop,1)
	 enddo
	 close(99)
	ENDDO
ccccccccccccccccccccc
	open(99,file='../D/INPUT/GALAXY_CABECERAS_IRTF',status='old') !NGSL
	do k=1,61
		read(99,'(A80)')cabk(k,4)
	enddo
	close(99)
ccccccccccccccccccccc
	DO NNN=1,nnIMF
	 ZMU=zi(NNN)+1.0d0
	 if(shape.eq.'b')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/Kbi'//chato !IRTF
	 elseif(shape.eq.'u')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/Kun'//chato !IRTF
	 elseif(shape.eq.'k')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/Kku'//chato !IRTF
	 elseif(shape.eq.'r')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/Kkb'//chato !IRTF
	 elseif(shape.eq.'c')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/Kch'//chato !IRTF
	 elseif(shape.eq.'f')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/KFe'//chato
	 elseif(shape.eq.'x')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/K/KFx'//chato
	 elseif(shape.eq.'l')then
	  int1=int(abs(ZMU1*10))
	  int2=int(abs(ZMU2*10))
	  if(ZMU1.ge.0.0d0)then
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/K/Kp0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/K/Kp'//chato
	   endif
	  else
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/K/Km0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/K/Km'//chato
	   endif
	  endif
	  if(ZMU2.ge.0.0d0)then
	   if(int2.lt.10)then
	    write(chato(1:1),'(I1)')int2
  	    write(chatom(13:15),'(A3)')'p0'//chato
	   else
	    write(chato(1:2),'(I2)')int2
  	    write(chatom(13:15),'(A3)')'p'//chato
	   endif
	  else
	   if(int2.lt.10)then
	    write(chato(1:1),'(I1)')int2
  	    write(chatom(13:15),'(A3)')'m0'//chato
	   else
	    write(chato(1:2),'(I2)')int2
  	    write(chatom(13:15),'(A3)')'m'//chato
	   endif
	  endif
	 endif
c LaBarbera 3-segments
            ZMU1=ZMU1+1.0d0
            ZMU2=ZMU2+1.0d0
            ZMU3=ZMU3+1.0d0
calculo coeficientes IMF bimodal
	 bb(1)=bicp**(-ZMU+1.0d0)
	 bb(2)=0.0d0
	 bb(3)=bich**(-ZMU+1.0d0)
	 bb(4)=(-ZMU+1.0d0)*bich**(-ZMU)
	 do iih=1,4
	  sson(iih)=0.0d0
	  do jjh=1,4
	   sson(iih)=sson(iih)+a(iih,jjh)*bb(jjh)
	  enddo
	 enddo
c*******COMIENZO DEL BUCLE DEL PORCENTAGE******************************
	 DO ialpha=1,nalpha      
	  alpha=zalpha(ialpha)
        if(alpha.eq.1.00d0)then
         ntyng=1
         ttyng(1)=ttold(1)
         nzyng=1
         izyng(1)=izold(1)
        endif
c*******COMIENZO DEL BUCLE DE LA SSPo**********************************
	  DO izo=1,nzold
	   if(pizold(izold(izo)).ge.0.0d0)then
	    write(chato(1:6),'(A,A,f4.2)')'Z','p',pizold(izold(izo))
	   else
	    write(chato(1:6),'(A,A,f4.2)')'Z','m',abs(pizold(izold(izo)))
	   endif
	   fehf=pizold(izold(izo)) ![Fe/H] que pasara a STU
  	   write(chatom(16:21),'(A6)')chato
	   DO iao=1,ntold
	    ageo=ttold(iao)
          call STU_IRTF(izold(izo),ageo,BETAo,realeo,ZMISOo,
     &wwffko,a9ko,f9mo,qumo,qu10o,qu05o,qu15o)
	    attao=realeo
	    if(attao.lt.10.0d0)then
	     write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	    else
	     write(chato(1:8),'(A,f7.4)')'T',attao
	    endif
  	    write(chatom(22:29),'(A8)')chato
c         name isoc+alpha
  	    write(chatom(30:37),'(A8)')cisoc
  	    write(chatom(38:44),'(A7)')cmgfe
	    chatom=chatom(1:44)
c*******COMIENZO DEL BUCLE DE LA SSPy**********************************
	    DO izy=1,nzyng
	     if(alpha.ne.1.0d0)then
  	      write(chatom(45:49),'(A,f4.2)')'a',alpha
	      if(pizyng(izyng(izy)).ge.0.0d0)then
	       write(chato(1:6),'(A,A,f4.2)')'Z','p',pizyng(izyng(izy))
	      else
	       write(chato(1:6),'(A,A,f4.2)')'Z','m',
     &abs(pizyng(izyng(izy)))
	      endif
  	      write(chatom(50:55),'(A6)')chato
	     endif
	     DO iay=1,ntyng
		agey=ttyng(iay)
	      if(agey.gt.ageo.and.alpha.ne.1.00d0)goto 2325
	      IF(alpha.ne.1.0d0)THEN
	       call STU_IRTF(izyng(izy),agey,BETAy,realey,ZMISOy,
     &wwffky,a9ky,f9my,qumy,qu10y,qu05y,qu15y)
	       attao=realey
	       if(attao.lt.10.0d0)then
	        write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	       else
	        write(chato(1:8),'(A,f7.4)')'T',attao
	       endif
  	       write(chatom(56:63),'(A8)')chato
	       chatom=chatom(1:63)
     	      ELSE
	       chatom=chatom(1:44)
     	       BETAy=0.0d0
     	       realey=0.0d0
     	       do ini=1,npxk !IRTF
     	  	  wwffky(ini,2)=0.0d0
      	 enddo
     	       f9my=0.0d0 !IRTF
     	       do ini=1,3
     	  	  a9ky(ini)=99.0d0
     	       enddo
     	      ENDIF
     	      do nt=1,3 !IRTF
	       a90k(nt)=(alpha*BETAo*a9ko(nt)+(1.0d0-alpha)*BETAy*
     &a9ky(nt))/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      enddo
c      write(38888,'(A60,1x,F7.0,1x,2(F7.3,1x))')chatom,(a90k(nt),nt=1,3)
c           Parametros calidad IRTF
     	      quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
            do io=1,npxk !IRTF
             wfk(io,2)=alpha*BETAo*wwffko(io,2)+(1.d0-alpha)*
     &BETAy*wwffky(io,2)
            enddo
cccccccc
c_NaFep06_aFep04_MgFeCORR	 
	 lchat=1+lnblnk(chatom)
	 if(iaba.eq.1)then
	  write(chatom(lchat:lchat+lctunn),'(A15)')ctunn
	 elseif(iaba.eq.2)then
	  write(chatom(lchat:lchat+lctunn),'(A24)')ctunn
	 endif
cccccccc
	      OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	      CLOSE(77,STATUS='DELETE')
8979	      OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	      do nonoo=1,61
	       write(77,'(A80)')cabk(nonoo,4)
	      enddo
	      do lup=1,npxk
		 write(77,*)(wfk(lup,kini),kini=1,2)
	      enddo
	      close(77)
	      write(*,*)chatom	      
	      call e_IRTF(alpha,izold(izo),realeo,izyng(izy),realey)
2325	      continue
	     ENDDO !Cerrado loop Tyng
	    ENDDO !Cerrado loop Zyng
	   ENDDO !Cerrado loop Told
	  ENDDO !Cerrado loop Zold
	 ENDDO !Cerrado loop alpha
	ENDDO !Cerrado loop IMF slope
	RETURN
6336	END




