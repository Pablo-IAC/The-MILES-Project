	SUBROUTINE a_NGSL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c NGSL
	PARAMETER (npxu=6192) !NGSL #pixels in the spectra
	PARAMETER (nstu=1100) !NGSL #stars in the library
	DIMENSION wfu(npxu,2),wwffuo(npxu,2),wwffuy(npxu,2) !NGSL
	DIMENSION wfus(npxu,2),wwffuso(npxu,2),wwffusy(npxu,2) !NGSL
	DIMENSION a9uo(3),a9uy(3),a90u(3) !NGSL
	DIMENSION cabu(150,4) !cabecera NGSL
c
	DIMENSION bb(4),sson(4)
	DIMENSION z(12,150000,16) !para TERAMO, Z00(15) common/eststu
	CHARACTER*7 cmgfe
	CHARACTER*8 cisoc
c NGSL
	CHARACTER*80 staru
	CHARACTER*80 cabu
c
	CHARACTER*80 chatom,chato,chatoms
	CHARACTER*1 shape
	CHARACTER*9 anom
	CHARACTER*80 jmiku,ctunn
	COMMON/aaaa/izold(15),izyng(15),ttold(99),ttyng(99),pizold(15)
     & ,pizyng(15),zi(150),zalpha(22),nnIMF,nalpha,nzold,ntold
     & ,nzyng,ntyng
c NGSL
	COMMON/aamcou/aau(nstu,4) !NGSL
      COMMON/hru1/nstaru !NGSL
      COMMON/hru2/staru(nstu) !NGSL
      COMMON/uilsta/starus(nstu,npxu,2) !NGSL STARS ARRAY: common routines:a,sigmam
      COMMON/escU/a90u !a,e
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
247   format(A6)
	OPEN(98,FILE='../D/INPUT/PARAM_NGSL',STATUS='OLD') !NGSL
      nstaru=0
	do k=1,99999
       read(98,*,end=25)anom,(aau(k,l),l=1,4)
c      Aplicamos escala GonzalezHernandezBonifacio09?
       if(nghb.eq.1)then
        if(aau(k,1).gt.3750.0d0.and.aau(k,1).lt.7500.0d0)then
         aau(k,1)=aau(k,1)+(-116.d0+0.0312d0*aau(k,1))
        elseif(aau(k,1).le.3750.0d0)then
         aau(k,1)=aau(k,1)+(-116.d0+0.0312d0*3750.0d0)
        elseif(aau(k,1).ge.7500.0d0)then
         aau(k,1)=aau(k,1)+(-116.d0+0.0312d0*7500.0d0)
        endif
       endif
c       staru(k)='./STARS_NGSL/'//anom !NGSL
       staru(k)=jmiku(7)
	 write(staru(k)(lenj(7)+1:lenj(7)+1+6),247)anom
       nstaru=nstaru+1
      enddo
25    CLOSE(98)
c LECTURA DE FICHEROS ASCII DE LIBRERIA ESPECTRAL
c SE CARGA EN MATRIZ:starus(nstaru,npxu,2)
	DO ios=1,nstaru !NGSL
	 open(99,file=staru(ios),status='old')
	 do iop=1,npxu
	   read(99,*)starus(ios,iop,1),starus(ios,iop,2)
	   if(ios.eq.1)wfu(iop,1)=starus(ios,iop,1)
	 enddo
	 close(99)
	ENDDO
ccccccccccccccccccccc
	open(99,file='../D/INPUT/GALAXY_CABECERAS_NGSL',status='old') !NGSL
	do k=1,61
		read(99,'(A80)')cabu(k,4)
	enddo
	close(99)
ccccccccccccccccccccc
	DO NNN=1,nnIMF
	 ZMU=zi(NNN)+1.0d0
	 if(shape.eq.'b')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/Ubi'//chato !IRTF
	 elseif(shape.eq.'u')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/Uun'//chato !IRTF
	 elseif(shape.eq.'k')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/Uku'//chato !IRTF
	 elseif(shape.eq.'r')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/Ukb'//chato !IRTF
	 elseif(shape.eq.'c')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/Uch'//chato !IRTF
	 elseif(shape.eq.'f')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/UFe'//chato
	 elseif(shape.eq.'x')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/U/UFx'//chato
	 elseif(shape.eq.'l')then
	  int1=int(abs(ZMU1*10))
	  int2=int(abs(ZMU2*10))
	  if(ZMU1.ge.0.0d0)then
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/U/Up0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/U/Up'//chato
	   endif
	  else
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/U/Um0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/U/Um'//chato
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
          call STU_NGSL(izold(izo),ageo,BETAo,realeo,ZMISOo,wwffuo,
     &wwffuso,a9uo,f9mo,qumo,qu10o,qu05o,qu15o)
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
	       call STU_NGSL(izyng(izy),agey,BETAy,realey,ZMISOy,wwffuy
     &,wwffusy,a9uy,f9my,qumy,qu10y,qu05y,qu15y)
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
     	       do ini=1,npxu !NGSL
     	  	  wwffuy(ini,2)=0.0d0
     	  	  wwffusy(ini,2)=0.0d0
      	       enddo
     	       f9my=0.0d0 !NGSL
     	       do ini=1,3
     	  	  a9uy(ini)=99.0d0
     	       enddo
     	      ENDIF
     	      do nt=1,3 !NGSL
	       a90u(nt)=(alpha*BETAo*a9uo(nt)+(1.0d0-alpha)*BETAy*
     &a9uy(nt))/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      enddo
c           Parametros calidad NGSL
     	      quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
            do io=1,npxu !NGSL
             wfu(io,2)=alpha*BETAo*wwffuo(io,2)+(1.d0-alpha)*
     &BETAy*wwffuy(io,2)
             wfus(io,2)=alpha*BETAo*wwffuso(io,2)+(1.d0-alpha)*
     &BETAy*wwffusy(io,2)
             wfus(io,1)=wfu(io,1)
            enddo
cccccccc
c_NaFep06_aFep04_MgFeCORR	 
	 lchat=1+lnblnk(chatom)
	 if(iaba.eq.1)then
	  write(chatom(lchat:lchat+lctunn),'(A15)')ctunn
	 elseif(iaba.eq.2)then
	  write(chatom(lchat:lchat+lctunn),'(A24)')ctunn
	 endif
         chatoms=chatom	    
	 lchat=1+lnblnk(chatoms)
c	 write(chatoms(lchat:lchat+4),'(A4)')'_sbf'
	 write(chatoms(lchat:lchat+4),'(A4)')'_var'
cccccccc
	      OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	      CLOSE(77,STATUS='DELETE')
8979	      OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	      do nonoo=1,61
	       write(77,'(A80)')cabu(nonoo,4)
	      enddo
	      do lup=1,npxu
		 write(77,*)(wfu(lup,kini),kini=1,2)
	      enddo
	      close(77)
	      OPEN(77,FILE=chatoms,STATUS='OLD',ERR=8989)
	      CLOSE(77,STATUS='DELETE')
8989	      OPEN(77,IOSTAT=IOS,FILE=chatoms,STATUS='NEW')
	      do nonoo=1,61
	       write(77,'(A80)')cabu(nonoo,4)
	      enddo
	      do lup=1,npxu
c		 write(77,*)wfus(lup,1),wfus(lup,2)/wfu(lup,2)
		 write(77,*)wfus(lup,1),wfus(lup,2)
	      enddo
	      close(77)
	      write(*,*)chatom,chatoms	      
	      call e_NGSL(alpha,izold(izo),realeo,izyng(izy),realey)
2325	      continue
	     ENDDO !Cerrado loop Tyng
	    ENDDO !Cerrado loop Zyng
	   ENDDO !Cerrado loop Told
	  ENDDO !Cerrado loop Zold
	 ENDDO !Cerrado loop alpha
	ENDDO !Cerrado loop IMF slope
	RETURN
6336	END




