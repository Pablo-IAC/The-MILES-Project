      SUBROUTINE a_XSL_276A
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c XSL
c      PARAMETER (npxi=106301) !XSL #pixels in the spectra
      PARAMETER (npxi=23578) !XSL #pixels in the spectra
      PARAMETER (nsti=1100) !XSL #stars in the library
      DIMENSION wfi(npxi,2),wwffio(npxi,2),wwffiy(npxi,2) !XSL
      DIMENSION wfis(npxi,2),wwffiso(npxi,2),wwffisy(npxi,2) !XSL
      DIMENSION a9io(3),a9iy(3),a90i(3) !XSL
      DIMENSION cabi(150,4) !cabecera XSL
      DIMENSION bb(4),sson(4)
      DIMENSION z(12,150000,16) !para TERAMO, Z00(15) common/eststu
      CHARACTER*7 cmgfe
      CHARACTER*8 cisoc
c XSL
      CHARACTER*80 stari
      CHARACTER*80 cabi
c
      CHARACTER*80 chatom,chato,chatoms
      CHARACTER*1 shape,cv(80),ansmin,ansmax !shape or SFR and galIMF output
      CHARACTER*80 chau
      CHARACTER*9 anom
      CHARACTER*80 jmiku,ctunn
      CHARACTER*25 fga_str
      CHARACTER*25 sfr_str
      CHARACTER*25 ZMLow_str
      CHARACTER*2 iru_str
      COMMON/aaaa/izold(15),izyng(15),ttold(99),ttyng(99),pizold(15)
     & ,pizyng(15),zi(150),zalpha(22),nnIMF,nalpha,nzold,ntold
     & ,nzyng,ntyng
c XSL
      COMMON/aamcoi/aai(nsti,4) !XSL
      COMMON/hri1/nstari !XSL
      COMMON/hri2/stari(nsti) !XSL
      COMMON/iilsta/staris(nsti,npxi,2) !XSL STARS ARRAY: common routines:a,sigmam
      COMMON/escI/a90i !a,e
c
      COMMON/labeli/cisoc,cmgfe !aaa,a
      COMMON/matina/a(4,4)
      COMMON/fezsol/fez(15)
      COMMON/tefsca/nghb
      COMMON/lmu0/bicl,bicp,bich
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
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
cc      COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !a,e
      COMMON/ab/iaba
      COMMON/jmik/jmiku(20),lenj(20) !nombre dirs estrellas
      COMMON/jmik2/ctunn,lctunn
      COMMON/IGIMF/yIGMFm(9999),yIGMFn(9999),igm !also in limf.f 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
245   format(A6)
      OPEN(98,FILE='../D/INPUT/PARAM_XSL',STATUS='OLD') !XSL
      nstari=0
      do k=1,99999
       read(98,*,end=25)anom,(aai(k,l),l=1,4)
c      Aplicamos escala GonzalezHernandezBonifacio09?
       if(nghb.eq.1)then
        if(aai(k,1).gt.3750.0d0.and.aai(k,1).lt.7500.0d0)then
         aai(k,1)=aai(k,1)+(-116.d0+0.0312d0*aai(k,1))
        elseif(aai(k,1).le.3750.0d0)then
         aai(k,1)=aai(k,1)+(-116.d0+0.0312d0*3750.0d0)
        elseif(aai(k,1).ge.7500.0d0)then
         aai(k,1)=aai(k,1)+(-116.d0+0.0312d0*7500.0d0)
        endif
       endif
       stari(k)=jmiku(5)
       write(stari(k)(lenj(5)+1:lenj(5)+1+6),245)anom
c       stari(k)='./STARS_XSL/'//anom !XSL
       nstari=nstari+1
      enddo
25    CLOSE(98)
c LECTURA DE FICHEROS ASCII DE LIBRERIA ESPECTRAL
c SE CARGA EN MATRIZ:staris(nstari,npxi,2)
	DO ios=1,nstari !XSL
	 open(99,file=stari(ios),status='old')
	 do iop=1,npxi
	   read(99,*)staris(ios,iop,1),staris(ios,iop,2)
	   if(ios.eq.1)wfi(iop,1)=staris(ios,iop,1)
	 enddo
	 close(99)
	ENDDO
ccccccccccccccccccccc
	open(99,file='../D/INPUT/GALAXY_CABECERAS_XSL_276A',status='old') !XSL
	do k=1,61
		read(99,'(A80)')cabi(k,4)
	enddo
	close(99)
ccccccccccccccccccccc
	DO NNN=1,nnIMF
	 if(shape.eq.'g')then
	     ZMU=zi(NNN)
	 else
	     ZMU=zi(NNN)+1.0d0
	 endif
	 if(shape.eq.'b')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xbi'//chato !XSL
	 elseif(shape.eq.'u')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xun'//chato !XSL
         elseif(shape.eq.'g')then
          if(ZMU.ge.0.d0)then
           write(chato(1:1),'(A1)')'p'
 	  else
           write(chato(1:1),'(A1)')'m'
	  endif
          write(chato(2:5),'(f4.2)')abs(ZMU)
          chatom='./OUT/X/Xg'//chato(1:5)
	 elseif(shape.eq.'t')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xba'//chato !MIUSC
	 elseif(shape.eq.'k')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xku'//chato !MIUSC
	 elseif(shape.eq.'r')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xkb'//chato !MIUSC
	 elseif(shape.eq.'c')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/Xch'//chato !MIUSC
	 elseif(shape.eq.'f')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/XFe'//chato
	 elseif(shape.eq.'x')then
	  write(chato(1:4),'(f4.2)')ZMU-1.0d0
	  chatom='./OUT/X/XFx'//chato
	 elseif(shape.eq.'l')then
	  int1=int(abs(ZMU1*10))
	  int2=int(abs(ZMU2*10))
	  if(ZMU1.ge.0.0d0)then
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/X/Xp0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/X/Xp'//chato
	   endif
	  else
	   if(int1.lt.10)then
	    write(chato(1:1),'(I1)')int1
	    chatom='./OUT/X/Xm0'//chato
	   else
	    write(chato(1:2),'(I2)')int1
	    chatom='./OUT/X/Xm'//chato
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
c     ==================== IGIMF OPTION ===========================
         if(shape.eq.'g')then
c	call system('./galIMF-master/Galaxy_wide_IMF.txt')
c       fga = metalicity [Fe/H] or [M/H] ?
	  fga=pizold(izold(izo))
c            sfr = 10**log(SFR)
          sfr=10.**ZMU
c            iru = resolution OSD (integer)
! 	  iru=0
	  write(sfr_str,'(f20.7)')sfr
	  write(fga_str,'(f20.7)')fga
	  write(ZMLow_str,'(f20.7)')ZMLow
! 	  write(iru_str,'(i0)')iru
	  call system
     &('python3 ./galIMF-master/igimf_calculator.py '//sfr_str//' '//fga_str//' d d '//ZMLow_str//' d')
      call system('rm ./OUT/IGMF_'//chatom(11:21)//'.txt')
      call system('cp ./Galaxy_wide_IMF.txt ./OUT/IGMF_'//chatom(11:21)//'.txt')
          OPEN(37,FILE='./Galaxy_wide_IMF.txt',
     &STATUS='OLD',IOSTAT=IER5)
          igm=0
          do nni=1,99999999
 	        read(37,'(A)',end=1910,ERR=1911)chau
		read(chau,'(80(A1))')(cv(lnni),lnni=1,80)
		ansmax=cv(1)
		ansmin=cv(1)
		do lnni=2,80
		 if(lgt(cv(lnni),ansmax)) ansmax=cv(lnni)
		 if(llt(cv(lnni),ansmax)) ansmin=cv(lnni)
		enddo
		if(ansmin.eq.ansmax.and.ansmin.eq.' ') goto 1911		 
	        read(chau,*,ERR=1911)xaton1,xaton2
 		igm=igm+1
                yIGMFm(igm)=xaton1
		yIGMFn(igm)=xaton2
1911		continue
          enddo
1910      close(37)
         ZMUS=yIGMFm(igm)
         endif
c     ==================== END IGIMF OPTION ===========================
	   DO iao=1,ntold
	    ageo=ttold(iao)
          call STU_XSL_276A(izold(izo),ageo,BETAo,realeo,ZMISOo,wwffio,
     &wwffiso,a9io,f9mo,qumo,qu10o,qu05o,qu15o)
cc          call STU_XSL(izold(izo),ageo,BETAo,realeo,ZMISOo,wwffio,
cc     &wwffiso,a9io,f9mo,qumo,qu10o,qu05o,qu15o,fxUmo,fxUio,fxUo,
cc     & fxUmso,fxUiso,fxUso,fxUis004o)
	    attao=realeo
cccc	    attao=ageo !lo he puesto para comprobar
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
cc	     call STU_XSL(izyng(izy),agey,BETAy,realey,ZMISOy,wwffiy,
cc     &wwffisy,a9iy,f9my,qumy,qu10y,qu05y,qu15y,fxUmy,fxUiy,fxUy,
cc     & fxUmsy,fxUisy,fxUsy,fxUis004y)
	     call STU_XSL_276A(izyng(izy),agey,BETAy,realey,ZMISOy,wwffiy,
     &wwffisy,a9iy,f9my,qumy,qu10y,qu05y,qu15y)
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
     	       do ini=1,npxi !XSL
     	  	  wwffiy(ini,2)=0.0d0
      	 enddo
     	       f9my=0.0d0 !XSL
     	       do ini=1,3
     	  	  a9iy(ini)=99.0d0
     	       enddo
     	      ENDIF
     	      do nt=1,3 !XSL
	       a90i(nt)=(alpha*BETAo*a9io(nt)+(1.0d0-alpha)*BETAy*
     &a9iy(nt))/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      enddo
c U fluxes
cc        fxUtm=(alpha*BETAo*fxUmo+(1.0d0-alpha)*BETAy*fxUmy)
cc        fxUti=(alpha*BETAo*fxUio+(1.0d0-alpha)*BETAy*fxUiy)
cc        fxUt=(alpha*BETAo*fxUo+(1.0d0-alpha)*BETAy*fxUy)
cc        fxUtms=(alpha*BETAo*fxUmso+(1.0d0-alpha)*BETAy*fxUmsy)
cc        fxUtis=(alpha*BETAo*fxUiso+(1.0d0-alpha)*BETAy*fxUisy)
cc        fxUtis004=(alpha*BETAo*fxUis004o+(1.0d0-alpha)*BETAy*fxUis004y)
cc        fxUts=(alpha*BETAo*fxUso+(1.0d0-alpha)*BETAy*fxUsy)
c           Parametros calidad XSL
     	      quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	      quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
            do io=1,npxi !MIUSC
             wfi(io,2)=alpha*BETAo*wwffio(io,2)+(1.d0-alpha)*
     &BETAy*wwffiy(io,2)
             wfis(io,2)=alpha*BETAo*wwffiso(io,2)+(1.d0-alpha)*
     &BETAy*wwffisy(io,2)
             wfis(io,1)=wfi(io,1)
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
	 write(chatoms(7:7),'(A1)')'x'
	 lchat=1+lnblnk(chatoms)
c	 write(chatoms(lchat:lchat+4),'(A4)')'_sbf'
	 write(chatoms(lchat:lchat+4),'(A4)')'_var'
cccccccc
	      OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	      CLOSE(77,STATUS='DELETE')
8979	      OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	      do nonoo=1,61
	       write(77,'(A80)')cabi(nonoo,4)
	      enddo
	      do lup=1,npxi
		 write(77,*)(wfi(lup,kini),kini=1,2)
	      enddo
	      close(77)
	      OPEN(77,FILE=chatoms,STATUS='OLD',ERR=8989)
	      CLOSE(77,STATUS='DELETE')
8989	      OPEN(77,IOSTAT=IOS,FILE=chatoms,STATUS='NEW')
	      do nonoo=1,61
	       write(77,'(A80)')cabi(nonoo,4)
	      enddo
	      do lup=1,npxi
c		 write(77,*)wfis(lup,1),wfis(lup,2)/wfi(lup,2)
		 write(77,*)wfis(lup,1),wfis(lup,2)
	      enddo
	      close(77)
	      write(*,*)chatom,chatoms	      
	      call e_XSL_276A(alpha,izold(izo),realeo,izyng(izy),realey) !XSL
2325	      continue
	     ENDDO !Cerrado loop Tyng
	    ENDDO !Cerrado loop Zyng
	   ENDDO !Cerrado loop Told
	  ENDDO !Cerrado loop Zold
	 ENDDO !Cerrado loop alpha
	ENDDO !Cerrado loop IMF slope
	RETURN
6336	END




