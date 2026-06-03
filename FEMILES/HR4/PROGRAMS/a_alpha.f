      SUBROUTINE a_alpha
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nstm=1999)
      PARAMETER (npxm=4300)
      DIMENSION starms(nstm,npxm,2) !MILES STARS ARRAY
      DIMENSION starcs(nstm,710,2) !CaT STARS ARRAY
      DIMENSION stars(nstm,2,1107,2) !Jones b&v STARS ARRAY
      DIMENSION aam(nstm,6)
      DIMENSION aamlo(nstm,6)
      DIMENSION aamhi(nstm,6)
      DIMENSION starm(nstm) !miles: added Mg/Fe column
      DIMENSION bb(4),sson(4)
      DIMENSION wfb(1120,2),wwffbo(1120,2),wwffby(1120,2) !jones
      DIMENSION wfv(1120,2),wwffvo(1120,2),wwffvy(1120,2) !jones
      DIMENSION wfc(710,2),wwffco(710,2),wwffcy(710,2) !cat
      DIMENSION wfm(npxm,2),wwffmo(npxm,2),wwffmy(npxm,2) !miles
      DIMENSION wfcs(710,2),wwffcso(710,2),wwffcsy(710,2) !cat
      DIMENSION wfms(npxm,2),wwffmso(npxm,2),wwffmsy(npxm,2) !miles
      DIMENSION a9o(3),a9y(3),a90(3) !jones
      DIMENSION a9co(3),a9cy(3),a90c(3) !cat
      DIMENSION a9mo(5),a9my(5),a90m(5) !miles
c	DIMENSION Q(150,15)
      DIMENSION cab(150,4) !cabeceras de jones (b,v), cat, miles
      DIMENSION ZZZ(8),ZZZo(8),ZZZy(8),SBF(50),ZZsbfo(50),DWGIT(2,8)
      DIMENSION DWGI(2,8),dwgio(2,8),dwgiy(2,8),ZZsbfy(50)
      DIMENSION ZISOLo(8),ZISOLy(8),ZISOL(8),VMo(8),VMy(8),VM(8)
      DIMENSION fZISOLo(8),fZISOLy(8),fZISOL(8)
      DIMENSION fhsto(4),fhsbfo(4),fhsty(4),fhsbfy(4),fhst(4),fhsbf(4)
      DIMENSION co(50)
      DIMENSION vsin(50),vsino(50),vsiny(50),stmo(50),stmy(50)
      DIMENSION z(12,150000,16) !para TERAMO, Z00(15) common/eststu
      DIMENSION ho(5,9),hy(5,9)
      DIMENSION abolo(2),aboly(2)
c      INTEGER npxm,nstm
      CHARACTER*4 gc          !no cambiar que sea de 4, fich sum*
      CHARACTER*1 shape,cv(80),ansmin,ansmax !shape or SFR and galIMF output
      CHARACTER*80 chau
      CHARACTER*10 co
      CHARACTER*80 cab,chatb,chator,chatoc,chatom,chatocs,chatoms,chato
      CHARACTER*80 ast,star,starc,starm
      CHARACTER*9 anom
      CHARACTER*7 cmgfe
      CHARACTER*8 cisoc
      CHARACTER*80 jmiku,ctunn
      CHARACTER*25 fga_str
      CHARACTER*25 sfr_str
      CHARACTER*25 ZMLow_str
      CHARACTER*2 iru_str
      COMMON/labeli/cisoc,cmgfe !aaa,a
      COMMON/aaaa/izold(15),izyng(15),ttold(99),ttyng(99),pizold(15),
     &pizyng(15),zi(150),zalpha(22),nnIMF,nalpha,nzold,ntold,nzyng,
     &ntyng
      COMMON/matina/a(4,4)
      COMMON/fezsol/fez(15)
      COMMON/tefsca/nghb
c	COMMON/ll/Q,nyr
      COMMON/lmu0/bicl,bicp,bich
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
      COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3
      COMMON/lshape/shape
      COMMON/lson/sson
      COMMON/stus/ZMASA,tiso(75),num(15),numtis    
      COMMON/zsolo/z
      COMMON/tgf/t222,g222,f222
      COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/chrbus/ast(650,2)!jones
      COMMON/hrbus/aa(650,6),nstar !jones
      COMMON/lstd/nsss1,nsss2 !jones
      COMMON/clstd/star(650,2) !jones
      COMMON/hrcat/aac(710,4),nstarc !cat
      COMMON/lstdc/nsssc !cat
      COMMON/clstdc/starc(710) !cat
      COMMON/mimgfe/aam,aamlo,aamhi !miles interpolador 4D
      COMMON/hrm1/nstarm !miles
      COMMON/hrm2/starm !miles
      COMMON/lstdm/nsssm !miles
      COMMON/esc/vsin,ZZZ,SBF,DWGI,a90,a90c,ZISOL,fZISOL,ncols,nobs,
     &fhst,fhsbf
      COMMON/TOAGB2/hh(5,9)
      COMMON/TOAGB/g(15,2000,7),ng(15)
      COMMON/TOAGBC/gc(15,2000)
      COMMON/KEYS/keypT(18)
      COMMON/REMNAN/re(15,2,99),nr(15)
      COMMON/CTFIN/CATFIN
      COMMON/NC/ncha !indice tipo isocronas
      COMMON/qua/quaqn,quaq10,quaq05,quaq15
cc      COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,vkvega,
     &zerovk
c	COMMON/FNOR/C2Code,FNORMM,FNORMC
      COMMON/enhanc/amgfe,DeltFeH ![Mg/Fe],[Fe/H] correc. value
      COMMON/enhand/damgfe !overlapping Mg/Fe value for 4D interpolation
      COMMON/enhanw/tamgfe !Largest deviation from [Mg/Fe] value
      COMMON/ab/iaba,iabaj
      COMMON/milsta/starms !MILES STARS ARRAY: common routines: a,sigmam
      COMMON/catsta/starcs !CaT STARS ARRAY: common routines: a,sigmac
      COMMON/jonsta/stars !Jones b&v STARS ARRAY:commons:a,lstdor(en lbusc.f)
      COMMON/FRB/frbol(15),frbolt !fraccion bolom fases evolutivas + total
      COMMON/jmik/jmiku(20),lenj(20) !nombre dirs estrellas
      COMMON/jmik2/ctunn,lctunn
      COMMON/bol/abola(3)
      COMMON/IGIMF/yIGMFm(9999),yIGMFn(9999),igm !also in limf.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      bolsol=4.70d0
      do i6=1,7
	 write(*,*)jmiku(i6)
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Default for scaled-solar element partition:
ccccccccccccccccccccccccccccccccccccccccccccc
c      [Z/H]=[Fe/H]+Gamma*[alpha/Fe] ; [Fe/H]=[Z/H]-Gamma*[alpha/Fe]
c      Dmico=gamma*[alpha/Fe]=gamma*amgfe
c Example for Coelho's library: gamma=0.75, amgfe=0.4 => DeltFeH=0.3
c      [Fe/H]=[Z/H]-0.3 
c       gamma=0.0d0
c       DeltFeH=0.0d0
c
cVER DONDE ESTA IMPLEMENTADO:"Fijamos Mg/Fe en funcion de las temperaturas" 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Introducimos valores criticos para interpolador sigmam_alpha.f:
c intervalo comun de Mg/Fe para *l y *h busqueda Mg/Fe
c               damgfe=Delta[Mg/Fe]
c rango dinamico maximo de Mg/Fe para la busqueda:
c               tamgfe*damgfe
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      damgfe=0.05d0 !ovlap Mg/Fe region
c      damgfe=0.025d0 !ovlap Mg/Fe region
c      damgfe=0.075d0 !ovlap Mg/Fe region
      tamgfe=0.25d0 !i.e. allowed Mg/Fe range: amgfe+-tamgfe
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      tzadi=4500.0d0 !temp mas baja distinguimos[Mg/Fe],dwarfs
      tzadf=8000.0d0 !temp mas alta distinguimos [Mg/Fe],dwarfs
      tzagi=4000.0d0 !temp mas baja distinguimos[Mg/Fe],giants
      tzagf=7000.0d0 !temp mas alta distinguimos [Mg/Fe],giants
c	numero de colores
	ncols=7
c	numero total de colores+indices
      nobs=39
      co(1)='U-V      ='
      co(2)='B-V      ='
      co(3)='V-R      ='
      co(4)='V-I      ='
      co(5)='V-J      ='
      co(6)='V-H      ='
      co(7)='V-K      ='
      co(8)='CN1 (mag)='
      co(9)='CN2 (mag)='
      co(10)='Ca4227   ='
      co(11)='G-band   ='
      co(12)='Fe4383   ='
      co(13)='Ca4455   ='
      co(14)='Fe4531   ='
      co(15)='Fe4668   ='
      co(16)='H-beta   ='
      co(17)='Fe5015   ='
      co(18)='Mg1 (mag)='
      co(19)='Mg2 (mag)='
      co(20)='Mgb      ='
      co(21)='Fe5270   ='
      co(22)='Fe5335   ='
      co(23)='Fe5406   ='
      co(24)='Fe5709   ='
      co(25)='Fe5782   ='
      co(26)='NaD      ='
      co(27)='TiO1(mag)='
      co(28)='TiO2(mag)='
c       co(29)='D4000    ='
      co(29)='HdeltaA  ='
      co(30)='HgammaA  ='
      co(31)='HdeltaF  ='
      co(32)='HgammaF  ='
      co(33)='D4000    ='
      co(34)='CaT*     ='
      co(35)='PaT      ='
      co(36)='CaT      ='
      co(37)='STio     ='
      co(38)='MgI      ='
      co(39)='CO       ='
c	co(29)='CaII(1)  ='
c	co(30)='CaII(2)  ='
c	co(32)='MgI	 ='
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	 write(*,*)'Metallicity scale: [Fe/H](1) or [O/H](0)?'
c	 write(*,*)'IGNORED: ADOPTED [Fe/H] scale'
c	 read(*,*)io01
      io01=1
      if(io01.ne.0)then
       io01=1 !base Fe/H (aam(k,3) en INPUT/PARAM_MILES
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Etiquetas nombre modelo:
      IF(iaba.eq.9)THEN
        if(amgfe.ge.0.0d0)then
	   write(cmgfe(1:7),'(A2,A1,f4.2)')'_M','p',amgfe
        else
	   write(cmgfe(1:7),'(A2,A1,f4.2)')'_M','m',abs(amgfe)
        endif
      ELSE
        if(io01.eq.1)then
	   cmgfe='_baseFe'
        else
	   cmgfe='_baseOH'
        endif
      ENDIF
cccccccccccJones
2333  format(A9,1x,f6.0,1x,f4.2,4(1x,f5.2))
      IF(iaba.eq.0.and.iabaj.eq.1)THEN
       OPEN(98,FILE='../D/INPUT/PARAM_STARS',STATUS='OLD')
       nstar=0
241    format(A9)
       do k=1,99999
         read(98,2333,end=23)anom,(aa(k,l),l=1,6)
c         ast(k,1)='./STARS_4/4'//anom
c         ast(k,2)='./STARS_5/5'//anom
         ast(k,1)=jmiku(1)
	 write(ast(k,1)(lenj(1)+1:lenj(1)+1),'(A1)')'4'
	 write(ast(k,1)(lenj(1)+2:lenj(1)+2+9),241)anom
         ast(k,2)=jmiku(2)
	 write(ast(k,2)(lenj(2)+1:lenj(2)+1),'(A1)')'5'
	 write(ast(k,2)(lenj(2)+2:lenj(2)+2+9),241)anom
c         ast(k,1)=jmiku(1)//anom
c         ast(k,2)=jmiku(2)//anom
         star(k,1)=ast(k,1)
         star(k,2)=ast(k,2)
         nstar=nstar+1
       enddo
23     CLOSE(98)
      ENDIF
cccccccccccJones
2334  FORMAT(A6,1x,f6.0,1x,f4.2,1x,f5.2,1x,f5.3)
c Es MUY IMPORTANTE no cambiar de formato PARAM_CAT!!!!
244   FORMAT(A6)
      OPEN(98,FILE='../D/INPUT/PARAM_CAT',STATUS='OLD')
      nstarc=0
      do k=1,99999
         read(98,2334,end=24)anom,(aac(k,l),l=1,4)
c         starc(k)='./STARS_C/'//anom
         starc(k)=jmiku(4)
	 write(starc(k)(lenj(4)+1:lenj(4)+1+6),244)anom
         nstarc=nstarc+1
      enddo
24    CLOSE(98)
243   FORMAT(A6)
      OPEN(98,FILE='../D/INPUT/PARAM_MILES',STATUS='OLD')
      nstarm=0
      DO k=1,99999
       if(iaba.eq.9.and.io01.eq.0)then !solo en caso baseOH
        read(98,*,end=25)anom,aam(k,1),aam(k,2),aam(k,6),aam(k,4),
     &  aam(k,5),aam(k,3)
       else
        read(98,*,end=25)anom,(aam(k,l9),l9=1,6)
       endif
c      starm(k)='./STARS_MILES/'//anom
       starm(k)=jmiku(3)
       write(starm(k)(lenj(3)+1:lenj(3)+1+6),243)anom
       nstarm=nstarm+1
c Aplicamos escala GonzalezHernandezBonifacio09?
       if(nghb.eq.1)then
        if(aam(k,1).gt.3750.0d0.and.aam(k,1).lt.7500.0d0)then
           aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*aam(k,1))
        elseif(aam(k,1).le.3750.0d0)then
           aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*3750.0d0)
        elseif(aam(k,1).ge.7500.0d0)then
           aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*7500.0d0)
        endif
       endif
c Fijamos [Mg/Fe] en funcion de las temperaturas
       IF(iaba.eq.9)THEN
	if(aam(k,2).gt.3.d0)then
            if(aam(k,1).lt.tzadi.or.aam(k,1).gt.tzadf)then
              aam(k,5)=amgfe
	    endif
	else
            if(aam(k,1).lt.tzagi.or.aam(k,1).gt.tzagf)then
              aam(k,5)=amgfe
	    endif
	endif
       ENDIF         
      ENDDO
25    CLOSE(98)
      nsss1=nstar !jones
      nsss2=nstar !jones
      nsssc=nstarc
      nsssm=nstarm
c Construimos las dos listas de espectros MILES para sigmam_alpha 3.5D
      DO k00=1,nstarm
	  do k09=1,6
	   aamlo(k00,k09)=aam(k00,k09)
	   aamhi(k00,k09)=aam(k00,k09)
	  enddo
	  if(aam(k00,5).lt.(amgfe+damgfe))then
	   aamlo(k00,5)=aam(k00,5)
	  else
	   aamlo(k00,5)=99.d0
	  endif
	  if(aam(k00,5).gt.(amgfe-damgfe))then
	   aamhi(k00,5)=aam(k00,5)
	  else
	   aamhi(k00,5)=99.d0
	  endif
      ENDDO
c LECTURA DE FICHEROS ASCII DE LIBRERIAS ESPECTRALES
c SE CARGA EN MATRIZ:starms(nstarm,npxm,2)
      DO ios=1,nstarm !MILES
	 open(99,file=starm(ios),status='old')
	 do iop=1,npxm
	   read(99,*)starms(ios,iop,1),starms(ios,iop,2)
	   if(ios.eq.1)wfm(iop,1)=starms(ios,iop,1)
	 enddo
	 close(99)
      ENDDO
      DO ios=1,nstarc !CaT
	 open(99,file=starc(ios),status='old')
	 do iop=1,710
	   read(99,*)starcs(ios,iop,1),starcs(ios,iop,2)
	   if(ios.eq.1)wfc(iop,1)=starcs(ios,iop,1)
	 enddo
	 close(99)
      ENDDO
      IF(iaba.eq.0.and.iabaj.eq.1)THEN
	DO ios=1,nstar !4 Jones
	  open(99,file=star(ios,1),status='old')
	  do iop=1,1107
	   read(99,*)stars(ios,1,iop,1),stars(ios,1,iop,2)
	   if(ios.eq.1)wfb(iop,1)=stars(ios,1,iop,1)
	  enddo
	  close(99)
	ENDDO
	DO ios=1,nstar !5 Jones
	  open(99,file=star(ios,2),status='old')
	  do iop=1,1107
	   read(99,*)stars(ios,2,iop,1),stars(ios,2,iop,2)
	   if(ios.eq.1)wfv(iop,1)=stars(ios,2,iop,1)
	  enddo
	  close(99)
	ENDDO
	open(99,file='../D/INPUT/GALAXY_CABECERAS_BLUE',status='old')
	do k=1,61
		read(99,'(A80)')cab(k,1)
	enddo
	close(99)
	open(99,file='../D/INPUT/GALAXY_CABECERAS_RED',status='old')
	do k=1,61
		read(99,'(A80)')cab(k,2)
	enddo
	close(99)
      ENDIF
cccccccccccMILES+CAT
      open(99,file='../D/INPUT/GALAXY_CABECERAS_CAT',status='old')
      do k=1,61
	read(99,'(A80)')cab(k,3)
      enddo
      close(99)
      open(99,file='../D/INPUT/GALAXY_CABECERAS_MILES',status='old')
      do k=1,61
	read(99,'(A80)')cab(k,4)
      enddo
      close(99)
c*******COMIENZO DEL BUCLE DE MU o SFR***************************************
C Abierto loop IMF slope
      DO NNN=1,nnIMF
	 if(shape.eq.'g')then
	     ZMU=zi(NNN)
	 else
	     ZMU=zi(NNN)+1.0d0
	 endif
       IF(shape.eq.'b')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bbi'//chato(1:4)
	  chator='./OUT/R/Rbi'//chato(1:4)
	 endif
         chatoc='./OUT/C/Cbi'//chato(1:4)
         chatom='./OUT/M/Mbi'//chato(1:4)
       ELSEIF(shape.eq.'u')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bun'//chato(1:4)
          chator='./OUT/R/Run'//chato(1:4)
	 endif
         chatoc='./OUT/C/Cun'//chato(1:4)
         chatom='./OUT/M/Mun'//chato(1:4)
       ELSEIF(shape.eq.'g')THEN
         if(ZMU.ge.0.d0)then
          write(chato(1:1),'(A1)')'p'
 	 else
          write(chato(1:1),'(A1)')'m'
	 endif
         write(chato(2:5),'(f4.2)')abs(ZMU)
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bg'//chato(1:5)
          chator='./OUT/R/Rg'//chato(1:5)
	 endif
         chatoc='./OUT/C/Cg'//chato(1:5)
         chatom='./OUT/M/Mg'//chato(1:5)
       ELSEIF(shape.eq.'t')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bba'//chato(1:4)
          chator='./OUT/R/Rba'//chato(1:4)
	 endif
         chatoc='./OUT/C/Cba'//chato(1:4)
         chatom='./OUT/M/Mba'//chato(1:4)
       ELSEIF(shape.eq.'k')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bku'//chato(1:4)
          chator='./OUT/R/Rku'//chato(1:4)
	 endif
         chatoc='./OUT/C/Cku'//chato(1:4)
         chatom='./OUT/M/Mku'//chato(1:4)
       ELSEIF(shape.eq.'r')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bkb'//chato(1:4)
	  chator='./OUT/R/Rkb'//chato(1:4)
	 endif
         chatoc='./OUT/C/Ckb'//chato(1:4)
         chatom='./OUT/M/Mkb'//chato(1:4)
       ELSEIF(shape.eq.'c')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/Bch'//chato(1:4)
	  chator='./OUT/R/Rch'//chato(1:4)
	 endif
         chatoc='./OUT/C/Cch'//chato(1:4)
         chatom='./OUT/M/Mch'//chato(1:4)
       ELSEIF(shape.eq.'f')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/BFe'//chato(1:4)
	  chator='./OUT/R/RFe'//chato(1:4)
	 endif
         chatoc='./OUT/C/CFe'//chato(1:4)
         chatom='./OUT/M/MFe'//chato(1:4)
       ELSEIF(shape.eq.'x')THEN
         write(chato(1:4),'(f4.2)')ZMU-1.0d0
         if(iaba.eq.0.and.iabaj.eq.1)then
	  chatb='./OUT/B/BFx'//chato(1:4)
	  chator='./OUT/R/RFx'//chato(1:4)
	 endif
         chatoc='./OUT/C/CFx'//chato(1:4)
         chatom='./OUT/M/MFx'//chato(1:4)
       ELSEIF(shape.eq.'l')THEN
         int1=int(abs(ZMU1*10))
         int2=int(abs(ZMU2*10))
         if(ZMU1.ge.0.0d0)then
           if(int1.lt.10)then
              write(chato(1:1),'(I1)')int1
              if(iaba.eq.0.and.iabaj.eq.1)then
	       chatb='./OUT/B/Bp0'//chato(1:1)
	       chator='./OUT/R/Rp0'//chato(1:1)
	      endif
              chatoc='./OUT/C/Cp0'//chato(1:1)
              chatom='./OUT/M/Mp0'//chato(1:1)
           else
              write(chato(1:2),'(I2)')int1
              if(iaba.eq.0.and.iabaj.eq.1)then
	       chatb='./OUT/B/Bp'//chato(1:2)
	       chator='./OUT/R/Rp'//chato(1:2)
	      endif
              chatoc='./OUT/C/Cp'//chato(1:2)
              chatom='./OUT/M/Mp'//chato(1:2)
           endif
         else
           if(int1.lt.10)then
              write(chato(1:1),'(I1)')int1
              if(iaba.eq.0.and.iabaj.eq.1)then
	       chatb='./OUT/B/Bm0'//chato(1:1)
	       chator='./OUT/R/Rm0'//chato(1:1)
	      endif
              chatoc='./OUT/C/Cm0'//chato(1:1)
              chatom='./OUT/M/Mm0'//chato(1:1)
           else
              write(chato(1:2),'(I2)')int1
              if(iaba.eq.0.and.iabaj.eq.1)then
	       chatb='./OUT/B/Bm'//chato(1:2)
	       chator='./OUT/R/Rm'//chato(1:2)
	      endif
              chatoc='./OUT/C/Cm'//chato(1:2)
              chatom='./OUT/M/Mm'//chato(1:2)
           endif
         endif
         if(ZMU2.ge.0.0d0)then
           if(int2.lt.10)then
              write(chato(1:1),'(I1)')int2
	      if(iaba.eq.0.and.iabaj.eq.1)then 
               write(chatb(13:15),'(A3)')'p0'//chato(1:1)
               write(chator(13:15),'(A3)')'p0'//chato(1:1)
	      endif
              write(chatoc(13:15),'(A3)')'p0'//chato(1:1)
              write(chatom(13:15),'(A3)')'p0'//chato(1:1)
           else
              write(chato(1:2),'(I2)')int2
	      if(iaba.eq.0.and.iabaj.eq.1)then 
               write(chatb(13:15),'(A3)')'p'//chato(1:2)
               write(chator(13:15),'(A3)')'p'//chato(1:2)
	      endif
              write(chatoc(13:15),'(A3)')'p'//chato(1:2)
              write(chatom(13:15),'(A3)')'p'//chato(1:2)
           endif
         else
           if(int2.lt.10)then
              write(chato(1:1),'(I1)')int2
 	      if(iaba.eq.0.and.iabaj.eq.1)then 
               write(chatb(13:15),'(A3)')'m0'//chato(1:1)
               write(chator(13:15),'(A3)')'m0'//chato(1:1)
	      endif
              write(chatoc(13:15),'(A3)')'m0'//chato(1:1)
              write(chatom(13:15),'(A3)')'m0'//chato(1:1)
           else
              write(chato(1:2),'(I2)')int2
	      if(iaba.eq.0.and.iabaj.eq.1)then 
               write(chatb(13:15),'(A3)')'m'//chato(1:2)
               write(chator(13:15),'(A3)')'m'//chato(1:2)
	      endif
              write(chatoc(13:15),'(A3)')'m'//chato(1:2)
              write(chatom(13:15),'(A3)')'m'//chato(1:2)
           endif
         endif
       ENDIF
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
c*******COMIENZO DEL BUCLE DEL PORCENTAGE*******************************
C Abierto loop alpha
       DO ialpha=1,nalpha      
	alpha=zalpha(ialpha)
        if(alpha.eq.1.0d0)then
          ntyng=1
          ttyng(1)=ttold(1)
          nzyng=1
          izyng(1)=izold(1)
        endif
c*******COMIENZO DEL BUCLE DE LA SSPo***********************************
C Abierto loop Zold
        DO izo=1,nzold
	 if(pizold(izold(izo)).ge.0.0d0)then
	  write(chato(1:6),'(A,A,f4.2)')'Z','p',pizold(izold(izo))
	 else
	  write(chato(1:6),'(A,A,f4.2)')'Z','m',abs(pizold(izold(izo)))
	 endif
         fehf=pizold(izold(izo)) ![Fe/H] que pasara a STU
  	 write(chatb(16:21),'(A6)')chato(1:6)
  	 write(chator(16:21),'(A6)')chato(1:6)
  	 write(chatoc(16:21),'(A6)')chato(1:6)
  	 write(chatom(16:21),'(A6)')chato(1:6)
c     ==================== IGIMF OPTION ===========================
         if(shape.eq.'g')then
c            fga = metalicity [Fe/H] or [M/H] ?
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
C Abierto loop Told
	 DO iao=1,ntold
	  ageo=ttold(iao)
c      write(57,*)' '
c      write(57,'(A21)')chatom
c      write(57,'(A2,1x,(F5.2,1x),2(F7.4,1x))')
c     &shape,ZMU-1.,fehf,ageo*0.001d0
c Fracciones bolometricas:
	  frbolt=0.0d0
          do iiii=1,15
            frbol(iiii)=0.0d0
          enddo
cc 	  call STU_alpha(izold(izo),ageo,BETAo,realeo,ZZZo,ZMISOo,
cc     &vsino,stmo,wwffbo,wwffvo,wwffco,wwffmo,wwffcso,wwffmso,a9o,
cc     &f9o,a9co,f9co,a9mo,f9mo,ZISOLo,fZISOLo,VMo,dwgio,ZZsbfo,fhsto,
cc     &fhsbfo,t222o,g222o,f222o,Wfcatto,ho,fluxxo,qumo,qu10o,qu05o,
cc     &qu15o,fxUmo,fxUio,fxUo,fxUmso,fxUiso,fxUso,fxUis004o,abolo)
 	  call STU_alpha(izold(izo),ageo,BETAo,realeo,ZZZo,ZMISOo,
     &vsino,stmo,wwffbo,wwffvo,wwffco,wwffmo,wwffcso,wwffmso,a9o,
     &f9o,a9co,f9co,a9mo,f9mo,ZISOLo,fZISOLo,VMo,dwgio,ZZsbfo,fhsto,
     &fhsbfo,t222o,g222o,f222o,Wfcatto,ho,fluxxo,qumo,qu10o,qu05o,
     &qu15o,abolo)
	  attao=realeo
	  if(attao.lt.10.0d0)then
	   write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	  else
	   write(chato(1:8),'(A,f7.4)')'T',attao
	  endif
	  if(iaba.eq.0.and.iabaj.eq.1)then 
	   write(chatb(22:29),'(A8)')chato(1:8)
  	   write(chator(22:29),'(A8)')chato(1:8)
	  endif
  	  write(chatoc(22:29),'(A8)')chato(1:8)
  	  write(chatom(22:29),'(A8)')chato(1:8)
c     name isoc+alpha
	  if(iaba.eq.0.and.iabaj.eq.1)then 
  	   write(chatb(30:37),'(A8)')cisoc
  	   write(chatb(38:44),'(A7)')cmgfe
  	   write(chator(30:37),'(A8)')cisoc
  	   write(chator(38:44),'(A7)')cmgfe
	   chatb=chatb(1:44)
	   chator=chator(1:44)
	  endif
  	  write(chatoc(30:37),'(A8)')cisoc
  	  write(chatoc(38:44),'(A7)')cmgfe
  	  write(chatom(30:37),'(A8)')cisoc
  	  write(chatom(38:44),'(A7)')cmgfe
	  chatoc=chatoc(1:44)
	  chatom=chatom(1:44)
c*******COMIENZO DEL BUCLE DE LA SSPy***********************************
C Abierto loop Zyng
	  DO izy=1,nzyng !Abierto loop Zyng
	   IF(alpha.ne.1.0d0)THEN
	    if(iaba.eq.0.and.iabaj.eq.1)then 
  	     write(chatb(45:49),'(A,f4.2)')'a',alpha
  	     write(chator(45:49),'(A,f4.2)')'a',alpha
	    endif
  	    write(chatoc(45:49),'(A,f4.2)')'a',alpha
  	    write(chatom(45:49),'(A,f4.2)')'a',alpha
	    if(pizyng(izyng(izy)).ge.0.0d0)then
	     write(chato(1:6),'(A,A,f4.2)')'Z','p',pizyng(izyng(izy))
	    else
	     write(chato(1:6),'(A,A,f4.2)')'Z','m',
     &abs(pizyng(izyng(izy)))
	    endif
	    if(iaba.eq.0.and.iabaj.eq.1)then 
  	     write(chatb(50:55),'(A6)')chato(1:6)
  	     write(chator(50:55),'(A6)')chato(1:6)
	    endif
  	    write(chatoc(50:55),'(A6)')chato(1:6)
  	    write(chatom(50:55),'(A6)')chato(1:6)
	   ENDIF
C Abierto loop Tyng
           DO iay=1,ntyng !Abierto loop Tyng
	    agey=ttyng(iay)
	    if(agey.gt.ageo.and.alpha.ne.1.00d0)goto 2325
	    IF(alpha.ne.1.0d0)THEN
cc	     call STU_alpha(izyng(izy),agey,BETAy,realey,ZZZy,ZMISOy,
cc     &vsiny,stmy,wwffby,wwffvy,wwffcy,wwffmy,wwffcsy,wwffmsy,a9y,f9y,
cc     &a9cy,f9cy,a9my,f9my,ZISOLy,fZISOLy,VMy,dwgiy,ZZsbfy,fhsty,fhsbfy,
cc     &t222y,g222y,f222y,Wfcatty,hy,fluxxy,qumy,qu10y,qu05y,qu15y,fxUmy,
cc     &fxUiy,fxUy,fxUmsy,fxUisy,fxUsy,fxUis004y,aboly)
	     call STU_alpha(izyng(izy),agey,BETAy,realey,ZZZy,ZMISOy,
     &vsiny,stmy,wwffby,wwffvy,wwffcy,wwffmy,wwffcsy,wwffmsy,a9y,f9y,
     &a9cy,f9cy,a9my,f9my,ZISOLy,fZISOLy,VMy,dwgiy,ZZsbfy,fhsty,fhsbfy,
     &t222y,g222y,f222y,Wfcatty,hy,fluxxy,qumy,qu10y,qu05y,qu15y,aboly)
             attao=realey
	     if(attao.lt.10.0d0)then
	      write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	     else
	      write(chato(1:8),'(A,f7.4)')'T',attao
	     endif
	     if(iaba.eq.0.and.iabaj.eq.1)then 
              write(chatb(56:63),'(A8)')chato(1:8)
  	      write(chator(56:63),'(A8)')chato(1:8)
	      chatb=chatb(1:63)
	      chator=chator(1:63)
	     endif
  	     write(chatoc(56:63),'(A8)')chato(1:8)
  	     write(chatom(56:63),'(A8)')chato(1:8)
	     chatoc=chatoc(1:63)
	     chatom=chatom(1:63)
     	    ELSE
	     if(iaba.eq.0.and.iabaj.eq.1)then 
	      chatb=chatb(1:44)
	      chator=chator(1:44)
	     endif
	     chatoc=chatoc(1:44)
	     chatom=chatom(1:44)
     	     BETAy=0.0d0
     	     realey=0.0d0
     	     do ini=1,8
     	  	ZZZy(ini)=0.
     	     enddo
     	     ZMISOy=0.0d0
     	     do ini=1,50
     	  	vsiny(ini)=0.0d0
     	  	stmy(ini)=0.0d0
     	     enddo
     	     do ini=1,1107 !jones
     	  	wwffby(ini,2)=0.0d0
     	  	wwffvy(ini,2)=0.0d0
     	     enddo
     	     do ini=1,710 !cat
     	  	wwffcy(ini,2)=0.0d0
     	  	wwffcsy(ini,2)=0.0d0
      	     enddo
     	     do ini=1,npxm !miles
     	  	wwffmy(ini,2)=0.0d0
      	  	wwffmsy(ini,2)=0.0d0
     	     enddo
     	     do ini=1,3 !jones
     	  	a9y(ini)=99.0d0
     	     enddo
     	     f9y=0.0d0 !cat
     	     do ini=1,3
     	  	a9cy(ini)=99.0d0
     	     enddo
     	     f9my=0.0d0 !miles
     	     do ini=1,5
     	  	a9my(ini)=99.0d0
     	     enddo
     	    ENDIF
c cat
     	    t222=(alpha*BETAo*t222o+(1.0d0-alpha)*BETAy*
     &t222y)/(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
     	    g222=(alpha*BETAo*g222o+(1.0d0-alpha)*BETAy*
     &g222y)/(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
     	    f222=(alpha*BETAo*f222o+(1.0d0-alpha)*BETAy*
     &f222y)/(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
c miles
c     	t222m=(alpha*BETAo*t222mo+(1.0d0-alpha)*BETAy*t222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c     	g222m=(alpha*BETAo*g222mo+(1.0d0-alpha)*BETAy*g222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c     	f222m=(alpha*BETAo*f222mo+(1.0d0-alpha)*BETAy*f222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
cccccccccccJones
	    IF(iaba.eq.0.and.iabaj.eq.1)THEN 
             do nt=1,3 !jones
	      a90(nt)=(alpha*BETAo*a9o(nt)+(1.0d0-alpha)*BETAy*
     &a9y(nt))/(alpha*BETAo*f9o+(1.0d0-alpha)*BETAy*f9y)
     	     enddo
     	    ENDIF
cccccccccccMILES+CAT
     	    do nt=1,3 !cat
	     a90c(nt)=(alpha*BETAo*a9co(nt)+(1.0d0-alpha)*BETAy*
     &a9cy(nt))/(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
     	    enddo
     	    do nt=1,5 !miles
	     a90m(nt)=(alpha*BETAo*a9mo(nt)+(1.0d0-alpha)*BETAy*
     &a9my(nt))/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	    enddo
c      write(38,'(A60,1x,F7.0,1x,4(F7.3,1x))')chatom,(a90m(nt),nt=1,5)
     	    do nt=1,4
	     fhst(nt)=(-2.50d0)*dlog10(alpha*BETAo*fhsto(nt)+
     &(1.-alpha)*BETAy*fhsty(nt))
	     fhsbf(nt)=(-2.50d0)*dlog10(
     &(alpha*BETAo*fhsbfo(nt)+(1.-alpha)*BETAy*fhsbfy(nt))/
     &(alpha*BETAo*fhsto(nt)+(1.-alpha)*BETAy*fhsty(nt)))
            enddo
c U fluxes
cc            fxUtm=(alpha*BETAo*fxUmo+(1.0d0-alpha)*BETAy*fxUmy)
cc            fxUti=(alpha*BETAo*fxUio+(1.0d0-alpha)*BETAy*fxUiy)
cc            fxUt=(alpha*BETAo*fxUo+(1.0d0-alpha)*BETAy*fxUy)
cc            fxUtms=(alpha*BETAo*fxUmso+(1.0d0-alpha)*BETAy*fxUmsy)
cc            fxUtis=(alpha*BETAo*fxUiso+(1.0d0-alpha)*BETAy*fxUisy)
cc            fxUtis004=(alpha*BETAo*fxUis004o+(1.0d0-alpha)*BETAy*
cc     &fxUis004y)
cc            fxUts=(alpha*BETAo*fxUso+(1.0d0-alpha)*BETAy*fxUsy)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(64,'(A44,16(1X,F7.4))')chatom,
c     &(100.0d0*(frbol(nt)/frbolt),nt=1,15),
c     &bolsol-2.5d0*dlog10(frbolt*BETAo)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	    abola(2)=bolsol-2.50d0*dlog10(alpha*BETAo*abolo(1)+
     &(1.0d0-alpha)*BETAy*aboly(1)) !M_bol
            abola(1)=BETAo*abolo(1) !solo se calcula para poblacione vieja
	    abola(3)=abolo(2) !solo se calcula para poblacione vieja
	    DO nt=1,8
	     ZZZ(nt)=(-2.50d0)*dlog10(alpha*BETAo*ZZZo(nt)+
     &(1.0d0-alpha)*BETAy*ZZZy(nt))
	     SBF(nt)=(-2.50d0)*dlog10(
     &(alpha*BETAo*ZZsbfo(nt)+(1.0d0-alpha)*BETAy*ZZsbfy(nt))/
     &(alpha*BETAo*ZZZo(nt)+(1.0d0-alpha)*BETAy*ZZZy(nt)))
             do ndo=1,5
              hh(ndo,nt)=100.0d0*(alpha*BETAo*ho(ndo,nt)+(1.0d0-alpha)
     &*BETAy*hy(ndo,nt))/
     &(alpha*BETAo*ZZZo(nt)+(1.0d0-alpha)*BETAy*ZZZy(nt))
             enddo
	     do ndo=1,2
	      DWGIT(ndo,nt)=alpha*BETAo*dwgio(ndo,nt)+
     &(1.0d0-alpha)*BETAy*dwgiy(ndo,nt)
             enddo
             DWGI(1,nt)=DWGIT(1,nt)*100./(alpha*BETAo*ZZZo(nt)+
     &(1.0d0-alpha)*BETAy*ZZZy(nt))
             DWGI(2,nt)=DWGIT(2,nt)*100.0d0/(alpha*BETAo*ZZZo(nt)+
     &(1.0d0-alpha)*BETAy*ZZZy(nt))
	    ENDDO
            DO nt=1,5
             hh(nt,9)=100.0d0*(alpha*BETAo*ho(nt,9)+(1.0d0-alpha)
     &*BETAy*hy(nt,9))/(alpha*BETAo*fluxxo+(1.0d0-alpha)*BETAy*fluxxy)
            ENDDO
	    DO nt=8,nobs
	     vsin(nt)=(alpha*BETAo*vsino(nt)+(1.0d0-alpha)*BETAy
     &*vsiny(nt))/(alpha*BETAo*stmo(nt)+(1.0d0-alpha)*BETAy*stmy(nt))
	    ENDDO
     	    WFFCAT=(alpha*BETAo*Wfcatto+(1.0d0-alpha)*BETAy*Wfcatty)
     &/(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
c   Parametros calidad MILES
     	    quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	    quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	    quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	    quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &/(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c        write(*,*)'a_alpha,WFFCAT',WFFCAT
            vsin(8)=(-2.50d0)*dlog10(1.0d0-(vsin(8)/35.0d0))
            vsin(9)=(-2.50d0)*dlog10(1.0d0-(vsin(9)/35.0d0))
            vsin(18)=(-2.50d0)*dlog10(1.0d0-(vsin(18)/65.0d0))
            vsin(19)=(-2.50d0)*dlog10(1.0d0-(vsin(19)/42.50d0))
            vsin(27)=(-2.50d0)*dlog10(1.0d0-(vsin(27)/57.50d0))
            vsin(28)=(-2.50d0)*dlog10(1.0d0-(vsin(28)/82.50d0))
cccccccccccMILES+CAT
            do io=1,710 !cat
             wfc(io,2)=alpha*BETAo*wwffco(io,2)+(1.d0-alpha)*BETAy*
     &wwffcy(io,2)
             wfcs(io,2)=alpha*BETAo*wwffcso(io,2)+(1.d0-alpha)*BETAy*
     &wwffcsy(io,2)
              wfcs(io,1)=wfc(io,1)
            enddo
            do io=1,npxm !miles
             wfm(io,2)=alpha*BETAo*wwffmo(io,2)+(1.d0-alpha)*BETAy*
     &wwffmy(io,2)
             wfms(io,2)=alpha*BETAo*wwffmso(io,2)+(1.d0-alpha)*BETAy*
     &wwffmsy(io,2)
             wfms(io,1)=wfm(io,1)
            enddo
cccccccccccJones
	    IF(iaba.eq.0.and.iabaj.eq.1)THEN 
             do io=1,1107 !jones
              wfb(io,2)=alpha*BETAo*wwffbo(io,2)+(1.d0-alpha)*BETAy*
     &wwffby(io,2)		
              wfv(io,2)=alpha*BETAo*wwffvo(io,2)+(1.d0-alpha)*BETAy*
     &wwffvy(io,2)
             enddo
	     OPEN(77,FILE=chatb,STATUS='OLD',ERR=8976)
	     CLOSE(77,STATUS='DELETE')
8976	     OPEN(77,IOSTAT=IOS,FILE=chatb,STATUS='NEW')
	     do nonoo=1,61
	      write(77,'(A80)')cab(nonoo,1)
	     enddo
	     do nonoo=1,1107
	      if(wfb(nonoo,1).gt.3855.4d0.and.wfb(nonoo,1).lt.4476.5d0)
     &then
		write(77,*)(wfb(nonoo,kini),kini=1,2)
	      endif
	     enddo
	     close(77)
	     write(*,*)chatb
c jones v
             OPEN(77,FILE=chator,STATUS='OLD',ERR=8977)
             CLOSE(77,STATUS='DELETE')
8977	     OPEN(77,IOSTAT=IOS,FILE=chator,STATUS='NEW')
             do nonoo=1,61
              write(77,'(A80)')cab(nonoo,2)
             enddo
             do nonoo=1,1107
	      if(wfv(nonoo,1).gt.4794.9d0.and.wfv(nonoo,1).lt.5465.1d0)
     &then
		write(77,*)(wfv(nonoo,kini),kini=1,2)
	      endif
	     enddo
	     close(77)
	     write(*,*)chator
	    ENDIF
c_NaFep06_aFep04_MgFeCORR	 
	    lchat=1+lnblnk(chatoc)
	    if(iaba.eq.1)then
	     write(chatoc(lchat:lchat+lctunn),'(A15)')ctunn
	    elseif(iaba.eq.2)then
	     write(chatoc(lchat:lchat+lctunn),'(A24)')ctunn
	    elseif(iaba.eq.3)then
	     write(chatoc(lchat:lchat+lctunn),'(A15)')ctunn
	    endif
	    lchat=1+lnblnk(chatom)
	    if(iaba.eq.1)then
	     write(chatom(lchat:lchat+lctunn),'(A15)')ctunn
	    elseif(iaba.eq.2)then
	     write(chatom(lchat:lchat+lctunn),'(A24)')ctunn
 	    elseif(iaba.eq.3)then
	     write(chatom(lchat:lchat+lctunn),'(A15)')ctunn
           endif
c nombre espectros sbf:
            chatocs=chatoc	    
	    write(chatocs(7:7),'(A1)')'c'
	    lchat=1+lnblnk(chatocs)
c	    write(chatocs(lchat:lchat+4),'(A4)')'_sbf'
	    write(chatocs(lchat:lchat+4),'(A4)')'_var'
c            chatocs='sbf_'//chatoc	   
            chatoms=chatom	    
	    write(chatoms(7:7),'(A1)')'m'
	    lchat=1+lnblnk(chatoms)
c	    write(chatoms(lchat:lchat+4),'(A4)')'_sbf'
	    write(chatoms(lchat:lchat+4),'(A4)')'_var'
c cat
	    OPEN(77,FILE=chatoc,STATUS='OLD',ERR=8978)
	    CLOSE(77,STATUS='DELETE')
8978	    OPEN(77,IOSTAT=IOS,FILE=chatoc,STATUS='NEW')
	    do nonoo=1,61
	     write(77,'(A80)')cab(nonoo,3)
	    enddo
	    do lup=1,710
             write(77,*)(wfc(lup,kini),kini=1,2)
	    enddo
	    close(77)
	    OPEN(77,FILE=chatocs,STATUS='OLD',ERR=8981)
	    CLOSE(77,STATUS='DELETE')
8981	    OPEN(77,IOSTAT=IOS,FILE=chatocs,STATUS='NEW')
	    do nonoo=1,61
	     write(77,'(A80)')cab(nonoo,3)
	    enddo
	    do lup=1,710
c             write(77,*)wfcs(lup,1),wfcs(lup,2)/wfc(lup,2)
             write(77,*)wfcs(lup,1),wfcs(lup,2)
	    enddo
	    close(77)
	    write(*,*)chatoc,chatocs
c miles
	    OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	    CLOSE(77,STATUS='DELETE')
8979	    OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	    do nonoo=1,61
	     write(77,'(A80)')cab(nonoo,4)
	    enddo
	    do lup=1,npxm
             write(77,*)(wfm(lup,kini),kini=1,2)
	    enddo
	    close(77)
	    OPEN(77,FILE=chatoms,STATUS='OLD',ERR=8980)
	    CLOSE(77,STATUS='DELETE')
8980	    OPEN(77,IOSTAT=IOS,FILE=chatoms,STATUS='NEW')
	    do nonoo=1,61
	     write(77,'(A80)')cab(nonoo,4)
	    enddo
	    do lup=1,npxm
c             write(77,*)wfms(lup,1),wfms(lup,2)/wfm(lup,2)
             write(77,*)wfms(lup,1),wfms(lup,2)
	    enddo
	    close(77)
	    write(*,*)chatom,chatoms
c    Medida del indice en el espectro sintetico del CaT
	    call NICO(wfc,CATFIN)
c-----------------------------------------------------------------------
c   ATENCION: solo se calcula la M/L para la pobl. vieja y no para 2pob.
	    do itton=1,8
              VM(itton)=VMo(itton)
              ZISOL(itton)=ZISOLo(itton)
              fZISOL(itton)=fZISOLo(itton)
	    enddo
	    if(alpha.eq.1.)then
	     do izsoly=1,8
	      ZISOLy(izsoly)=0.0d0
	     enddo
	    endif
c-----------------------------------------------------------------------
	    call e_alpha(alpha,izold(izo),realeo,izyng(izy),realey,
     &a90m)
2325	    continue
           ENDDO !Cerrado loop Tyng
          ENDDO !Cerrado loop Zyng
         ENDDO !Cerrado loop Told
        ENDDO !Cerrado loop Zold
       ENDDO !Cerrado loop alpha
      ENDDO !Cerrado loop IMF slope
      RETURN
      END




