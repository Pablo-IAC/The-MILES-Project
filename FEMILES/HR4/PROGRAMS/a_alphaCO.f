       	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (npxm=4300) !miles
	PARAMETER (nstm=1999) !coelho
	DIMENSION a(4,4),bb(4),sson(4)
	DIMENSION wfm(npxm,2),wwffmo(npxm,2),wwffmy(npxm,2) !miles
	DIMENSION a9mo(3),a9my(3),a90m(3) !miles
	DIMENSION izold(15),izyng(15),ttold(99),ttyng(99)
	DIMENSION pizold(15),pizyng(15),yter(12),ygir(7)
	DIMENSION cab(150) !cabecera de miles
c	DIMENSION ZZZ(8),ZZZo(8),ZZZy(8),DWGIT(2,8),DWGI(2,8),
c     &dwgio(2,8),dwgiy(2,8),tiso(75)
	DIMENSION ZZZ(8),ZZZo(8),ZZZy(8),SBF(50),ZZsbfo(50),
     &  ZZsbfy(50),DWGIT(2,8),DWGI(2,8),dwgio(2,8),dwgiy(2,8),tiso(75)
        DIMENSION ZISOLo(8),ZISOLy(8),ZISOL(8),VMo(8),VMy(8),VM(8)
        DIMENSION fZISOLo(8),fZISOLy(8),fZISOL(8)
        DIMENSION fhsto(4),fhsbfo(4),fhsty(4),fhsbfy(4),fhst(4),fhsbf(4)
	DIMENSION ze(75),zi(150),num(15)
	DIMENSION zalpha(22),stmo(50),stmy(50)
	DIMENSION z(12,150000,16),Z00s(15) !para TERAMO, Z00(15) common/eststu
	DIMENSION ho(5,9),hy(5,9)
        DIMENSION aam(nstm,5),starm(nstm) !miles
        DIMENSION starms(nstm,npxm,2) !COELHO STELLAR SPECTRA ARRAY
	character*8 dirout
	character*43 mkout
	integer staout,mkout2,system
        CHARACTER*4 gc          !no cambiar que sea de 4, fich sum*
	CHARACTER*1 shape,lowcut
	CHARACTER*80 cab,chatom,chato
	CHARACTER*20 ast,starm
	CHARACTER*9 anom
      COMMON/milsta/starms !COELHO STELLAR SPECTRA ARRAY: common routines: a,sigmam        
	COMMON/mic/mico
	COMMON/fezsol/fez(15)
      COMMON/lshape/shape
	COMMON/lson/sson
	COMMON/stus/ZMASA,tiso,num,numtis    
	COMMON/zsolo/z
	COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/rat/ra(3,50)
      COMMON/hrm1/aam,nstarm !miles
      COMMON/hrm2/starm !miles
	COMMON/lstdm/nsssm !miles
	common/zalpCO/ZZZ,ZISOL,fZISOL,SBF,fhst,fhsbf !common a.f, e.f
      COMMON/esc/vsin,DWGI,a90,a90c,ncols,nobs
      COMMON/TOAGB2/hh(5,9)
      COMMON/TOAGB/g(15,2000,7),ng(15)
      COMMON/TOAGBC/gc(15,2000)
	COMMON/KEYS/keypT(18)
	COMMON/REMNAN/re(15,2,99),nr(15)
      COMMON/NC/ncha !indice tipo isocronas
	COMMON/qua/quaqn,quaq10,quaq05,quaq15
	COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !common a.f
c	COMMON/lmu/ZMU,bicl,bicp,bich
	COMMON/lmu0/bicl,bicp,bich !!!!!!!!!!!!!
	COMMON/lmu1/ZMU !!!!!!!!!!!!!
	COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3 !!!!!!!!!!!!!
	COMMON/turnp/xtpoint !!!!!!!!!!!!!
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk !!!!!!!!!!!!!
      C2Code=0.958d0 !calculada para BC=-0.07, pero no importa
c para calcular valores [Fe/H]
	DATA ygir/.23d0,.23d0,.23d0,.24d0,.25d0,.273d0,.30d0/
c	DATA yter/.245d0,.245d0,.246d0,.246d0,.248d0,.251d0,.256d0,
c     & .259d0,.2734d0,.288d0,.303d0/
	DATA yter/.245d0,.245d0,.246d0,.246d0,.248d0,.251d0,.256d0,
     & .259d0,.2734d0,.279d0,.288d0,.303d0/
      CALL FNORM
	zxlogG=dlog10(0.019d0/(1.d0-0.019d0-0.273d0))
	zxlogT=dlog10(0.0244d0)
c Key points Teramo:
c	DATA keypT/1,200,260,320,390,760,790,1190,1200,1350,1450,
c     & 1550,1630,1710,1850,2000,2010,2250/
	keypT(1)=1
	keypT(2)=200
        keypT(3)=260
        keypT(4)=320
        keypT(5)=390
        keypT(6)=760
        keypT(7)=790
        keypT(8)=1190
        keypT(9)=1200
        keypT(10)=1350
        keypT(11)=1450
        keypT(12)=1550
        keypT(13)=1630
        keypT(14)=1710
        keypT(15)=1850
        keypT(16)=2000
        keypT(17)=2010
        keypT(18)=2250
	dirout='ls ./OUT'
	mkout='mkdir ./OUT ./OUT/A' !coelho
        nm05T=15 !15 lineas estrellas de masa < 0.5 no presentes en isocronas
        do i=1,18
	  keypT(i)=keypT(i)+nm05T
	enddo
c ISOC TERAMO: SI MASA MINIMA NO TAN PEQUEÑA, MODIFICAR keypT
c ESTO SE HACE EN STU.f
c
c keypT(3)=  Maximum in Teff along the Main Sequence - TURN OFF POINT
c keypT(5)=***Minimum in logL for high-mass or Base of the RGB for low-mass stars
c keypT(8)=  Tip of the RGB
c keypT(9)=***Start quiescent central He-burning phase
c keypT(15)=***Central abundance of He equal to 0.00
c keypT(16)=*** When the energy produced by the CNO cycle is larger than that provided by the He burning during the AGB (Lcno > L3alpha)
c keypT(17)=  The maximum luminosity before the first Thermal Pulse
c keypT(18)=  The AGB termination
c Creamos el subdirectorio OUT si no existiese
c	dirout='ls ./OUT' 
c	mkout='mkdir ./OUT ./OUT/B ./OUT/C ./OUT/R ./OUT/M' 
	staout=system( dirout )
	if(staout.ne.0) then
      write(*,*)'Creating "./OUT" subdirectory to include output files'
      mkout2=system( mkout )
	else
      write(*,*)'Output files will be included in "./OUT" subdirectory'	 
	endif
c  si cambio orden arrays en TOAGB produce alignment warning: 
c  el orden debe ser real,natural
31963	write(*,'(A)')'The options for the IMF shape are:'
        write(*,'(A)')'UNIMODAL=u,BIMODAL=b,KROUPA=k,KROUPA BIN=r' 
        write(*,'(A)')'u/b/k/r/c ?'
	read(*,'(A1)')shape
	if(shape.eq.'b') then
	 write(*,'(A)')'---BIMODAL IMF shape---'
	elseif(shape.eq.'u') then
	 write(*,'(A)')'---UNIMODAL IMF shape---'
	elseif(shape.eq.'k') then
	 write(*,'(A)')'---KROUPA IMF shape---'
	elseif(shape.eq.'r') then
       write(*,'(A)')'---KROUPA BINARIES IMF shape (binaries corr.)---'
	elseif(shape.eq.'c') then
       write(*,'(A)')'---CHABRIER IMF shape---'
	elseif(shape.eq.'f') then
       write(*,'(A)')'--FERRERAS(turning point =0.5(limf.f))?'
       read(*,*)xtpoint
       write(*,'(A,F6.3)')'Turning Point=',xtpoint
	elseif(shape.eq.'x') then
       write(*,'(A)')'--FERRERAS X-SHAPED (turning point=0.5)(limf.f(xtpoint))--'
       write(*,'(A)')'--FERRERAS(turning point =0.5(limf.f))?'
       read(*,*)xtpoint
	elseif(shape.eq.'l') then
       write(*,'(A)')'--LaBarbera 3-segment(m1=.4,m2=.7)(limf.f(xlplow,xlpmed))-'
       write(*,'(A)')'xlplow(0.4d0)?,xlpmed(0.7d0)?'
       read(*,*)xlplow,xlpmed
c       xlplow=0.3d0 !still can be modified when asked about
c       xlpmed=0.7d0 !still can be modified when asked about
       write(*,*)xlplow,xlpmed
	 ZMU3=1.30d0
       write(*,'(A)')'ADOPTED log.slope for the third segment: 1.3'
       write(*,'(A)')'First and second segments log.slope (Salpeter=1.3)?'
       read(*,*)ZMU1,ZMU2
	else
	 goto 31963
	endif
      write(*,'(A)')'Set low mass cutoff diff.than 0.1Mo?(y/n)'
       write(*,*)'IGNORED THIS OPTION'
c	read(*,*)lowcut
        lowcut='N'
	if(lowcut.eq.'y'.or.lowcut.eq.'Y')then
      		write(*,*)'Lower mass cutoff (Msun)?'
		read(*,*)ZMlwlw
	endif
       write(*,*)'Teff scale: MILES(0),GonzalezHernandezBonifacio(1)?'
       write(*,*)'IGNORED THIS OPTION'
c      read(*,*)nghb
        nghb=0	
        if(nghb.ne.1)nghb=0
      write(*,'(A)')'Isochrone ?'
      write(*,'(A)')'Teramo(ae)    (-1)'
      write(*,'(A)')'Teramo(ss)     (0)'
      write(*,'(A)')'Padova_00      (1)'
      write(*,'(A)')'Padova_94      (2)'
      write(*,'(A)')'Salasnich_ae   (3)'
      write(*,'(A)')'Salasnich      (4)'
      read(*,*)ncha
	if(ncha.eq.0)then
	 write(*,*)'Teramo isochrones (scaled-solar)'
	elseif(ncha.eq.-1)then
	 write(*,*)'Teramo isochrones (alpha-enhanced=0.4)'
	elseif(ncha.eq.4)then
	 write(*,*)'Padova scaled-solar isochrones of Salasnich etal'
	elseif(ncha.eq.3)then
	 write(*,*)'Padova alpha-enhanced isochrones of Salasnich etal'
	elseif(ncha.eq.2)then
	 write(*,*)'OLD Padova isochrones of Bertelli etal'
	else
	 write(*,*)'NEW Padova isochrones of Girardi etal'
	endif
      WRITE(*,'(A14)')'Star spectra ?'
      WRITE(*,'(A14)')'MILES     (-1)'
      WRITE(*,'(A14)')'COELHO ss  (0)'
      WRITE(*,'(A14)')'COELHO aa  (4)'
      READ(*,*)mico
      if(mico.eq.-1)write(*,*)'Stellar library: MILES'
      if(mico.eq.0)write(*,*)'Stellar library: COELHO scaled-solar'
      if(mico.eq.4)write(*,*)'Stellar library: COELHO alpha-enhanced'
      	bicl=0.2
	bicp=0.4
	bich=0.6
	a(1,1)=bicl*bicl*bicl
	a(2,1)=3.0d0*bicl*bicl
	a(3,1)=bich*bich*bich
	a(4,1)=3.0d0*bich*bich
	a(1,2)=bicl*bicl
	a(2,2)=2.0d0*bicl
	a(3,2)=bich*bich
	a(4,2)=2.0d0*bich
	a(1,3)=bicl
	a(2,3)=1.0d0
	a(3,3)=bich
	a(4,3)=1.0d0
	a(1,4)=1.0d0
	a(2,4)=0.0d0
	a(3,4)=1.0d0
	a(4,4)=0.0d0
	call matinx(a)
c APERTURA FICHEROS PARAMETROS ESTELARES, NOMBRE ESTRELLAS Y GALAXY_CABECERAS
      nstarm=0
      IF(mico.eq.-1)THEN
	OPEN(99,file='./INPUT/GALAXY_CABECERAS_MILES',status='old')
	 do k=1,61
		read(99,'(A80)')cab(k)
	 enddo
	CLOSE(99)
c       Es probable que no sea necesario formatear PARAM_MILES
	OPEN(98,FILE='./INPUT/PARAM_MILES_MGFE_FeH',STATUS='OLD')
	 do k=1,99999
          read(98,*,end=25)anom,(aam(k,l),l=1,5)
c         Aplicamos escala GonzalezHernandezBonifacio09?
          if(nghb.eq.1)then
           if(aam(k,1).gt.3750.0d0.and.aam(k,1).lt.7500.0d0)then
            aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*aam(k,1))
           elseif(aam(k,1).le.3750.0d0)then
            aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*3750.0d0)
           elseif(aam(k,1).ge.7500.0d0)then
            aam(k,1)=aam(k,1)+(-116.d0+0.0312d0*7500.0d0)
           endif
          endif
          starm(k)='./STARS_MILES/'//anom
          nstarm=nstarm+1
         enddo
25      CLOSE(98)
      ELSEIF(mico.eq.0)THEN
	OPEN(99,file='./INPUT/GALAXY_CABECERAS_CO_ss',status='old')
	 do k=1,61
		read(99,'(A80)')cab(k)
	 enddo
	CLOSE(99)
	OPEN(98,FILE='./INPUT/PARAM_CO_ss',STATUS='OLD')
	 do k=1,99999
          read(98,*,end=26)anom,(aam(k,l),l=1,5)
          starm(k)='./STARS_CO_ss/'//anom
          nstarm=nstarm+1
         enddo
26      CLOSE(98)
      ELSEIF(mico.eq.4)THEN
	OPEN(99,file='./INPUT/GALAXY_CABECERAS_CO_aa',status='old')
	 do k=1,61
		read(99,'(A80)')cab(k)
	 enddo
	CLOSE(99)
	OPEN(98,FILE='./INPUT/PARAM_CO_aa',STATUS='OLD')
	 do k=1,99999
          read(98,*,end=27)anom,(aam(k,l),l=1,5)
          starm(k)='./STARS_CO_aa/'//anom
          nstarm=nstarm+1
         enddo
27      CLOSE(98)
      ELSE
        write(*,*)'run it again and include the value of mico: -1,0,4'
	STOP
	goto 6336
      ENDIF
      nsssm=nstarm
cccccccccccccccccccccccccccccc
c LECTURA DE ESPECTROS ASCII COELHO
c SE CARGA EN MATRIZ:starms(nstarm,npxm,2)
	DO ios=1,nstarm
	 open(99,file=starm(ios),status='old')
	 do iop=1,npxm
	   read(99,*)starms(ios,iop,1),starms(ios,iop,2)
	   if(ios.eq.1) wfm(iop,1)=starms(ios,iop,1) !escribe valores lambda
	 enddo
	 close(99)
	ENDDO
cccccccccccccccccccccccccccccc
**LECTURA DEL FICHERO dat_MASS_HR, dat_MAZT_HR y ratios.dat****************
       if(ncha.eq.0.or.ncha.eq.-1)then
          ZML=0.10d0
       elseif(ncha.eq.2)then
          ZML=0.0992d0
       else
          ZML=0.15d0
       endif
       OPEN(22,FILE='./INPUT/dat_MASS_HR',STATUS='OLD',IOSTAT=IER5)
	  READ(22,*) ZMASA
	  READ(22,*) ZMLow
	  READ(22,*) ZMUS
       CLOSE(22)
       if(lowcut.eq.'Y'.or.lowcut.eq.'y')then
		ZMlow=ZMlwlw
		if(ZMlwlw.gt.ZML) ZML=ZMlwlw
       endif
c ISOC TERAMO: SI LA MASA MINIMA NO ES TAN PEQUEÑA, los keypT(se modifican)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ASIGNAR METALICIDADES ISOCRONAS
	if(ncha.eq.0.or.ncha.eq.-1)then !TERAMO
	Z00(1)=0.0001d0
	Z00s(1)=0.00020d0
	Z00(2)=0.0003d0
	Z00s(2)=0.00045d0
	Z00(3)=0.0006d0
	Z00s(3)=0.00080d0
	Z00(4)=0.0010d0
	Z00s(4)=0.00150d0
	Z00(5)=0.0020d0
	Z00s(5)=0.00300d0
	Z00(6)=0.0040d0
	Z00s(6)=0.00600d0
	Z00(7)=0.0080d0
	Z00s(7)=0.00900d0
	Z00(8)=0.0100d0
	Z00s(8)=0.01490d0
	Z00(9)=0.0198d0
	Z00s(9)=0.02190d0
	Z00(10)=0.02400d0
	Z00s(10)=0.02700d0
	Z00(11)=0.0300d0
	Z00s(11)=0.03500d0
	Z00(12)=0.0400d0
	do nft=1,12
         fez(nft)=dlog10(Z00(nft)/(1.d0-Z00(nft)-yter(nft)))-zxlogT
	enddo
	elseif(ncha.eq.1)then !PADOVA 00
	Z00(1)=0.00010d0
	Z00s(1)=0.00025d0
	Z00(2)=0.00040d0
	Z00s(2)=0.00070d0
	Z00(3)=0.0010d0
	Z00s(3)=0.00250d0
	Z00(4)=0.0040d0
	Z00s(4)=0.00600d0
	Z00(5)=0.0080d0
	Z00s(5)=0.01400d0
        Z00(6)=0.0190d0
	Z00s(6)=0.02450d0
	Z00(7)=0.030d0
	do nft=1,7
         fez(nft)=dlog10(Z00(nft)/(1.d0-Z00(nft)-ygir(nft)))-zxlogG
c	 write(*,*)fez(nft),zxlogG
	enddo
	elseif(ncha.eq.2)then !PADOVA 94
	Z00(1)=0.00010d0
	Z00s(1)=0.00025d0
	Z00(2)=0.00040d0
	Z00s(2)=0.00070d0
	Z00(3)=0.0010d0
	Z00s(3)=0.00250d0
	Z00(4)=0.0040d0
	Z00s(4)=0.00600d0
	Z00(5)=0.0080d0
	Z00s(5)=0.01400d0
        Z00(6)=0.020d0
	Z00s(6)=0.03500d0
	Z00(7)=0.050d0
	do nft=1,7
         fez(nft)=dlog10(Z00(nft)/0.0190d0)
	enddo
	elseif(ncha.eq.3.or.ncha.eq.4)then !Salasnich 00
	Z00(1)=0.0080d0
	Z00s(1)=0.01400d0
        Z00(2)=0.019d0
	Z00s(2)=0.0300d0
	Z00(3)=0.040d0
	do nft=1,3
         fez(nft)=dlog10(Z00(nft)/0.0190d0)
	enddo
	else
	write(*,*)'isochrone metalicity and type to be specified'
	goto 6336
	endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do ii=1,17
		ze(ii)=dble(ii*1000)
	enddo
	OPEN(22,FILE='dat_MAZT_HR',STATUS='OLD',IOSTAT=IER5)
        READ(22,'(A25)')porque
        READ(22,'(A25)')porque
	nnIMF=0
	nalpha=0
	nzold=0
	ntold=0
	nzyng=0
	ntyng=0
        DO i=3,100
         read(22,*,end=322)qmu,qalpha,qqqq1,qqqq2,qqqq3,qqqq4
         if(qmu.lt.99)then
         	nnIMF=nnIMF+1
                zi(nnIMF)=qmu
         endif
         if(qalpha.lt.99)then
         	nalpha=nalpha+1
                zalpha(nalpha)=qalpha
         endif
	 if(qqqq1.le.0.0060d0.or.qqqq3.le.0.0060d0)then
	  if(ncha.eq.3.or.ncha.eq.4)then
	   write(*,*)'Salasnich etal isochrones require [M/H]>0.006'
	   goto 6336
	  endif
	 endif
c Valor de Zsun:
         if(ncha.eq.0.or.ncha.eq.-1)then
		Zsun=0.0198d0
	 elseif(ncha.eq.2)then
		Zsun=0.020d0
	 else
		Zsun=0.0190d0
	 endif
         IF(qqqq1.lt.99)THEN
          nzold=nzold+1
	  if(ncha.eq.0.or.ncha.eq.-1)then
         	if(qqqq1.le.Z00s(1))then
                	izold(nzold)=1
		elseif(qqqq1.gt.Z00s(1).and.qqqq1.le.Z00s(2))then
                	izold(nzold)=2 
		elseif(qqqq1.gt.Z00s(2).and.qqqq1.le.Z00s(3))then
                	izold(nzold)=3 
		elseif(qqqq1.gt.Z00s(3).and.qqqq1.le.Z00s(4))then
                	izold(nzold)=4 
		elseif(qqqq1.gt.Z00s(4).and.qqqq1.le.Z00s(5))then
                	izold(nzold)=5 
		elseif(qqqq1.gt.Z00s(5).and.qqqq1.le.Z00s(6))then
                	izold(nzold)=6 
		elseif(qqqq1.gt.Z00s(6).and.qqqq1.le.Z00s(7))then
                	izold(nzold)=7 
		elseif(qqqq1.gt.Z00s(7).and.qqqq1.le.Z00s(8))then
                	izold(nzold)=8 
		elseif(qqqq1.gt.Z00s(8).and.qqqq1.le.Z00s(9))then
                	izold(nzold)=9 
		elseif(qqqq1.gt.Z00s(9).and.qqqq1.le.Z00s(10))then
                	izold(nzold)=10 
		elseif(qqqq1.gt.Z00s(10).and.qqqq1.le.Z00s(11))then
                	izold(nzold)=11 
		elseif(qqqq1.gt.Z00s(11))then
                	izold(nzold)=12 
		endif
	  elseif(ncha.eq.1)then
         	if(qqqq1.le.Z00s(1))then
                	izold(nzold)=1
		elseif(qqqq1.gt.Z00s(1).and.qqqq1.le.Z00s(2))then
                	izold(nzold)=2 
		elseif(qqqq1.gt.Z00s(2).and.qqqq1.le.Z00s(3))then
                	izold(nzold)=3 
		elseif(qqqq1.gt.Z00s(3).and.qqqq1.le.Z00s(4))then
                	izold(nzold)=4 
		elseif(qqqq1.gt.Z00s(4).and.qqqq1.le.Z00s(5))then
                	izold(nzold)=5 
		elseif(qqqq1.gt.Z00s(5).and.qqqq1.le.Z00s(6))then
                	izold(nzold)=6 
		elseif(qqqq1.gt.Z00s(6))then
                	izold(nzold)=7 
	  	endif
	  elseif(ncha.eq.2)then
         	if(qqqq1.le.Z00s(1))then
                	izold(nzold)=1
		elseif(qqqq1.gt.Z00s(1).and.qqqq1.le.Z00s(2))then
                	izold(nzold)=2 
		elseif(qqqq1.gt.Z00s(2).and.qqqq1.le.Z00s(3))then
                	izold(nzold)=3 
		elseif(qqqq1.gt.Z00s(3).and.qqqq1.le.Z00s(4))then
                	izold(nzold)=4 
		elseif(qqqq1.gt.Z00s(4).and.qqqq1.le.Z00s(5))then
                	izold(nzold)=5 
		elseif(qqqq1.gt.Z00s(5).and.qqqq1.le.Z00s(6))then
                	izold(nzold)=6 
		elseif(qqqq1.gt.Z00s(6))then
                	izold(nzold)=7 
	  	endif
	  elseif(ncha.eq.3.or.ncha.eq.4)then
         	if(qqqq1.le.Z00s(1))then
                	izold(nzold)=1
		elseif(qqqq1.gt.Z00s(1).and.qqqq1.le.Z00s(2))then
                	izold(nzold)=2 
		elseif(qqqq1.gt.Z00s(2).and.qqqq1.le.Z00s(3))then
                	izold(nzold)=3 
	  	endif
          endif
          pizold(izold(nzold))=fez(izold(nzold))
		write(*,*)izold(nzold),pizold(izold(nzold))
	 ENDIF
         IF(qqqq2.lt.99)THEN
          ntold=ntold+1
          ttold(ntold)=qqqq2*1000.0d0
         ENDIF
         IF(qqqq3.lt.99)THEN
          nzyng=nzyng+1
	  if(ncha.eq.0.or.ncha.eq.-1)then
         	if(qqqq3.le.Z00s(1))then
                	izyng(nzyng)=1
		elseif(qqqq3.gt.Z00s(1).and.qqqq3.le.Z00s(2))then
                	izyng(nzyng)=2 
		elseif(qqqq3.gt.Z00s(2).and.qqqq3.le.Z00s(3))then
                	izyng(nzyng)=3 
		elseif(qqqq3.gt.Z00s(3).and.qqqq3.le.Z00s(4))then
                	izyng(nzyng)=4 
		elseif(qqqq3.gt.Z00s(4).and.qqqq3.le.Z00s(5))then
                	izyng(nzyng)=5 
		elseif(qqqq3.gt.Z00s(5).and.qqqq3.le.Z00s(6))then
                	izyng(nzyng)=6 
		elseif(qqqq3.gt.Z00s(6).and.qqqq3.le.Z00s(7))then
                	izyng(nzyng)=7 
		elseif(qqqq3.gt.Z00s(7).and.qqqq3.le.Z00s(8))then
                	izyng(nzyng)=8 
		elseif(qqqq3.gt.Z00s(8).and.qqqq3.le.Z00s(9))then
                	izyng(nzyng)=9 
		elseif(qqqq3.gt.Z00s(9).and.qqqq3.le.Z00s(10))then
                	izyng(nzyng)=10 
		elseif(qqqq3.gt.Z00s(10).and.qqqq3.le.Z00s(11))then
                	izyng(nzyng)=11 
		elseif(qqqq3.gt.Z00s(11))then
                	izyng(nzyng)=12 
		endif
	  elseif(ncha.eq.1)then
         	if(qqqq3.le.Z00s(1))then
                	izyng(nzyng)=1
		elseif(qqqq3.gt.Z00s(1).and.qqqq3.le.Z00s(2))then
                	izyng(nzyng)=2 
		elseif(qqqq3.gt.Z00s(2).and.qqqq3.le.Z00s(3))then
                	izyng(nzyng)=3 
		elseif(qqqq3.gt.Z00s(3).and.qqqq3.le.Z00s(4))then
                	izyng(nzyng)=4 
		elseif(qqqq3.gt.Z00s(4).and.qqqq3.le.Z00s(5))then
                	izyng(nzyng)=5 
		elseif(qqqq3.gt.Z00s(5).and.qqqq3.le.Z00s(6))then
                	izyng(nzyng)=6 
		elseif(qqqq3.gt.Z00s(6))then
                	izyng(nzyng)=7 
	        endif
	  elseif(ncha.eq.2)then	 
         	if(qqqq3.le.Z00s(1))then
                	izyng(nzyng)=1
		elseif(qqqq3.gt.Z00s(1).and.qqqq3.le.Z00s(2))then
                	izyng(nzyng)=2 
		elseif(qqqq3.gt.Z00s(2).and.qqqq3.le.Z00s(3))then
                	izyng(nzyng)=3 
		elseif(qqqq3.gt.Z00s(3).and.qqqq3.le.Z00s(4))then
                	izyng(nzyng)=4 
		elseif(qqqq3.gt.Z00s(4).and.qqqq3.le.Z00s(5))then
                	izyng(nzyng)=5 
		elseif(qqqq3.gt.Z00s(5).and.qqqq3.le.Z00s(6))then
                	izyng(nzyng)=6 
		elseif(qqqq3.gt.Z00s(6))then
                	izyng(nzyng)=7 
	        endif
	  elseif(ncha.eq.3.or.ncha.eq.4)then
         	if(qqqq3.le.Z00s(1))then
                	izyng(nzyng)=1
		elseif(qqqq3.gt.Z00s(1).and.qqqq3.le.Z00s(2))then
                	izyng(nzyng)=2 
		elseif(qqqq3.gt.Z00s(2).and.qqqq3.le.Z00s(3))then
                	izyng(nzyng)=3 
	  	endif
          endif
          pizyng(izyng(nzyng))=fez(izyng(nzyng))
	 ENDIF
         IF(qqqq4.lt.99)THEN
          ntyng=ntyng+1
          ttyng(ntyng)=qqqq4*1000.0d0
         ENDIF
        ENDDO
322     CLOSE(22)
        if(zalpha(1).eq.zalpha(2).and.zalpha(1).eq.1.0d0)then
        	ntyng=1
        	ttyng(1)=ttold(1)
        	nzyng=1
        	izyng(1)=izold(1)
        endif
c ahora nos aseguramos que si la IMF es Kroupa solo hay
c un solo calculo
	IF(shape.eq.'k'.or.shape.eq.'r')THEN
         nnIMF=1
         zi(nnIMF)=1.30 !de todas formas no se usara en el calculo
      ENDIF
c*******LECTURA DE ISOCRONAS***********************************************
       DO 1500 inii=1,nzold
	m=izold(inii)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IF(ncha.eq.0)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0003T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.12)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.-1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0003T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.12)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0010_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0040_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0080_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0190_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0300_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
        ELSEIF(ncha.eq.2)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z001NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z004NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z008NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z02NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z05NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
	ELSE
	 if(m.eq.1)then !Z=0.008
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z008_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z008_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.2)then !Z=0.019 or Z=0.020
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z019_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z019_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.3)then !Z=0.04
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z040_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z040_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.3)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif		
	ENDIF		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO I=1,300000
	   READ(99,*,end=11111)(z(m,I,J),J=1,16)
	   num(m)=I
	ENDDO
11111	CLOSE (99)
	if(ncha.gt.0)then
	 DO I=1,30000
       	   READ(23,*,end=11107)(g(m,I,J),J=1,7),tin,tin,tin,gc(m,I)
       	   ng(m)=I
	 ENDDO
	endif
11107	CLOSE (23)
	DO I=1,300
       	   READ(24,*,end=11108)tin,re(m,1,I),re(m,2,I)
       	   nr(m)=I
	ENDDO
11108	CLOSE (24)
1500   ENDDO
       DO 1501 inii=1,nzyng
	nttico=0
	do j=1,nzold
	 if(izyng(inii).eq.izold(j))then
		nttico=1
	 endif
	enddo
       IF(nttico.ne.1)THEN
	m=izyng(inii)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IF(ncha.eq.0)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0003T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T_ss',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.-1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0003T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T_aa',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0010_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0040_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0080_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0190_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0300_G',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
        ELSEIF(ncha.eq.2)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z001NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z004NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z008NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z02NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z05NNN',STATUS='OLD')
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
	ELSE
	 if(m.eq.1)then !Z=0.008
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z008_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z008_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum008',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.2)then !Z=0.019 or Z=0.020
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z019_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z019_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum019',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.3)then !Z=0.04
		 if(ncha.eq.3)then
		  OPEN(99,FILE='./INPUT/Z040_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='./INPUT/Z040_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='./INPUT/sum030',STATUS='OLD')
		 OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.3)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif		
	ENDIF		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO I=1,300000
	  READ(99,*,end=11112)(z(m,I,J),J=1,16)
	  num(m)=I
	ENDDO
11112	CLOSE (99)
	if(ncha.gt.0)then
	 DO I=1,30000
       	   READ(23,*,end=11106)(g(m,I,J),J=1,7),tin,tin,tin,gc(m,I)
       	   ng(m)=I
	 ENDDO
	endif
11106	CLOSE (23)
	DO I=1,300
       	  READ(24,*,end=11105)tin,re(m,1,I),re(m,2,I)
       	  nr(m)=I
	ENDDO
11105	CLOSE (24)
       ENDIF
1501   ENDDO
c********LECTURA DE LOS FICHEROS 'tiempos'**************************
       if(ncha.eq.0.or.ncha.eq.-1)then
        OPEN(99,FILE='./INPUT/tiemposT',STATUS='OLD')
       else
        OPEN(99,FILE='./INPUT/tiemposP',STATUS='OLD')
       endif
       numtis=0
       DO I=1,30000
	 READ(99,*,end=12345) tiso(I)
		numtis=I
       ENDDO
12345  CLOSE (99)
c*******APERTURA DE LOS FICHEROS out_**************************************
	OPEN(29,FILE='./OUT/out_u',STATUS='OLD',ERR=6329)
	CLOSE(29,STATUS='DELETE')
6329	OPEN(29,IOSTAT=IOS,FILE='./OUT/out_u',STATUS='NEW')
c	OPEN(30,FILE='./OUT/out_li',STATUS='OLD',ERR=6330)
c	CLOSE(30,STATUS='DELETE')
c6330	OPEN(30,IOSTAT=IOS,FILE='./OUT/out_li',STATUS='NEW')
c	OPEN(31,FILE='./OUT/out_ef',STATUS='OLD',ERR=6331)
c	CLOSE(31,STATUS='DELETE')
c6331	OPEN(31,IOSTAT=IOS,FILE='./OUT/out_ef',STATUS='NEW')
c	OPEN(32,FILE='./OUT/out_dg',STATUS='OLD',ERR=6332)
c	CLOSE(32,STATUS='DELETE')
c6332	OPEN(32,IOSTAT=IOS,FILE='./OUT/out_dg',STATUS='NEW')
	OPEN(33,FILE='./OUT/out_sbf',STATUS='OLD',ERR=6333)
	CLOSE(33,STATUS='DELETE')
6333	OPEN(33,IOSTAT=IOS,FILE='./OUT/out_sbf',STATUS='NEW')
	OPEN(34,FILE='./OUT/out_phot',STATUS='OLD',ERR=6334)
	CLOSE(34,STATUS='DELETE')
6334	OPEN(34,IOSTAT=IOS,FILE='./OUT/out_phot',STATUS='NEW')
	OPEN(35,FILE='./OUT/out_contr',STATUS='OLD',ERR=6335)
	CLOSE(35,STATUS='DELETE')
6335	OPEN(35,IOSTAT=IOS,FILE='./OUT/out_contr',STATUS='NEW')
	OPEN(40,FILE='./OUT/out_mass',STATUS='OLD',ERR=6339)
	CLOSE(40,STATUS='DELETE')
6339	OPEN(40,IOSTAT=IOS,FILE='./OUT/out_mass',STATUS='NEW')
c*******COMIENZO DEL BUCLE DE MU*******************************************
	DO NNN=1,nnIMF
c	    ZMU=zi(NNN)
		ZMU=zi(NNN)+1.0d0
		write(chato(1:4),'(f4.2)')ZMU-1.0d0
		if(shape.eq.'b')then
		 chatom='./OUT/A/Abi'//chato
		elseif(shape.eq.'u')then
		 chatom='./OUT/A/Aun'//chato
		elseif(shape.eq.'k')then
		 chatom='./OUT/A/Aku'//chato
		elseif(shape.eq.'r')then
		 chatom='./OUT/A/Akb'//chato
		elseif(shape.eq.'c')then
		 chatom='./OUT/A/Ach'//chato
		endif
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
		DO ialpha=1,nalpha      
	alpha=zalpha(ialpha)
        if(alpha.eq.1.00d0)then
        	ntyng=1
        	ttyng(1)=ttold(1)
        	nzyng=1
        	izyng(1)=izold(1)
        endif
c*******COMIENZO DEL BUCLE DE LA SSPo***********************************
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
      write(57,*)' '
      write(57,'(A21)')chatom
      write(57,'(A2,1x,(F5.2,1x),2(F7.4,1x))')
     &shape,ZMU-1.,fehf,ageo*0.001d0
	call STU_alphaCO(npxm,izold(izo),ageo,BETAo,realeo,ZZZo,ZMISOo,vsino,
     & stmo,wwffbo,wwffvo,wwffco,wwffmo,a9o,f9o,a9co,f9co,a9mo,f9mo,
     & ZISOLo,fZISOLo,VMo,dwgio,ZZsbfo,fhsto,fhsbfo,t222o,g222o,f222o,
     & Wfcatto,ho,fluxxo,qumo,qu10o,qu05o,qu15o,fxUmo,fxUio,fxUo,
     & fxUmso,fxUiso,fxUso,fxUis004o)
	attao=realeo
	if(attao.lt.10.0d0)then
	 write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	else
	 write(chato(1:8),'(A,f7.4)')'T',attao
	endif
  	write(chatom(22:29),'(A8)')chato
	if(ncha.eq.0)then
	 write(chato(1:8),'(A8)')'_iTp0.00'
	elseif(ncha.eq.-1)then
	 write(chato(1:8),'(A8)')'_iTp0.40'
	elseif(ncha.eq.1)then
	 write(chato(1:8),'(A8)')'_iPp0.00'
	endif
  	write(chatom(30:37),'(A8)')chato
	if(mico.eq.4)then
	 write(chato(1:7),'(A7)')'_Ap0.40'
	else
	 write(chato(1:7),'(A7)')'_Ap0.00'
	endif
  	write(chatom(38:44),'(A7)')chato
	write(*,*)chatom
c*******COMIENZO DEL BUCLE DE LA SSPy***********************************
				  DO izy=1,nzyng
	if(alpha.ne.1.0d0)then
  	 write(chatom(45:49),'(A,f4.2)')'a',alpha
	  if(pizyng(izyng(izy)).ge.0.0d0)then
	 write(chato(1:6),'(A,A,f4.2)')'Z','p',pizyng(izyng(izy))
	  else
	 write(chato(1:6),'(A,A,f4.2)')'Z','m',abs(pizyng(izyng(izy)))
	  endif
  	 write(chatom(50:55),'(A6)')chato
	endif
				    DO iay=1,ntyng
					agey=ttyng(iay)
	if(agey.gt.ageo.and.alpha.ne.1.0d0)goto 2325
	IF(alpha.ne.1.0d0)THEN
	call STU_alphaCO(npxm,izyng(izy),agey,BETAy,realey,ZZZy,ZMISOy,vsiny,
     & stmy,wwffby,wwffvy,wwffcy,wwffmy,a9y,f9y,a9cy,f9cy,a9my,f9my,
     & ZISOLy,fZISOLy,VMy,dwgiy,ZZsbfy,fhsty,fhsbfy,t222y,g222y,f222y,
     & Wfcatty,hy,fluxxy,qumy,qu10y,qu05y,qu15y,fxUmy,fxUiy,fxUy,
     & fxUmsy,fxUisy,fxUsy,fxUis004y)
c	call STU_alphaCO(npxm,izyng(izy),agey,BETAy,realey,ZZZy,ZMISOy,
c     &vsiny,stmy,wwffmy,a9my,f9my,ZISOLy,fZISOLy,VMy,dwgiy,ZZsbfy,fhsty
c     &,fhsbfy,hy,fluxxy,qumy,qu10y,qu05y,qu15y,fxUmy,fxUy,fxUmsy,fxUsy)
c
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
     	  do ini=1,8
     	  	ZZZy(ini)=0.0d0
     	  enddo
     	  ZMISOy=0.0d0
     	  do ini=1,npxm !miles
     	  	wwffmy(ini,2)=0.0d0
      	  enddo
     	  f9my=0.0d0 !miles
     	  do ini=1,3
     	  	a9my(ini)=99.0d0
     	  enddo
     	ENDIF
c	 write(*,'(A51)')chatom		
c miles
     	do nt=1,3 !miles
	a90m(nt)=(alpha*BETAo*a9mo(nt)+(1.0d0-alpha)*BETAy*a9my(nt))/
     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	enddo
     	do nt=1,4
	 fhst(nt)=(-2.50d0)*dlog10(alpha*BETAo*fhsto(nt)+
     &   (1.-alpha)*BETAy*fhsty(nt))
	 fhsbf(nt)=(-2.50d0)*dlog10(
     &   (alpha*BETAo*fhsbfo(nt)+(1.-alpha)*BETAy*fhsbfy(nt))/
     &   (alpha*BETAo*fhsto(nt)+(1.-alpha)*BETAy*fhsty(nt)))
        enddo
c U fluxes
        fxUtm=(alpha*BETAo*fxUmo+(1.0d0-alpha)*BETAy*fxUmy)
        fxUt=(alpha*BETAo*fxUo+(1.0d0-alpha)*BETAy*fxUy)
        fxUtms=(alpha*BETAo*fxUmso+(1.0d0-alpha)*BETAy*fxUmsy)
        fxUts=(alpha*BETAo*fxUso+(1.0d0-alpha)*BETAy*fxUsy)
	DO nt=1,8
	 ZZZ(nt)=(-2.50d0)*dlog10(alpha*BETAo*ZZZo(nt)+
     &   (1.0d0-alpha)*BETAy*ZZZy(nt))
	 SBF(nt)=(-2.50d0)*dlog10(
     &   (alpha*BETAo*ZZsbfo(nt)+(1.0d0-alpha)*BETAy*ZZsbfy(nt))/
     &   (alpha*BETAo*ZZZo(nt)+(1.0d0-alpha)*BETAy*ZZZy(nt)))
         do ndo=1,5
          hh(ndo,nt)=100.0d0*(alpha*BETAo*ho(ndo,nt)+(1.0d0-alpha)*
     &    BETAy*hy(ndo,nt))/
     &    (alpha*BETAo*ZZZo(nt)+(1.0d0-alpha)*BETAy*ZZZy(nt))
         enddo
	 do ndo=1,2
	  DWGIT(ndo,nt)=alpha*BETAo*dwgio(ndo,nt)+
     &    (1.0d0-alpha)*BETAy*dwgiy(ndo,nt)
         enddo
         DWGI(1,nt)=DWGIT(1,nt)*100.d0/(alpha*BETAo*ZZZo(nt)+
     &   (1.0d0-alpha)*BETAy*ZZZy(nt))
         DWGI(2,nt)=DWGIT(2,nt)*100.0d0/(alpha*BETAo*ZZZo(nt)+
     &   (1.0d0-alpha)*BETAy*ZZZy(nt))
	ENDDO
        do nt=1,5
         hh(nt,9)=100.0d0*(alpha*BETAo*ho(nt,9)+(1.0d0-alpha)*BETAy*
     &   hy(nt,9))/(alpha*BETAo*fluxxo+(1.0d0-alpha)*BETAy*fluxxy)
        enddo
c   Parametros calidad MILES
     	quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
       do io=1,npxm !miles
        wfm(io,2)=alpha*BETAo*wwffmo(io,2)+(1.d0-alpha)*
     &  BETAy*wwffmy(io,2)
       enddo
c miles
	OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	CLOSE(77,STATUS='DELETE')
8979	OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	do nonoo=1,61
	  write(77,'(A80)')cab(nonoo)
	enddo
	do lup=1,npxm
		write(77,*)(wfm(lup,kini),kini=1,2)
	enddo
	close(77)
c   ATENCION: solo se calcula la M/L para la pobl. vieja y no para 2pob.
	do itton=1,8
		VM(itton)=VMo(itton)
		ZISOL(itton)=ZISOLo(itton)
		fZISOL(itton)=fZISOLo(itton)
	enddo
	if(alpha.eq.1.) then
	 do izsoly=1,8
	 ZISOLy(izsoly)=0.0d0
	 enddo
	endif
c        write(*,*)'ZZZ(3),ZISOL(3),fZISOL(3)',ZZZ(3),ZISOL(3),fZISOL(3)
c	write(*,*)'realeo escrito desde a',realeo
	call e_alphaCO(alpha,izold(izo),realeo,izyng(izy),realey)
c        write(*,*)'ZZZ(3),ZISOL(3),fZISOL(3)',ZZZ(3),ZISOL(3),fZISOL(3)
2325	continue
				    ENDDO !Cerrado loop Tyng
				  ENDDO !Cerrado loop Zyng
				ENDDO !Cerrado loop Told
			ENDDO !Cerrado loop Zold
		ENDDO !Cerrado loop alpha
	ENDDO !Cerrado loop IMF slope
        close(29)
c	close(30)
c	close(31)
c	close(32)
	close(33)
	close(34)
	close(35)
	close(40)
	STOP
6336	END




