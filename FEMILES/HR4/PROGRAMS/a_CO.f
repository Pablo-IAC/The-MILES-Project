       	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c	DOUBLE PRECISION a(4,4),bb(4),sson(4)
	PARAMETER (npxm=4300) !miles
c	PARAMETER (nstm=1100) !miles
	PARAMETER (nstm=1999) !coelho
	DIMENSION a(4,4),bb(4),sson(4)
cc	DIMENSION wfb(1120,2),wwffbo(1120,2),wwffby(1120,2) !jones
cc	DIMENSION wfv(1120,2),wwffvo(1120,2),wwffvy(1120,2) !jones
cc	DIMENSION wfc(710,2),wwffco(710,2),wwffcy(710,2) !cat
	DIMENSION wfm(npxm,2),wwffmo(npxm,2),wwffmy(npxm,2) !miles
cc	DIMENSION a9o(3),a9y(3),a90(3) !jones
cc	DIMENSION a9co(3),a9cy(3),a90c(3) !cat
	DIMENSION a9mo(3),a9my(3),a90m(3) !miles
	DIMENSION izold(15),izyng(15),ttold(99),ttyng(99)
c	DIMENSION pizold(15),pizyng(15),yter(11),ygir(7)
	DIMENSION pizold(15),pizyng(15),yter(12),ygir(7)
c	DIMENSION Q(150,15)
	DIMENSION cab(150,4) !cabeceras de jones (b,v), cat, miles
	DIMENSION ZZZ(8),ZZZo(8),ZZZy(8),SBF(50),ZZsbfo(50),
     &  ZZsbfy(50),DWGIT(2,8),DWGI(2,8),dwgio(2,8),dwgiy(2,8),tiso(75)
        DIMENSION ZISOLo(8),ZISOLy(8),ZISOL(8),VMo(8),VMy(8),VM(8)
        DIMENSION fZISOLo(8),fZISOLy(8),fZISOL(8)
        DIMENSION fhsto(4),fhsbfo(4),fhsty(4),fhsbfy(4),fhst(4),fhsbf(4)
	DIMENSION ze(75),zi(150),num(15),co(50)
	DIMENSION zalpha(22),vsin(50),vsino(50),vsiny(50),stmo(50),stmy(50)
c	DIMENSION z(7,14250,15)
	DIMENSION z(12,150000,16),Z00s(15) !para TERAMO, Z00(15) common/eststu
	DIMENSION ho(5,9),hy(5,9)
        DIMENSION aam(nstm,4),starm(nstm) !miles
c	DIMENSION keypT(18)
c	integer keypT
c	character*8 dirout / 'ls ./OUT' /
c	character*43 mkout / 'mkdir ./OUT ./OUT/B ./OUT/C ./OUT/R ./OUT/M' /
	character*8 dirout
	character*43 mkout
	integer staout,mkout2,system
        CHARACTER*4 gc          !no cambiar que sea de 4, fich sum*
	CHARACTER*1 shape,lowcut
	CHARACTER*10 co
	CHARACTER*80 cab,chatb,chator,chatoc,chatom,chato
	CHARACTER*20 ast,star,starc,starm
	CHARACTER*9 anom
c	COMMON/ll/Q,nyr
	COMMON/fezsol/fez(15)
	COMMON/lmu/ZMU,bicl,bicp,bich
	COMMON/lshape/shape
	COMMON/lson/sson
	COMMON/stus/ZMASA,tiso,num,numtis    
	COMMON/zsolo/z
	COMMON/tgf/t222,g222,f222
	COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
        COMMON/rat/ra(3,50)
        COMMON/chrbus/ast(650,2)!jones
        COMMON/hrbus/aa(650,6),nstar !jones
	COMMON/lstd/nsss1,nsss2 !jones
	COMMON/clstd/star(650,2) !jones
        COMMON/hrcat/aac(710,4),nstarc !cat
	COMMON/lstdc/nsssc !cat
	COMMON/clstdc/starc(710) !cat
        COMMON/hrm1/aam,nstarm !miles
        COMMON/hrm2/starm !miles
	COMMON/lstdm/nsssm !miles
      COMMON/esc/vsin,ZZZ,SBF,DWGI,a90,a90c,ZISOL,fZISOL,ncols,nobs,
     & fhst,fhsbf
        COMMON/TOAGB2/hh(5,9)
        COMMON/TOAGB/g(15,2000,7),ng(15)
        COMMON/TOAGBC/gc(15,2000)
	COMMON/KEYS/keypT(18)
	COMMON/REMNAN/re(15,2,99),nr(15)
        COMMON/CTFIN/CATFIN
        COMMON/NC/ncha !indice tipo isocronas
	COMMON/qua/quaqn,quaq10,quaq05,quaq15
	COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004
        COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu
c	COMMON/FNOR/C2Code,FNORMM,FNORMC
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
c	mkout='mkdir ./OUT ./OUT/B ./OUT/C ./OUT/R ./OUT/M'
	mkout='mkdir ./OUT ./OUT/0 ./OUT/4 ./OUT/M ./OUT/A' !coelho
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
c	****************IMPORTANTEIMPORTANTE*************
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
	co(14)='Fe4531	 ='
	co(15)='Fe4668	 ='
	co(16)='H-beta	 ='
	co(17)='Fe5015	 ='
	co(18)='Mg1 (mag)='
	co(19)='Mg2 (mag)='
	co(20)='Mgb	 ='
	co(21)='Fe5270	 ='
	co(22)='Fe5335	 ='
	co(23)='Fe5406	 ='
	co(24)='Fe5709	 ='
	co(25)='Fe5782	 ='
	co(26)='NaD	 ='
	co(27)='TiO1(mag)='
	co(28)='TiO2(mag)='
c	co(29)='D4000    ='
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
c
c	co(29)='CaII(1)  ='
c	co(30)='CaII(2)  ='
c	co(32)='MgI	 ='
31963	write(*,'(A)')'The options for the IMF shape are:'
        write(*,'(A)')'UNIMODAL=u,BIMODAL=b,KROUPA=k,KROUPA BIN=r' 
        write(*,'(A)')'u/b/k/r ?'
	read(*,'(A1)')shape
c	shape='u'
	if(shape.eq.'b') then
	 write(*,'(A)')'---BIMODAL IMF shape---'
	elseif(shape.eq.'u') then
	 write(*,'(A)')'---UNIMODAL IMF shape---'
	elseif(shape.eq.'k') then
	 write(*,'(A)')'---KROUPA IMF shape---'
	elseif(shape.eq.'r') then
       write(*,'(A)')'---KROUPA BINARIES IMF shape (binaries corr.)---'
	else
	 goto 31963
	endif
      write(*,'(A)')'Set low mass cutoff diff.than 0.1Mo?(y/n)'
	read(*,*)lowcut
c	lowcut='n'
	if(lowcut.eq.'y'.or.lowcut.eq.'Y')then
      		write(*,*)'Lower mass cutoff (Msun)?'
		read(*,*)ZMlwlw
	endif
c
       write(*,*)'MILES (0) or GonzalezHernandezBonifacio (1) Teff scale?'
       read(*,*)nghb
       if(nghb.ne.1)nghb=0 	
      write(*,'(A)')'Teramo(ss)(0),Padova_00(1),Padova_94(2),Salasnich_alpha(3),Salasnich(4)? (1)'
	read(*,*)ncha
	if(ncha.eq.0)then
c	 ncha=0
	 write(*,*)'Teramo isochrones of Pietrinferni etal (scaled-solar)'
c	ncha=1
	elseif(ncha.eq.4)then
c	 ncha=4
	 write(*,*)'Padova scaled-solar isochrones of Salasnich etal'
	elseif(ncha.eq.3)then
c	 ncha=3
	 write(*,*)'Padova alpha-enhanced isochrones of Salasnich etal'
	elseif(ncha.eq.2)then
c	 ncha=2
	 write(*,*)'OLD Padova isochrones of Bertelli etal'
	else
c	 ncha=1
	 write(*,*)'NEW Padova isochrones of Girardi etal'
	endif
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
2333	format(A9,1x,f6.0,1x,f4.2,4(1x,f5.2))
	OPEN(98,FILE='./INPUT/PARAM_STARS',STATUS='OLD')
        nstar=0
	do k=1,99999
         read(98,2333,end=23)anom,(aa(k,l),l=1,6)
         ast(k,1)='./STARS_4/4'//anom
         ast(k,2)='./STARS_5/5'//anom
         star(k,1)=ast(k,1)
         star(k,2)=ast(k,2)
         nstar=nstar+1
        enddo
23      CLOSE(98)
2334	format(A6,1x,f6.0,1x,f4.2,1x,f5.2,1x,f5.3)
c Es MUY IMPORTANTE no cambiar de formato PARAM_CAT!!!!
	OPEN(98,FILE='./INPUT/PARAM_CAT',STATUS='OLD')
        nstarc=0
	do k=1,99999
         read(98,2334,end=24)anom,(aac(k,l),l=1,4)
         starc(k)='./STARS_C/'//anom
         nstarc=nstarc+1
        enddo
24      CLOSE(98)
c2335	format(A6,1x,f6.0,2(1x,f5.2),1x,f6.2)
c Es probable que no sea necesario formatear PARAM_MILES
	OPEN(98,FILE='./INPUT/PARAM_MILES',STATUS='OLD')
        nstarm=0
	do k=1,99999
c         read(98,2335,end=25)anom,(aam(k,l),l=1,4)
         read(98,*,end=25)anom,(aam(k,l),l=1,4)
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
c         
         starm(k)='./STARS_MILES/'//anom
         nstarm=nstarm+1
        enddo
25      CLOSE(98)
         nsss1=nstar
         nsss2=nstar
         nsssc=nstarc
         nsssm=nstarm
c	nsss1=0
c	open(98,file='LISTA_4WHD',status='old')
c	do nonoo=1,1000
c		read(98,'(A10)',end=56)star(nonoo,1)
c		nsss1=nsss1+1
c	enddo
c56	close(98)
c	nsss2=0
c	open(99,file='LISTA_5WHD',status='old')
c	do nonoo=1,1000
c		read(99,'(A10)',end=57)star(nonoo,2)
c		nsss2=nsss2+1
c	enddo
c57	close(99)
c	Escribir los valores de las longitudes de onda
	open(99,file=star(1,1),status='old')
	do io=1,1107
		read(99,*)wfb(io,1)
	enddo
	close(99)
	open(99,file=star(1,2),status='old')
	do io=1,1107
		read(99,*)wfv(io,1)
	enddo
	close(99)
c
	open(99,file=starc(1),status='old')
	do io=1,710
		read(99,*)wfc(io,1)
	enddo
	close(99)
c
	open(99,file=starm(1),status='old')
	do io=1,npxm
		read(99,*)wfm(io,1)
	enddo
	close(99)
cccc
	open(99,file='./INPUT/GALAXY_CABECERAS_MILES',status='old')
	do k=1,61
		read(99,'(A80)')cab(k,1)
	enddo
	close(99)
	open(99,file='./INPUT/GALAXY_CABECERAS_MILESstatus='old')
	do k=1,61
		read(99,'(A80)')cab(k,2)
	enddo
	close(99)
	open(99,file='./INPUT/GALAXY_CABECERAS_MILES',status='old')
	do k=1,61
		read(99,'(A80)')cab(k,3)
	enddo
	close(99)
	open(99,file='./INPUT/GALAXY_CABECERAS_MILES',status='old')
	do k=1,61
		read(99,'(A80)')cab(k,4)
	enddo
	close(99)
c	elmejc=1.e9
c	elmejl=1.e9
c	elmejt=1.e9
c	indice=0
c	OPEN(19,FILE='obs_dat',STATUS='OLD')
c	DO I=1,100
c		READ(19,*,end=10004) (FILT(I,J),J=1,3)
c		FILT(I,2)=abs(FILT(I,2))
c		nobs=I
c	ENDDO
c10004	CLOSE(19)
c	do i1=1,nobs
c		vobs(i1)=FILT(i1,1)
c	enddo
c	NUMERO DE COLORES A USAR:
c	ncols=7
c	nobs=
c	FILT(10,1)=-2.5*log10(abs(1.-(FILT(10,1)/12.5)))
c	FILT(10,2)=abs((2.5/log(10.))*(FILT(10,2)/(12.5-FILT(10,2))))
c	FILT(11,1)=-2.5*log10(abs(1.-(FILT(11,1)/35.0)))
c	FILT(11,2)=abs((2.5/log(10.))*(FILT(11,2)/(35.0-FILT(11,2))))
c	FILT(12,1)=-2.5*log10(abs(1.-(FILT(12,1)/51.25)))
c	FILT(12,2)=abs((2.5/log(10.))*(FILT(12,2)/(51.25-FILT(12,2))))
c	FILT(13,1)=-2.5*log10(abs(1.-(FILT(13,1)/22.5)))
c	FILT(13,2)=abs((2.5/log(10.))*(FILT(13,2)/(22.5-FILT(13,2))))
c	FILT(14,1)=-2.5*log10(abs(1.-(FILT(14,1)/45.0)))
c	FILT(14,2)=abs((2.5/log(10.))*(FILT(14,2)/(45.0-FILT(14,2))))
c	FILT(15,1)=-2.5*log10(abs(1.-(FILT(15,1)/86.25)))
c	FILT(15,2)=abs((2.5/log(10.))*(FILT(15,2)/(86.25-FILT(15,2))))
c	FILT(16,1)=-2.5*log10(abs(1.-(FILT(16,1)/28.75)))
c	FILT(16,2)=abs((2.5/log(10.))*(FILT(16,2)/(28.75-FILT(16,2))))
c	FILT(17,1)=-2.5*log10(abs(1.-(FILT(17,1)/76.25)))
c	FILT(17,2)=abs((2.5/log(10.))*(FILT(17,2)/(76.25-FILT(17,2))))
c	FILT(20,1)=-2.5*log10(abs(1.-(FILT(20,1)/32.5)))
c	FILT(20,2)=abs((2.5/log(10.))*(FILT(20,2)/(32.5-FILT(20,2))))
c	FILT(21,1)=-2.5*log10(abs(1.-(FILT(21,1)/40.0)))
c	FILT(21,2)=abs((2.5/log(10.))*(FILT(21,2)/(40.0-FILT(21,2))))
c	FILT(22,1)=-2.5*log10(abs(1.-(FILT(22,1)/40.0)))
c	FILT(22,2)=abs((2.5/log(10.))*(FILT(22,2)/(40.0-FILT(22,2))))
c	FILT(23,1)=-2.5*log10(abs(1.-(FILT(23,1)/27.5)))
c	FILT(23,2)=abs((2.5/log(10.))*(FILT(23,2)/(27.5-FILT(23,2))))
c	FILT(24,1)=-2.5*log10(abs(1.-(FILT(24,1)/23.75)))
c	FILT(24,2)=abs((2.5/log(10.))*(FILT(24,2)/(23.75-FILT(24,2))))
c	FILT(25,1)=-2.5*log10(abs(1.-(FILT(25,1)/20.0)))
c	FILT(25,2)=abs((2.5/log(10.))*(FILT(25,2)/(20.0-FILT(25,2))))
c	FILT(26,1)=-2.5*log10(abs(1.-(FILT(26,1)/32.5)))
c	FILT(26,2)=abs((2.5/log(10.))*(FILT(26,2)/(32.5-FILT(26,2))))
c	FILT(29,1)=-2.5*log10(abs(1.-(FILT(29,1)/30.0)))
c	FILT(29,2)=abs((2.5/log(10.))*(FILT(29,2)/(30.0-FILT(29,2))))
c	FILT(30,1)=-2.5*log10(abs(1.-(FILT(30,1)/30.0)))
c	FILT(30,2)=abs((2.5/log(10.))*(FILT(30,2)/(30.0-FILT(30,2))))
c	FILT(31,1)=-2.5*log10(abs(1.-(FILT(31,1)/30.0)))
c	FILT(31,2)=abs((2.5/log(10.))*(FILT(31,2)/(30.0-FILT(31,2))))
c	FILT(32,1)=-2.5*log10(abs(1.-(FILT(32,1)/15.0)))
c	FILT(32,2)=abs((2.5/log(10.))*(FILT(32,2)/(15.0-FILT(32,2))))
c	FILT(33,1)=-2.5*log10(abs(1.-(FILT(33,1)/38.75)))
c	FILT(33,2)=abs((2.5/log(10.))*(FILT(33,2)/(38.75-FILT(33,2))))
c	FILT(34,1)=-2.5*log10(abs(1.-(FILT(34,1)/43.75)))
c	FILT(34,2)=abs((2.5/log(10.))*(FILT(34,2)/(43.75-FILT(34,2))))
c	FILT(35,1)=-2.5*log10(abs(1.-(FILT(35,1)/21.25)))
c	FILT(35,2)=abs((2.5/log(10.))*(FILT(35,2)/(21.25-FILT(35,2))))
c	FILT(36,1)=-2.5*log10(abs(1.-(FILT(36,1)/21.0)))
c	FILT(36,2)=abs((2.5/log(10.))*(FILT(36,2)/(21.0-FILT(36,2))))
c	do ix=1,nobs
c		if(vobs(ix).gt.98.)then
c			FILT(ix,1)=100.
c			FILT(ix,2)=0.
c		endif
c	enddo
c	CALCULO DEL MERITO ACEPTABLE***************************************
c	calculo del merito de los colores
c	rr=0.
c	rmc=0.
c	nada1=0
c	do na=1,ncols
c		if(FILT(na,1).ge.98.)then
c			rr=0.0
c			nada1=nada1+1
c		else
c			rr=FILT(na,3)
c		endif
c		rmc=rmc+rr
c	enddo
c	divc=dble(ncols-nada1)
c	calculo del merito de las lineas espectrales
c	rr=0.
c	rml=0.
c	nada2=0
c	do na=ncols+1,nobs
c		if(FILT(na,1).ge.98.)then
c			rr=0.0
c			nada2=nada2+1
c		else
c			rr=FILT(na,3)
c		endif
c		rml=rml+rr
c	enddo
c	divl=dble(nobs-(ncols+nada2))
c	calculo del merito total
c	rmt=rmc+rml*(rmc/rml)
c**LECTURA DEL FICHERO dat_MASS_HR, dat_MAZT_HR y ratios.dat****************
	OPEN(22,FILE='./INPUT/ratios.dat',STATUS='OLD',IOSTAT=IER5)
        READ(22,'(A150)')porque
        READ(22,'(A150)')porque
        do i=3,5
           read(22,*)(ra(i-2,j),j=1,20)
        enddo
        CLOSE(22)
       if(ncha.eq.0)then
          ZML=0.10d0
       elseif(ncha.eq.2)then
          ZML=0.0992d0
       else
          ZML=0.15d0
       endif
	OPEN(22,FILE='./INPUT/dat_MASS_HR',STATUS='OLD',IOSTAT=IER5)
	READ(22,*) ZMASA
	READ(22,*) ZMLow,ZMUS
c	READ(22,*) ZML,ZMUS
	CLOSE(22)
	if(lowcut.eq.'Y'.or.lowcut.eq.'y')then
		ZMlow=ZMlwlw
		if(ZMlwlw.gt.ZML) ZML=ZMlwlw
	endif
c	write(*,*)'ZMlow,ZML,ZMlwlw',ZMlow,ZML,ZMlwlw
c ISOC TERAMO: SI LA MASA MINIMA NO ES TAN PEQUEÑA, los keypT(se modifican)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ASIGNAR METALICIDADES ISOCRONAS
	if(ncha.eq.0)then !TERAMO
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
	do nft=1,11
         fez(nft)=dlog10(Z00(nft)/(1.d0-Z00(nft)-yter(nft)))-zxlogT
c	 write(*,*)fez(nft)
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
	 if(ncha.eq.0)then
		Zsun=0.0198d0
	 elseif(ncha.eq.2)then
		Zsun=0.020d0
	 else
		Zsun=0.0190d0
	 endif
         IF(qqqq1.lt.99)THEN
          nzold=nzold+1
	  if(ncha.eq.0)then
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
	  if(ncha.eq.0)then
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
c	do i=1,nzold
c                write(*,*)izold(i)
c        enddo
c	do i=1,ntold
c                write(*,*)ttold(i)
c        enddo
c	do i=1,nzyng
c                write(*,*)izyng(i)
c        enddo
c	do i=1,ntyng
c                write(*,*)ttyng(i)
c        enddo
c	write(*,*)nzold,ntold,nzyng,ntyng
c*******LECTURA DE ISOCRONAS***********************************************
       DO 1500 inii=1,nzold
	m=izold(inii)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IF(ncha.eq.0)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='./INPUT/Z0001T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T',STATUS='OLD')		  
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
		OPEN(99,FILE='./INPUT/Z0001T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='./INPUT/Z0004T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='./INPUT/Z0006T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='./INPUT/Z0010T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='./INPUT/Z0020T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='./INPUT/Z0040T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='./INPUT/Z0080T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='./INPUT/Z0100T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='./INPUT/Z0198T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='./INPUT/Z0240T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='./INPUT/Z0300T',STATUS='OLD')		  
		OPEN(24,FILE='./INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='./INPUT/Z0400T',STATUS='OLD')		  
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
       if(ncha.eq.0)then
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
c********LECTURA DEL FICHEROS 'yr'**************************
c	OPEN(99,FILE='yr',STATUS='OLD')
c	DO l=1,1000
c	  READ(99,*,end=11131) (Q(l,J),J=1,15)
c	  nyr=l
c	ENDDO
c11131	CLOSE (99)
c	if((ZMUS.lt.Q(nyr,1)).and.(ZMUS.gt.Q(1,1)))then
c	  call llim(1,ZMUS,k2,yy1,yr1)
c	  nyr=k2  
c	  do n=1,7
c	     call llim(n,ZMUS,k2,yy1,yr1)
c	     Q(nyr,2*n)=yr1
c	     Q(nyr,2*n+1)=yy1
c	  enddo
c	  Q(nyr,1)=ZMUS
c	elseif(ZMUS.le.Q(1,1))then
c	   nyr=1
c	elseif(ZMUS.ge.Q(nyr,1))then
c	   ZMUS=Q(nyr,1)	   
c	endif
c*******CALCULO DE LOS VALORES DE alpha y MU SEGUN LOS RANGOS DADOS********
c	nalpha=int((alpha2-alpha1)/dalpha+1.1)
c	salpha=alpha1-dalpha
c	do Nni=1,nalpha      
c		salpha=salpha+dalpha
c		zalpha(Nni)=salpha
c	enddo
c	nnIMF=int((slopeIMF2-slopeIMF1)/dmu+1.1)
c	ssIMF=slopeIMF1-dmu
c	do N=1,nnIMF
c		ssIMF=ssIMF+dmu
c		zi(N)=ssIMF
c	enddo
c*******INCLUSION LEY DE EXTINCION*****************************************
c	npol=int((xxpol2-xxpol1)*10.0+1.1)
c	xpolvo=xxpol1-0.10
c	do k=1,npol
c		xpolvo=xpolvo+0.1
c		zd(k)=xpolvo
c		AD(3,k)=xpolvo/1.324
c		AD(2,k)=xpolvo
c		AD(1,k)=AD(3,k)*1.531
c 			AD(4,k)=AD(3,k)*0.849
c 			AD(5,k)=AD(3,k)*0.645
c 		AD(6,k)=AD(3,k)*0.282
c 		AD(7,k)=AD(3,k)*0.175
c 		AD(8,k)=AD(3,k)*0.112
c	enddo
c*******APERTURA DE LOS FICHEROS out_**************************************
c	OPEN(29,FILE='./OUT/out_u',STATUS='OLD',ERR=6329)
cc	CLOSE(29,STATUS='DELETE')
c6329	OPEN(29,IOSTAT=IOS,FILE='./OUT/out_u',STATUS='NEW')
c	OPEN(30,FILE='./OUT/out_li',STATUS='OLD',ERR=6330)
c	CLOSE(30,STATUS='DELETE')
c6330	OPEN(30,IOSTAT=IOS,FILE='./OUT/out_li',STATUS='NEW')
c	OPEN(31,FILE='./OUT/out_ef',STATUS='OLD',ERR=6331)
c	CLOSE(31,STATUS='DELETE')
c6331	OPEN(31,IOSTAT=IOS,FILE='./OUT/out_ef',STATUS='NEW')
	OPEN(32,FILE='./OUT/out_dg',STATUS='OLD',ERR=6332)
	CLOSE(32,STATUS='DELETE')
6332	OPEN(32,IOSTAT=IOS,FILE='./OUT/out_dg',STATUS='NEW')
c	OPEN(33,FILE='./OUT/out_sbf',STATUS='OLD',ERR=6333)
c	CLOSE(33,STATUS='DELETE')
c6333	OPEN(33,IOSTAT=IOS,FILE='./OUT/out_sbf',STATUS='NEW')
c	OPEN(34,FILE='./OUT/out_phot',STATUS='OLD',ERR=6334)
c	CLOSE(34,STATUS='DELETE')
c6334	OPEN(34,IOSTAT=IOS,FILE='./OUT/out_phot',STATUS='NEW')
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
		write(chato(1:4),'(f4.2)')ZMU-1.
		if(shape.eq.'b')then
		 chatb='./OUT/0/0bi'//chato
		 chator='./OUT/4/4bi'//chato
		 chatoc='./OUT/M/Mbi'//chato
		 chatom='./OUT/A/Abi'//chato
		elseif(shape.eq.'u')then
		 chatb='./OUT/0/0un'//chato
		 chator='./OUT/4/4un'//chato
		 chatoc='./OUT/M/Mun'//chato
		 chatom='./OUT/A/Aun'//chato
		elseif(shape.eq.'k')then
		 chatb='./OUT/0/0ku'//chato
		 chator='./OUT/4/4ku'//chato
		 chatoc='./OUT/M/Mku'//chato
		 chatom='./OUT/A/Aku'//chato
		elseif(shape.eq.'r')then
		 chatb='./OUT/0/0kb'//chato
		 chator='./OUT/4/4kb'//chato
		 chatoc='./OUT/M/Mkb'//chato
		 chatom='./OUT/A/Akb'//chato
		endif
calculo de los coeficientes de la IMF bimodal
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
  	write(chatb(16:21),'(A6)')chato
  	write(chator(16:21),'(A6)')chato
  	write(chatoc(16:21),'(A6)')chato
  	write(chatom(16:21),'(A6)')chato
				DO iao=1,ntold
					ageo=ttold(iao)
	call STU(npxm,izold(izo),ageo,BETAo,realeo,ZZZo,ZMISOo,vsino,
     & stmo,wwffbo,wwffvo,wwffco,wwffmo,a9o,f9o,a9co,f9co,a9mo,f9mo,
     & ZISOLo,fZISOLo,VMo,dwgio,ZZsbfo,fhsto,fhsbfo,t222o,g222o,f222o,
     & Wfcatto,ho,fluxxo,qumo,qu10o,qu05o,qu15o,fxUmo,fxUio,fxUo,
     & fxUmso,fxUiso,fxUso,fxUis004o)
c
	attao=realeo
	if(attao.lt.10.0d0)then
	 write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	else
	 write(chato(1:8),'(A,f7.4)')'T',attao
	endif
	write(chatb(22:29),'(A8)')chato
  	write(chator(22:29),'(A8)')chato
  	write(chatoc(22:29),'(A8)')chato
  	write(chatom(22:29),'(A8)')chato
c*******COMIENZO DEL BUCLE DE LA SSPy***********************************
				  DO izy=1,nzyng
	if(alpha.ne.1.0d0)then
  	  write(chatb(30:34),'(A,f4.2)')'a',alpha
  	  write(chator(30:34),'(A,f4.2)')'a',alpha
  	  write(chatoc(30:34),'(A,f4.2)')'a',alpha
  	  write(chatom(30:34),'(A,f4.2)')'a',alpha
c
	  if(pizyng(izyng(izy)).ge.0.0d0)then
	    write(chato(1:6),'(A,A,f4.2)')'Z','p',pizyng(izyng(izy))
	  else
	    write(chato(1:6),'(A,A,f4.2)')'Z','m',abs(pizyng(izyng(izy)))
	  endif
  	  write(chatb(35:40),'(A6)')chato
  	  write(chator(35:40),'(A6)')chato
  	  write(chatoc(35:40),'(A6)')chato
  	  write(chatom(35:40),'(A6)')chato
	endif
				    DO iay=1,ntyng
					agey=ttyng(iay)
	if(agey.gt.ageo.and.alpha.ne.1.00d0)goto 2325
	IF(alpha.ne.1.0d0)THEN
	call STU(npxm,izyng(izy),agey,BETAy,realey,ZZZy,ZMISOy,vsiny,
     & stmy,wwffby,wwffvy,wwffcy,wwffmy,a9y,f9y,a9cy,f9cy,a9my,f9my,
     & ZISOLy,fZISOLy,VMy,dwgiy,ZZsbfy,fhsty,fhsbfy,t222y,g222y,f222y,
     & Wfcatty,hy,fluxxy,qumy,qu10y,qu05y,qu15y,fxUmy,fxUiy,fxUy,
     & fxUmsy,fxUisy,fxUsy,fxUis004y)
c
	  attao=realey
	  if(attao.lt.10.0d0)then
	   write(chato(1:8),'(A,A,f6.4)')'T','0',attao
	  else
	   write(chato(1:8),'(A,f7.4)')'T',attao
	  endif
 	  write(chatb(41:48),'(A8)')chato
  	  write(chator(41:48),'(A8)')chato
  	  write(chatoc(41:48),'(A8)')chato
  	  write(chatom(41:48),'(A8)')chato
	  chatb=chatb(1:48)
	  chator=chator(1:48)
	  chatoc=chatoc(1:48)
	  chatom=chatom(1:48)
     	ELSE
	  chatb=chatb(1:29)
	  chator=chator(1:29)
	  chatoc=chatoc(1:29)
	  chatom=chatom(1:29)
c
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
      	  enddo
     	  do ini=1,npxm !miles
     	  	wwffmy(ini,2)=0.0d0
      	  enddo
     	  do ini=1,3 !jones
     	  	a9y(ini)=99.0d0
     	  enddo
     	  f9y=0.0d0 !cat
     	  do ini=1,3
     	  	a9cy(ini)=99.0d0
     	  enddo
     	  f9my=0.0d0 !miles
     	  do ini=1,3
     	  	a9my(ini)=99.0d0
     	  enddo
     	ENDIF
	 write(*,'(4(A50))')chatb,chator,chatoc,chatom		
c cat
c	t222o=5040./t222o
c	t222y=5040./t222y
     	t222=(alpha*BETAo*t222o+(1.0d0-alpha)*BETAy*t222y)/
     &  (alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
c        t222=5040./t222
     	g222=(alpha*BETAo*g222o+(1.0d0-alpha)*BETAy*g222y)/
     &  (alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
     	f222=(alpha*BETAo*f222o+(1.0d0-alpha)*BETAy*f222y)/
     &  (alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
c miles
c     	t222m=(alpha*BETAo*t222mo+(1.0d0-alpha)*BETAy*t222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c     	g222m=(alpha*BETAo*g222mo+(1.0d0-alpha)*BETAy*g222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c     	f222m=(alpha*BETAo*f222mo+(1.0d0-alpha)*BETAy*f222my)/
c     &  (alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c
      	do nt=1,3 !jones
	a90(nt)=(alpha*BETAo*a9o(nt)+(1.0d0-alpha)*BETAy*a9y(nt))/
     &  (alpha*BETAo*f9o+(1.0d0-alpha)*BETAy*f9y)
     	enddo
     	do nt=1,3 !cat
	a90c(nt)=(alpha*BETAo*a9co(nt)+(1.0d0-alpha)*BETAy*a9cy(nt))/
     &  (alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
     	enddo
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
c        fxUtm=(alpha*BETAo*fxUmo+(1.0d0-alpha)*BETAy*fxUmy)
c        fxUti=(alpha*BETAo*fxUio+(1.0d0-alpha)*BETAy*fxUiy)
c        fxUt=(alpha*BETAo*fxUo+(1.0d0-alpha)*BETAy*fxUy)
c        fxUtms=(alpha*BETAo*fxUmso+(1.0d0-alpha)*BETAy*fxUmsy)
c        fxUtis=(alpha*BETAo*fxUiso+(1.0d0-alpha)*BETAy*fxUisy)
c       fxUtis004=(alpha*BETAo*fxUis004o+(1.0d0-alpha)*BETAy*fxUis004y)
c        fxUts=(alpha*BETAo*fxUso+(1.0d0-alpha)*BETAy*fxUsy)
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
         DWGI(1,nt)=DWGIT(1,nt)*100./(alpha*BETAo*ZZZo(nt)+
     &   (1.0d0-alpha)*BETAy*ZZZy(nt))
         DWGI(2,nt)=DWGIT(2,nt)*100.0d0/(alpha*BETAo*ZZZo(nt)+
     &   (1.0d0-alpha)*BETAy*ZZZy(nt))
	ENDDO
        do nt=1,5
         hh(nt,9)=100.0d0*(alpha*BETAo*ho(nt,9)+(1.0d0-alpha)*BETAy*
     &   hy(nt,9))/(alpha*BETAo*fluxxo+(1.0d0-alpha)*BETAy*fluxxy)
        enddo
	do nt=8,nobs
	vsin(nt)=(alpha*BETAo*vsino(nt)+(1.0d0-alpha)*BETAy*vsiny(nt))
     &  /(alpha*BETAo*stmo(nt)+(1.0d0-alpha)*BETAy*stmy(nt))
	enddo
     	WFFCAT=(alpha*BETAo*Wfcatto+(1.0d0-alpha)*BETAy*Wfcatty)
     &  /(alpha*BETAo*f9co+(1.0d0-alpha)*BETAy*f9cy)
c   Parametros calidad MILES
     	quaqn=(alpha*BETAo*qumo+(1.0d0-alpha)*BETAy*qumy)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq10=(alpha*BETAo*qu10o+(1.0d0-alpha)*BETAy*qu10y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq05=(alpha*BETAo*qu05o+(1.0d0-alpha)*BETAy*qu05y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
     	quaq15=(alpha*BETAo*qu15o+(1.0d0-alpha)*BETAy*qu15y)
     &  /(alpha*BETAo*f9mo+(1.0d0-alpha)*BETAy*f9my)
c        write(*,*)WFFCAT
	vsin(8)=(-2.50d0)*dlog10(1.0d0-(vsin(8)/35.0d0))
	vsin(9)=(-2.50d0)*dlog10(1.0d0-(vsin(9)/35.0d0))
	vsin(18)=(-2.50d0)*dlog10(1.0d0-(vsin(18)/65.0d0))
	vsin(19)=(-2.50d0)*dlog10(1.0d0-(vsin(19)/42.50d0))
	vsin(27)=(-2.50d0)*dlog10(1.0d0-(vsin(27)/57.50d0))
	vsin(28)=(-2.50d0)*dlog10(1.0d0-(vsin(28)/82.50d0))
       do io=1,mpxm !jones
      wfb(io,2)=alpha*BETAo*wwffbo(io,2)+(1.d0-alpha)*BETAy*wwffby(io,2)		
      wfv(io,2)=alpha*BETAo*wwffvo(io,2)+(1.d0-alpha)*BETAy*wwffvy(io,2)
       enddo
       do io=1,npxm !cat
      wfc(io,2)=alpha*BETAo*wwffco(io,2)+(1.d0-alpha)*BETAy*wwffcy(io,2)
       enddo
       do io=1,npxm !miles
      wfm(io,2)=alpha*BETAo*wwffmo(io,2)+(1.d0-alpha)*BETAy*wwffmy(io,2)
       enddo
c jones b
	OPEN(77,FILE=chatb,STATUS='OLD',ERR=8976)
	CLOSE(77,STATUS='DELETE')
8976	OPEN(77,IOSTAT=IOS,FILE=chatb,STATUS='NEW')
	do nonoo=1,61
	  write(77,'(A80)')cab(nonoo,1)
	enddo
	do nonoo=1,1107
	  if(wfb(nonoo,1).gt.3855.4d0.and.wfb(nonoo,1).lt.4476.5d0)then
		write(77,*)(wfb(nonoo,kini),kini=1,2)
	  endif
	enddo
	close(77)
c jones v
	OPEN(77,FILE=chator,STATUS='OLD',ERR=8977)
	CLOSE(77,STATUS='DELETE')
8977	OPEN(77,IOSTAT=IOS,FILE=chator,STATUS='NEW')
	do nonoo=1,61
	  write(77,'(A80)')cab(nonoo,2)
	enddo
	do nonoo=1,1107
	if(wfv(nonoo,1).gt.4794.9d0.and.wfv(nonoo,1).lt.5465.1d0)then
		write(77,*)(wfv(nonoo,kini),kini=1,2)
	endif
	enddo
	close(77)
c cat
	OPEN(77,FILE=chatoc,STATUS='OLD',ERR=8978)
	CLOSE(77,STATUS='DELETE')
8978	OPEN(77,IOSTAT=IOS,FILE=chatoc,STATUS='NEW')
	do nonoo=1,61
	  write(77,'(A80)')cab(nonoo,3)
	enddo
	do lup=1,710
		write(77,*)(wfc(lup,kini),kini=1,2)
	enddo
	close(77)
c miles
	OPEN(77,FILE=chatom,STATUS='OLD',ERR=8979)
	CLOSE(77,STATUS='DELETE')
8979	OPEN(77,IOSTAT=IOS,FILE=chatom,STATUS='NEW')
	do nonoo=1,61
	  write(77,'(A80)')cab(nonoo,4)
	enddo
	do lup=1,npxm
		write(77,*)(wfm(lup,kini),kini=1,2)
	enddo
	close(77)
c    Medida del indice en el espectro sintetico del CaT
cc	call NICO(wfc,CATFIN)
c	write(*,*)'cat_espectro=',CATFIN
c
c	write(79,'(A5,f10.4)')'MEAS=',Wcatt	
c
c	VMy=-2.5*log10(BETAy*ZZZy(3))
c	ZISOL=(ZMISOo*alpha+ZMISOy*(1.-alpha))/
c     &  (alpha*10**(-.4*(VMo-4.79))+(1.-alpha)*10**(-.4*(VMy-4.79)))
c	VMy=0.
c	ZISOL=0.
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
	call e(alpha,izold(izo),realeo,izyng(izy),realey)
c	write(*,*)alpha,izold(izo),realeo,izyng(izy),realey,quaqum,quaqu0
2325	continue
				    ENDDO !Cerrado loop Tyng
				  ENDDO !Cerrado loop Zyng
				ENDDO !Cerrado loop Told
			ENDDO !Cerrado loop Zold
		ENDDO !Cerrado loop alpha
	ENDDO !Cerrado loop IMF slope
	close(30)
	close(31)
	close(32)
	close(33)
	close(34)
	close(35)
	STOP
6336	END




