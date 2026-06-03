      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION z(12,150000,16),Z00s(15),yter(12),ygir(7),ze(75)
      CHARACTER*1 miuk,shape,lowcut
      INTEGER miukn
      CHARACTER*8 dirout,cisoc
      CHARACTER*80 mkout
      CHARACTER*4 gc !no cambiar que sea de 4, fich sum*
      CHARACTER*7 cmgfe
      CHARACTER*80 jmiku
      INTEGER staout,mkout2,system
      COMMON/labeli/cisoc,cmgfe !aaa,a
      COMMON/aaaa/izold(15),izyng(15),ttold(99),ttyng(99),pizold(15),
     &pizyng(15),zi(150),zalpha(22),nnIMF,nalpha,nzold,ntold,nzyng,
     &ntyng
      COMMON/stus/ZMASA,tiso(75),num(15),numtis    
      COMMON/zsolo/z
      COMMON/NC/ncha !indice tipo isocronas
      COMMON/eststu/ZMUS,ZML,ZMLow,Z00(15)
      COMMON/rat/ra(3,50)
      COMMON/lmu0/bicl,bicp,bich
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
      COMMON/lshape/shape
      COMMON/turnp/xtpoint
      COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3
      COMMON/matina/a(4,4)
      COMMON/tefsca/nghb
      COMMON/fezsol/fez(15)
      COMMON/REMNAN/re(15,2,99),nr(15)
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,vkvega,
     &zerovk
c      COMMON/ypate/yter(12),ygir(7)
      COMMON/KEYS/keypT(18)
      COMMON/TOAGB/g(15,2000,7),ng(15)
      COMMON/TOAGBC/gc(15,2000)
      COMMON/part0/ipartial
      COMMON/part1/pmassL,pmassH
      COMMON/ab/iaba,iabaj
      COMMON/jmik/jmiku(20),lenj(20) !nombre dirs estrellas
      COMMON/enhanc/amgfe,DeltFeH ![Mg/Fe],[Fe/H] correc. value
c para calcular valores [Fe/H]
      DATA ygir/.23d0,.23d0,.23d0,.24d0,.25d0,.273d0,.30d0/
      DATA yter/.245d0,.245d0,.246d0,.246d0,.248d0,.251d0,.256d0,
     &.259d0,.2734d0,.279d0,.288d0,.303d0/
      ipartial=0
      write(*,*)'Compute SSP (0) or partial SSP (1) (default=0)?'
      read(*,*)ipartial
c      write(*,*)'IGNORED, COMPUTING SSPs'
      if(ipartial.eq.1)then
       write(*,*)'Low and High Mass cutoff PartialSSP(def. 0., 1.e9)?'
       read(*,*)pmassL,pmassH
       write(*,*)'Will compute Partial SSPs within',pmassL,pmassH
      else
       ipartial=0
       pmassL=0.
       pmassH=1.e9
       write(*,*)'Will compute SSPs'
      endif
      zxlogG=dlog10(0.019d0/(1.d0-0.019d0-0.273d0))
      zxlogT=dlog10(0.0244d0)
c      COMMON/FNOR/C2Code,FNORMM,FNORMC
c      C2Code=0.958d0 !calculada para BC=-0.07, pero no importa
cccccccccccccc Abundance ratios modelling approach ccccccccccccccccc
      iaba=0
      iabaj=0
      WRITE(*,*)'-----------------------------------------------------'
      WRITE(*,*)'0 Base modeling (default)?'
      WRITE(*,*)'1 [a/Fe] (stars based)?'
      WRITE(*,*)'2 [a/Fe] + correc [Mg/Fe](stars-based)?'
      WRITE(*,*)'3 [a/Fe] + correc [Mg/Fe] Adam (stars-based)?'
      WRITE(*,*)'9 [Mg/Fe] & fully empirical(MILES (SSP-based))?'
      WRITE(*,*)'-----------------------------------------------------'
      READ(*,*)iaba
      DeltFeH=0.0d0
      IF(iaba.eq.0)THEN
       write(*,*)'JONES AND LICK LIBRARY BASED MODELS(no(0),yes(1))?'
       read(*,*)iabaj
      ENDIF
      IF(iaba.eq.1.or.iaba.eq.2)THEN
       write(*,*)'Adopting default base modeling + varying abundance'
       WRITE(*,*)'[Na/Fe](-0.3,0.0,0.3,0.6,0.9,1.2)?'
       READ(*,*)aaNa
       WRITE(*,*)'[a/Fe](-0.2,0.0,0.2,0.4,0.6)?'
       READ(*,*)aaFe
       gamma=0.0d0
c      [Z/H]=[Fe/H]+0.57166d0*aaFe+0.24443d0*aaFe*aaFe
       DeltFeH=0.57166d0*aaFe+0.24443d0*aaFe*aaFe
       amgfe=99.99999999d0
      ELSEIF(iaba.eq.3)THEN
       write(*,*)'Adopting default base modeling + varying [a/Fe]'
       WRITE(*,*)'[a/Fe](-0.2,0.0,0.2,0.4,0.6)?'
       READ(*,*)aaFe
       gamma=0.0d0
       DeltFeH=0.64706d0*aaFe+0.19110d0*aaFe*aaFe
       amgfe=99.99999999d0
       WRITE(*,*)'[C/Fe](-0.25,0.0,+0.25)?'
c       READ(*,*)CFe
       write(*,*)'cancelled, for now [C/Fe]=0.0'
       CFe=0.0d0
       aaNa=0.0d0      
      ELSEIF(iaba.eq.9)THEN
478    FORMAT(A30,2(F7.3,1X))
       WRITE(*,*)'Adopting Vazdekis+15 modeling...'
       WRITE(*,*)'Which [Mg/Fe] value?'
       READ(*,*)amgfe
       ncoap=0
       IF(amgfe.ne.0.0d0)THEN
        write(*,*)'Allende-Prieto(0), Coelho(1), Adam(2) gamma factor?'
	read(*,*)ncoap
        if(ncoap.eq.1.and.amgfe.gt.0.35.and.amgfe.lt.0.45)then
           gamma=0.75d0
           WRITE(*,*)'[M/H]=[Fe/H]+Gamma*[alpha/Fe]'
	   WRITE(*,*)'(Gamma=0.75;[a/Fe]=0.4)=>[Fe/H]=[M/H]-0.3'
	   WRITE(*,*)'Strictly valid for [a/Fe]=0.4'
	   DeltFeH=gamma*amgfe
        elseif(ncoap.eq.0)then
	   DeltFeH=0.57166d0*amgfe+0.24443d0*amgfe*amgfe
        else
           DeltFeH=0.64706d0*amgfe+0.19110d0*amgfe*amgfe
c ADAM [Z/H]=0.64706*[alpha/Fe]+0.19110*[alpha/Fe]*[alpha/Fe]
cOur equation coefficients are different from that of La Barbera because we included Ne in the alpha elements and the alpha abundance range in our models is larger (-0.25 to +0.75).

cWe attach a plot to show the differences in the equations. We plot [Z/H] as a function of [alpha/Fe] with Ne included (black points) using a routine from Carlos. This routine calculates [Z/H] for our range of [alpha/Fe] abundances (for Asplund 2005 solar abundances). We also show our best fit to the functional form of [Z/H]=a*[alpha/Fe]+b*[alpha/Fe]*[alpha/Fe] (Knowles 2018 relation), as well as the relation used in La Barbera et al. 2017 work. Our fit shows good agreement to the [Z/H] calculations in the alpha abundance range of our models.        
        endif
        WRITE(*,478)'[Mg/Fe],Delta[Fe/H]=',amgfe,-DeltFeH
       ENDIF
      ELSE
       write(*,*)'Adopting default base modeling...'
       gamma=0.0d0
       DeltFeH=0.0d0
       amgfe=99.99999999d0
      ENDIF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
5478  write(*,*)'MILES(m),XSL(x),NGSL(u),IRTF(k),ALL(a)?'
      read(*,*)miuk
      IF(miuk.eq.'m'.or.miuk.eq.'M')THEN
       miuk='m'
       dirout='ls ./OUT'
       IF(iaba.eq.0)THEN
        dirout='ls ./OUT'
        if(iabaj.eq.1)then
         write(*,*)'IN ADDITION IT WILL BE COMPUTED: B,R & C'
      mkout='mkdir ./OUT && cd OUT && mkdir B R M m C c'
	else
         write(*,*)'IN ADDITION IT WILL BE COMPUTED: CaT(C) models'
         mkout='mkdir ./OUT && cd OUT && mkdir M m C c'
        endif
       ELSE
        write(*,*)'IN ADDITION IT WILL BE COMPUTED: CaT(C) models'
        mkout='mkdir ./OUT && cd OUT && mkdir M m C c'
       ENDIF
      ELSEIF(miuk.eq.'x'.or.miuk.eq.'X')THEN
       miuk='x'
       dirout='ls ./OUT'
       mkout='mkdir ./OUT && cd OUT && mkdir X x'
      ELSEIF(miuk.eq.'u'.or.miuk.eq.'U')THEN
       miuk='u'
       dirout='ls ./OUT'
       mkout='mkdir ./OUT && cd OUT && mkdir U u'
       ELSEIF(miuk.eq.'k'.or.miuk.eq.'K')THEN
       miuk='k'
       dirout='ls ./OUT'
       mkout='mkdir ./OUT && cd OUT && mkdir K k'
      ELSEIF(miuk.eq.'a'.or.miuk.eq.'A')THEN
       miuk='a'
       dirout='ls ./OUT'
       if(iaba.eq.0.and.iabaj.eq.1)then
      mkout='mkdir ./OUT && cd OUT && mkdir B R M C X U K m c x u k' 
       else
      mkout='mkdir ./OUT && cd OUT && mkdir M C X U K m c x u k'
       endif
       write(*,*)mkout
      ELSE
       goto 5478
      ENDIF
         write(*,*)mkout
c
      IF(miuk.ne.'m'.and.iaba.ge.3)THEN !includes option 9
       write(*,*)'NOT YET POSSIBLE, ONLY AVAILABLE FOR THE MILES RANGE'
       goto 6336
      ENDIF
c Creamos el subdirectorio OUT si no existiese
      staout=system( dirout )
      if(staout.ne.0) then
       write(*,*)'Creating "./OUT" directory to include outputs'
       mkout2=system( mkout )
       write(*,*)'Outputs to be included in "./OUT" subdirectory'	 
      else
c       write(*,*)'Outputs to be included in "./OUT" subdirectory'	 
       write(*,*)'________________________________________'
       write(*,*)'Please, remove or rename ./OUT directory'
       write(*,*)'________________________________________'
       goto 6336	 
      endif
c  si cambio orden arrays en TOAGB produce alignment warning: 
c  el orden debe ser real,natural
      giadif=0.0d0
      write(*,*)'PLEASE VARY TEMPERATURE OF GIANTS (e.g. -100.0d0)'
      write(*,*)'IGNORED, TEMPERATURE OF GIANTS NOT VARIED'
c	read(*,*)giadif
      if(abs(giadif).lt.1.0d0) giadif=0.0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL FNORM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      nm05T=15 !15 lineas estrellas de masa < 0.5 no presentes en isocronas
      do ikeyk=1,18
	 keypT(ikeyk)=keypT(ikeyk)+nm05T
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
c
31963 write(*,'(A)')'The options for the IMF shape are:'
      write(*,'(A)')'UNIMODAL  = u' 
      write(*,'(A)')'BIMODAL   = b' 
      write(*,'(A)')'KROUPA	 = k' 
      write(*,'(A)')'KROUPA BIN= r' 
      write(*,'(A)')'CHABRIER  = c'
      write(*,'(A)')'galactic IMF  = g'
      write(*,'(A)')'FERRERAS  = f' 
      write(*,'(A)')'FERRERAS X= x' 
      write(*,'(A)')'LaBARBERA = l' 
      write(*,'(A)')'u/b/k/r/c/g/f/x/l ?'
      read(*,'(A1)')shape
      if(shape.eq.'b') then
	 write(*,'(A)')'---BIMODAL IMF shape---'
c	 bicl=0.2d0
c	 bicp=0.4d0
c	 bich=0.6d0
c	 bicl=0.3d0
c	 bicp=0.45d0
c	 bich=0.6d0
	 bicl=0.4d0
	 bicp=0.55d0
	 bich=0.7d0
c For Mg
          bicl=0.6d0
          bicp=0.65d0
          bich=0.7d0
c
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
      elseif(shape.eq.'u') then
       write(*,'(A)')'UNIMODAL IMF shape'
      elseif(shape.eq.'k') then
       write(*,'(A)')'KROUPA IMF shape'
      elseif(shape.eq.'r') then
       write(*,'(A)')'KROUPA BINARIES IMF shape (binaries corr.)'
      elseif(shape.eq.'c') then
       write(*,'(A)')'CHABRIER shape'
      elseif(shape.eq.'g') then
       write(*,'(A)')'galactic IMF (depends on [M/H] and SFR)'
      elseif(shape.eq.'f') then
       write(*,'(A)')'FERRERAS(turning point =0.5(limf.f))?'
       read(*,*)xtpoint
       write(*,'(A,F6.3)')'Turning Point=',xtpoint
      elseif(shape.eq.'x') then
       write(*,'(A)')'FERRERAS X-SHAPED(turnpoint=0.5)(limf(xtpoint))'
       write(*,'(A)')'FERRERAS(turnpoint=0.5(limf))?'
       read(*,*)xtpoint
      elseif(shape.eq.'l') then
       write(*,'(A)')'LaBarbera 3segm(m1=.4,m2=.7)limf(xlplow,xlpmed)'
       write(*,'(A)')'xlplow(0.4d0)?,xlpmed(0.7d0)?'
       read(*,*)xlplow,xlpmed
c       xlplow=0.3d0 !still can be modified when asked about
c       xlpmed=0.7d0 !still can be modified when asked about
       write(*,*)xlplow,xlpmed
       ZMU3=1.30d0
       write(*,'(A)')'ADOPTED log.slope for the third segment: 1.3'
       write(*,'(A)')'1st and 2nd segments log.slope(Salpeter=1.3)?'
       read(*,*)ZMU1,ZMU2
      else
       goto 31963
      endif
c      write(*,'(A)')'Set low mass cutoff diff.than 0.1Mo?(y/n)'
      write(*,'(A)')'Set low mass cutoff diff.than 0.1Mo?(y/n);IGNORED'
c      read(*,*)lowcut
      lowcut='n'
      write(*,*)'IGNORED: ADOPTED 0.1Mo'
      if(lowcut.eq.'y'.or.lowcut.eq.'Y')then
         write(*,*)'Lower mass cutoff (Msun)?'
	 read(*,*)ZMlwlw
      endif
      write(*,*)'MILES(0) or GonzalezHernandezBonifacio(1) Teff scale?'
      write(*,*)'IGNORED: ADOPTED MILES SCALE'
c       read(*,*)nghb
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
c	ncha=1
c      write(*,'(A)')'ADOPTED Padova_00      (1)'
      if(ncha.eq.0)then
	 write(*,*)'Teramo isochrones (scaled-solar)'
	 write(cisoc(1:8),'(A8)')'_iTp0.00'
      elseif(ncha.eq.-1)then
	 write(*,*)'Teramo isochrones (alpha-enhanced=0.4)'
	 write(cisoc(1:8),'(A8)')'_iTp0.40'
      elseif(ncha.eq.4)then
	 write(*,*)'Padova scaled-solar isochrones of Salasnich etal'
	 write(cisoc(1:8),'(A8)')'_iSp0.00'
      elseif(ncha.eq.3)then
	 write(*,*)'Padova alpha-enhanced isochrones of Salasnich etal'
	 write(cisoc(1:8),'(A8)')'_iSp0.40'
      elseif(ncha.eq.2)then
	 write(*,*)'OLD Padova isochrones of Bertelli etal'
	 write(cisoc(1:8),'(A8)')'_iBp0.00'
      else
	 write(*,*)'NEW Padova isochrones of Girardi etal'
	 write(cisoc(1:8),'(A8)')'_iPp0.00'
      endif
c Etiqueta por defecto modelos base
      cmgfe='_baseFe'
c**LECTURA DEL FICHERO dat_MASS_HR, dat_MAZT_HR y ratios.dat****************
      OPEN(22,FILE='../D/INPUT/ratios.dat',STATUS='OLD',IOSTAT=IER5)
      READ(22,'(A150)')porque
      READ(22,'(A150)')porque
      do i=3,5
         read(22,*)(ra(i-2,j),j=1,20)
      enddo
      CLOSE(22)
      if(ncha.eq.0)then
          ZML=0.10d0
      elseif(ncha.eq.-1)then
          ZML=0.10d0
      elseif(ncha.eq.2)then
          ZML=0.0992d0
      else
          ZML=0.15d0
      endif
      OPEN(22,FILE='../D/INPUT/dat_MASS_HR',STATUS='OLD',IOSTAT=IER5)
      READ(22,*) ZMASA
      READ(22,*) ZMLow,ZMUS
c	READ(22,*) ZML,ZMUS
      CLOSE(22)
      if(lowcut.eq.'Y'.or.lowcut.eq.'y')then
	 ZMlow=ZMlwlw
	 if(ZMlwlw.gt.ZML) ZML=ZMlwlw
      endif
c      if(iaba.ge.1.and.iaba.le.3)then
       CALL dirs(aaNa,CFe,aaFe)
c      endif
c	write(*,*)ia
c	do inii=1,7
c	 write(*,*)jmiku(inii)
c	enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ISOC TERAMO: SI LA MASA MINIMA NO ES TAN PEQUEÑA, los keypT(se modifican)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ASIGNAR METALICIDADES ISOCRONAS
      IF(ncha.eq.0.or.ncha.eq.-1)THEN !TERAMO
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
	   write(*,*)fez(nft)
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
c         fez(nft)=dlog10(Z00(nft)/0.0200d0)
         fez(nft)=1.024d0*dlog10(Z00(nft))+1.739d0
c         fez(nft)=dlog10(Z00(nft)/0.0190d0)
	enddo
      ELSEIF(ncha.eq.3.or.ncha.eq.4)THEN !Salasnich 00
	Z00(1)=0.0080d0
	Z00s(1)=0.01400d0
        Z00(2)=0.019d0
	Z00s(2)=0.0300d0
	Z00(3)=0.040d0
	do nft=1,3
         fez(nft)=dlog10(Z00(nft)/0.0190d0)
c         fez(nft)=dlog10(Z00(nft)/0.0200d0)
c         fez(nft)=1.024d0*dlog10(Z00(nft)+1.739
	enddo
      ELSE
	write(*,*)'isochrone metalicity and type to be specified'
	goto 6336
      ENDIF
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
         read(22,*,end=322)qmusfr,qalpha,qqqq1,qqqq2,qqqq3,qqqq4
         if(qmusfr.lt.99.999.and.qmusfr.gt.-99.999)then
         	nnIMF=nnIMF+1
                zi(nnIMF)=qmusfr
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
322   CLOSE(22)
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
         zi(nnIMF)=1.30d0 !de todas formas no se usara en el calculo
      ENDIF
c*******LECTURA DE ISOCRONAS***********************************************
      DO 1500 inii=1,nzold
	m=izold(inii)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	IF(ncha.eq.0)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0003T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0006T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0010T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0020T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0040T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0080T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='../D/INPUT/Z0100T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='../D/INPUT/Z0198T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='../D/INPUT/Z0240T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='../D/INPUT/Z0300T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='../D/INPUT/Z0400T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.12)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.-1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0003T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0006T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0010T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0020T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0040T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0080T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='../D/INPUT/Z0100T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='../D/INPUT/Z0198T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='../D/INPUT/Z0240T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='../D/INPUT/Z0300T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='../D/INPUT/Z0400T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.12)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0004_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0010_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0040_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0080_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0190_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0300_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
        ELSEIF(ncha.eq.2)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0004NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z001NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z004NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z008NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z02NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z05NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
        ELSE
	 if(m.eq.1)then !Z=0.008
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z008_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z008_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.2)then !Z=0.019 or Z=0.020
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z019_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z019_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.3)then !Z=0.04
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z040_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z040_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.3)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif		
        ENDIF		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	DO I=1,300000
	 READ(99,*,end=11111)(z(m,I,J),J=1,16)
c esto solo para cambiar teff giants:
         if(giadif.ne.0.0d0)then
	  if(z(m,I,4).lt.3.0d0)then
	   tusto=0.d0
	   tusto=giadif+(10.**(z(m,I,3)))
	   z(m,I,3)=dlog10(tusto)
	  endif
	 endif
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
1500  ENDDO
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
		OPEN(99,FILE='../D/INPUT/Z0001T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0003T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0006T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0010T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0020T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0040T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0080T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='../D/INPUT/Z0100T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='../D/INPUT/Z0198T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='../D/INPUT/Z0240T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='../D/INPUT/Z0300T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='../D/INPUT/Z0400T_ss',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.-1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0003T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0006T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0010T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0020T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0040T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0080T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.8)then
		OPEN(99,FILE='../D/INPUT/Z0100T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.9)then
		OPEN(99,FILE='../D/INPUT/Z0198T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.10)then
		OPEN(99,FILE='../D/INPUT/Z0240T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.11)then
		OPEN(99,FILE='../D/INPUT/Z0300T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.12)then
		OPEN(99,FILE='../D/INPUT/Z0400T_aa',STATUS='OLD')		  
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif
        ELSEIF(ncha.eq.1)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0004_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z0010_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z0040_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z0080_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z0190_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z0300_G',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
	 else
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
        ELSEIF(ncha.eq.2)THEN
	 if(m.eq.1)then
		OPEN(99,FILE='../D/INPUT/Z0001_G',STATUS='OLD')
	 	OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
		OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD') !FALTA sum0001
	 elseif(m.eq.2)then
		OPEN(99,FILE='../D/INPUT/Z0004NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum0004',STATUS='OLD')
	 elseif(m.eq.3)then
		OPEN(99,FILE='../D/INPUT/Z001NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum001',STATUS='OLD')
	 elseif(m.eq.4)then
		OPEN(99,FILE='../D/INPUT/Z004NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_004',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum004',STATUS='OLD')
	 elseif(m.eq.5)then
		OPEN(99,FILE='../D/INPUT/Z008NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
	 elseif(m.eq.6)then
		OPEN(99,FILE='../D/INPUT/Z02NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
	 elseif(m.eq.7)then
		OPEN(99,FILE='../D/INPUT/Z05NNN',STATUS='OLD')
		OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 	OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
	 elseif(m.gt.7)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif	
	ELSE
	 if(m.eq.1)then !Z=0.008
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z008_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z008_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum008',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_008',STATUS='OLD')
	 elseif(m.eq.2)then !Z=0.019 or Z=0.020
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z019_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z019_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum019',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.eq.3)then !Z=0.04
		 if(ncha.eq.3)then
		  OPEN(99,FILE='../D/INPUT/Z040_SA',STATUS='OLD')
		 elseif(ncha.eq.4)then
		  OPEN(99,FILE='../D/INPUT/Z040_SS',STATUS='OLD')
		 endif
		 OPEN(23,FILE='../D/INPUT/sum030',STATUS='OLD')
		 OPEN(24,FILE='../D/INPUT/REMNANTS_019',STATUS='OLD')
	 elseif(m.gt.3)then
	  write(*,*)'Not available such high metallicity isochrones'
		  goto 6336		
	 endif		
	ENDIF		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        DO I=1,300000
	 READ(99,*,end=11112)(z(m,I,J),J=1,16)
c esto solo para cambiar teff giants:
         if(giadif.ne.0.0d0)then
	   if(z(m,I,4).lt.3.0d0)then
	    tusto=0.d0
	    tusto=giadif+(10.**(z(m,I,3)))
	    z(m,I,3)=dlog10(tusto)
	   endif
	 endif
	 num(m)=I
        ENDDO
11112   CLOSE (99)
        if(ncha.gt.0)then
         DO I=1,30000
       	   READ(23,*,end=11106)(g(m,I,J),J=1,7),tin,tin,tin,gc(m,I)
       	   ng(m)=I
         ENDDO
        endif
11106   CLOSE (23)
        DO I=1,300
       	  READ(24,*,end=11105)tin,re(m,1,I),re(m,2,I)
       	  nr(m)=I
        ENDDO
11105   CLOSE (24)
       ENDIF
1501  ENDDO
c********LECTURA DE LOS FICHEROS 'tiempos'**************************
      IF(ncha.eq.0.or.ncha.eq.-1)THEN
        OPEN(99,FILE='../D/INPUT/tiemposT',STATUS='OLD')
      ELSE
        OPEN(99,FILE='../D/INPUT/tiemposP',STATUS='OLD')
      ENDIF
      numtis=0
      DO I=1,30000
	 READ(99,*,end=12345) tiso(I)
	 numtis=I
      ENDDO
12345 CLOSE (99)
cccccccccccccccccccccccccccccc
      IF(miuk.eq.'m'.or.miuk.eq.'a')THEN
	 OPEN(29,FILE='./OUT/out_u',STATUS='OLD',ERR=6329)
	 CLOSE(29,STATUS='DELETE')
6329	 OPEN(29,IOSTAT=IOS,FILE='./OUT/out_u',STATUS='NEW')
         if(iaba.eq.0.and.iabaj.eq.1)then
           OPEN(30,FILE='./OUT/out_li',STATUS='OLD',ERR=6330)
           CLOSE(30,STATUS='DELETE')
6330       OPEN(30,IOSTAT=IOS,FILE='./OUT/out_li',STATUS='NEW')
         endif
	 OPEN(31,FILE='./OUT/out_ef',STATUS='OLD',ERR=6331)
	 CLOSE(31,STATUS='DELETE')
6331	 OPEN(31,IOSTAT=IOS,FILE='./OUT/out_ef',STATUS='NEW')
	 OPEN(32,FILE='./OUT/out_dg',STATUS='OLD',ERR=6332)
	 CLOSE(32,STATUS='DELETE')
6332	 OPEN(32,IOSTAT=IOS,FILE='./OUT/out_dg',STATUS='NEW')
	 OPEN(33,FILE='./OUT/out_sbf',STATUS='OLD',ERR=6333)
	 CLOSE(33,STATUS='DELETE')
6333	 OPEN(33,IOSTAT=IOS,FILE='./OUT/out_sbf',STATUS='NEW')
	 OPEN(34,FILE='./OUT/out_phot',STATUS='OLD',ERR=6334)
	 CLOSE(34,STATUS='DELETE')
6334	 OPEN(34,IOSTAT=IOS,FILE='./OUT/out_phot',STATUS='NEW')
	 OPEN(35,FILE='./OUT/out_contr',STATUS='OLD',ERR=6335)
	 CLOSE(35,STATUS='DELETE')
6335	 OPEN(35,IOSTAT=IOS,FILE='./OUT/out_contr',STATUS='NEW')
	 OPEN(39,FILE='./OUT/out_qua',STATUS='OLD',ERR=6337)
	 CLOSE(39,STATUS='DELETE')
6337	 OPEN(39,IOSTAT=IOS,FILE='./OUT/out_qua',STATUS='NEW')
	 OPEN(40,FILE='./OUT/out_mass',STATUS='OLD',ERR=6339)
	 CLOSE(40,STATUS='DELETE')
6339	 OPEN(40,IOSTAT=IOS,FILE='./OUT/out_mass',STATUS='NEW')
c
         CALL a_alpha
c
	 CLOSE(29)
         if(iaba.eq.0.and.iabaj.eq.1) CLOSE(30)
	 CLOSE(31)
	 CLOSE(32)
	 CLOSE(33)
	 CLOSE(34)
	 CLOSE(35)
	 CLOSE(40)
      ENDIF
      IF(miuk.eq.'x'.or.miuk.eq.'a')THEN
cc	 OPEN(49,FILE='./OUT/out_u_X',STATUS='OLD',ERR=6648)
cc	 CLOSE(49,STATUS='DELETE')
cc6648	 OPEN(49,IOSTAT=IOS,FILE='./OUT/out_u_X',STATUS='NEW')
	 OPEN(50,FILE='./OUT/out_ef_X',STATUS='OLD',ERR=6649)
	 CLOSE(50,STATUS='DELETE')
6649	 OPEN(50,IOSTAT=IOS,FILE='./OUT/out_ef_X',STATUS='NEW')
	 OPEN(399,FILE='./OUT/out_qua_X',STATUS='OLD',ERR=6650)
	 CLOSE(399,STATUS='DELETE')
6650	 OPEN(399,IOSTAT=IOS,FILE='./OUT/out_qua_X',STATUS='NEW')
c
         CALL a_XSL
c
cc	 CLOSE(49)
	 CLOSE(50)
      ENDIF
      IF(miuk.eq.'u'.or.miuk.eq.'a')THEN
	 OPEN(51,FILE='./OUT/out_ef_U',STATUS='OLD',ERR=6450)
	 CLOSE(51,STATUS='DELETE')
6450	 OPEN(51,IOSTAT=IOS,FILE='./OUT/out_ef_U',STATUS='NEW')
	 OPEN(3999,FILE='./OUT/out_qua_U',STATUS='OLD',ERR=6451)
	 CLOSE(3999,STATUS='DELETE')
6451	 OPEN(3999,IOSTAT=IOS,FILE='./OUT/out_qua_U',STATUS='NEW')
c
	 CALL a_NGSL
c
	 CLOSE(51)
      ENDIF
      IF(miuk.eq.'k'.or.miuk.eq.'a')THEN
	 OPEN(52,FILE='./OUT/out_ef_K',STATUS='OLD',ERR=6551)
	 CLOSE(52,STATUS='DELETE')
6551	 OPEN(52,IOSTAT=IOS,FILE='./OUT/out_ef_K',STATUS='NEW')
	 OPEN(39999,FILE='./OUT/out_qua_K',STATUS='OLD',ERR=6552)
	 CLOSE(39999,STATUS='DELETE')
6552	 OPEN(39999,IOSTAT=IOS,FILE='./OUT/out_qua_K',STATUS='NEW')
c
	 CALL a_IRTF
c
	 CLOSE(52)
      ENDIF
4785  FORMAT(A8,F6.3,1X,A17,F6.3)
      write(*,*)'================================================='
      write(*,'(A)')'PLEASE BE AWARE THAT:      '
c      WRITE(*,478)'[Mg/Fe],Delta[Fe/H]=',amgfe,-DeltFeH
      if(iaba.eq.0)then
       amgfe=0.0d0
       DeltFeH=0.0d0
      endif
      write(*,4785)'[Mg/Fe]=',amgfe,' ; [Fe/H]=[M/H] -',DeltFeH
      write(*,*)'================================================='
c478    FORMAT(A30,2(F7.3,1X))
6336  END	
