c Determinacion constantes escalado espectros
c Llamada desde a.f; 
c COMMON FNOR: a.f,hrsl.f,STU.f,xflabs.f
      SUBROUTINE FNORM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (npick2=2299) !alphaLyrae sed
c      PARAMETER (npick2=1221) !A0V_KURUCZ_92.SED
      DIMENSION aaa(npick2,2),aaa1(npick2),aaa2(npick2)
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Para calcular M/L_3.6 ver explicacion STU_alpha.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Flujo referencia en 5556A (Fukugita+95,PASP,107,945) (erg/cm2.s.A)
      F5556A=3.44e-09 
c Luminosidad Sol (erg/s):
      xxLsun=3.826e33
c Flujo sol en superficie Tierra desde 10pc (erg/s.A)
      pi=3.1415926536d0
      fsol=xxLsun/(4.0d0*pi*(10.0d0*3.0856e18)**2.d0)
c      fsol=3.1978335158079E-7
c limites lambdas filtros U V I
      UBUS1=3050.0d0
      UBUS2=4200.0d0
      V1=4750.d0
      V2=7400.d0
c      XK1=19400.0d0 !IRTF
c      XK2=24800.0d0 !IRTF
c      XK1=18500.0d0 !IRTF lambda inicial filtro K donde supera 0.
c	XK2=25500.0d0 !IRTF lambda final filtro K donde supera 0.
      XK1=18000.0d0 !IRTF lambda inicial filtro K donde supera 0.
      XK2=26000.0d0 !IRTF lambda final filtro K donde supera 0.
c      XIJ1=6800.0d0 
c      XIJ2=12000.0d0
      XIC1=7000.0d0 
      XIC2=9300.0d0
c
      bcsolV=-0.12d0
      bolsol=4.70d0
      Tsol=5770.d0
      gsol=4.4d0
      fesol=0.0d0
      VsolV=4.82d0
      UVsol=0.78d0
      BVsol=0.621d0
      VRsolc=0.361d0
      VIsolc=0.672d0
      VJsol=1.109d0
      VHsol=1.428d0
      VKsol=1.486d0
c  ESCALAMOS VEGA SI FUESE NECESARIO
      OPEN(32,file='../D/INPUT/vega_sed/alpha_lyr_stis_008.SED',status='OLD')
c      OPEN(32,file='INPUT/vega_sed/alpha_lyr_005.sed',status='OLD')
c      OPEN(32,file='INPUT/vega_sed/A0V_KURUCZ_92.SED',status='OLD')
      do l=1,npick2
       read(32,*)aaa(l,1),aaa(l,2)
       aaa1(l)=aaa(l,1)
       aaa2(l)=aaa(l,2)
      enddo
      CLOSE(32)
      call hunt(aaa1,npick2,5556.0d0,k1)
      call polint(aaa1(k1-1),aaa2(k1-1),3,5556.0d0,rpol,dy)
c dividimos por flujo(5556A) y escalamos a F5556A.
      zerove=F5556A/rpol
c      zerove=1.0d0
c DETERMINAMOS CONSTANTES PARA CALCULAR FLUJO ABSOLUTO
      fvl=0.0d0
      fil=0.0d0
      ful=0.0d0
      fkl=0.0d0 !IRTF
      rfV=0.0d0
      rfi=0.0d0
      rfU=0.0d0
      rfK=0.0d0 !IRTF
      DO l=1,npick2
       aaa(l,2)=aaa(l,2)*zerove
       aaa2(l)=aaa(l,2)
       if(l.eq.npick2)then
        dl=abs(aaa(l-1,1)-aaa(l,1))
       else
        dl=abs(aaa(l+1,1)-aaa(l,1))
       endif
       rfV=0.d0
       rfi=0.d0
       rfU=0.d0
       rfK=0.d0
       if(aaa(l,1).ge.V1.and.aaa(l,1).le.V2)then
        call respv(aaa(l,1),rfV)
        fvl=fvl+rfV*aaa(l,2)*dl
       endif
       if(aaa(l,1).ge.XIC1.and.aaa(l,1).le.XIC2)then
c	if(aaa(l,1).ge.XIJ1.and.aaa(l,1).le.XIJ2)then
c	 call respij(aaa(l,1),rfi)
        call respic(aaa(l,1),rfi)
        fil=fil+rfi*aaa(l,2)*dl
       endif
       if(aaa(l,1).ge.UBUS1.and.aaa(l,1).le.UBUS2)then
        call respu(aaa(l,1),rfU)
        ful=ful+rfU*aaa(l,2)*dl
       endif
       if(aaa(l,1).ge.XK1.and.aaa(l,1).le.XK2)then
        call respk(aaa(l,1),rfK)
        fkl=fkl+rfK*aaa(l,2)*dl
       endif
c	 fvl=fvl+rfV*aaa(l,2)*dl
c	 fil=fil+rfi*aaa(l,2)*dl
c	 ful=ful+rfU*aaa(l,2)*dl
c	 fkl=fkl+rfK*aaa(l,2)*dl
      ENDDO
      zerovv=-2.5d0*dlog10(fvl)
      zerovi=-2.5d0*dlog10(fil)
      zerovu=-2.5d0*dlog10(ful)
      zerovk=-2.5d0*dlog10(fkl)
c      write(*,*)'ZP_V,ZP_I,ZP_U,ZP_K',zerovu,zerovv,zerovi,zerovk
c      zerovv=13.762802d0 !Vega Hayes
c      zerovi=14.383596d0 !Vega Hayes
c      zerovu=14.119219d0 !Vega alphaLyr but very similar to Hayes.
c Vega: alpha_lyr_stis_008.SED
c The directory CURRENT_CALSPEC includes all the latest HST Calibration
c Spectra, listed in the webpage 
c http://www.stsci.edu/hst/observatory/cdbs/calspec.html
      zerovu=14.116818875 !U_Buser
c            13.0117035057 B2_Buser
c            13.0390112888 B3_Buser
      zerovv=13.7588072654 !V_Buser
c            13.66441674 R_Cousins
      zerovi=14.3768163283 !I_Cousins
c            14.9779199337 J_Johnson
c            16.2030103398 H_Bessell
      zerovk=16.5257323033 !K_Johnson
c            18.3461782842 [3.6]_Spitzer
c            19.0853423617 [4.5]_Spitzer
c      vijvega=-0.040d0 !adoptado tras varias transf. y mas proximo a Makeiso.f
c      vkvega=-0.075d0 !Alonso+95(A&A,297,197 (Table 2) COLORS.f ?
c      vkvega=0.011d0 !Alonso+95(A&A,297,197 (Table 2) !AFINAR
c      vkvega=0.059d0 !Alonso+95(A&A,297,197 (Table 2) filtros TCS
c      write(*,*)'ADOPTED VALUES     ',zerovu,zerovv,zerovi,zerovk
      uvvega=0.0d0 !Alonso+95(A&A,297,197 (Table 2)
      vijvega=-0.00774500d0 !Tabla2Alonso+95->Eq.1Alonso+95->Transf.Cousin Fernie
      vkvega=0.01d0 !Alonso+95(A&A,297,197 (Table 2) (me sale: V-Kj=.009054)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Valores que se deducen de las eqs. de Alonso+94 y Tabla 2 de Alonso+95
c Tab.2:U=B=V=0.03,V-Rj=-0.04,V-Ij=-0.07,V-Jtcs=0.043,V-Htcs=0.035,V-Ktcs=0.059
c=>V-Jj=.014471,V-Hj=.001159,V-Kj=.009054,Jj-Kj=.0066,Jj-Hj=-.005417,Hj-Kj=.007895
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RETURN
      END
c
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
cElena: 
czerov=-2.5*log( int(Fvega*Tv) )
cEstoy normalizando a F(5556)=3.44e-9
c
c Calculado (alphaLyr) V,Ic,K   13.761882276840447 14.382294627813838 16.573655756543499
c
c Adoptado (Hayes) V,Ic,K  13.762802000000001 14.383596000000001 16.573655756543499 
c-------------------------------------------------------------------------------
c Elena's colours (Johnson) of Vega computed with the spectrum:
c './INPUT/vega_sed/alpha_lyr_005.sed'
c V-I=-0.02248
c V-R=-0.004
c R-I=-0.01848       
c      vijvega=-0.039419996757721d0 !from MAKEISO.f with inv. transf. Fernie        
c      vijvega=-0.010291002397314d0
c      vijvega=-0.0365d0
c      vijvega=-0.035633545779014d0
c      vijvega=-0.02248d0
c      vijvega=0.d0
c vijvega I Cousin:
c      vijvega=0.006d0
c      vijvega=-0.025473d0
c      vijvega=-0.04455d0 !(Alonso+ transf.inv. Fernie)
c
c      FNORMM=10.**(-0.4d0*(VsolV-(bolsol-C2Code)))!MILES,V99b,V99r SSPs
c      solvfo=10.**(0.4d0*(-0.12d0-C2Code))     
c      FNORMM=10.**(0.4d0*(bcsolV-C2Code))!MILES,V99b,V99r SSPs
c      RIsolc=VIsolc-VRsolc
c      VIsolj=-0.005d0+1.273d0*VIsolc
c      fnCaT=0.198935997d0-0.00599133915d0*dlog10(Tsol)!CaT range fraction of Ij
c      FNORMC=(10.**(0.4d0*VIsolj))*FNORMM*fnCaT!CaT SSPs
c
c	open(32,file='g2v.ascii',status='old')
c	do l=1,npick
c	 read(32,*)aaa(l,1),aaa(l,2)
c	 aaa(l,1)=aaa(l,1)*10.0d0
c	enddo
c	close(32)
c	 fcat=0.0d0
c	 fcat1=0.0d0
c	 fij=0.0d0
c	 rfC=0.0d0
c	 do l=1,npick
c	  if(l.eq.npick)then
c	   dl=abs(aaa(l-1,1)-aaa(l,1))
c	  else
c	   dl=abs(aaa(l+1,1)-aaa(l,1))
c	  endif
c	  call respij(aaa(l,1),rfC)
c	  fij=fij+rfC*aaa(l,2)*dl
c	  if(aaa(l,1).ge.C1.and.aaa(l,1).le.C2)then
c	    fcat=fcat+rfC*aaa(l,2)*dl
c	    fcat1=fcat1+aaa(l,2)*dl
c	  endif
c	 enddo
c	 write(*,*)'COCIENTECAT',fcat/fij,fcat1/fij
c	 open(32,file='m0iii.ascii',status='old')
c	 do l=1,npick
c	  read(32,*)aaa(l,1),aaa(l,2)
c	  aaa(l,1)=aaa(l,1)*10.0d0
c	 enddo
c	 close(32)
c	 fcat=0.0d0
cc	 fcat1=0.0d0
c	 fij=0.0d0
c	 rfC=0.0d0
c	 do l=1,npick
c	  if(l.eq.npick)then
c	   dl=abs(aaa(l-1,1)-aaa(l,1))
c	  else
c	   dl=abs(aaa(l+1,1)-aaa(l,1))
c	  endif
c	  call respij(aaa(l,1),rfC)
c	  fij=fij+rfC*aaa(l,2)*dl
c	  if(aaa(l,1).ge.C1.and.aaa(l,1).le.C2)then
c	    fcat=fcat+rfC*aaa(l,2)*dl
cc	    fcat1=fcat1+aaa(l,2)*dl
c	  endif
c	 enddo
c	 write(*,*)'COCIENTECAT',fcat/fij,fcat1/fij
