c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c*******SUBRUTINA CAT + Near-IR 888********************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c Esto es lo que voy a llamar yo:
c         call nIR(tetas,ggrav,zfeff,ZCATS,ecats,ZPAT,epat,ZCAT,
c     +   ecat,ZSTIO,estio,ZMGI,emgi,iflag,iflag2,iflag3,iflag4)
c
      subroutine nIR(t,g,z,findex_cats,eindex_cats,
     +                           findex_pat,eindex_pat,
     +                           findex_cat,eindex_cat,
     +                           findex_stio,eindex_stio,
     +                           findex_mgi,eindex_mgi,
     +                           iflag,iflag2,iflag3,iflag4)
c Version 09/11/01
C------------------------------------------------------------------------------
c This program predicts the indices <CaT*>, <CaT> and <PaT> using the empirical 
c fitting functions from Cenarro et al. (Paper III, 2000, MNRAS, accepted),
c as well as the indices sTiO (continuum slope) and MgI (Mg8807) using fitting
c functions which will be published in a forthcoming paper... (I hope so!) 
C-----------------------------------------------------------------------------
c INPUT:
c       t = effective temperature, in K                      (REAL)
c       g = logarithm (base 10) of the surface gravity       (REAL)
c       z = metallicity [Fe/H]                               (REAL)
c  (A value higher or equal to 99 in g or z in input means that this 
c    parameter is unknown. In same cases, an index value can still be 
c    computed) 
c  
c
c OUTPUT:
c       findex_cats = CaT* index value                                 (REAL)
c       eindex_cats = CaT* error                                       (REAL)
c       findex_pat = PaT index value                                   (REAL)
c       eindex_pat = PaT error                                         (REAL)
c       findex_cat = CaT index value                                   (REAL)
c       eindex_cat = CaT error                                         (REAL)
c       findex_cat = sTiO index value                                  (REAL)
c       eindex_cat = sTiO error (flux cal. errors not included!!)      (REAL)
c       findex_cat = MgI index value                                   (REAL)
c       eindex_cat = MgI error (flux cal. errors not included!!)       (REAL)
c       iflag = indicates sucess of the functions:                  (INTEGER)
c              0 -> OK
c              1 -> extrapolation to high temperatures
c              2 -> extrapolation to low temperatures
c              3 -> extrapolation in g (uncertain value)
c             -1 -> Error (no Teff)
c             -2 -> Error (z needed)
c             -3 -> Error (g needed)
c             -5 -> Error
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer iflag,iflag2,iflag3,iflag4
      double precision t,g,z,t0,cdef
      double precision findex_cats,eindex_cats
      double precision findex_pat,eindex_pat
      double precision findex_cat,eindex_cat
      double precision findex_stio,eindex_stio
      double precision findex_mgi,eindex_mgi
 
c LA VARIABLE TETAS YA ESTA PASADA A 5040/t, ASI QUE COMENTAMOS ESTA LINEA
c      t0=5040./t !Converts to theta=5040/Teff 
	t0=t

      call finfitCatpat(t0,g,z,findex_cats,eindex_cats,iflag)
      call finfitpat(t0,g,z,findex_pat,eindex_pat,iflag2)
      cdef=0.93d0
      findex_cat=findex_cats+cdef*findex_pat
      eindex_cat=sqrt(eindex_cats**2+(cdef*eindex_pat)**2)
      call finfitsTiO(t0,g,z,findex_stio,eindex_stio,iflag3)
      call finfitMgI(t0,g,z,findex_mgi,eindex_mgi,iflag4)

      return
      end

c-----------------------------------------------------------------------------------------------
c ################## INTERPOLACION DE FUNCIONES Y COEFICIENTES: <CaT*> #########################
c-----------------------------------------------------------------------------------------------

      subroutine finfitCatpat(t,g,z,findex,eindex,iflag)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
      integer iflag
      double precision t,g,z,findex,eindex,pi
      
      double precision findex1,findex2,findex3
      double precision eindex1,eindex2,eindex3
      double precision thet,x
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      double precision theta,geta,zeta
      double precision xf(20)
      logical nog,noz
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      common/blkx/xf

C Regions:
C 1 = hot dwarfs (hd)
C 2 = hot giants (hg)
C 3 = warm stars (w)
C 4 = cool stars (c)
C 5 = cold dwarfs (cd)
C 6 = cold giants (cg)
C 7 = intermediate giants (i)  
C 8 = metal poor warm giants (lm)
c Coefficients and variance-covariance matrices are stored in this subroutine
      call readcoefcatpat

      pi=3.14159d0

      nog=.false.
      noz=.false.
      if(g.ge.99.d0) nog=.true.
      if(z.ge.99.d0) noz=.true.

      iflag=0
      findex=0.d0
      eindex=0.d0
      if(t.le.0.d0) then
         iflag=-1
         return
      end if
      if(nog) then
         iflag=-3
         return
      end if

c Program works in theta (=5040/Teff)
      theta=dble(t)
      thet=t
      geta=dble(g)
      zeta=dble(z)
      xf(1)=1.D0
      xf(2)=theta
      xf(3)=geta
      xf(4)=zeta
      xf(5)=theta*zeta
      xf(6)=theta*geta
      xf(7)=theta*theta
      xf(8)=geta*geta
      xf(9)=zeta*zeta
      xf(10)=geta*zeta
      xf(11)=theta*theta*geta
      xf(12)=theta*theta*theta
      xf(13)=geta*geta*geta
      xf(14)=zeta*zeta*zeta
      xf(15)=geta*geta*zeta
      xf(16)=theta*theta*zeta
      xf(17)=theta*zeta*zeta
      xf(18)=theta*geta*geta
      xf(19)=geta*zeta*zeta
      xf(20)=theta*geta*zeta
 
c Hot stars
      if(thet.le.0.70d0) then
         if(thet.lt.0.13d0) iflag=1
         if(g.ge.2.80d0) then
            call compindex(1,findex1,eindex1)
         end if
         if(g.le.3.0d0) then
            call compindex(2,findex2,eindex2)
         end if
         if(g.le.2.80d0) then
            findex3=findex2
            eindex3=eindex2
         else
            if(g.ge.3.0d0) then
               findex3=findex1
               eindex3=eindex1
            else
               x=(g-2.80d0)/0.20d0
               findex3=(1.d0-x)*findex2+x*findex1
               eindex3=(1.d0-x)*eindex2+x*eindex1
            end if
         end if

c If really hot or no z values, it keeps this index

         if(thet.le.0.50d0.or.noz) then
            findex=findex3
            eindex=eindex3
            if(g.lt.1.d0.or.g.gt.4.2d0) iflag=3
            return
         end if
      end if

c Warm stars

c      if(thet.le.0.90) then
c         if(noz) then
c            iflag=-2
c            return
c         end if
c         call compindex(3,findex1,eindex1)
c         if(thet.le.0.70) then
c            x=cos(pi/2.*(0.70-thet)/0.20)
c            findex=(1.-x)*findex3+x*findex1
c            eindex=(1.-x)*eindex3+x*eindex1
c            if(g.lt.0.5.or.g.gt.4.6) iflag=3
c            return
c         end if
c      end if
      
      if(thet.le.1.00d0) then
         if(noz) then
            iflag=-2
            return
         end if
         call compindex(3,findex1,eindex1)
         call compindex(8,findex2,eindex2) ! solo lo usaremos para gigantes poco metalicas.
         if(g.lt.2.6d0.and.z.lt.-0.25d0)then !todas menos las rojas y naranjas (quedan igual). OJO, 2.6 por HB
            if(z.lt.-0.25d0.and.thet.le.0.65d0) then !azules, verdes y amarillas hasta theta=0.65
               x=cos(pi/2.*(0.65d0-thet)/0.15d0)
               findex=(1.d0-x)*findex3+x*findex2
               eindex=(1.d0-x)*eindex3+x*eindex2
               if(g.lt.0.5d0) iflag=3
               return
            else
               findex1=findex2
               eindex1=eindex2
            endif   
         else !gigantes naranjas y rojas y todas las enanas 
            if(thet.le.0.70d0) then   
               x=cos(pi/2.*(0.70d0-thet)/0.20d0)
               findex=(1.d0-x)*findex3+x*findex1
               eindex=(1.d0-x)*eindex3+x*eindex1
               if(g.gt.4.6d0) iflag=3
               return
            end if
         endif   
      end if
      
c Cool stars
      if(thet.lt.1.50d0) then
         if(noz) then
            if((g.ge.3.d0.and.thet.lt.1.06d0).or.
     c         (g.lt.3.d0.and.thet.lt.1.3d0)) then
              iflag=-2
              return
            else
               goto 10
            end if
         end if
         call compindex(4,findex2,eindex2)
         call compindex(7,findex3,eindex3) !intermediate zone, only giants, not dwarfs

         if(g.lt.2.6d0)then !solo para gigantes OJO 2.6 para excluir las HB
            if(z.lt.-1.25d0)then ! azules desde theta 0.65 a 0.95
               if(thet.le.0.95d0) then
                  x=cos(pi/2.d0*(thet-0.65d0)/0.30d0)
                  findex=(1.d0-x)*findex2+x*findex1
                  eindex=(1.d0-x)*eindex2+x*eindex1
                  if(g.lt.0.d0) iflag=3
                  return
               endif   
            endif   
            if(z.ge.-1.25d0.and.z.lt.-0.75d0)then ! verdes desde theta 0.65 a 0.90
               if(thet.le.0.90d0) then
                  x=cos(pi/2.*(thet-0.65d0)/0.25d0) 
                  findex=(1.d0-x)*findex2+x*findex1
                  eindex=(1.d0-x)*eindex2+x*eindex1
                  if(g.lt.0.d0) iflag=3
                  return
               endif   
            endif   
            if(z.ge.-1.25d0.and.z.lt.-0.25d0)then ! amarillas desde theta 0.65 a 0.85
               if(thet.le.0.85d0) then
                  x=cos(pi/2.d0*(thet-0.65d0)/0.20d0) 
                  findex=(1.d0-x)*findex2+x*findex1
                  eindex=(1.d0-x)*eindex2+x*eindex1
                  if(g.lt.0.d0) iflag=3
                  return
               endif   
            endif   
            if(z.ge.-0.25d0) then !naranjas y rojas desde 0.70, como siempre
               if(thet.le.0.90d0) then
                  x=cos(pi/2.d0*(thet-0.70d0)/0.20d0)
                  findex=(1.d0-x)*findex2+x*findex1
                  eindex=(1.d0-x)*eindex2+x*eindex1
                  if(g.lt.0.d0) iflag=3
                  return
               endif   
            endif   
         else !enanas   
            if(thet.le.0.90d0) then
               x=cos(pi/2.d0*(thet-0.70d0)/0.20d0)
               findex=(1.d0-x)*findex2+x*findex1
               eindex=(1.d0-x)*eindex2+x*eindex1
               if(g.lt.0.d0.or.g.gt.5.d0) iflag=3
               return
            endif   
         end if
         if(g.ge.3.0d0) then 
            if(thet.lt.1.0d0)then
               findex=findex2
               eindex=eindex2
               return
            endif   
         else 
            if(thet.lt.1.20d0) then
               if(thet.ge.1.05d0) then
                  x=cos(pi/2.d0*(thet-1.05d0)/0.15d0)
                  findex=(1.d0-x)*findex3+x*findex2
                  eindex=(1.d0-x)*eindex3+x*eindex2
                  if(g.lt.0.d0.or.g.gt.5.d0) iflag=3
                  return
               else   
                  findex=findex2
                  eindex=eindex2
                  return
               endif   
            endif  
         endif   
      endif   
c         if(thet.ge.1..and.g.ge.3.) then
c            findex3=findex2
c            eindex3=eindex2
c         else
c            if(thet.le.1.2) then
c              findex=findex2
c              eindex=eindex2
c              if(g.lt.-0.1.or.g.gt.5.) iflag=3
c              return
c            else
c              findex3=findex2
c              eindex3=eindex2
c            end if
c         end if
c      end if

 10   if(g.ge.3.0) then 
         if(g.lt.4.4d0.or.g.gt.5.2d0) iflag=3
         call compindex(5,findex1,eindex1)
         if(thet.ge.1.25d0.or.noz) then
            findex=findex1
            eindex=eindex1
            if(thet.gt.1.89d0) iflag=2
            return
         else
            x=cos(pi/2.*(thet-1.0d0)/0.25d0)
            findex=(1.d0-x)*findex1+x*findex2
            eindex=(1.d0-x)*eindex1+x*eindex2
            return
         end if 
      else
         if(g.gt.1.8d0.or.g.lt.-0.1d0) iflag=3
         call compindex(6,findex1,eindex1)
         if(thet.ge.1.50d0.or.noz) then
            findex=findex1
            eindex=eindex1
            if(thet.gt.1.74d0) iflag=2
            return
         else
            x=cos(pi/2.d0*(thet-1.20d0)/0.30d0)
            findex=(1.d0-x)*findex1+x*findex3
            eindex=(1.d0-x)*eindex1+x*eindex3
            return
         end if
      end if

c Error: No value was computed
      iflag=-5
c      write(*,*)t,g,z,findex
c      write(*,'(A)') 'ERROR: NO index value was computed "finfitCatpat"'
c      stop

      end

c-------------------------------------------------------------------------------------------
      
      subroutine readcoefcatpat
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer i,j
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm

c-Initialization----------------------------------------------------------
      do i=1,25
         fhd(i)=0.D0
         fhg(i)=0.D0
         fw(i)=0.D0
         fc(i)=0.D0
         fcd(i)=0.D0
         fcg(i)=0.D0
         fi(i)=0.D0
         flm(i)=0.D0
         lhd(i)=.false.
         lhg(i)=.false.
         lw(i)=.false.
         lc(i)=.false.
         lcd(i)=.false.
         lcg(i)=.false.
         li(i)=.false.
         llm(i)=.false.
         do j=1,25
            chd(i,j)=0.D0
            chg(i,j)=0.D0
            cw(i,j)=0.D0
            cc(i,j)=0.D0
            ccd(i,j)=0.D0
            ccg(i,j)=0.D0
            ci(i,j)=0.D0
            clm(i,j)=0.D0
         end do
      end do
      lehd=.false.
      lehg=.false.
      lew=.false.
      lec=.false.
      lecd=.false.
      lecg=.false.
      lei=.false.
      lelm=.false.
      ctehd=0.D0
      ctehg=0.D0
      ctew=0.D0
      ctec=0.D0
      ctecd=0.D0
      ctecg=0.D0
      ctei=0.D0
      ctelm=0.D0


c For hot dwarfs
c N=44 sigma_res=0.64117 sigma_typ=0.19320 r**2=0.852
Coef  1 :         0.1010234840772093E+01        0.3328633433839919E+00
Coef  7 :        -0.4079786758081508E+02        0.5242578829816516E+01
Coef 12 :         0.6712183060874379E+02        0.7292626079044976E+01

      fhd(1)=  0.1010234840772093E+01
      fhd(7)= -0.4079786758081508E+02
      fhd(12)= 0.6712183060874379E+02
      lhd(1)=.true.
      lhd(7)=.true.
      lhd(12)=.true.
      srhd=0.6411670277E+00
      chd( 1, 1)=  0.2695191207614669E+00
      chd( 1, 2)= -0.3472429151576591E+01
      chd( 1, 3)=  0.4398603795783181E+01
      chd( 2, 2)=  0.6685710666456619E+02
      chd( 2, 3)= -0.9210177382583753E+02
      chd( 3, 3)=  0.1293676030302836E+03

c For hot giants
c N=26 sigma_res=0.61231 sigma_typ=0.14823 r**2=0.925
Coef  1 :        -0.4492070013738538E+00        0.2930124101618676E+00
Coef 12 :         0.2046337453401181E+02        0.1642371038782907E+01

      fhg(1)= -0.4492070013738538E+00
      fhg(12)= 0.2046337453401181E+02
      lhg(1)=.true.
      lhg(12)=.true.
      srhg=0.6123062197E+00
      chg( 1, 1)=  0.2289995808942432E+00
      chg( 1, 2)= -0.1073025401417251E+01
      chg( 2, 2)=  0.7194576162110260E+01

c For warm stars
c N=193 sigma_res=0.64901 sigma_typ=0.19186 r**2=0.936
Coef  1 :        -0.2712443296702683E+02        0.5618020964894731E+01
Coef  3 :         0.1827325591342781E+02        0.5694843876505995E+01
Coef  5 :         0.3607983477791699E+02        0.6158257142615004E+01
Coef  6 :        -0.2944080627123500E+02        0.7250298365863318E+01
Coef  7 :         0.1277746531094536E+03        0.1720088801471742E+02
Coef  8 :        -0.2931841909454243E+01        0.9825535571598784E+00
Coef  9 :         0.3611452207610669E+01        0.1258214871917101E+01
Coef 10 :        -0.5038345594941118E+01        0.1050915869665071E+01
Coef 12 :        -0.7318048016203497E+02        0.1086107084405936E+02
Coef 15 :         0.4445755626664886E+00        0.1701235786382027E+00
Coef 16 :        -0.2050782869070491E+02        0.4693153012041808E+01
Coef 18 :         0.4405810536616976E+01        0.1239147296851348E+01
Coef 19 :        -0.8167748381113530E+00        0.3178820717500236E+00

      fw(1)= -0.2712443296702683E+02
      fw(3)=  0.1827325591342781E+02
      fw(5)=  0.3607983477791699E+02
      fw(6)= -0.2944080627123500E+02
      fw(7)=  0.1277746531094536E+03
      fw(8)= -0.2931841909454243E+01
      fw(9)=  0.3611452207610669E+01
      fw(10)=-0.5038345594941118E+01
      fw(12)=-0.7318048016203497E+02
      fw(15)= 0.4445755626664886E+00
      fw(16)=-0.2050782869070491E+02
      fw(18)= 0.4405810536616976E+01
      fw(19)=-0.8167748381113530E+00
      lw(1)=.true. 
c      lw(2)=.true.
      lw(3)=.true.
      lw(5)=.true.
      lw(6)=.true.
      lw(7)=.true.
      lw(8)=.true.
      lw(9)=.true.
      lw(10)=.true.
      lw(12)=.true. 
      lw(15)=.true. 
      lw(16)=.true. 
      lw(18)=.true. 
      lw(19)=.true. 
      srw=0.6490083081E+00
      cw( 1, 1)=  0.7493180572153959E+02
      cw( 1, 2)= -0.7287507911189091E+02
      cw( 1, 3)= -0.1984589195722098E+02
      cw( 1, 4)=  0.9121968185807818E+02
      cw( 1, 5)= -0.2160595724326513E+03
      cw( 1, 6)=  0.1226696365087670E+02
      cw( 1, 7)= -0.1660015577629536E+01
      cw( 1, 8)=  0.1891251509959769E+01
      cw( 1, 9)=  0.1221490452723156E+03
      cw( 1,10)= -0.9943263591572606E-01
      cw( 1,11)=  0.1567638580172093E+02
      cw( 1,12)= -0.1531442505240781E+02
      cw( 1,13)=  0.4322638527511040E+00
      cw( 2, 2)=  0.7699510795787947E+02
      cw( 2, 3)=  0.1684546212631264E+02
      cw( 2, 4)= -0.9650706289802018E+02
      cw( 2, 5)=  0.1998828461196264E+03
      cw( 2, 6)= -0.1316764851971396E+02
      cw( 2, 7)=  0.9412684841106710E+00
      cw( 2, 8)= -0.1403854896373619E+01
      cw( 2, 9)= -0.1059654471717934E+03
      cw( 2,10)=  0.6910988317110785E-01
      cw( 2,11)= -0.1384864724298732E+02
      cw( 2,12)=  0.1646220576720622E+02
      cw( 2,13)= -0.2251612482713273E+00
      cw( 3, 3)=  0.9003577886546677E+02
      cw( 3, 4)= -0.2439501848701901E+02
      cw( 3, 5)=  0.7199671862668183E+02
      cw( 3, 6)= -0.2703391024164132E+01
      cw( 3, 7)=  0.1339286752160108E+02
      cw( 3, 8)= -0.1257399946089560E+02
      cw( 3, 9)= -0.4649605991087689E+02
      cw( 3,10)=  0.1291889718575365E+01
      cw( 3,11)= -0.6642024885544814E+02
      cw( 3,12)=  0.4010379464615109E+01
      cw( 3,13)= -0.3431775583898482E+01
      cw( 4, 4)=  0.1247990402880702E+03
      cw( 4, 5)= -0.2614088544841073E+03
      cw( 4, 6)=  0.1639702354755679E+02
      cw( 4, 7)= -0.1537951898278491E+00
      cw( 4, 8)=  0.3289336296381872E+01
      cw( 4, 9)=  0.1385407648373453E+03
      cw( 4,10)= -0.4096744426911602E+00
      cw( 4,11)=  0.2004927198755444E+02
      cw( 4,12)= -0.2112379453531347E+02
      cw( 4,13)=  0.2383872746858495E-01
      cw( 5, 5)=  0.7024270444800600E+03
      cw( 5, 6)= -0.3376177605985115E+02
      cw( 5, 7)=  0.3383059349997920E+01
      cw( 5, 8)= -0.1061818752471341E+02
      cw( 5, 9)= -0.4269987402612845E+03
      cw( 5,10)=  0.1220651264863564E+01
      cw( 5,11)= -0.5611470498615714E+02
      cw( 5,12)=  0.4397583312489982E+02
      cw( 5,13)= -0.9658971669110616E+00
      cw( 6, 6)=  0.2291985954551302E+01
      cw( 6, 7)= -0.2184310526902095E+00
      cw( 6, 8)=  0.1875477384629954E+00
      cw( 6, 9)=  0.1824278351439266E+02
      cw( 6,10)= -0.2116926876614674E-03
      cw( 6,11)=  0.2169199742837241E+01
      cw( 6,12)= -0.2849442334008708E+01
      cw( 6,13)=  0.5118126291750101E-01
      cw( 7, 7)=  0.3758452930960190E+01
      cw( 7, 8)= -0.1349498262868165E+01
      cw( 7, 9)= -0.3882998377698184E+01
      cw( 7,10)=  0.7249574867212937E-02
      cw( 7,11)= -0.8394086331322780E+01
      cw( 7,12)=  0.1352359891012749E+00
      cw( 7,13)= -0.9408527584351993E+00
      cw( 8, 8)=  0.2622016304427186E+01
      cw( 8, 9)=  0.6635412812112852E+01
      cw( 8,10)= -0.3780683277159227E+00
      cw( 8,11)=  0.8904280546187730E+01
      cw( 8,12)= -0.4964202159894149E+00
      cw( 8,13)=  0.3525170178945272E+00
      cw( 9, 9)=  0.2800559347499195E+03
      cw( 9,10)= -0.6682475117434016E+00
      cw( 9,11)=  0.3518627356247458E+02
      cw( 9,12)= -0.2377323930466710E+02
      cw( 9,13)=  0.1085128547215229E+01
      cw(10,10)=  0.6871135403140803E-01
      cw(10,11)= -0.9844161870710365E+00
      cw(10,12)=  0.5381099927939942E-01
      cw(10,13)= -0.3528489284187723E-02
      cw(11,11)=  0.5229123693616295E+02
      cw(11,12)= -0.3209946113892272E+01
      cw(11,13)=  0.2190581169027858E+01
      cw(12,12)=  0.3645401391485986E+01
      cw(12,13)= -0.2896623579261971E-01
      cw(13,13)=  0.2399007230858955E+00

c For cool stars
c N=551 sigma_res=0.53953 sigma_typ=0.19062 r**2=0.948
Coef  1 :        -0.7087350516888279E+02        0.1456751070393353E+02
Coef  2 :         0.3022258042530659E+03        0.4401446452747689E+02
Coef  3 :        -0.2058493937239117E+02        0.2112468865975549E+01
Coef  4 :         0.5760244588427595E+02        0.1119946292558866E+02
Coef  5 :        -0.8081379882056881E+02        0.2067339004209639E+02
Coef  6 :         0.1277313725542240E+02        0.1881607989553389E+01
Coef  7 :        -0.3124706353011856E+03        0.4421579356818125E+02
Coef  8 :         0.3514332599689722E+01        0.4534293525888189E+00
Coef  9 :         0.2314309647301204E+01        0.8884964695303000E+00
Coef 10 :        -0.5789383407667002E+01        0.9708723394443542E+00
Coef 12 :         0.9974744972708994E+02        0.1450852324449692E+02
Coef 13 :        -0.1520157613785365E+00        0.3226902585814677E-01
Coef 15 :         0.1314551535472353E+00        0.7382020411302488E-01
Coef 16 :         0.2958219688333445E+02        0.9289226750612295E+01
Coef 17 :        -0.2103162000202632E+01        0.9157040312900007E+00
Coef 18 :        -0.1674306650092225E+01        0.3290789631448154E+00
Coef 20 :         0.4070640493841211E+01        0.8855185761972609E+00

      fc(1)= -0.7087350516888279E+02
      fc(2)=  0.3022258042530659E+03
      fc(3)= -0.2058493937239117E+02
      fc(4)=  0.5760244588427595E+02
      fc(5)= -0.8081379882056881E+02
      fc(6)=  0.1277313725542240E+02
      fc(7)= -0.3124706353011856E+03
      fc(8)=  0.3514332599689722E+01
      fc(9)=  0.2314309647301204E+01
      fc(10)=-0.5789383407667002E+01
      fc(12)= 0.9974744972708994E+02
      fc(13)=-0.1520157613785365E+00
      fc(15)= 0.1314551535472353E+00
      fc(16)= 0.2958219688333445E+02
      fc(17)=-0.2103162000202632E+01
      fc(18)=-0.1674306650092225E+01
      fc(20)= 0.4070640493841211E+01
      lc(1)=.true. 
      lc(2)=.true. 
      lc(3)=.true. 
      lc(4)=.true. 
      lc(5)=.true. 
      lc(6)=.true. 
      lc(7)=.true. 
      lc(8)=.true. 
      lc(9)=.true. 
      lc(10)=.true.
      lc(12)=.true.
      lc(13)=.true.
      lc(15)=.true.
      lc(16)=.true.
      lc(17)=.true.
      lc(18)=.true.
      lc(20)=.true.
      src=0.5395257782E+00
      cc( 1, 1)=  0.7290314971722066E+03
      cc( 1, 2)= -0.2172104856318294E+04
      cc( 1, 3)= -0.3991269238222725E+01
      cc( 1, 4)= -0.3832230725091786E+01
      cc( 1, 5)=  0.6034823384371128E+01
      cc( 1, 6)= -0.1032427306006564E+01
      cc( 1, 7)=  0.2130066794854820E+04
      cc( 1, 8)=  0.9299463185206085E+00
      cc( 1, 9)= -0.2661622133854840E+01
      cc( 1,10)= -0.5438637831244631E+01
      cc( 1,11)= -0.6847882703823188E+03
      cc( 1,12)= -0.1076148352736095E+00
      cc( 1,13)=  0.4853278960493291E+00
      cc( 1,14)= -0.9373500459267928E+00
      cc( 1,15)=  0.2101304740459673E+01
      cc( 1,16)=  0.6056108623939594E+00
      cc( 1,17)=  0.2999948020083868E+01
      cc( 2, 2)=  0.6655281744875102E+04
      cc( 2, 3)= -0.3160303044936469E+02
      cc( 2, 4)=  0.6808424545376316E+02
      cc( 2, 5)= -0.9875528336950508E+02
      cc( 2, 6)=  0.4144171594275024E+02
      cc( 2, 7)= -0.6645175961884298E+04
      cc( 2, 8)=  0.4803377928345587E+01
      cc( 2, 9)=  0.5873665952514834E+01
      cc( 2,10)=  0.6243631314242021E+01
      cc( 2,11)=  0.2158070158228403E+04
      cc( 2,12)=  0.6493201868523413E-01
      cc( 2,13)= -0.1089167416458036E+01
      cc( 2,14)=  0.3045965264595750E+02
      cc( 2,15)= -0.3868235331572442E+01
      cc( 2,16)= -0.7239877123989358E+01
      cc( 2,17)= -0.1327370769575413E+01
      cc( 3, 3)=  0.1533049698883655E+02
      cc( 3, 4)= -0.8282383443126394E+01
      cc( 3, 5)=  0.7186077537715893E+01
      cc( 3, 6)= -0.1321074017256200E+02
      cc( 3, 7)=  0.5509760538943265E+02
      cc( 3, 8)= -0.3076185252053324E+01
      cc( 3, 9)=  0.3887313169493996E+00
      cc( 3,10)=  0.2445716807814792E+01
      cc( 3,11)= -0.2073128867484598E+02
      cc( 3,12)=  0.1024575885802276E+00
      cc( 3,13)= -0.1207333081240571E+00
      cc( 3,14)= -0.3364500918157235E+00
      cc( 3,15)= -0.5507069139529219E+00
      cc( 3,16)=  0.2211909490653990E+01
      cc( 3,17)= -0.1617297112945044E+01
      cc( 4, 4)=  0.4308935498956192E+03
      cc( 4, 5)= -0.7890936257790466E+03
      cc( 4, 6)=  0.5299865238599524E+01
      cc( 4, 7)= -0.1117728661871987E+03
      cc( 4, 8)=  0.1766774154600654E+01
      cc( 4, 9)=  0.1475695117493788E+02
      cc( 4,10)= -0.2385233885262151E+02
      cc( 4,11)=  0.4836340288312837E+02
      cc( 4,12)= -0.1861737067497626E+00
      cc( 4,13)= -0.8149479436991713E+00
      cc( 4,14)=  0.3479695222577717E+03
      cc( 4,15)= -0.1610342028075634E+02
      cc( 4,16)= -0.3228123822500148E+00
      cc( 4,17)=  0.2891574349601288E+02
      cc( 5, 5)=  0.1468246577825917E+04
      cc( 5, 6)= -0.2925126743194039E+01
      cc( 5, 7)=  0.1696427683082147E+03
      cc( 5, 8)= -0.1829875356847343E+01
      cc( 5, 9)= -0.2596289105891584E+02
      cc( 5,10)=  0.3860713115511245E+02
      cc( 5,11)= -0.7723206085675164E+02
      cc( 5,12)=  0.3098459892839541E+00
      cc( 5,13)=  0.1879560368971104E+01
      cc( 5,14)= -0.6558101695604281E+03
      cc( 5,15)=  0.2881841866195740E+02
      cc( 5,16)= -0.4617814168566072E+00
      cc( 5,17)= -0.5046204880636127E+02
      cc( 6, 6)=  0.1216280929206568E+02
      cc( 6, 7)= -0.6246831509446773E+02
      cc( 6, 8)=  0.2458081358787715E+01
      cc( 6, 9)= -0.3373130602629392E+00
      cc( 6,10)= -0.1944060377927810E+01
      cc( 6,11)=  0.2267540496435518E+02
      cc( 6,12)= -0.4685965323663489E-01
      cc( 6,13)=  0.7855979697867212E-01
      cc( 6,14)= -0.1322924240488459E+01
      cc( 6,15)=  0.4665164628818430E+00
      cc( 6,16)= -0.2051357025710202E+01
      cc( 6,17)=  0.1383466499391804E+01
      cc( 7, 7)=  0.6716305590702032E+04
      cc( 7, 8)= -0.8710515102271312E+01
      cc( 7, 9)= -0.4336349515525471E+01
      cc( 7,10)= -0.2445951738698011E+00
      cc( 7,11)= -0.2197679612580000E+04
      cc( 7,12)=  0.8275433670233230E-01
      cc( 7,13)=  0.9407167836222756E+00
      cc( 7,14)= -0.5862285546963736E+02
      cc( 7,15)=  0.2247986162758540E+01
      cc( 7,16)=  0.9968018825821179E+01
      cc( 7,17)= -0.3772775101951983E+01
      cc( 8, 8)=  0.7063091972675920E+00
      cc( 8, 9)= -0.6067250613026813E-01
      cc( 8,10)= -0.4001403527442108E+00
      cc( 8,11)=  0.3291101744050278E+01
      cc( 8,12)= -0.3431108419515742E-01
      cc( 8,13)=  0.2159897548794922E-01
      cc( 8,14)=  0.3449373501081838E+00
      cc( 8,15)=  0.9807077089353204E-01
      cc( 8,16)= -0.4329968245655792E+00
      cc( 8,17)=  0.2427984704306224E+00
      cc( 9, 9)=  0.2711983314572205E+01
      cc( 9,10)= -0.1662427393851166E+00
      cc( 9,11)=  0.1051753762560688E+01
      cc( 9,12)=  0.3728585697246172E-04
      cc( 9,13)= -0.8855491225140164E-01
      cc( 9,14)=  0.1031797975354038E+02
      cc( 9,15)= -0.2773225387290768E+01
      cc( 9,16)=  0.5431145719107300E-01
      cc( 9,17)=  0.7339088467678652E+00
      cc(10,10)=  0.3238171576303523E+01
      cc(10,11)= -0.9447417885760687E+00
      cc(10,12)=  0.1482338381901688E-01
      cc(10,13)= -0.1031390696524925E+00
      cc(10,14)= -0.1548270780800556E+02
      cc(10,15)=  0.1741371447171051E+00
      cc(10,16)=  0.2370575076206407E+00
      cc(10,17)= -0.2635309905277247E+01
      cc(11,11)=  0.7231393924206724E+03
      cc(11,12)= -0.6280115856724940E-01
      cc(11,13)= -0.2974296378327303E+00
      cc(11,14)=  0.2897743046062516E+02
      cc(11,15)= -0.3814117664916540E+00
      cc(11,16)= -0.3435377779502425E+01
      cc(11,17)=  0.2226192366172138E+01
      cc(12,12)=  0.3577233674945322E-02
      cc(12,13)= -0.7315694518960863E-03
      cc(12,14)= -0.1326801233578677E+00
      cc(12,15)= -0.1654556420146421E-02
      cc(12,16)=  0.7192511116071472E-02
      cc(12,17)= -0.9201287887715884E-02
      cc(13,13)=  0.1872087242144528E-01
      cc(13,14)= -0.9147325413481322E+00
      cc(13,15)=  0.9881930230880608E-01
      cc(13,16)= -0.1034698921469510E-01
      cc(13,17)= -0.5574061694479755E-02
      cc(14,14)=  0.2964385829874507E+03
      cc(14,15)= -0.1164239537380561E+02
      cc(14,16)=  0.6094449275264173E+00
      cc(14,17)=  0.2129962118161965E+02
      cc(15,15)=  0.2880619210000385E+01
      cc(15,16)= -0.7836857072505710E-01
      cc(15,17)= -0.8134827770706332E+00
      cc(16,16)=  0.3720281827597321E+00
      cc(16,17)= -0.1537386114468435E+00
      cc(17,17)=  0.2693834760214337E+01

c For cold dwarfs
c N=23 sigma_res=0.28989 sigma_typ=0.20047 r**2=0.985
Coef  1 :     -0.1046722017025981E+03       0.2376398710269291E+02
Coef  2 :      0.2451506517215240E+03       0.5105058125125843E+02
Coef  7 :     -0.1731898012961590E+03       0.3603653336964771E+02
Coef 12 :      0.3871691678519028E+02       0.8349521982443495E+01

      fcd(1)= -0.1046722017025981E+03
      fcd(2)=  0.2451506517215240E+03
      fcd(7)= -0.1731898012961590E+03
      fcd(12)= 0.3871691678519028E+02
      lcd(1)=.true.
      lcd(2)=.true.
      lcd(7)=.true.
      lcd(12)=.true.
      srcd=0.2898909992E+00
      ccd( 1, 1)=  0.6719998153312355E+04
      ccd( 1, 2)= -0.1441994567400875E+05
      ccd( 1, 3)=  0.1014643546988935E+05
      ccd( 1, 4)= -0.2339353083081731E+04
      ccd( 2, 2)=  0.3101215316149643E+05
      ccd( 2, 3)= -0.2186841785637790E+05
      ccd( 2, 4)=  0.5052029967243325E+04
      ccd( 3, 3)=  0.1545313327258727E+05
      ccd( 3, 4)= -0.3577042152414013E+04
      ccd( 4, 4)=  0.8295713838435843E+03

c For cold giants
c N=27 sigma_res=0.88810 sigma_typ=0.12661 r**2=0.831
Coef  1 :        -0.2965821442013717E+02        0.1428390444337743E+02
Coef  2 :         0.4810243560543612E+02        0.1975588735859654E+02
Coef  7 :        -0.1821468380242129E+02        0.6803748507505433E+01

      lecg=.true.
      ctecg=2.0D0
      fcg(1)=-0.2965821442013717E+02
      fcg(2)= 0.4810243560543612E+02
      fcg(7)=-0.1821468380242129E+02
      lcg(1)=.true.
      lcg(2)=.true.
      lcg(7)=.true.
      srcg=0.2481292052E+00
      ccg( 1, 1)=  0.3313890113648655E+04
      ccg( 1, 2)= -0.4580491053703806E+04
      ccg( 1, 3)=  0.1574286679515214E+04
      ccg( 2, 2)=  0.6339241743050667E+04
      ccg( 2, 3)= -0.2181645620594679E+04
      ccg( 3, 3)=  0.7518665004043432E+03

c      return
c      end

c For intermediate giants (1<theta<1.4) to perform a smooth iterpolation.
Coef  1 :          0.3688936997817394E+03        0.1463868381411887E+03
Coef  2 :        -0.8849125022585878E+03        0.3829445497969912E+03
Coef  3 :        -0.7259751833879855E+01        0.2101212933112567E+01
Coef  4 :          0.1040247204760071E+02        0.2958125388217792E+01
Coef  5 :        -0.5575005421343860E+01        0.2476710759482117E+01
Coef  7 :          0.7427265519059204E+03        0.3334113477968885E+03
Coef  8 :          0.6285107277827846E+00        0.1689087628230818E+00
Coef 10 :       -0.8552241053414129E+00        0.2304699850569047E+00
Coef 11 :         0.2389086442786061E+01        0.1311606747790137E+01
Coef 12 :       -0.2098484867656973E+03        0.9640982132362731E+02

      sri=0.5286706830E+00

      fi(1)=  0.3688936997817394E+03
      fi(2)=-0.8849125022585878E+03
      fi(3)=-0.7259751833879855E+01
      fi(4)=  0.1040247204760071E+02 
      fi(5)=-0.5575005421343860E+01
      fi(7)=  0.7427265519059204E+03
      fi(8)=  0.6285107277827846E+00 
      fi(10)=-0.8552241053414129E+00
      fi(11)=  0.2389086442786061E+01
      fi(12)=-0.2098484867656973E+03 
      li(1)=.true.
      li(2)=.true.
      li(3)=.true.
      li(4)=.true.
      li(5)=.true.
      li(7)=.true.
      li(8)=.true.
      li(10)=.true.
      li(11)=.true.
      li(12)=.true.

      ci(1, 1)= 0.7667143746461114E+05
      ci(1, 2)=-0.2001099035603601E+06
      ci(1, 3)=-0.2872849777070574E+02
      ci(1, 4)=-0.3879825564626826E+03
      ci(1, 5)= 0.3290332614047143E+03
      ci(1, 6)= 0.1733346402965799E+06
      ci(1, 7)= 0.3242800709798209E+01
      ci(1, 8)= 0.2159120536943892E+02
      ci(1, 9)= 0.1839679942197472E+02
      ci(1,10)=-0.4982423498395699E+05
      ci(2, 2)= 0.5246882402779374E+06
      ci(2, 3)=-0.9592150707315798E+02
      ci(2, 4)= 0.1146996601079436E+04
      ci(2, 5)=-0.9767133854697819E+03
      ci(2, 6)=-0.4561498068342509E+06
      ci(2, 7)=-0.4767673552600417E-01
      ci(2, 8)=-0.6016201700320376E+02
      ci(2, 9)= 0.6439136344848664E+02
      ci(2,10)= 0.1314847657452301E+06
      ci(3, 3)= 0.1579682021109497E+02
      ci(3, 4)=-0.8323779191593324E+01
      ci(3, 5)= 0.7101960450909760E+01
      ci(3, 6)= 0.1942139219443618E+03
      ci(3, 7)=-0.1047640338468818E+01
      ci(3, 8)= 0.3941339760057139E+00
      ci(3, 9)=-0.9584485072063480E+01
      ci(3,10)=-0.7857930087584137E+02
      ci(4, 4)= 0.3130853182829086E+02
      ci(4, 5)=-0.2607334466576119E+02
      ci(4, 6)=-0.1086753738693339E+04
      ci(4, 7)= 0.2807996568585902E+00
      ci(4, 8)=-0.1693267151715232E+01
      ci(4, 9)= 0.5706437853263377E+01
      ci(4,10)= 0.3327788191188060E+03
      ci(5, 5)= 0.2194725079878623E+02
      ci(5, 6)= 0.9287366363294755E+03
      ci(5, 7)=-0.2174211001233173E+00
      ci(5, 8)= 0.1279502072092572E+01
      ci(5, 9)=-0.4998827619411261E+01
      ci(5,10)=-0.2852027961324281E+03
      ci(6, 6)= 0.3977317848135094E+06
      ci(6, 7)=-0.4989243490900344E+01
      ci(6, 8)= 0.5431084175364924E+02
      ci(6, 9)=-0.1306873850764194E+03
      ci(6,10)=-0.1149063941236902E+06
      ci(7, 7)= 0.1020784123363358E+00
      ci(7, 8)=-0.2997626178797578E-01
      ci(7, 9)= 0.5402107743601779E+00
      ci(7,10)= 0.2392930264654306E+01
      ci(8, 8)= 0.1900458069913316E+00
      ci(8, 9)=-0.1932802921352775E+00
      ci(8,10)=-0.1601769461815670E+02
      ci(9, 9)= 0.6155124324006462E+01
      ci(9,10)= 0.5320222975531255E+02
      ci(10,10)=0.3325615998721759E+05

c For warm giants (0<theta<1.2) of low metallicity ([Fe/H]<-0.29) 
c to perform a smooth iterpolation with the hot stars
Coef  1 :              0.1049380017437647E+01
Coef  7 :              0.4438353082387340E+01
Coef  9 :              0.4209784606407982E+00
Coef 11 :              0.2460264291903057E+00
Coef 12 :              0.3060608465022252E+01
Coef 14 :              0.1609101065613140E+00

c      srlm=0.5193621786E+00
c      flm(1)=   -0.2337677246572492E+01 
c      flm(7)=    0.2442633756775564E+02 
c      flm(9)=   -0.1992840282196657E+01 
c      flm(11)=  -0.5404998643511945E+00 
c      flm(12)=  -0.1348398781369407E+02 
c      flm(14)=  -0.4481361165768376E+00 
c      llm(1)=.true.
c      llm(7)=.true.
c      llm(9)=.true.
c      llm(11)=.true.
c      llm(12)=.true.
c      llm(14)=.true.
c      clm(1,1)=  0.4082488219526129E+01  
c      clm(1,2)= -0.1372421626331690E+02  
c      clm(1,3)=  0.1342981601471380E+00  
c      clm(1,4)=  0.1725020113847485E+00  
c      clm(1,5)=  0.8977069124583632E+01  
c      clm(1,6)=  0.3404275981734390E-01  
c      clm(2,2)=  0.7303029538630648E+02  
c      clm(2,3)= -0.2985367554377732E+01  
c      clm(2,4)= -0.2529146513521098E+01  
c      clm(2,5)= -0.5002755253098195E+02  
c      clm(2,6)= -0.9180755202012417E+00  
c      clm(3,3)=  0.6570207893832742E+00  
c      clm(3,4)=  0.1827839708526419E+00  
c      clm(3,5)=  0.2087659377453877E+01  
c      clm(3,6)=  0.2454467346787623E+00  
c      clm(4,4)=  0.2244000177323338E+00  
c      clm(4,5)=  0.1670269381737149E+01  
c      clm(4,6)=  0.5214813295371387E-01  
c      clm(5,5)=  0.3472761118113564E+02  
c      clm(5,6)=  0.6493498024138842E+00  
c      clm(6,6)=  0.9599000296899736E-01  

c For warm giants (0<theta<1.2) of low metallicity ([Fe/H]<-0.29) 
c to perform a smooth iterpolation with the hot stars

Coef  1 :                0.1006654245098970E+01
Coef  5 :                0.6694668393825272E+00
Coef  7 :                0.4561546509208224E+01
Coef 11 :                0.2790099963599746E+00
Coef 12 :                0.3573057702352454E+01
Coef 20 :                0.2929488265968303E+00

c      srlm=0.8261527952E+00
c      flm(1)= -0.3829214280950368E+01
c      flm(5)=  0.4322377752693303E+01
c      flm(7)=  0.4242936241454660E+02
c      flm(11)= -0.1872579130809023E+01
c      flm(12)=-0.2538853919413639E+02
c      flm(20)=-0.6320995997011418E+00
c      llm(1)=.true.
c      llm(5)=.true.
c      llm(7)=.true.
c      llm(11)=.true.
c      llm(12)=.true.
c      llm(20)=.true.
c      clm(1,1)= 0.1484704091936107E+01
c      clm(1,2)= 0.1701929847086155E-02
c      clm(1,3)=-0.6089109243851533E+01
c      clm(1,4)= 0.1271179481287679E-01
c      clm(1,5)= 0.4489577414239366E+01
c      clm(1,6)=-0.1400455919709659E-01
c      clm(2,2)= 0.6566551987104466E+00
c      clm(2,3)= 0.7643441136956576E-01
c      clm(2,4)=-0.1590851730316463E+00
c      clm(2,5)= 0.3241183815456172E+00
c      clm(2,6)=-0.2717232091785816E+00
c      clm(3,3)= 0.3048621171888249E+02
c      clm(3,4)=-0.1949087172871688E+00
c      clm(3,5)=-0.2335047722315808E+02
c      clm(3,6)= 0.7342065207950883E-01
c      clm(4,4)= 0.1140561673265751E+00
c      clm(4,5)=-0.1085947543660751E+00
c      clm(4,6)= 0.6524894558268005E-01
c      clm(5,5)= 0.1870506865052425E+02
c      clm(5,6)=-0.2127527584160448E+00
c      clm(6,6)= 0.1257369068490164E+00

c For warm giants (0<theta<1.1) of low metallicity ([Fe/H]<-0.25) 
c to perform a smooth iterpolation with the hot stars
Coef  1 :        -0.7939379746366495E+00        0.2028073086067468E+01
Coef  2 :         0.8593839071143471E+01        0.1992861782565105E+01
Coef  5 :         0.2258308001031946E+01        0.2706416819460465E+00

c      srlm=0.6044705908E+00
c      flm(1)= -0.7939379746366495E+00
c      flm(2)=  0.8593839071143471E+01
c      flm(5)=  0.2258308001031946E+01
c      llm(1)=.true.
c      llm(2)=.true.
c      llm(5)=.true.
c      clm(1,1)= 0.1125684928102006E+02
c      clm(1,2)=-0.1096218961408537E+02
c      clm(1,3)= 0.6823816786104085E-01
c      clm(2,2)= 0.1086936080676830E+02
c      clm(2,3)= 0.9112026507806012E-01
c      clm(3,3)= 0.2004652110148580E+00

c For warm giants (0<theta<1.1) of low metallicity ([Fe/H]<-0.25) including <g> term
c to perform a smooth iterpolation with the hot stars
Coef  1 :        -0.1212936084752947E+00        0.2277141998307367E+01
Coef  2 :         0.8780042585941301E+01        0.2055739609606027E+01
Coef  3 :        -0.3118236523224460E+00        0.4344705528762989E+00
Coef  5 :         0.2395626257426897E+01        0.3366112233614504E+00

      srlm=0.6185569535E+00
      flm(1)=-0.1212936084752947E+00
      flm(2)= 0.8780042585941301E+01
      flm(3)=-0.3118236523224460E+00
      flm(5)= 0.2395626257426897E+01
      llm(1)=.true.
      llm(2)=.true.
      llm(3)=.true.
      llm(5)=.true.
      clm(1,1)= 0.1355254442944315E+02
      clm(1,2)=-0.1032668811129193E+02
      clm(1,3)=-0.1064235548251992E+01
      clm(1,4)= 0.5368971961004119E+00
      clm(2,2)= 0.1104528231809331E+02
      clm(2,3)=-0.2946050091645182E+00
      clm(2,4)= 0.2208559213773388E+00
      clm(3,3)= 0.4933570134262479E+00
      clm(3,4)=-0.2172603789332134E+00
      clm(4,4)= 0.2961404947741711E+00

      return
      end

c-----------------------------------------------------------------------------------------------
c ################## INTERPOLACION DE FUNCIONES Y COEFICIENTES: <PaT> ##########################
c-----------------------------------------------------------------------------------------------

      subroutine finfitpat(t,g,z,findex,eindex,iflag)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
      integer iflag
      double precision t,g,z,findex,eindex,pi
      
      double precision findex1,findex2,findex3
      double precision eindex1,eindex2,eindex3
      double precision thet,x
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c   fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      double precision theta,geta,zeta
      double precision xf(20)
      logical nog,noz
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      common/blkx/xf

C Regions:
C 1 = hot dwarfs (hd)
C 2 = hot giants (hg)
C 3 = warm stars (w)
C 4 = cool stars (c)
C 5 = cold dwarfs (cd)
C 6 = cold giants (cg)
C 7 = intermediate (i)

c Coefficients and variance-covariance matrices are stored in this subroutine
      call readcoefpat

c-----------------------------------------------------------------------
      pi=3.14159d0

      nog=.false.
      noz=.false.
      if(g.ge.99.d0) nog=.true.
      if(z.ge.99.d0) noz=.true.

      iflag=0
      findex=0.d0
      eindex=0.d0
      if(t.le.0.d0) then
         iflag=-1
         return
      end if
      if(nog) then
         iflag=-3
         return
      end if

c Program works in theta (=5040/Teff)
      theta=dble(t)
      thet=t
      geta=dble(g)
      zeta=dble(z)
      xf(1)=1.D0
      xf(2)=theta
      xf(3)=geta
      xf(4)=zeta
      xf(5)=theta*zeta
      xf(6)=theta*geta
      xf(7)=theta*theta
      xf(8)=geta*geta
      xf(9)=zeta*zeta
      xf(10)=geta*zeta
      xf(11)=theta*theta*geta
      xf(12)=theta*theta*theta
      xf(13)=geta*geta*geta
      xf(14)=zeta*zeta*zeta
      xf(15)=geta*geta*zeta
      xf(16)=theta*theta*zeta
      xf(17)=theta*zeta*zeta
      xf(18)=theta*geta*geta
      xf(19)=geta*zeta*zeta
      xf(20)=theta*geta*zeta
 

c Hot stars

      if(thet.le.0.70d0) then
         if(thet.lt.0.13d0) iflag=1
         if(g.ge.2.80d0) then
            call compindex(1,findex1,eindex1)
         end if
         if(g.le.3.0d0) then
            call compindex(2,findex2,eindex2)
         end if
         if(g.le.2.80d0) then
            findex3=findex2
            eindex3=eindex2
         else
            if(g.ge.3.0d0) then
               findex3=findex1
               eindex3=eindex1
            else
               x=(g-2.80d0)/0.20d0
               findex3=(1.d0-x)*findex2+x*findex1
               eindex3=(1.d0-x)*eindex2+x*eindex1
            end if
         end if

c If really hot or no z values, it keeps this index

         if(thet.le.0.45d0.or.noz) then
            findex=findex3
            eindex=eindex3
            if(g.lt.1.d0.or.g.gt.4.2d0) iflag=3
            return
         end if
      end if

c Warm stars

      if(thet.le.0.90d0) then
         if(noz) then
            iflag=-2
            return
         end if
         call compindex(3,findex1,eindex1)
         if(thet.le.0.70d0) then
            x=cos(pi/2.d0*(0.70d0-thet)/0.25d0)
            findex=(1.d0-x)*findex3+x*findex1
            eindex=(1.d0-x)*eindex3+x*eindex1
            if(g.lt.0.5d0.or.g.gt.4.6d0) iflag=3
            return
         end if
      end if

c Cool stars
      if(thet.le.1.35d0) then
         if(noz) then
            if((g.ge.3.d0.and.thet.lt.1.06d0).or.
     c         (g.lt.3.d0.and.thet.lt.1.3d0)) then
              iflag=-2
              return
            else
               goto 10
            end if
         end if
         call compindex(4,findex2,eindex2)
         if(thet.le.0.90d0) then
            x=cos(pi/2.d0*(thet-0.70d0)/0.20d0)
            findex=(1.d0-x)*findex2+x*findex1
            eindex=(1.d0-x)*eindex2+x*eindex1
            if(g.lt.0.d0.or.g.gt.5.d0) iflag=3
            return
         end if
         if(thet.ge.1.05d0) then
            findex3=findex2
            eindex3=eindex2
         else
            findex=findex2
            eindex=eindex2
            if(g.lt.-0.1d0.or.g.gt.5.d0) iflag=3
            return
         end if
      end if

c intermediate stars:
      if(thet.lt.1.3) then
         call compindex(7,findex2,eindex2)
         if(thet.le.1.1d0) then
            x=cos(pi/2.d0*(thet-1.05d0)/0.05d0)
            findex=(1.d0-x)*findex2+x*findex3
            eindex=(1.d0-x)*eindex2+x*eindex3
            if(g.lt.0.d0.or.g.gt.5.d0) iflag=3
            return
         end if
         if(thet.le.1.2d0) then
            findex=findex2
            eindex=eindex2
            if(g.lt.-0.1d0.or.g.gt.5.d0) iflag=3
            return
         else
            findex3=findex2
            eindex3=eindex2
         end if
      end if

 10   if(g.ge.3.d0) then                      
         if(g.lt.4.4d0.or.g.gt.5.2d0) iflag=3   
         call compindex(5,findex1,eindex1)  
         if(thet.ge.1.3d0.or.noz) then        
            findex=findex1                  
            eindex=eindex1                  
            if(thet.gt.1.89d0) iflag=2        
            return                          
         else                               
            x=cos(pi/2.d0*(thet-1.2d0)/0.1d0)     
            findex=(1.d0-x)*findex1+x*findex3 
            eindex=(1.d0-x)*eindex1+x*eindex3 
            return                          
         end if                             
      else
         if(g.gt.1.8d0.or.g.lt.-0.1d0) iflag=3  
         call compindex(6,findex1,eindex1)  
         if(thet.ge.1.3d0.or.noz) then        
            findex=findex1                  
            eindex=eindex1                  
            if(thet.gt.1.74d0) iflag=2        
            return                          
         else                               
            x=cos(pi/2.d0*(thet-1.2d0)/0.1d0)     
            findex=(1.d0-x)*findex1+x*findex3 
            eindex=(1.d0-x)*eindex1+x*eindex3 
            return                          
         end if                             
      end if


c Error: No value was computed
      iflag=-5
c      write(*,*)t,g,z,findex
c      write(*,'(A)') 'ERROR: NO index value was computed "finfitpat"'
c      stop

      end


      subroutine readcoefpat
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer i,j
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm

c-Initialization----------------------------------------------------------
      do i=1,25
         fhd(i)=0.D0
         fhg(i)=0.D0
         fw(i)=0.D0
         fc(i)=0.D0
         fcd(i)=0.D0
         fcg(i)=0.D0
         fi(i)=0.D0
         lhd(i)=.false.
         lhg(i)=.false.
         lw(i)=.false.
         lc(i)=.false.
         lcd(i)=.false.
         lcg(i)=.false.
         li(i)=.false.
         do j=1,25
            chd(i,j)=0.D0
            chg(i,j)=0.D0
            cw(i,j)=0.D0
            cc(i,j)=0.D0
            ccd(i,j)=0.D0
            ccg(i,j)=0.D0
            ci(i,j)=0.D0
         end do
      end do
      lehd=.false.
      lehg=.false.
      lew=.false.
      lec=.false.
      lecd=.false.
      lecg=.false.
      lei=.false.
      ctehd=0.D0
      ctehg=0.D0
      ctew=0.D0
      ctec=0.D0
      ctecd=0.D0
      ctecg=0.D0
      ctei=0.D0

c For hot dwarfs
c N=46 sigma_res=0.91055 sigma_typ=0.16547 r**2=0.907
Coef  1 :        -0.6414842041998517E+01        0.7915653920861390E+00
Coef  2 :         0.5471170690564588E+02        0.3391977387347575E+01
Coef 11 :        -0.4283154729334703E+01        0.1091360469602459E+01
Coef 12 :        -0.5543770369501858E+02        0.5660695864045900E+01

      srhd=0.9105525189E+00
      fhd(1)= -0.6414842041998517E+01
      fhd(2)=  0.5471170690564588E+02
      fhd(11)=-0.4283154729334703E+01
      fhd(12)=-0.5543770369501858E+02
      lhd(1)=.true.
      lhd(2)=.true.
      lhd(11)=.true.
      lhd(12)=.true.
      chd(1, 1)=  0.7557246640808089E+00
      chd(1, 2)= -0.2960558549963638E+01
      chd(1, 3)=  0.2973113221130755E+00
      chd(1, 4)=  0.2613220064363448E+01
      chd(2, 2)=  0.1387700984856905E+02
      chd(2, 3)= -0.2296138146016453E+01
      chd(2, 4)= -0.9211030832773542E+01
      chd(3, 3)=  0.1436568826092399E+01
      chd(3, 4)= -0.4120317405826246E+01
      chd(4, 4)=  0.3864823307247444E+02

c For hot giants
c N=29 sigma_res=1.31028 sigma_typ=0.12337 r**2=0.638
Coef  1 :         0.1600578784001784E+01        0.1512014713844868E+01
Coef  7 :         0.4481441967825036E+02        0.2141114317773037E+02
Coef 12 :        -0.4688020793569266E+02        0.2920405282121500E+02     

      srhg=0.1310275642E+01
      fhg(1)=  0.1600578784001784E+01
      fhg(7)=  0.4481441967825036E+02
      fhg(12)=-0.4688020793569266E+02
      lhg(1)=.true.
      lhg(7)=.true.
      lhg(12)=.true.
      chg(1, 1)=  0.1331639593680920E+01
      chg(1, 2)= -0.1658956419068594E+02
      chg(1, 3)=  0.2104426335542746E+02
      chg(2, 2)=  0.2670265077686691E+03
      chg(2, 3)= -0.3610334950577474E+03
      chg(3, 3)=  0.4967763534664019E+03

c For warm stars ***** Type (1) *****
c N=193 sigma_res=0.55323 sigma_typ=0.16987 r**2=0.962
Coef  1 :        -0.1490223712122600E+03        0.2168265245152611E+02
Coef  2 :         0.6511855098857886E+03        0.8841202255706503E+02
Coef  3 :         0.1873879029352709E+02        0.3834616413010322E+01
Coef  5 :        -0.4651097177403340E+01        0.1765983253306492E+01
Coef  6 :        -0.5639648093112812E+02        0.1091799036371259E+02
Coef  7 :        -0.8409103925848722E+03        0.1239715601668017E+03
Coef 11 :         0.3925573822172255E+02        0.7616416142122161E+01
Coef 12 :         0.3379593962162174E+03        0.6000671265103185E+02
Coef 16 :         0.5879643313356195E+01        0.2139014787958550E+01

      srw=0.5532316326E+00
      fw(1)= -0.1490223712122600E+03
      fw(2)=  0.6511855098857886E+03
      fw(3)=  0.1873879029352709E+02 
      fw(5)= -0.4651097177403340E+01
      fw(6)= -0.5639648093112812E+02
      fw(7)= -0.8409103925848722E+03
      fw(11)= 0.3925573822172255E+02
      fw(12)= 0.3379593962162174E+03
      fw(16)= 0.5879643313356195E+01
      lw(1)=.true.
      lw(2)=.true.
      lw(3)=.true.
      lw(5)=.true.
      lw(6)=.true.
      lw(7)=.true.
      lw(11)=.true.
      lw(12)=.true.
      lw(16)=.true.
      cw(1, 1)=  0.1536069298095563E+04
      cw(1, 2)= -0.6113941338131892E+04
      cw(1, 3)= -0.9838261497798115E+02
      cw(1, 4)=  0.4315108233296499E+01
      cw(1, 5)=  0.2775752401684439E+03
      cw(1, 6)=  0.7911205942323512E+04
      cw(1, 7)= -0.1914596276075634E+03
      cw(1, 8)= -0.3339055523758496E+04
      cw(1, 9)= -0.5220826143284483E+01
      cw(2, 2)=  0.2553927963201097E+05
      cw(2, 3)=  0.1756104434707228E+03
      cw(2, 4)= -0.4075077816530224E+02
      cw(2, 5)= -0.4959889322844957E+03
      cw(2, 6)= -0.3485098808087657E+05
      cw(2, 7)=  0.3430475111595326E+03
      cw(2, 8)=  0.1558689143415327E+05
      cw(2, 9)=  0.4915444909333141E+02
      cw(3, 3)=  0.4804296975245864E+02
      cw(3, 4)=  0.4671557332845460E+01
      cw(3, 5)= -0.1363090983195444E+03
      cw(3, 6)=  0.1024441096845876E+03
      cw(3, 7)=  0.9432484123741678E+02
      cw(3, 8)= -0.2056177197184068E+03
      cw(3, 9)= -0.5549482959094648E+01
      cw(4, 4)=  0.1018964733756857E+02
      cw(4, 5)= -0.1252506295225635E+02
      cw(4, 6)=  0.8740347728984709E+02
      cw(4, 7)=  0.8274948831219591E+01
      cw(4, 8)= -0.5363112203927325E+02
      cw(4, 9)= -0.1229567671835213E+02
      cw(5, 5)=  0.3894676633222779E+03
      cw(5, 6)= -0.3018624340130482E+03
      cw(5, 7)= -0.2710389410536912E+03
      cw(5, 8)=  0.5969555314583603E+03
      cw(5, 9)=  0.1487074706098748E+02
      cw(6, 6)=  0.5021461360440202E+05
      cw(6, 7)=  0.2135570781143478E+03
      cw(6, 8)= -0.2371047150490940E+05
      cw(6, 9)= -0.1054489928829617E+03
      cw(7, 7)=  0.1895340842327231E+03
      cw(7, 8)= -0.4201831391928452E+03
      cw(7, 9)= -0.9820046876808831E+01
      cw(8, 8)=  0.1176483017537042E+05
      cw(8, 9)=  0.6478126764065532E+02
      cw(9, 9)=  0.1494904901082926E+02


c For cool stars
c N=551 sigma_res=0.29565 sigma_typ=0.17091 r**2=0.888
Coef  1 :         0.1771057843811058E+03        0.1195332457234410E+02
Coef  2 :        -0.4454748100851141E+03        0.3370782940513188E+02
Coef  3 :        -0.1599897246926430E+02        0.1377647465160583E+01
Coef  6 :         0.2920219627575549E+02        0.2821488576566488E+01
Coef  7 :         0.3729055830763336E+03        0.3121541987725868E+02
Coef  9 :        -0.3610242222085497E+00        0.1032968208245468E+00
Coef 11 :        -0.1333366141867366E+02        0.1430307087690397E+01
Coef 12 :        -0.1032931599724490E+03        0.9500331108348696E+01
Coef 14 :        -0.1079659157514062E+00        0.4527165813017063E-01

      src=0.2956477348E+00
      fc( 1)=  0.1771057843811058E+03 
      fc( 2)= -0.4454748100851141E+03 
      fc( 3)= -0.1599897246926430E+02 
      fc( 6)=  0.2920219627575549E+02 
      fc( 7)=  0.3729055830763336E+03 
      fc( 9)= -0.3610242222085497E+00 
      fc(11)= -0.1333366141867366E+02
      fc(12)= -0.1032931599724490E+03
      fc(14)= -0.1079659157514062E+00
      lc(1)=.true.
      lc(2)=.true.
      lc(3)=.true.
      lc(6)=.true.
      lc(7)=.true.
      lc(9)=.true.
      lc(11)=.true.
      lc(12)=.true.
      lc(14)=.true.
      cc(1, 1)=  0.1634663301604898E+04
      cc(1, 2)= -0.4589732160033514E+04
      cc(1, 3)= -0.1479740563543551E+03
      cc(1, 4)=  0.3072551161921372E+03
      cc(1, 5)=  0.4197640777779317E+04
      cc(1, 6)= -0.1247111415769479E+01
      cc(1, 7)= -0.1563496400350022E+03
      cc(1, 8)= -0.1251643685716961E+04
      cc(1, 9)= -0.4875891669393784E+00
      cc(2, 2)=  0.1299907540348569E+05
      cc(2, 3)=  0.3932208160508007E+03
      cc(2, 4)= -0.8234027051034848E+03
      cc(2, 5)= -0.1198856836465734E+05
      cc(2, 6)=  0.3453038599776237E+01
      cc(2, 7)=  0.4221717881692833E+03
      cc(2, 8)=  0.3604238171042190E+04
      cc(2, 9)=  0.1315754297423843E+01
      cc(3, 3)=  0.2171336251996711E+02
      cc(3, 4)= -0.4431423711620204E+02
      cc(3, 5)= -0.3381508048304661E+03
      cc(3, 6)=  0.2522183912346759E+00
      cc(3, 7)=  0.2224107135380521E+02
      cc(3, 8)=  0.9384187481022262E+02
      cc(3, 9)=  0.1011531394874571E+00
      cc(4, 4)=  0.9107674080205160E+02
      cc(4, 5)=  0.7139904016076794E+03
      cc(4, 6)= -0.5431739023388669E+00
      cc(4, 7)= -0.4601612682004582E+02
      cc(4, 8)= -0.1998010077250525E+03
      cc(4, 9)= -0.2148041083530422E+00
      cc(5, 5)=  0.1114780209970000E+05
      cc(5, 6)= -0.3100381412286014E+01
      cc(5, 7)= -0.3688377998722546E+03
      cc(5, 8)= -0.3378919598804511E+04
      cc(5, 9)= -0.1145648758759504E+01
      cc(6, 6)=  0.1220744564406719E+00
      cc(6, 7)=  0.2868608204945824E+00
      cc(6, 8)=  0.8877884372850690E+00
      cc(6, 9)=  0.5201724887153460E-01
      cc(7, 7)=  0.2340504442713192E+02
      cc(7, 8)=  0.1039870465021625E+03
      cc(7, 9)=  0.1113685710828198E+00
      cc(8, 8)=  0.1032591086435052E+04
      cc(8, 9)=  0.3171871255005078E+00
      cc(9, 9)=  0.2344788584461537E-01

c For cold dwarfs
c N=29 sigma_res=0.37429 sigma_typ=0.17734 r**2=0.626
Coef  1 :         0.1200491649229753E+01        0.2075604807638792E+00
Coef 12 :        -0.3093617468826813E+00        0.6790070081313161E-01

      srcd=0.3742928812E+00
      fcd(1)=  0.1200491649229753E+01
      fcd(12)=-0.3093617468826813E+00
      lcd(1)=.true.
      lcd(12)=.true.
      ccd(1, 1)=  0.3075149269772798E+00
      ccd(1, 2)= -0.8849937686531108E-01
      ccd(2, 2)=  0.3290981031178655E-01

c For cold giants
c N=44 sigma_res=0.35057 sigma_typ=0.11180 r**2=0.787
Coef  1 :         0.3442965807069610E+03        0.1837755822433069E+03
Coef  2 :        -0.7402356340413226E+03        0.3745100729770819E+03
Coef  7 :         0.5270228706861221E+03        0.2530399365580612E+03
Coef 12 :        -0.1238108039872834E+03        0.5669263182847755E+02

      srcg=0.3505687328E+00
      fcg(1)=  0.3442965807069610E+03
      fcg(2)= -0.7402356340413226E+03
      fcg(7)=  0.5270228706861221E+03
      fcg(12)=-0.1238108039872834E+03
      lcg(1)=.true.
      lcg(2)=.true.
      lcg(7)=.true.
      lcg(12)=.true.
      ccg(1, 1)=  0.2748079276703048E+06
      ccg(1, 2)= -0.5598849094337926E+06
      ccg(1, 3)=  0.3780019645695414E+06
      ccg(1, 4)= -0.8457987522765558E+05
      ccg(2, 2)=  0.1141249627230004E+07
      ccg(2, 3)= -0.7708925668275768E+06
      ccg(2, 4)=  0.1725786040683132E+06
      ccg(3, 3)=  0.5209928730908013E+06
      ccg(3, 4)= -0.1166951925485475E+06
      ccg(4, 4)=  0.2615211874974796E+05

c For intermediate stars (1<theta<1.3; z>-1)

      sri=0.2506505188E+00
      fi(1)=0.7600288487198743E+00
      fi(4)=0.1743380405206003E+00
      li(1)=.true.
      li(4)=.true.
      ci(1, 1)=0.1867809416616686E-01
      ci(1, 2)=0.1533157262353484E-01
      ci(2, 2)=0.4809370444541371E-01
      return
      end

c-----------------------------------------------------------------------------------------------
c ################## INTERPOLACION DE FUNCIONES Y COEFICIENTES: <sTiO> ##########################
c-----------------------------------------------------------------------------------------------

      subroutine finfitsTiO(t,g,z,findex,eindex,iflag)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
      integer iflag
      double precision t,g,z,findex,eindex,pi
      
      double precision findex1,findex2,findex3
      double precision eindex1,eindex2,eindex3
      double precision thet,x
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c   fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      double precision theta,geta,zeta
      double precision xf(20)
      logical nog,noz
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      common/blkx/xf

C Regions:
C 1 = hot dwarfs (hd)
C 2 = hot giants (hg)
C 3 = intermediate (i)
C 4 = cold dwarfs (cd)
C 5 = cold giants (c)
C 6 = very cold giants (cg)


c Coefficients and variance-covariance matrices are stored in this subroutine
      call readcoefsTiO

c-----------------------------------------------------------------------
      pi=3.14159d0

      nog=.false.
      noz=.false.
      if(g.ge.99.d0) nog=.true.
      if(z.ge.99.d0) noz=.true.

      iflag=0
      findex=0.d0
      eindex=0.d0
      if(t.le.0.d0) then
         iflag=-1
         return
      end if
      if(nog) then
         iflag=-3
         return
      end if

c Program works in theta (=5040/Teff)
      theta=dble(t)
      thet=t
      geta=dble(g)
      zeta=dble(z)
      xf(1)=1.D0
      xf(2)=theta
      xf(3)=geta
      xf(4)=zeta
      xf(5)=theta*zeta
      xf(6)=theta*geta
      xf(7)=theta*theta
      xf(8)=geta*geta
      xf(9)=zeta*zeta
      xf(10)=geta*zeta
      xf(11)=theta*theta*geta
      xf(12)=theta*theta*theta
      xf(13)=geta*geta*geta
      xf(14)=zeta*zeta*zeta
      xf(15)=geta*geta*zeta
      xf(16)=theta*theta*zeta
      xf(17)=theta*zeta*zeta
      xf(18)=theta*geta*geta
      xf(19)=geta*zeta*zeta
      xf(20)=theta*geta*zeta
 

c Hot stars

      if(thet.le.0.85d0) then
         if(thet.lt.0.13d0) iflag=1
         if(g.ge.2.80d0) then
            call compindex(1,findex1,eindex1)
         end if
         if(g.le.3.0d0) then
            call compindex(2,findex2,eindex2)
         end if
         if(g.le.2.80d0) then
            findex3=findex2
            eindex3=eindex2
         else
            if(g.ge.3.0d0) then
               findex3=findex1
               eindex3=eindex1
            else
               x=(g-2.80d0)/0.20d0 !Interp en logg (2.8-3.0) hdw-hgi
               findex3=(1.d0-x)*findex2+x*findex1
               eindex3=(1.d0-x)*eindex2+x*eindex1
            end if
         end if

         if(thet.le.0.45d0.or.noz) then
            findex=findex3
            eindex=eindex3
            if(g.lt.1.d0.or.g.gt.4.2d0) iflag=3
            return
         end if
      end if

c Intermediate stars

      if(thet.le.1.35d0) then
         if(noz) then
            iflag=-2
            return
         end if
         call compindex(7,findex1,eindex1)
         if(thet.le.0.60d0)then
            findex=findex3
            eindex=eindex3
            return            
c         elseif(thet.le.0.70d0) then
         elseif(thet.le.0.80d0) then ! Interp en theta (0.80-0.60), hdw-i, hgi-i 
            x=cos(pi/2.d0*(0.80d0-thet)/0.20d0)
c            x=cos(pi/2.d0*(0.70d0-thet)/0.10d0)
            findex=(1.d0-x)*findex3+x*findex1
            eindex=(1.d0-x)*eindex3+x*eindex1
            if(g.lt.0.5d0.or.g.gt.4.6d0) iflag=3
            return
         end if
      end if

c 10   if(g.ge.3.d0) then   !10 label not used                   
      if(g.ge.3.d0) then                      
         if(g.lt.4.4d0.or.g.gt.5.2d0) iflag=3   
         call compindex(5,findex2,eindex2)  
         if(thet.le.1.07d0) then        
            findex=findex1                  
            eindex=eindex1                  
            return                          
         endif   
         if(thet.ge.1.27d0.or.noz) then        
            findex=findex2                  
            eindex=eindex2                  
            if(thet.gt.1.89d0) iflag=2        
            return                          
         else                               
            x=cos(pi/2.d0*(thet-1.07d0)/0.2d0) ! Interp en theta (1.27-1.07), i-cdw
            findex=(1.d0-x)*findex2+x*findex1 
            eindex=(1.d0-x)*eindex2+x*eindex1 
            return                          
         end if                             
      else
         if(g.gt.1.8d0.or.g.lt.-0.1d0) iflag=3  
         call compindex(4,findex2,eindex2)  
         call compindex(6,findex3,eindex3)  
         if(thet.le.1.28d0)then
            findex=findex1
            eindex=eindex1
            return
         else  ! mayor que theta 1.28d0
            if(thet.le.1.40d0)then
               if(thet.ge.1.35d0)then
                  findex=findex2                  
                  eindex=eindex2                  
                  return                          
               else
                  x=cos(pi/2.d0*(thet-1.28d0)/0.07d0) ! Interp en theta (1.35-1.28), i-cgi
                  findex=(1.d0-x)*findex2+x*findex1 
                  eindex=(1.d0-x)*eindex2+x*eindex1
                  return
               endif
            else  ! mayor que theta 1.40
               if(thet.ge.1.47d0) then        
                  findex=findex3                  
                  eindex=eindex3                  
                  if(thet.ge.1.70d0)then
                     findex=2.6327d0
                     eindex=0.d0
                  endif   
                  if(thet.gt.1.74d0) iflag=2        
                  return                          
               else
                  x=cos(pi/2.*(thet-1.40d0)/0.07d0) !Interp en theta (1.47-1.40),cgi-vcgi
                  findex=(1.d0-x)*findex3+x*findex2 
                  eindex=(1.d0-x)*eindex3+x*eindex2
                  return
               endif
            end if
         endif
      endif   
c Error: No value was computed
      iflag=-5
c      write(*,*)t,g,z,findex
c      write(*,'(A)') 'ERROR: NO index value was computed "finfitsTiO"'
c      stop

      end

c-------------------------------------------------------------------------------------------

      subroutine readcoefsTiO
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer i,j
      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm

c-Initialization----------------------------------------------------------
      do i=1,25
         fhd(i)=0.D0
         fhg(i)=0.D0
         fw(i)=0.D0
         fc(i)=0.D0
         fcd(i)=0.D0
         fcg(i)=0.D0
         fi(i)=0.D0
         lhd(i)=.false.
         lhg(i)=.false.
         lw(i)=.false.
         lc(i)=.false.
         lcd(i)=.false.
         lcg(i)=.false.
         li(i)=.false.
         do j=1,25
            chd(i,j)=0.D0
            chg(i,j)=0.D0
            cw(i,j)=0.D0
            cc(i,j)=0.D0
            ccd(i,j)=0.D0
            ccg(i,j)=0.D0
            ci(i,j)=0.D0
         end do
      end do
      lehd=.false.
      lehg=.false.
      lew=.false.
      lec=.false.
      lecd=.false.
      lecg=.false.
      lei=.false.
      ctehd=0.D0
      ctehg=0.D0
      ctew=0.D0
      ctec=0.D0
      ctecd=0.D0
      ctecg=0.D0
      ctei=0.D0
C########################################################################
c For hot dwarfs WEIGHTED
c N=35 sigma_res=0.01331  sigma_typ=0.00404 r**2=0.778 
C EXPONENCIAL
Constant:         0.78
Coef  1 :        -0.2261233080668275E+01        0.5934562331055201E-01
Coef  7 :         0.8061797366742040E+01        0.9155868953042732E+00
Coef 12 :        -0.1106751606048752E+02        0.1310667967660797E+01

      lehd=.true.
      ctehd=0.78

      srhd=0.01331
      fhd(1)= -0.2261233080668275E+01
      fhd(7)=  0.8061797366742040E+01
      fhd(12)=-0.1106751606048752E+02
      lhd(1)=.true.
      lhd(7)=.true.
      lhd(12)=.true.
      chd(1, 1)=  0.4833522211254657E+00
      chd(1, 2)= -0.6214702562058652E+01
      chd(1, 3)=  0.8062019936497267E+01
      chd(2, 2)=  0.1150496928220778E+03
      chd(2, 3)= -0.1629387461759924E+03
      chd(3, 3)=  0.2357608553282230E+03

c For hot dwarfs UNWEIGHTED !!
c N=35 sigma_res=0.01472  sigma_typ=0.00404 r**2=0.770 
C EXPONENCIAL
Constant:         0.78
Coef  1 :        -0.2333337208746753E+01        0.5788387993413698E-01
Coef  7 :         0.9237687253950606E+01        0.9101051733908466E+00
Coef 12 :        -0.1272060661697236E+02        0.1303611580425314E+01

C      lehd=.true.
C      ctehd=0.78
 
C      srhd=0.01472
C      fhd(1)= -0.2333337208746753E+01
C      fhd(7)=  0.9237687253950606E+01
C      fhd(12)=-0.1272060661697236E+02
C      lhd(1)=.true.
C      lhd(7)=.true.
C      lhd(12)=.true.
C      chd(1, 1)=  0.3529945983179774E+00
C      chd(1, 2)= -0.4626986781075922E+01
C      chd(1, 3)=  0.6028053723386295E+01
C      chd(2, 2)=  0.8726416909007017E+02
C      chd(2, 3)= -0.1237358666408604E+03
C      chd(3, 3)=  0.1790396463056568E+03

C########################################################################
c For hot giants WEIGHTED
c N=35 sigma_res=0.2994189736E-01 sigma_typ=0.00336 r**2=0
Coef  1 :         0.9550751486998874E+00        0.6603954246153728E-02

      srhg=0.2994189736E-01
      fhg(1)= 0.9550751486998874E+00
      lhg(1)=.true.
      chg(1, 1)=  0.4864626227774354E-01


c For hot giants UNWEIGHTED
c N=35 sigma_res=0.2933440419E-01 sigma_typ=0.00336 r**2=0
Coef  1 :         0.9461750015616417E+00        0.4638176552230728E-02

C      srhg=0.2933440419E-01
C      fhg(1)= 0.9461750015616417E+00
C      lhg(1)=.true.
C      chg(1, 1)=  0.2500000000000000E-01

C########################################################################
c For intermediate stars WEIGHTED!!
c N=569 sigma_res=0.1538049923E-01 sigma_typ=0.00386 r**2=0.754
Coef  1 :         0.1854752154516079E+01        0.2232798701662782E+00
Coef  2 :        -0.2617652104385464E+01        0.7003789023707058E+00
Coef  3 :        -0.7342855910128371E-01        0.1508437149588484E-01
Coef  4 :         0.3045724703471940E-01        0.5105520588759884E-02
Coef  7 :         0.2690895177483312E+01        0.7266804630860749E+00
Coef  8 :         0.1767549185343616E-01        0.6612570175225384E-02
Coef 10 :        -0.8611503421385915E-02        0.1667622777979342E-02
Coef 12 :        -0.8812754426127665E+00        0.2470761282616948E+00
Coef 13 :        -0.1650805226102687E-02        0.8573226599441148E-03

      sri=0.1538049923E-01
      fi(1)=  0.1854752154516079E+01
      fi(2)= -0.2617652104385464E+01
      fi(3)= -0.7342855910128371E-01
      fi(4)=  0.3045724703471940E-01
      fi(7)=  0.2690895177483312E+01
      fi(8)=  0.1767549185343616E-01
      fi(10)=-0.8611503421385915E-02
      fi(12)=-0.8812754426127665E+00
      fi(13)=-0.1650805226102687E-02 
      li(1)=.true.
      li(2)=.true.
      li(3)=.true.
      li(4)=.true.
      li(7)=.true.
      li(8)=.true.
      li(10)=.true.
      li(12)=.true.
      li(13)=.true.
      ci(1, 1)= 0.2107454842543804E+03
      ci(1, 2)=-0.6588945700436742E+03
      ci(1, 3)=-0.3298183230841380E+01
      ci(1, 4)= 0.6263280533023270E+00
      ci(1, 5)= 0.6790493049526253E+03
      ci(1, 6)= 0.1054066896597554E+01
      ci(1, 7)=-0.2301155291295764E+00
      ci(1, 8)=-0.2283720614350136E+03
      ci(1, 9)=-0.9034019994684314E-01
      ci(2, 2)= 0.2073601251172390E+04
      ci(2, 3)= 0.9033667313186347E+01
      ci(2, 4)=-0.1655771657237430E+01
      ci(2, 5)=-0.2146190071238452E+04
      ci(2, 6)=-0.2808596072355688E+01
      ci(2, 7)= 0.6027529973992993E+00
      ci(2, 8)= 0.7245335523013720E+03
      ci(2, 9)= 0.2219142078452992E+00
      ci(3, 3)= 0.9618637881226449E+00
      ci(3, 4)=-0.9552260493021560E-01
      ci(3, 5)=-0.9623400842706879E+01
      ci(3, 6)=-0.4046099946518070E+00
      ci(3, 7)= 0.2947402135085784E-01
      ci(3, 8)= 0.3264826726229500E+01
      ci(3, 9)= 0.4950713608274269E-01
      ci(4, 4)= 0.1101892429926264E+00
      ci(4, 5)= 0.1779091870780121E+01
      ci(4, 6)= 0.2579600318986813E-01
      ci(4, 7)=-0.3196982747820084E-01
      ci(4, 8)=-0.6259556886914209E+00
      ci(4, 9)=-0.2465226183675684E-02
      ci(5, 5)= 0.2232266820977214E+04
      ci(5, 6)= 0.2907929897031515E+01
      ci(5, 7)=-0.6089107396187290E+00
      ci(5, 8)=-0.7571695251574179E+03
      ci(5, 9)=-0.2212401857602275E+00
      ci(6, 6)= 0.1848416019035949E+00
      ci(6, 7)=-0.8668595062954606E-02
      ci(6, 8)=-0.9382284629956974E+00
      ci(6, 9)=-0.2365173928056231E-01
      ci(7, 7)= 0.1175586994060908E-01
      ci(7, 8)= 0.2022173178895492E+00
      ci(7, 9)= 0.9225570922270887E-03
      ci(8, 8)= 0.2580600102108933E+03
      ci(8, 9)= 0.6615423854313809E-01
      ci(9, 9)= 0.3107046415595799E-02


c For intermediate stars UNWEIGHTED!!
c N=569 sigma_res=0.1932650741E-01 sigma_typ=0.00386 r**2=0.633

Coef  1 :         0.2152788996975835E+01        0.3083736123365859E+00
Coef  2 :        -0.3494024108112372E+01        0.9272405111847989E+00
Coef  3 :        -0.7939999562902490E-01        0.1329595655874592E-01
Coef  4 :         0.9695894064332151E-01        0.1074152555956326E-01
Coef  7 :         0.3541275766556135E+01        0.9238616743747999E+00
Coef  8 :         0.1925137214729037E-01        0.5292716197242333E-02
Coef  9 :         0.2945354068348422E-01        0.5279149640234604E-02
Coef 10 :        -0.3903380122822667E-01        0.6607633608267910E-02
Coef 12 :        -0.1151599154492152E+01        0.3029075881682156E+00
Coef 13 :        -0.1759863520338429E-02        0.6462658760319349E-03
Coef 14 :         0.3084819003541224E-02        0.1355345137780355E-02
Coef 15 :         0.3648909493787331E-02        0.1036204931572107E-02
Coef 19 :        -0.5080082168192546E-02        0.1157415615007583E-02

c      sri=0.1932650741E-01
c      fi(1)=   0.2152788996975835E+01
c      fi(2)=  -0.3494024108112372E+01
c      fi(3)=  -0.7939999562902490E-01
c      fi(4)=   0.9695894064332151E-01
c      fi(7)=   0.3541275766556135E+01
c      fi(8)=   0.1925137214729037E-01
c      fi(9)=   0.2945354068348422E-01
c      fi(10)= -0.3903380122822667E-01
c      fi(12)= -0.1151599154492152E+01
c      fi(13)= -0.1759863520338429E-02
c      fi(14)=  0.3084819003541224E-02
c      fi(15)=  0.3648909493787331E-02
c      fi(19)= -0.5080082168192546E-02
c      li(1)=.true.
c      li(2)=.true.
c      li(3)=.true.
c      li(4)=.true.
c      li(7)=.true.
c      li(8)=.true.
c      li(9)=.true.
c      li(10)=.true.
c      li(12)=.true.
c      li(13)=.true.
c      li(14)=.true.
c      li(15)=.true.
c      li(19)=.true.
c      ci(1, 1)=  0.2545937051084736E+03
c      ci(1, 2)= -0.7644182595755190E+03
c      ci(1, 3)= -0.2824245756321125E+01
c      ci(1, 4)=  0.1724211864847978E+01
c      ci(1, 5)=  0.7590647147831573E+03
c      ci(1, 6)=  0.9699015085103523E+00
c      ci(1, 7)=  0.1832215880897246E+00
c      ci(1, 8)= -0.1341612806929590E+01
c      ci(1, 9)= -0.2474472165547955E+03
c      ci(1,10)= -0.9443421230568591E-01
c      ci(1,11)= -0.1717011813603724E-01
c      ci(1,12)=  0.1876256470793063E+00
c      ci(1,13)= -0.1068685504367764E+00
c      ci(2, 2)=  0.2301855411613938E+04
c      ci(2, 3)=  0.7931117103346655E+01
c      ci(2, 4)= -0.4820768361459741E+01
c      ci(2, 5)= -0.2290557121803216E+04
c      ci(2, 6)= -0.2717785818343277E+01
c      ci(2, 7)= -0.5102339459590599E+00
c      ci(2, 8)=  0.3782108219272097E+01
c      ci(2, 9)=  0.7481797522160847E+03
c      ci(2,10)=  0.2608656372627980E+00
c      ci(2,11)=  0.3936057861425905E-01
c      ci(2,12)= -0.5315703377825873E+00
c      ci(2,13)=  0.2899551334660131E+00
c      ci(3, 3)=  0.4732955487056909E+00
c      ci(3, 4)= -0.1584293808923516E+00
c      ci(3, 5)= -0.8210316302257914E+01
c      ci(3, 6)= -0.1817724958860764E+00
c      ci(3, 7)= -0.1551165657846210E-01
c      ci(3, 8)=  0.1158956486955563E+00
c      ci(3, 9)=  0.2753734761947737E+01
c      ci(3,10)=  0.2081675589173117E-01
c      ci(3,11)=  0.1087878091923635E-02
c      ci(3,12)= -0.1762248523762673E-01
c      ci(3,13)=  0.6414791445359164E-02
c      ci(4, 4)=  0.3089051703182414E+00
c      ci(4, 5)=  0.4781537641227523E+01
c      ci(4, 6)=  0.5043661644158664E-01
c      ci(4, 7)=  0.9941372029940587E-01
c      ci(4, 8)= -0.1632995640380851E+00
c      ci(4, 9)= -0.1539598834332882E+01
c      ci(4,10)= -0.4926606804482511E-02
c      ci(4,11)=  0.5701418764362049E-02
c      ci(4,12)=  0.2043954901658532E-01
c      ci(4,13)= -0.2260055554701840E-01
c      ci(5, 5)=  0.2285110191701045E+04
c      ci(5, 6)=  0.2788601726411324E+01
c      ci(5, 7)=  0.4684580898300869E+00
c      ci(5, 8)= -0.3746437480948889E+01
c      ci(5, 9)= -0.7482668853646976E+03
c      ci(5,10)= -0.2657218654146882E+00
c      ci(5,11)= -0.3746338824941547E-01
c      ci(5,12)=  0.5301652090034690E+00
c      ci(5,13)= -0.2693606816725717E+00
c      ci(6, 6)=  0.7499813422257430E-01
c      ci(6, 7)=  0.7507137657306033E-02
c      ci(6, 8)= -0.3659316239792770E-01
c      ci(6, 9)= -0.9219157893345988E+00
c      ci(6,10)= -0.9009103186517901E-02
c      ci(6,11)=  0.2694153814687723E-03
c      ci(6,12)=  0.5585535974912365E-02
c      ci(6,13)= -0.2181467102042634E-02
c      ci(7, 7)=  0.7461414898139920E-01
c      ci(7, 8)= -0.2602438923211626E-01
c      ci(7, 9)= -0.1344879299896889E+00
c      ci(7,10)= -0.1008290142149144E-02
c      ci(7,11)=  0.1402378792319312E-01
c      ci(7,12)=  0.2004469301156310E-02
c      ci(7,13)= -0.6579006512343740E-02
c      ci(8, 8)=  0.1168920975740739E+00
c      ci(8, 9)=  0.1205652468868790E+01
c      ci(8,10)=  0.3433898844373122E-02
c      ci(8,11)=  0.2773455747011985E-02
c      ci(8,12)= -0.1706405870235887E-01
c      ci(8,13)=  0.1146692269000566E-01
c      ci(9, 9)=  0.2456481801403160E+03
c      ci(9,10)=  0.8676006357427010E-01
c      ci(9,11)=  0.1227555035673118E-01
c      ci(9,12)= -0.1715279722113083E+00
c      ci(9,13)=  0.8082700250398131E-01
c      ci(10,10)= 0.1118190234339597E-02
c      ci(10,11)=-0.8584500347243533E-04
c      ci(10,12)=-0.5123535107978734E-03
c      ci(10,13)= 0.2365672117758820E-03
c      ci(11,11)= 0.4918051239881125E-02
c      ci(11,12)=-0.2124041760490929E-03
c      ci(11,13)= 0.1163693494778856E-02
c      ci(12,12)= 0.2874647217254517E-02
c      ci(12,13)=-0.8255924468323892E-03
c      ci(13,13)= 0.3586508954736126E-02

C#######################################################################
C For cold dwarfs WEIGHTED!!
c N=21 sigma_res=0.1707349865E-01 sigma_typ=0.00410 r**2=0.982
Coef  1 :         0.2454211706120712E+01        0.2958641819585296E+00
Coef  2 :        -0.2547328043299685E+01        0.4176511731406848E+00
Coef  7 :         0.1070011927384788E+01        0.1447020601722790E+00

      srcd=0.1707349865E-01
      fcd(1)= 0.2454211706120712E+01
      fcd(2)=-0.2547328043299685E+01
      fcd(7)= 0.1070011927384788E+01
      lcd(1)=.true.
      lcd(2)=.true.
      lcd(7)=.true.
      ccd(1, 1)=  0.3002892144416614E+03
      ccd(1, 2)= -0.4229166715330222E+03
      ccd(1, 3)=  0.1453787197936784E+03
      ccd(2, 2)=  0.5983872921392599E+03
      ccd(2, 3)= -0.2067251366622285E+03
      ccd(3, 3)=  0.7182975404674985E+02


C For cold dwarfs UNWEIGHTED!!
c     N=21 sigma_res=0.02388 sigma_typ=0.00410 r**2=0.969 UNWEIGHTED!!
Coef  1 :         0.1968153183997134E+01        0.2454741107084623E+00
Coef  2 :        -0.1859798835553561E+01        0.3413639088402062E+00
Coef  7 :         0.8333784313867646E+00        0.1167383821855962E+00

c      srcd=0.2387929993E-01
c      fcd(1)= 0.1968153183997134E+01
c      fcd(2)=-0.1859798835553561E+01
c      fcd(7)= 0.8333784313867646E+00
c      lcd(1)=.true.
c      lcd(2)=.true.
c      lcd(7)=.true.
c      ccd(1, 1)= 0.1056740153955372E+03 
c      ccd(1, 2)=-0.1465524968068300E+03 
c      ccd(1, 3)= 0.4972028179838071E+02 
c      ccd(2, 2)= 0.2043581794130855E+03 
c      ccd(2, 3)=-0.6969631448509452E+02 
c      ccd(3, 3)= 0.2389924382507220E+02 

C########################################################################
c     For cold giants WEIGHTED
c N=25 sigma_res=0.2357431340E-01 sigma_typ=0.00327 r**2=0.946
Coef  1 :         0.5996788996485434E+01        0.1924627780398984E+01
Coef  7 :        -0.9245555811031059E+01        0.3118596465408741E+01
Coef 12 :         0.4836582212691041E+01        0.1524354601184322E+01

      src=0.2357431340E-01
      fc(1)=  0.5996788996485434E+01
      fc(7)= -0.9245555811031059E+01
      fc(12)= 0.4836582212691041E+01
      lc(1)=.true.
      lc(7)=.true.
      lc(12)=.true.
      cc(1, 1)=  0.6665233903101890E+04
      cc(1, 2)= -0.1079374279359434E+05
      cc(1, 3)=  0.5271717343928391E+04
      cc(2, 2)=  0.1750008906571762E+05
      cc(2, 3)= -0.8552506505761716E+04
      cc(3, 3)=  0.4181132266938564E+04

c For cold giants UNWEIGHTED
c N=25 sigma_res=0.02890 sigma_typ=0.00327 r**2=0.925
Coef  1 :         0.7887374130908312E+01        0.1531281672755719E+01
Coef  7 :        -0.1227120331563940E+02        0.2475265541254331E+01
Coef 12 :         0.6304003980471593E+01        0.1208296460510297E+01

C      src=0.2889557795E-01
C      fc(1)=  0.7887374130908312E+01
C      fc(7)= -0.1227120331563940E+02
C      fc(12)= 0.6304003980471593E+01
C      lc(1)=.true.
C      lc(7)=.true.
C      lc(12)=.true.
C      cc(1, 1)=  0.2808325210506040E+04
C      cc(1, 2)= -0.4537319066245092E+04
C      cc(1, 3)=  0.2213436596485317E+04
C      cc(2, 2)=  0.7338052612642523E+04
C      cc(2, 3)= -0.3581567006516093E+04
C      cc(3, 3)=  0.1748574883596293E+04

C########################################################################
c For very cold giants WEIGHTED
c N=14 sigma_res=0.6465804641E-01 sigma_typ=0.00578 r**2=0.990
Coef  1 :        -0.3979267524044312E+02        0.5448400753716350E+01
Coef  2 :         0.3743776077340574E+02        0.5236015429403985E+01
Coef 12 :        -0.4318914544663060E+01        0.7093003230222535E+00

      srcg=0.6465804641E-01
      fcg(1)= -0.3979267524044312E+02
      fcg(2)=  0.3743776077340574E+02
      fcg(12)=-0.4318914544663060E+01
      lcg(1)=.true.
      lcg(2)=.true.
      lcg(12)=.true.
      ccg(1, 1)=  0.7100565392186860E+04
      ccg(1, 2)= -0.6821735995994344E+04
      ccg(1, 3)=  0.9216843432359631E+03
      ccg(2, 2)=  0.6557777510319650E+04
      ccg(2, 3)= -0.8871391975287052E+03
      ccg(3, 3)=  0.1203414272630102E+03

c For very cold giants UNWEIGHTED
c N=14 sigma_res=0.5818126097E-01 sigma_typ=0.00578 r**2=0.989
Coef  1 :        -0.1964314951261601E+02        0.1980546381488807E+01
Coef  7 :         0.2322684330482054E+02        0.2403107861974504E+01
Coef 12 :        -0.9129099147478389E+01        0.1013330888163838E+01

C      srcg=0.5818126097E-01
C      fcg(1)= -0.1964314951261601E+02
C      fcg(7)=  0.2322684330482054E+02
C      fcg(12)=-0.9129099147478389E+01
C      lcg(1)=.true.
C      lcg(7)=.true.
C      lcg(12)=.true.
C      ccg(1, 1)=  0.1158787430513617E+04
C      ccg(1, 2)= -0.1404788726111892E+04
C      ccg(1, 3)=  0.5917048566078530E+03
C      ccg(2, 2)=  0.1706004881357891E+04
C      ccg(2, 3)= -0.7192142261127576E+03
C      ccg(3, 3)=  0.3033446241373636E+03

      return
      end

c-----------------------------------------------------------------------------------------------
c ################## INTERPOLACION DE FUNCIONES Y COEFICIENTES: <MgI> ##########################
c-----------------------------------------------------------------------------------------------


      subroutine finfitMgI(t,g,z,findex,eindex,iflag)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
      integer iflag
      double precision t,g,z,findex,eindex,pi
      
      double precision findex1,findex2,findex3
      double precision eindex1,eindex2,eindex3
      double precision thet,x

      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      double precision theta,geta,zeta
      double precision xf(20)
      logical nog,noz
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm

C Regions:
C 1 = hot dwarfs (hd)
C 2 = hot giants (hg)
C 3 = cool stars (c)
C 4 = cold dwarfs (cd)
C 5 = cold giants (cg)
C 6 = intermediate region (i)
c Coefficients and variance-covariance matrices are stored in this subroutine
      call readcoefMgI

c-----------------------------------------------------------------------
      pi=3.14159d0

      nog=.false.
      noz=.false.
      if(g.ge.99.d0) nog=.true.
      if(z.ge.99.d0) noz=.true.

      iflag=0
      findex=0.d0
      eindex=0.d0
      if(t.le.0.d0) then
         iflag=-1
         return
      end if
      if(nog) then
         iflag=-3
         return
      end if

c Program works in theta (=5040/Teff)
      theta=dble(t)
      thet=t
      geta=dble(g)
      zeta=dble(z)
      xf(1)=1.D0
      xf(2)=theta
      xf(3)=geta
      xf(4)=zeta
      xf(5)=theta*zeta
      xf(6)=theta*geta
      xf(7)=theta*theta
      xf(8)=geta*geta
      xf(9)=zeta*zeta
      xf(10)=geta*zeta
      xf(11)=theta*theta*geta
      xf(12)=theta*theta*theta
      xf(13)=geta*geta*geta
      xf(14)=zeta*zeta*zeta
      xf(15)=geta*geta*zeta
      xf(16)=theta*theta*zeta
      xf(17)=theta*zeta*zeta
      xf(18)=theta*geta*geta
      xf(19)=geta*zeta*zeta
      xf(20)=theta*geta*zeta
 
c Hot stars
      if(thet.le.0.70d0) then
         if(thet.lt.0.13d0) iflag=1
         if(g.ge.2.80d0.and.thet.le.0.70d0) then
            call compindex(1,findex1,eindex1)
         end if
c         if(g.le.3.0d0.and.thet.le.0.6d0) then
         if(g.le.3.0d0.and.thet.le.0.7d0) then
            call compindex(2,findex2,eindex2)
         end if
c         if(g.le.2.80d0.and.thet.le.0.6d0) then
         if(g.le.2.80d0.and.thet.le.0.7d0) then
            findex3=findex2
            eindex3=eindex2
c         elseif(g.gt.3.0) then
         elseif(g.gt.3.0d0.and.thet.le.0.70d0) then
            findex3=findex1
            eindex3=eindex1
         elseif(g.gt.2.8d0.and.g.le.3.0d0.and.thet.le.0.6d0)then
            x=(g-2.80d0)/0.20d0
            findex3=(1.d0-x)*findex2+x*findex1
            eindex3=(1.d0-x)*eindex2+x*eindex1
         end if
c      end if
c If really hot or no z values, it keeps this index
c         if((g.ge.3.0.and.thet.le.0.4).or.(g.lt.3.0.and.thet.le.0.4)
c     >                                                .or.noz) then
         if((g.ge.3.0d0.and.thet.le.0.45d0).or.(g.lt.3.0d0.and.
     >      thet.le.0.45d0).or.noz) then
            findex=findex3
            eindex=eindex3
            if(g.lt.1.d0.or.g.gt.4.2d0) iflag=3
            return
         end if
      end if


c Cool stars
c      if(thet.le.1.40d0) then
      if(thet.le.1.60d0) then
         if(noz) then
            iflag=-2
            return
         end if
         call compindex(4,findex1,eindex1)
c         if(thet.le.0.75d0.and.thet.gt.0.4d0.and.g.ge.3.0d0)then
         if(thet.le.0.75d0.and.thet.gt.0.45d0.and.g.ge.3.0d0)then
c            x=cos(pi/2.d0*(0.75d0-thet)/0.35d0)
            x=cos(pi/2.d0*(0.75d0-thet)/0.3d0)
            findex=(1.d0-x)*findex3+x*findex1
            eindex=(1.d0-x)*eindex3+x*eindex1
            if(g.lt.0.0d0.or.g.gt.4.85d0) iflag=3
            return
c         elseif(thet.le.0.75.and.thet.gt.0.4.and.g.lt.3.0)then   
         elseif(thet.le.0.75d0.and.thet.gt.0.45d0.and.g.lt.3.0d0)then   
c            x=cos(pi/2.*(0.75-thet)/0.35)
            x=cos(pi/2.d0*(0.75d0-thet)/0.3d0)
            findex=(1.d0-x)*findex3+x*findex1
            eindex=(1.d0-x)*eindex3+x*eindex1
            if(g.lt.0.0d0.or.g.gt.4.85d0) iflag=3
            return
         end if
      end if

c Cold stars
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(thet.lt.1.95d0)then
         if(thet.le.1.60d0)call compindex(7,findex3,eindex3) !!
         if(thet.ge.1.20d0.and.g.le.3.d0)then
            call compindex(6,findex2,eindex2)
            if(g.gt.1.65d0) iflag=3
            if(thet.gt.1.60d0)then !OJO
               findex=findex2
               eindex=eindex2
               return
            else
               x=cos(pi/2.d0*(1.60d0-thet)/0.40d0)
               findex=(1.d0-x)*findex3+x*findex2
               eindex=(1.d0-x)*eindex3+x*eindex2
               return
            endif   
         elseif(thet.ge.1.10d0.and.g.gt.3.d0)then
            call compindex(5,findex2,eindex2)
            if(g.lt.4.4d0) iflag=3
            if(thet.gt.1.50d0)then
               findex=findex2
               eindex=eindex2
               return
            else
               x=cos(pi/2.d0*(1.50d0-thet)/0.4d0)
               findex=(1.d0-x)*findex3+x*findex2
               eindex=(1.d0-x)*eindex3+x*eindex2
               return
            endif   
         elseif(thet.ge.1.0d0)then
            if(g.gt.3.d0)then
               x=cos(pi/2.d0*(1.10d0-thet)/0.10d0)
               findex=(1.d0-x)*findex1+x*findex3
               eindex=(1.d0-x)*eindex1+x*eindex3
               return
            else
               x=cos(pi/2.d0*(1.20d0-thet)/0.20d0)
               findex=(1.d0-x)*findex1+x*findex3
               eindex=(1.d0-x)*eindex1+x*eindex3
               return
            endif   
         else      
            findex=findex1
            eindex=eindex1
            return
         endif   
      endif    

      iflag=-5
c      write(*,'(A)') 'ERROR: NO index value was computed "finfitMgI"'
c      write(*,*)t,g,z,findex
c      stop

      end
      
c-----------------------------------------------------------------------------------
      subroutine readcoefMgI
c     ----------------------
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer i,j

      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm

c-Initialization----------------------------------------------------------
      do i=1,25
         fhd(i)=0.D0
         fhg(i)=0.D0
c         fw(i)=0.D0
         fc(i)=0.D0
         fi(i)=0.D0
         fcd(i)=0.D0
         fcg(i)=0.D0
         lhd(i)=.false.
         lhg(i)=.false.
c         lw(i)=.false.
         lc(i)=.false.
         li(i)=.false.
         lcd(i)=.false.
         lcg(i)=.false.
         do j=1,25
            chd(i,j)=0.D0
            chg(i,j)=0.D0
c            cw(i,j)=0.D0
            cc(i,j)=0.D0
            ci(i,j)=0.D0
            ccd(i,j)=0.D0
            ccg(i,j)=0.D0
         end do
      end do
      lehd=.false.
      lehg=.false.
c      lew=.false.
      lec=.false.
      lei=.false.
      lecd=.false.
      lecg=.false.
      ctehd=0.D0
      ctehg=0.D0
c      ctew=0.D0
      ctec=0.D0
      ctei=0.D0
      ctecd=0.D0
      ctecg=0.D0
c---------------
c For hot dwarfs
c---------------
      lehd=.true.
      ctehd=-1.5D0

      fhd(1)=0.4283691026649478E+00
      fhd(7)=-0.4833159426443133E+01
      fhd(12)=0.7475618317255569E+01
      lhd(1)=.true.
      lhd(7)=.true.
      lhd(12)=.true.
c OJO!! LA srhd NO ES CORRECTA, ES LA EXPONENCIAL
      srhd=0.9350953826E-01
c
      chd( 1, 1)=  0.1001756208054846E+01
      chd( 1, 2)= -0.1102672772245672E+02
      chd( 1, 3)=  0.1337259875294098E+02
      chd( 2, 2)=  0.1477618491062455E+03
      chd( 2, 3)= -0.1907285416946117E+03
      chd( 3, 3)=  0.2518970841061548E+03
c---------------
c For hot giants
c---------------
      fhg(1)=-0.1201917746652227E+00
      fhg(12)=0.6253823261104153E+00
      lhg(1)=.true.
      lhg(12)=.true.
      srhg=0.5567235031E-01
      chg( 1, 1)=  0.4417233795926754E+00
      chg( 1, 2)= -0.2631472255668900E+01
      chg( 2, 2)=  0.2129317719407705E+02
c---------------
c For cool stars
c---------------

Coef  1 :         0.1845554436074967E+01        0.4230046702279171E+00
Coef  2 :        -0.7960169355385649E+01        0.1649539709370115E+01
Coef  3 :        -0.1998660066764755E+00        0.3489405155349934E-01
Coef  5 :         0.3079053263251620E+00        0.3096350620445296E-01
Coef  7 :         0.1167968392995355E+02        0.2023801463181931E+01
Coef  9 :         0.1861894114599228E+00        0.6001181272385286E-01
Coef 11 :        -0.6817320165691089E-01        0.3318732810675449E-01
Coef 12 :        -0.4570981793624912E+01        0.7597292136727547E+00
Coef 17 :        -0.1624176015999102E+00        0.6046245985145578E-01
Coef 18 :         0.5072432272157244E-01        0.5677477144819303E-02

      fc(1)=  0.1845554436074967E+01
      fc(2)= -0.7960169355385649E+01
      fc(3)= -0.1998660066764755E+00
      fc(5)=  0.3079053263251620E+00
      fc(7)=  0.1167968392995355E+02
      fc(9)=  0.1861894114599228E+00
      fc(11)=-0.6817320165691089E-01
      fc(12)=-0.4570981793624912E+01
      fc(17)=-0.1624176015999102E+00
      fc(18)= 0.5072432272157244E-01
      lc(1)=.true. 
      lc(2)=.true. 
      lc(3)=.true. 
      lc(5)=.true. 
      lc(7)=.true. 
      lc(9)=.true. 
      lc(11)=.true.
      lc(12)=.true.
      lc(17)=.true.
      lc(18)=.true.
      src=0.9410767921E-01
      cc( 1, 1)=  0.2020413202110936E+02
      cc( 1, 2)= -0.7458261524987259E+02
      cc( 1, 3)=  0.4231479350000749E-01
      cc( 1, 4)=  0.8440331779995956E-01
      cc( 1, 5)=  0.8613751800580114E+02
      cc( 1, 6)=  0.3941306503294438E+00
      cc( 1, 7)= -0.5650646973667980E+00
      cc( 1, 8)= -0.3107014920603885E+02
      cc( 1, 9)= -0.4016444495721797E+00
      cc( 1,10)=  0.9821433410656737E-01
      cc( 2, 2)=  0.3072383490040916E+03
      cc( 2, 3)= -0.1872867798093312E+01
      cc( 2, 4)= -0.4123531493213390E+00
      cc( 2, 5)= -0.3724108063470849E+03
      cc( 2, 6)= -0.2005001651823398E+01
      cc( 2, 7)=  0.3525213025453493E+01
      cc( 2, 8)=  0.1373461461584094E+03
      cc( 2, 9)=  0.2020532480200087E+01
      cc( 2,10)= -0.2985776092312173E+00
      cc( 3, 3)=  0.1374841616842507E+00
      cc( 3, 4)=  0.6785867186776715E-03
      cc( 3, 5)=  0.2865265956906230E+01
      cc( 3, 6)= -0.4210064894573011E-02
      cc( 3, 7)= -0.7269466577958775E-01
      cc( 3, 8)= -0.1114655987519009E+01
      cc( 3, 9)=  0.5049378388627425E-02
      cc( 3,10)= -0.1209627251953279E-01
      cc( 4, 4)=  0.1082555420305153E+00
      cc( 4, 5)=  0.5459673518614844E+00
      cc( 4, 6)=  0.1025432364889745E-01
      cc( 4, 7)= -0.7965719238705232E-02
      cc( 4, 8)= -0.1996921223205734E+00
      cc( 4, 9)=  0.4097971357613309E-01
      cc( 4,10)=  0.1776924361834383E-02
      cc( 5, 5)=  0.4624722560008947E+03
      cc( 5, 6)=  0.2762042901343269E+01
      cc( 5, 7)= -0.4952032353258934E+01
      cc( 5, 8)= -0.1728315913454728E+03
      cc( 5, 9)= -0.2798722362625522E+01
      cc( 5,10)=  0.3522043578097100E+00
      cc( 6, 6)=  0.4066524224544624E+00
      cc( 6, 7)= -0.5722985234051747E-01
      cc( 6, 8)= -0.1059879973013045E+01
      cc( 6, 9)= -0.3940367306839343E+00
      cc( 6,10)=  0.9417051216469793E-02
      cc( 7, 7)=  0.1243639338695414E+00
      cc( 7, 8)=  0.1911563045043996E+01
      cc( 7, 9)=  0.6015611311220458E-01
      cc( 7,10)= -0.7590951194486823E-02
      cc( 8, 8)=  0.6517297202848211E+02
      cc( 8, 9)=  0.1083406026981396E+01
      cc( 8,10)= -0.1278626857606091E+00
      cc( 9, 9)=  0.4127827092546594E+00
      cc( 9,10)= -0.9657917635096921E-02
      cc(10,10)=  0.3639658714066504E-02

c For intermediate stars (1.05<theta<1.35) to perform a smooth iterpolation.
Coef  1 :        -0.9004806321043938E+01        0.3562934132919736E+01
Coef  2 :         0.1675975620113735E+02        0.6002674363568172E+01
Coef  5 :         0.1710887478998283E+01        0.6407403839066386E+00
Coef  6 :        -0.1261752802519696E+00        0.3999140313888341E-01
Coef  7 :        -0.6940718645584071E+01        0.2515696314481113E+01
Coef  9 :         0.8445398657025406E+00        0.4527176953396914E+00
Coef 13 :         0.7010998269749679E-02        0.1746505460553927E-02
Coef 16 :        -0.1246620331728085E+01        0.5363575661682110E+00
Coef 17 :        -0.7490548803177109E+00        0.3923877673952887E+00

      sri=0.9931482624E-01

      fi(1)= -0.9004806321043938E+01
      fi(2)=  0.1675975620113735E+02
      fi(5)=  0.1710887478998283E+01
      fi(6)= -0.1261752802519696E+00
      fi(7)= -0.6940718645584071E+01
      fi(9)=  0.8445398657025406E+00
      fi(13)= 0.7010998269749679E-02
      fi(16)=-0.1246620331728085E+01
      fi(17)=-0.7490548803177109E+00
      li(1)=.true.
      li(2)=.true.
      li(5)=.true.
      li(6)=.true.
      li(7)=.true.
      li(9)=.true.
      li(13)=.true.
      li(16)=.true.
      li(17)=.true.

      ci(1, 1)= 0.1287026275045267E+04
      ci(1, 2)=-0.2166195432519061E+04
      ci(1, 3)= 0.5153159395523907E+02
      ci(1, 4)= 0.5168665686169808E+00
      ci(1, 5)= 0.9062541780189847E+03
      ci(1, 6)= 0.1854125834223293E+02
      ci(1, 7)=-0.1133942183849994E+00
      ci(1, 8)=-0.4570886992625751E+02
      ci(1, 9)=-0.1631314219714071E+02
      ci(2, 2)= 0.3653098598004281E+04
      ci(2, 3)=-0.8287670315135324E+02
      ci(2, 4)=-0.1673946264286967E+01
      ci(2, 5)=-0.1530239396360016E+04
      ci(2, 6)=-0.2989062092981035E+02
      ci(2, 7)= 0.2196763136873048E+00
      ci(2, 8)= 0.7385019426242383E+02
      ci(2, 9)= 0.2626626508662325E+02
      ci(3, 3)= 0.4162325311507458E+02
      ci(3, 4)=-0.1161218763729956E+00
      ci(3, 5)= 0.3340329315575608E+02
      ci(3, 6)= 0.2370091819190702E+02
      ci(3, 7)=-0.3221771335988157E-02
      ci(3, 8)=-0.3473059252967840E+02
      ci(3, 9)=-0.1998361658368487E+02
      ci(4, 4)= 0.1621455782739226E+00
      ci(4, 5)= 0.8277376029942489E+00
      ci(4, 6)=-0.1389428133124717E+00
      ci(4, 7)=-0.6473951002104493E-02
      ci(4, 8)= 0.8901936656314548E-01
      ci(4, 9)= 0.1488106660712493E+00
      ci(5, 5)= 0.6416353057519885E+03
      ci(5, 6)= 0.1226906839998515E+02
      ci(5, 7)=-0.9603637112809334E-01
      ci(5, 8)=-0.2987430743147666E+02
      ci(5, 9)=-0.1079619039717306E+02
      ci(6, 6)= 0.2077910156800531E+02
      ci(6, 7)=-0.1230632786784978E-02
      ci(6, 8)=-0.2006568299609760E+02
      ci(6, 9)=-0.1795888235498498E+02
      ci(7, 7)= 0.3092514335069387E-03
      ci(7, 8)= 0.3730070999790371E-02
      ci(7, 9)= 0.1798440263167686E-03
      ci(8, 8)= 0.2916625366403541E+02
      ci(8, 9)= 0.1702820027649301E+02
      ci(9, 9)= 0.1560999433862792E+02

c----------------
c For cold dwarfs
c----------------

Coef  1 :         0.2638412678769896E+01        0.4694311889302759E+00
Coef  7 :        -0.2098307284274805E+01        0.6781919821977767E+00
Coef 12 :         0.8046143463274884E+00        0.3055404148450539E+00

      fcd(1)=  0.2638412678769896E+01
      fcd(7)= -0.2098307284274805E+01
      fcd(12)= 0.8046143463274884E+00
      lcd(1)=.true.
      lcd(7)=.true.
      lcd(12)=.true.
      srcd=0.8466422136E-01
      ccd( 1, 1)=  0.3074284420232750E+02
      ccd( 1, 2)= -0.4397261885034377E+02
      ccd( 1, 3)=  0.1960902850096256E+02
      ccd( 2, 2)=  0.6416607359052465E+02
      ccd( 2, 3)= -0.2885059522316201E+02
      ccd( 3, 3)=  0.1302379317386308E+02

c----------------
c For cold giants
c----------------

Coef  1 :         0.1896394087207532E+01        0.2540999469750470E+00
Coef  7 :        -0.5761852181356737E+00        0.1161344667339780E+00

      fcg(1)= 0.1896394087207532E+01
      fcg(7)=-0.5761852181356737E+00
      lcg(1)=.true.
      lcg(7)=.true.
      srcg=0.1571020172E+00
      ccg( 1, 1)= 0.2616047504753017E+01
      ccg( 1, 2)=-0.1174075465163657E+01
      ccg( 2, 2)= 0.5464604524762594E+00

      return
      end


c-----------------------------------------------------------------------------------------------
c ################################# SUBRUTINAS COMUNES #########################################
c-----------------------------------------------------------------------------------------------

      subroutine compindex(it,fin,ein)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer i,j,it,nter
      double precision fin,ein

      double precision fhd(25),fhg(25),fw(25),fc(25),fcd(25),fcg(25),
     c  fi(25),flm(25)
      double precision chd(25,25),chg(25,25),cw(25,25),cc(25,25)
      double precision ccd(25,25),ccg(25,25),ci(25,25),clm(25,25)
      double precision srhd,srhg,srw,src,srcd,srcg,sri,srlm
      double precision ctehd,ctehg,ctew,ctec,ctecd,ctecg,cte,ctei,ctelm
      double precision xf(20),x(25)
      double precision f(25),c(25,25)
      double precision sr,indt,derr
      logical lhd(25),lhg(25),lw(25),lc(25),lcd(25),lcg(25),li(25),
     c  llm(25)
      logical l(25)
      logical lehd,lehg,lew,lec,lecd,lecg,le,lei,lelm
      common/blkl/lhd,lhg,lw,lc,lcd,lcg,li,llm
      common/blkle/lehd,lehg,lew,lec,lecd,lecg,lei,lelm
      common/blkf/fhd,fhg,fw,fc,fcd,fcg,fi,flm
      common/blkc/chd,chg,cw,cc,ccd,ccg,ci,clm
      common/blksr/srhd,srhg,srw,src,srcd,srcg,sri,srlm
      common/blkcte/ctehd,ctehg,ctew,ctec,ctecd,ctecg,ctei,ctelm
      common/blkx/xf

      if(it.eq.1) then
         sr=srhd
         le=lehd
         cte=ctehd
         do i=1,25
            f(i)=fhd(i)
            l(i)=lhd(i)
            do j=1,25
               c(i,j)=chd(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.2) then
         sr=srhg
         le=lehg
         cte=ctehg
         do i=1,25
            f(i)=fhg(i)
            l(i)=lhg(i)
            do j=1,25
               c(i,j)=chg(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.3) then
         sr=srw
         le=lew
         cte=ctew
         do i=1,25
            f(i)=fw(i)
            l(i)=lw(i)
            do j=1,25
               c(i,j)=cw(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.4) then
         sr=src
         le=lec
         cte=ctec
         do i=1,25
            f(i)=fc(i)
            l(i)=lc(i)
            do j=1,25
               c(i,j)=cc(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.5) then
         sr=srcd
         le=lecd
         cte=ctecd
         do i=1,25
            f(i)=fcd(i)
            l(i)=lcd(i)
            do j=1,25
               c(i,j)=ccd(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.6) then
         sr=srcg
         le=lecg
         cte=ctecg
         do i=1,25
            f(i)=fcg(i)
            l(i)=lcg(i)
            do j=1,25
               c(i,j)=ccg(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.7) then
         sr=sri
         le=lei
         cte=ctei
         do i=1,25
            f(i)=fi(i)
            l(i)=li(i)
            do j=1,25
               c(i,j)=ci(i,j)
            end do
         end do
         goto 10
      end if
      if(it.eq.8) then
         sr=srlm
         le=lelm
         cte=ctelm
         do i=1,25
            f(i)=flm(i)
            l(i)=llm(i)
            do j=1,25
               c(i,j)=clm(i,j)
            end do
         end do
         goto 10
      end if

c Now it computes the index
 10   indt=0.D0
      j=0
      do i=1,25
         if(l(i)) then
            indt=indt+f(i)*xf(i)
            j=j+1
            x(j)=xf(i)
         end if
      end do
      nter=j
      if(.not.le) then
         fin=dble(indt)
      else
         fin=dble(cte+dexp(indt))
      end if
c and the error
      call comperr(nter,x,c,derr)
      ein=dble(sr*dsqrt(derr))
      if(le) ein=ein*(fin-dble(cte))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comperr(n,x,c,v)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer n
      double precision x(25),v,c(25,25)

      integer i,j
      double precision x1(25,25),x2(25,25),temp(25,25),fin(25,25)

c First it reconstructs the whole v-c matrix
      
      do i=2,n
         do j=1,i-1
            c(i,j)=c(j,i)
         end do
      end do

      do j=1,n
         x1(1,j)=x(j)
         x2(j,1)=x(j)
      end do

      call multmatrix(x1,c,temp,1,n,n,n)

      call multmatrix(temp,x2,fin,1,n,n,1)
      v=fin(1,1)

      return
      end

c Subroutine to multiply matrices
c Input: a -> matrix n1 x n2
c        b -> matrix n3 x n4
c n2 must be equal to n3
c Input : n1,n2,n3,n4 (integers, dimensions of matrices)
c Output: c -> matrix n1 x n4
c In the main program matrices must be defined as double precision with 
c   dimensions 25x25
 
      subroutine multmatrix(a,b,c,n1,n2,n3,n4)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none

      integer n1,n2,n3,n4
      integer i,j,k
      double precision a(25,25),b(25,25),c(25,25)

      if(n2.ne.n3) then
         write(6,*) 'Error in matrix dimensions'
         stop
      end if
      
      do i=1,n1
         do j=1,n4
            c(i,j)=0.D0
            do k=1,n2
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do

      return
      end

