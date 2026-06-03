      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nstm=1999)
      PARAMETER (npxm=4300)
      DIMENSION x(npxm,2) !MILES STARS ARRAY
      DIMENSION s(nstm),a(nstm) !miles
      CHARACTER*5 s
      CHARACTER*6 a
      CHARACTER*1 vname
      V1=4750.0d0
c      V1=4875.0d0 ;10%
      V2=7400.0d0
c      V2=6375.0d0 ;10%
      vname='V'
c LECTURA DE FICHEROS ASCII DE LIBRERIAS ESPECTRALES
c SE CARGA EN MATRIZ:starms(nstarm,npxm,2)
       OPEN(98,FILE='lista_to_V',STATUS='OLD')
       n=0
241    format(A5)
       do k=1,nstm
         read(98,241,end=23)s(k)
         a(k)=s(k)//vname
         n=n+1
       enddo
23     CLOSE(98)
c
188    format(F6.1,1X,F20.16)
       DO i=1,n !MILES
	 open(99,file=s(i),status='old')
         fvmil=0.0d0
      	 rfV=0.0d0
	 do l=1,npxm
	       read(99,*)x(l,1),x(l,2)
      	       rfV=0.0d0
      	       if(x(l,1).ge.V1.and.x(l,1).le.V2)then
      	         if(l.eq.1.or.l.eq.npxm)then
		   dl=0.9d0
		 else     	         
		   dl=abs(x(l,1)-x(l-1,1))
		 endif
      	         call respv(x(l,1),rfV)
                 fvmil=fvmil+rfV*x(l,2)*dl
c		 write(*,*)x(l,1),x(l,2),dl,rfV,fvmil
      	       endif
         enddo
	 close(99)
	 open(97,file=a(i),status='new')
         do l=1,npxm
       	    write(97,188)x(l,1),x(l,2)/fvmil
          enddo
	 close(97)
      ENDDO
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE respv(xl,s)
c Subroutine to calculate the response of the V Johnson filter
c It calls the "hunt" and "polint" routines of the numerical recipes 
c Filter in use: Buser & Kurucz 1978 (A&A,70,555)
c WARNING: the provided wavelength should be double-precision
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nf=54)
      DIMENSION XV(nf),SV(nf)
      data XV/4750.d0,4800.d0,4850.d0,4900.d0,4950.d0,5000.d0,
     &5050.d0,5100.d0,5150.d0,5200.d0,5250.d0,5300.d0,5350.d0,5400.d0,
     &5450.d0,5500.d0,5550.d0,5600.d0,5650.d0,5700.d0,5750.d0,5800.d0,
     &5850.d0,5900.d0,5950.d0,6000.d0,6050.d0,6100.d0,6150.d0,6200.d0,
     &6250.d0,6300.d0,6350.d0,6400.d0,6450.d0,6500.d0,6550.d0,6600.d0,
     &6650.d0,6700.d0,6750.d0,6800.d0,6850.d0,6900.d0,6950.d0,7000.d0,
     &7050.d0,7100.d0,7150.d0,7200.d0,7250.d0,7300.d0,7350.d0,7400.d0/
      data SV/0.0d0,0.0299999993,0.0839999989,0.163000003,0.300999999,
     &0.458000004,0.629999995,0.779999971,0.894999981,0.967000008,
     &0.996999979,1.000000000,0.987999976,0.958000004,0.919000030,
     &0.876999974,0.819000006,0.764999986,0.711000025,0.657000005,
     &0.601999998,0.545000017,0.488000005,0.433999985,0.386000007,
     &0.331000000,0.289000005,0.250000000,0.214000002,0.180999994,
     &0.150999993,0.119999997,0.0930000022,0.0689999983,0.050999999,
     &0.0359999985,0.0270000007,0.0209999997,0.0179999992,
     &0.0160000008,0.0140000004,0.0120000001,0.0109999999,
     &0.00999999978,0.00899999961,0.00800000038,0.00700000022,
     &0.00600000005,0.00499999989,0.00400000019,0.00300000003,
     &0.00200000009,0.00100000005,0.0d0/      
      xli=XV(1)
      xlf=XV(nf)
      call hunt(XV,nf,xl,k1)
      if(xl.gt.xli.and.xl.lt.xlf)then
	   call polint(XV(k1-1),SV(k1-1),3,xl,s,dy)
	   if(s.gt.1.0d0) s=1.0d0
	   if(s.lt.0.0d0) s=0.0d0
      else
	   s=0.0d0
      endif
      RETURN
      END
c      data XV/4700.,4750.,4800.,4850.,4900.,4950.,5000.,5050.,5100.,
c     &5150.,5200.,5250.,5300.,5350.,5400.,5450.,5500.,5550.,5600.,
c     &5650.,5700.,5750.,5800.,5850.,5900.,5950.,6000.,6050.,6100.,
c     &6150.,6200.,6250.,6300.,6350.,6400.,6450.,6500.,6550.,6600.,
c     &6650.,6700.,6750.,6800.,6850.,6900.,6950.,7000.,7050.,7100.,
c     &7150.,7200.,7250.,7300.,7350.,7400.,7450./
c      data SV/.0,.0,.03,.084,.163,.301,.458,.630,.780,.895,.967,.997,1.
c     &,.988,.958,.919,.877,.819,.765,.711,.657,.602,.545,.488,.434,.386
c     &,.331,.289,.250,.214,.181,.151,.120,.093,.069,.051,.036,.027,.021
c     &,.018,.016,.014,.012,.011,.010,.009,.008,.007,.006,.005,.004,.003
c     &,.002,.001,.0,.0/
c     xli=4750.0d0
c	xlf=7400.0d0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hunt(xx,n,x,jlo)
C subroutine from the Numerical Recipes
      INTEGER jlo,n
      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE polint(xa,ya,n,x,y,dy)
C    Numerical Recipes subroutine
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)STOP 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
