*******SUBRUTINA D4000************************************************
c        SUBROUTINE d4000(nra,i,W,stars)
      SUBROUTINE d4000(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION feh(15)
      COMMON/spec/R(12,2999,15),jotai(15)
      COMMON/fezsol/fez(15)
      COMMON/rat/ra(3,50)
      common/cd4000/fra4(2999),fbj(2999)
      COMMON/MAMA/R140(12,99000),imam(12)
      feh(i)=fez(i)
c        feh(i)=fez(i)+dlog10(ra(nra,15)/ra(1,15))
      spend=(4000.d0-3600.d0)/(4400.d0-3600.d0)
      W=0.0
      stars=0.0
      if(jotai(i).eq.0)goto 40497
      DO j=1,jotai(i)
	at=10**R(i,j,3)
c	 tetas=5040.d0/at
	ggrav=R(i,j,4)
	zfeff=feh(i)
        call d4000e(at,ggrav,zfeff,ZMAG,eindx,iflage)
	W=W+ZMAG*R(i,j,14)*fra4(j)*fbj(j)
	stars=stars+R(i,j,14)*fra4(j)*fbj(j)
        if(j.eq.1)then
	  if(imam(i).gt.0)then
	   do mma=1,imam(i)
	     W=W+ZMAG*R140(i,mma)*fra4(j)*fbj(j)
	     stars=stars+R140(i,mma)*fra4(j)*fbj(j)
	   enddo
	  endif
	 W=W+ZMAG*R(i,j,14)*fra4(j)*fbj(j)
	 stars=stars+R(i,j,14)*fra4(j)*fbj(j)
	else
	 W=W+ZMAG*R(i,j,14)*fra4(j)*fbj(j)
	 stars=stars+R(i,j,14)*fra4(j)*fbj(j)
	endif
      ENDDO
40497 continue
      RETURN
      END

c*******SUBRUTINA del D4000e********************************************
      SUBROUTINE d4000e(t,g,z,fxex,eixex,iflag)

c Version 11/12/98
C-----------------------------------------------------------------------------
c Program to evaluate the D4000 index using the empirical fitting functions
c from Gorgas, Cardiel, Pedraz, Gonzalez (1999, A&A Suppl., in preparation)
c
c INPUT:
c       t = effective temperature, in K                      (REAL)
c       g = logarithm (base 10) of the surface gravity (cgs) (REAL)
c       z = metallicity ([Fe/H])                             (REAL)
c  (A value higher or equal to 99 in g or z in input means that this 
c    parameter is unknown. In some cases, an index value can still be 
c    computed. T<=0 means that temperature is unknown) 
c  
c
c OUTPUT:
c       fxex = index value                                 (REAL)
c       eixex = error in the index value                    (REAL)
c       iflag = indicates sucess of the functions:           (INTEGER)
c              0 -> OK
c              1 -> extrapolation to high temperatures
c              2 -> extrapolation to low temperatures
c              3 -> extrapolation in g (uncertain value)
c             -1 -> Error (no Teff)
c             -2 -> Error (z needed)
c             -3 -> Error (g needed)
c             -5 -> Error
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
c      real t,g,z,fxex,eixex
      INTEGER i,j,iflag
c      real thet,x,fxex0,eixex0
c      real ctcg,ctcd,ctd,errt,errw
      DIMENSION xh(25),ccg(25,25),ccd(25,25),cw(25,25),ch(25,25)
      DIMENSION cg(25),cd(25,25),xh2(25),fcg(25),fcd(25),fw(25)
      DIMENSION fh(25),fg(25),fd(25),fi(25),ci(25,25)
c      double precision theta,thetaz,theta2,thet2z,xh(25)
c      double precision ccg(25,25),ccd(25,25),cw(25,25),ch(25,25),cg(25)
c      double precision cd(25,25),dumt,xh2(25)
c      double precision fcg(25),fcd(25),fw(25),fh(25),fg(25),fd(25)
c      double precision srcg,srcd,vh,srw,srh,xfh,dx,srg,xfg,srfg,srd,fw1
c      double precision sri,fi(25),ci(25,25)
      LOGICAL nog,noz
c-----------------------------------------------------------------------
c Coefficients and Variance-covariance matrices:

c For hot stars
      srh=0.4615183238D-01
      fh(1)=0.6548051838D+00
      fh(2)=0.1340106583D+01
      ch(1,1)=0.2313283739D+00
      ch(1,2)=-0.4761853888D+00
      ch(2,2)=0.1142256111D+01

c For warm stars
      srw=0.3422608269D-01
      fw(1)=0.1822983515D+01
      fw(2)=-0.5068388040D+00
      cw(1,1)=0.1663359277D+02
      cw(1,2)=-0.2414081677D+02
      cw(2,2)=0.3519089299D+02
      errw=0.159526188D-01
      fw1= 1.4428544D+00
c For hot supergiants
      sri=0.7879815732D-01
      fi(1)=0.8613232907D+00
      fi(2)=0.1848690037D+01
      ci(1,1)=0.1731426778D+00
      ci(1,2)=-0.8470102471D+00
      ci(2,2)=0.7135431404D+01

c For cool giants:
      ctcg=0.9
      srcg=0.1309941549D+00
      fcg(1)=-0.5664970674D+01
      fcg(2)=0.9278637462D+01
      fcg(4)=-0.3273398709D+01
      fcg(5)=0.7322437356D+01
      fcg(7)=-0.3080110968D+01
      fcg(16)=-0.3694125492D+01
      data ((ccg(i,j),j=1,6),i=1,6)/0.1424049628D+03,-0.2516045719D+03,
     c 0.8616301985D+02,-0.1533506540D+03,0.1102731351D+03,
     c 0.6779636477D+02,0.D0,0.4460123942D+03,-0.1534799668D+03,
     c 0.2742738641D+03,-0.1961177449D+03,-0.1217490895D+03,0.D0,0.D0,
     c 0.3013629553D+03,-0.5639211958D+03,0.6795204379D+02,
     c 0.2608696279D+03,0.D0,0.D0,0.D0,0.1060365074D+04,
     c -0.1219393397D+03,-0.4926576419D+03,0.D0,0.D0,0.D0,0.D0,
     c 0.8652160543D+02,0.5435617680D+02,0.D0,0.D0,0.D0,0.D0,0.D0,
     c 0.2298069839D+03/

c For cool dwarfs:
      ctcd=0.9
      srcd=0.1647942839D+00
      fcd(1)=-0.8153624600D+01
      fcd(2)=0.1344994016D+02
      fcd(4)=-0.1706301430D+02
      fcd(5)=0.3811728640D+02
      fcd(7)=-0.4768742248D+01
      fcd(16)=-0.2068889239D+02
      data ((ccd(i,j),j=1,6),i=1,6)/0.1482840358D+03,-0.3158204285D+03,
     c 0.1389514149D+02,-0.2538831761D+02,0.1668009227D+03,
     c 0.1078900852D+02,0.D0,0.6745000058D+03,-0.3087237636D+02,
     c 0.5813276884D+02,-0.3571841966D+03,-0.2581987859D+02,0.D0,0.D0,
     c 0.1950227161D+04,-0.4230766254D+04,0.1667832891D+02,
     c 0.2278759127D+04,0.D0,0.D0,0.D0,0.9188955680D+04,
     c -0.3219517176D+02,-0.4955197523D+04,0.D0,0.D0,0.D0,0.D0,
     c 0.1896456546D+03,0.1480236474D+02,0.D0,0.D0,0.D0,0.D0,0.D0,
     c 0.2675388494D+04/

c For cold giants
      srg=0.2107104366D+00
      fg(1)=0.9524702185D+01
      fg(2)=-0.4093622777D+01
      cg(1)=0.1861588935D+01
      cg(2)=0.2986724620D+01
      xfg=0.1301180D+01
      srfg=0.1310919D+00

c For cold dwarfs
      ctd=1.15
      srd=0.1868498444D+00
      fd(1)=-0.2908460064D+01
      fd(2)= 0.6534652942D+01
      fd(7)=-0.2975879712D+01
      data ((cd(i,j),j=1,3),i=1,3)/0.2490535985D+03,-0.3573457710D+03,
     c 0.1258060822D+03,0.D0,0.5155536330D+03,-0.1824972759D+03,0.D0,
     c 0.D0,0.6497403497D+02/

c-----------------------------------------------------------------------

      nog=.false.
      noz=.false.
      if(g.ge.99.d0) nog=.true.
      if(z.ge.99.d0) noz=.true.

      iflag=0
      fxex=0.d0
      eixex=0.d0
      if(t.le.0.d0) then
         iflag=-1
         return
      end if

c Program works in theta (=5040/Teff)
      theta=5040.D0/dble(t)
      thet=dble(theta)
c Now, it prepares coefficients for error evaluation
      thetaz=theta*dble(z)
      theta2=theta*theta
      thet2z=theta2*dble(z)
      xh(1)=1.D0
      xh(2)=theta
      xh(3)=dble(z)
      xh(4)=thetaz
      xh(5)=theta2
      xh(6)=thet2z
      dumt=0.8D+00
      xh2(1)=1.D0
      xh2(2)=dumt
      xh2(3)=dble(z)
      xh2(4)=dumt*dble(z)
      xh2(5)=dumt*dumt
      xh2(6)=dumt*dumt*dble(z)
 
c Hot supergiants
      if(thet.le.0.7d0) then
         if(g.lt.(3.456933d0-1.221595d0*thet)) then
            if(t.lt.0.19d0) iflag=1
            fxex=dble(fi(1)+fi(2)*theta*theta*theta)
            xh(2)=theta*theta*theta
            call copr4(2,xh,ci,vh)
            eixex=dble(sri*dsqrt(vh))
            return
         end if
      end if

c Theta=(0.127,0.6326)
c Not known dependence on g or z
c Even is theta is lower, the program extrapolates.
      if(thet.le.0.6326d0) then
         if(t.lt.0.127d0) iflag=1
         fxex=dble(fh(1)+fh(2)*theta)
         call copr4(2,xh,ch,vh)
         eixex=dble(srh*dsqrt(vh))
         return
      end if
      
c Theta=(0.6326,0.75)
      if(thet.le.0.75d0) then
         fxex=dble(fw(1)+fw(2)*theta)
         call copr4(2,xh,cw,vh)
         eixex=dble(srw*dsqrt(vh))
         return
      end if

      if(nog) then
         if(thet.le.0.8d0) then
            fxex=dble(fw(1)+fw(2)*theta)
            call copr4(2,xh,cw,vh)
            eixex=dble(srw*dsqrt(vh))
            return
         else
            iflag=-3
            return
         end if
      end if

c Giants (logg < 3.5)
      if(g.le.3.5d0) then
c Theta=(0.75,1.3)
         if(thet.lt.dble(xfg)) then
c If there is no value of z,g it cannot compute the index, except for theta 
c below 0.8
            if(noz) then
               if(thet.le.0.8d0) then
                  fxex=dble(fw(1)+fw(2)*theta)
                  call copr4(2,xh,cw,vh)
                  eixex=dble(srw*dsqrt(vh))
                  return
               else
                  iflag=-2
                  return
               end if
            end if
            if(thet.ge.0.8d0) then
               fxex=ctcg+dble(dexp(fcg(1)+fcg(2)*theta+fcg(4)*dble(z)+
     c           fcg(5)*thetaz+fcg(7)*theta2+fcg(16)*thet2z))
               call copr4(6,xh,ccg,vh)
               errt=dble(srcg*dsqrt(vh))
c This is the error in the mean index for these atm. par. If we wanted the 
c expected error for one measurement we should compute:
c               errt=sigres*sqrt(1+vh)
               eixex=errt*(fxex-ctcg)
c If theta is between 0.75 and 0.8 we interpolate between this value and the 
c constant value of warm stars at theta=0.75 (0.
            else
               dumt=0.8D+00
               fxex=ctcg+dble(dexp(fcg(1)+fcg(2)*dumt+fcg(4)*dble(z)+
     c           fcg(5)*dumt*dble(z)+fcg(7)*dumt*dumt+fcg(16)*dumt*dumt*
     c           dble(z)))
               call copr4(6,xh2,ccg,vh)
               errt=dble(srcg*dsqrt(vh))
               eixex=errt*(fxex-ctcg)
               x=(thet-0.75d0)/0.05d0
               fxex=dble(fw1)*(1.-x)+fxex*x
               eixex=errw*(1.-x)+eixex*x
            end if
         else
c Theta=(1.3,1.72) for giants
            fxex=dble(fg(1)+fg(2)*theta)
            dx=theta-xfg
            eixex=dble(dsqrt(srg*srg*(dx**2)*cg(1)+
     c        srfg*srfg*(1.D0-dx*cg(2))**2))
            if(g.gt.1.5d0) iflag=3
            if(thet.gt.1.73d0) iflag=2
         end if
         if(g.le.3.d0) then
            return
         else
            fxex0=fxex
            eixex0=eixex
         end if
      end if

c Dwarfs (logg > 3.)
      if(g.gt.3.d0) then
c Theta=(0.75,1.07)
c Cambiar lo siguiente:
         if(thet.le.1.0757d0) then
c If there is no value of z it cannot compute the index, except for theta 
c below 0.8
            if(noz) then
               if(thet.le.0.8d0) then
                  fxex=dble(fw(1)+fw(2)*theta)
                  call copr4(2,xh,cw,vh)
                  eixex=dble(srw*dsqrt(vh))
                  return
               else
                  iflag=-2
                  return
               end if
            end if
            if(thet.ge.0.8d0) then
               fxex=ctcd+dble(dexp(fcd(1)+fcd(2)*theta+fcd(4)*dble(z)+
     c           fcd(5)*thetaz+fcd(7)*theta2+fcd(16)*thet2z))
               call copr4(6,xh,ccd,vh)
               errt=dble(srcd*dsqrt(vh))
               eixex=errt*(fxex-ctcd)
c If theta is between 0.75 and 0.8 we interpolate between this value and the 
c constant value (for 0.6-0.75) 
            else
               dumt=0.8D+00
               fxex=ctcd+dble(dexp(fcd(1)+fcd(2)*dumt+fcd(4)*dble(z)+
     c           fcd(5)*dumt*dble(z)+fcd(7)*dumt*dumt+fcd(16)*dumt*dumt*
     c           dble(z)))
               call copr4(6,xh2,ccd,vh)
               errt=dble(srcd*dsqrt(vh))
               eixex=errt*(fxex-ctcd)
               x=(thet-0.75d0)/0.05d0
               fxex=dble(fw1)*(1.-x)+fxex*x
               eixex=errw*(1.d0-x)+eixex*x
            end if
         else
c Theta=(1.08,1.84) for dwarfs
            fxex=ctd+dble(dexp(fd(1)+fd(2)*theta+fd(7)*theta2))
            xh(3)=xh(5)
            call copr4(3,xh,cd,vh)
            errt=dble(srd*dsqrt(vh))
            eixex=errt*(fxex-ctd)
            if(g.lt.4.d0) iflag=3
            if(thet.gt.1.84d0) iflag=2
         end if
         if(g.ge.3.5d0) then
            return
         else
c Here we interpolate between this value and the one computed for giants 
c (this applies when log g is between 3.0 and 3.5):
            x=(g-3.d0)/0.5d0
            fxex=fxex0*(1.d0-x)+fxex*x
            eixex=eixex0*(1.d0-x)+eixex*x
            return
         end if
      end if
c Error: No value was computed
      iflag=-5
      return

      end
c
c
c
      SUBROUTINE copr4(n,x,c,v)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
      INTEGER n
c      double precision x(25),v,c(25,25)
      DIMENSION x(25),c(25,25)
      DIMENSION x1(25,25),x2(25,25),temp(25,25),fin(25,25)
c      INTEGER i,j,n1,n2,n3,n4
c      double precision x1(25,25),x2(25,25),temp(25,25),fin(25,25)
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

      call MMATX4(x1,c,temp,1,n,n,n)

      call MMATX4(temp,x2,fin,1,n,n,1)
      v=fin(1,1)

      return
      end

c SUBROUTINE to multiply matrices
c Input: a -> matrix n1 x n2
c        b -> matrix n3 x n4
c n2 must be equal to n3
c Input : n1,n2,n3,n4 (integers, DIMENSIONs of matrices)
c Output: c -> matrix n1 x n4
c In the main program matrices must be defined as double precision with 
c   DIMENSIONs 25x25
 
      SUBROUTINE MMATX4(a,b,c,n1,n2,n3,n4)
c      SUBROUTINE MMTX4(a,b,c,n1,n2,n3,n4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      implicit none
c      INTEGER n1,n2,n3,n4
      INTEGER i,j,k
      DIMENSION a(25,25),b(25,25),c(25,25)
c      double precision a(25,25),b(25,25),c(25,25)

      if(n2.ne.n3) then
         write(6,*) 'Error in matrix DIMENSIONs'
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


