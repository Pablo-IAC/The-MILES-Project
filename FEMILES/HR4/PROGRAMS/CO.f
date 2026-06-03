*******SUBRUTINA CO**************************************************
c        SUBROUTINE CO(nra,i,W,stars)
        SUBROUTINE CO(i,W,stars)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION feh(15)
	logical agb
	COMMON/spec/R(12,2999,15),jotai(15)
	COMMON/fezsol/fez(15)
        COMMON/rat/ra(3,50)
	COMMON/KEYS/keypT(18)
        COMMON/AGBpad/gsi(5,2),Rm0(12,2999,1)
        COMMON/NC/ncha
	feh(i)=fez(i)
c	common/cd4000/fra4(14250),fbj(14250)
c        feh(i)=fez(i)+dlog10(ra(nra,15)/ra(1,15))
c	spend=(4000.d0-3600.d0)/(4400.d0-3600.d0)
	W=0.0d0
	stars=0.0d0
	if(jotai(i).eq.0)goto 5924
	DO j=1,jotai(i)
	 at=10**R(i,j,3)
c	 tetas=5040.d0/at
         at=5040.D0/at  !pues finfit_CO trabaja con theta 
	 ggrav=R(i,j,4)
	 zfeff=feh(i)
	 if(ncha.eq.0)then
	  if(j.gt.keypT(15)) then
	   agb=.true.
	  else
	   agb=.false.
	  endif
	 else
	  if(Rm0(i,j,1).gt.gsi(5,1)-0.000025d0) then !aprox redondeo en sum0*
	   agb=.true.
	  else
	   agb=.false.
	  endif
	 endif
c	   agb=.false.
         call finfit_CO(at,ggrav,zfeff,agb,ZMAG,eindx,iflage)
c	 W=W+ZMAG*R(i,j,14)*fra4(j)*fbj(j)
c	 stars=stars+R(i,j,14)*fra4(j)*fbj(j)
	 flujv=10**((-0.4d0)*((-2.5d0)*dlog10(R(i,j,7))-3.762d0))
	 vk=(-2.5d0)*dlog10(R(i,j,7)/R(i,j,12))
	 flujk=(flujv*(10**(+0.4d0*(vk))))
	 if(ggrav.ge.3.5d0)then
	  if(R(i,j,3).lt.3.81d0)then
	   frxCO=-0.8455d0+0.7653d0*R(i,j,3)
	  else
	   frxCO=0.0207
	  endif
	 else
	  if(R(i,j,3).lt.3.86d0)then
	   frxCO=0.9641d0+0.2865d0*R(i,j,3)
	  else
	   frxCO=0.0207
	  endif
	 endif
	 flujCO=frxCO*flujk
	 W=W+ZMAG*R(i,j,14)*flujCO !flujo banda K (mejor FK*fn(spec_C0)*fn
	 stars=stars+R(i,j,14)*flujCO
c	 write(25,*)fra4(j),fbj(j),(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
c     & ,(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)/(fra4(j)*fbj(j))
c	 W=W+ZMAG*R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
c	 stars=stars+R(i,j,14)*(R(i,j,5)+(R(i,j,6)-R(i,j,5))*spend)
c	 write(*,*)'agb=',agb,at,ggrav,zfeff
c	 write(36,*)'agb=',agb,Rm0(i,j,1),gsi(5,1),R(i,j,3),ggrav,ZMAG
	ENDDO
5924   continue
	return
	END

c Version Jun 2008
C-----------------------------------------------------------------------------
c INPUT:
c       t = effective temperature, in theta (5040/Teff)      (REAL)
c       g = logarithm (base 10) of the surface gravity       (REAL)
c       z = metallicity [Fe/H]                               (REAL)
c     agb = AGB star                                         (logical)
c  (A value higher or equal to 99 in g or z in input means that this 
c    parameter is unknown. In same cases, an index value can still be 
c    computed) 
c  
c
c OUTPUT:
c       findex = index value                                 (REAL)
c       iflag = indicates sucess of the functions:           (integer)
c              0 -> OK
c              1 -> extrapolation to high temperatures
c              2 -> extrapolation to low temperatures
c              3 -> extrapolation in g (uncertain value)
c             -1 -> Error (no Teff)
c             -2 -> Error (z needed)
c             -3 -> Error (g needed)
c             -5 -> Error
c******************************************************************************
c INDEX DEFINITION: 
c The CO break is defined as the ratio flux_red/flux_blue, where flux_red is the
c averaged flux per Angstrom in the region of the CO absorption at 2.29 microns
c (namely in the range between 22880.0 and 23010.0 Angstroms), and flux_blue is
c the averaged flux per Angstrom in a pseudo-continuum region defined with two
c separated bands (first band: 22460.0-22550.0, second band: 22710.0-22770.0).
c For more details, see Marmol-Queralto et al. (2008).
c******************************************************************************
      subroutine finfit_CO(t,g,z,agb,findex,eindex,iflag)
      implicit none
      
      integer ncomax
      parameter(ncomax=20)
 
      integer iflag
      double precision t,g,z,findex,eindex,pi
c      double precision tt
      
      double precision findex1,findex2,findex3,findex4,findex5,findex6,findex7
      double precision eindex1,eindex2,eindex3,eindex4,eindex5,eindex6,eindex7
      double precision findex8,eindex8
      double precision thet,x
      double precision lim_t_g1,lim_t_g2,lim_t_g3,lim_g
      double precision lim_t_dc_d,lim_t_iz_d
      common/blklim/lim_t_dc_d,lim_t_iz_d,lim_g
      common/blklim_g/lim_t_g1,lim_t_g2,lim_t_g3
      
      double precision theta,geta,zeta
      double precision xf(ncomax)
      logical nog,noz
      
      common/blkx/xf

      logical agb

c-----------------------------------------------------------------------
      pi=3.14159d0

      nog=.false.
      noz=.false.
      if(g.ge.99.) nog=.true.
      if(z.ge.99.) noz=.true.

      iflag=0
      findex=0.
      eindex=0.
      if(t.le.0.) then
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

C Gravity limit 
      lim_g=3.50
C theta limits for dwarfs
      lim_t_iz_d=0.80
      lim_t_dc_d=0.90
C theta limits for giants
      lim_t_g1=0.90
      lim_t_g2=1.10
      lim_t_g3=1.55
      lim_t_g3=1.60

      if(agb)then
        if(thet.lt.lim_t_g1) then
          call compindexCO(3,findex3,eindex3)
          findex=findex3
          eindex=eindex3
          iflag=0
          if(thet.lt.5040./12136.) iflag=1
          return
        endif  
        if(thet.ge.0.90.and.thet.lt.1.55) then ! caja agb
          call compindexCO(8,findex8,eindex8)
          findex=findex8
          eindex=eindex8
          iflag=0
          return
        end if
        if(thet.ge.1.55) then ! caja 6
          findex=1.311
          eindex=0.010
          iflag=2
          return
        end if
      end if

      if(.not.agb)then
c Dwarfs (g > 3.5)
      if(g.ge.lim_g)then
c 1.- Hot dwarfs      
        if(thet.lt.lim_t_dc_d) then
          call compindexCO(1,findex1,eindex1)
          if(thet.lt.lim_t_iz_d) then
            findex=findex1
            eindex=eindex1
            iflag=0
            if(thet.lt.5040./13397.) iflag=1
            return
          end if
        end if  
c 2.- Cold dwarfs      
        if(thet.ge.lim_t_iz_d.and.thet.le.1.50) then
          if(noz) then 
            iflag=-2
            return
          end if
          call compindexCO(2,findex2,eindex2)
          if(thet.gt.lim_t_dc_d.and.thet.le.1.45) then
            findex=findex2
            eindex=eindex2
            iflag=0
            return
          end if
        end if
c Intersection
        if((thet.ge.lim_t_iz_d).and.(thet.le.lim_t_dc_d)) then
          x=(pi/2)*(thet-lim_t_iz_d)/(lim_t_dc_d-lim_t_iz_d)
          findex=findex1*cos(x)+findex2*(1.0-cos(x))
          eindex=eindex1*cos(x)+eindex2*(1.0-cos(x))
          iflag=0
          return
        end if
c 7.- Cool dwarfs      
        if((thet.gt.1.45)) then
          call compindexCO(7,findex7,eindex7)
          if(thet.gt.1.50) then
            findex=findex7
            eindex=eindex7
            iflag=0
            if(thet.gt.5040./2799.) iflag=2
            return
          endif
        endif
c Intersection
        if((thet.ge.1.45).and.(thet.le.1.50)) then
          x=(pi/2)*(thet-1.45)/(1.50-1.45)
          findex=findex2*cos(x)+findex7*(1.0-cos(x))
          eindex=eindex2*cos(x)+eindex7*(1.0-cos(x))
          iflag=0
          return
        end if

C Giants (g < 3.5)
      elseif(g.lt.lim_g)then
c 3.- Hot giants
        if(thet.lt.lim_t_g1+0.03) then 
          if(thet.lt.0.293) iflag=1
          call compindexCO(3,findex3,eindex3)
          if(thet.lt.lim_t_g1) then 
            findex=findex3
            eindex=eindex3
            iflag=0
            if(thet.lt.5040./12136.) iflag=1
            return
          endif
        endif  
c 4.- Warm giants
        if(thet.ge.lim_t_g1.and.thet.lt.lim_t_g2+0.03) then
          if(noz) then
            iflag=-2
            return
          end if
          call compindexCO(4,findex4,eindex4)
          if(thet.ge.lim_t_g1+0.03.and.thet.lt.lim_t_g2-0.00) then
            findex=findex4
            eindex=eindex4
            iflag=0
          return
          endif
        endif
c Intersection
        if((thet.ge.lim_t_g1).and.(thet.lt.lim_t_g1+0.03)) then
          x=(pi/2)*(thet-lim_t_g1)/(lim_t_g1+0.03-lim_t_g1)
          findex=findex3*cos(x)+findex4*(1.0-cos(x))
          eindex=eindex3*cos(x)+eindex4*(1.0-cos(x))
          iflag=0
          return
        endif
c 5.- Cool giants
        if(thet.ge.lim_t_g2-0.01.and.thet.lt.lim_t_g3) then 
          if(noz) then
            iflag=-2
            return
          end if
          call compindexCO(5,findex5,eindex5)
          if(thet.ge.lim_t_g2+0.01.and.thet.lt.lim_t_g3) then 
            findex=findex5
            eindex=eindex5
            iflag=0
            return
          end if
        end if  
c Intersection
        if((thet.ge.lim_t_g2-0.01).and.(thet.lt.lim_t_g2+0.01)) then
          x=(thet-lim_t_g2-0.01)/(lim_t_g2+0.01-lim_t_g2+0.01)
          findex=(findex4+findex5)/2
          eindex=(eindex4+eindex5)/2
          iflag=0
          return
        endif
c 6.- Cold giants        
        if(thet.ge.lim_t_g3) then ! caja 6
          call compindexCO(6,findex6,eindex6)
          if(thet.ge.lim_t_g3) then
            findex=findex6
            eindex=eindex6
            iflag=0
            if(thet.gt.5040./2485.) iflag=2
            return
          end if
        end if

      end if   
      end if ! end agb

c If it is not in any box:
c Error: No value was computed
      iflag=-5
      write(*,'(A)') 'ERROR: NO index value was computed'
      stop

      end
C
C******************************************************************************
C
      subroutine compindexCO(it,fin,ein)
      implicit none

      integer ncomax
      parameter(ncomax=20)
 
      integer i,j,it,nter
      double precision fin,ein
      
      double precision xf(ncomax),x(ncomax)  
      double precision indt,derr

      common/blkx/xf
      
      double precision theta_min,theta_max
      double precision g_min,g_max
      double precision feh_min,feh_max
      double precision coefbox(ncomax),cte,sr
      double precision covar(ncomax,ncomax)
      logical flpolybox,flweigbox
      logical flcoebox(ncomax)
      
      common/blkrange/theta_min,theta_max,g_min,g_max,feh_min,feh_max  
      common/blkbox2/coefbox,cte,sr      
      common/blkcovar/covar

      common/blkfitbox/flweigbox,flpolybox,flcoebox
      
      double precision lim_t_g1,lim_t_g2,lim_t_g3,lim_g
      double precision lim_t_dc_d,lim_t_iz_d
      common/blklim/lim_t_dc_d,lim_t_iz_d,lim_g
      common/blklim_g/lim_t_g1,lim_t_g2,lim_t_g3

c Boxes
c 1.- Hot dwarfs
      if(it.eq.1) then
        call read_coeffCO(1)
      end if
c 2.- Cool dwarfs
      if(it.eq.2) then
        call read_coeffCO(2)
      end if
c 7.- Cold dwarfs
      if(it.eq.7) then
        call read_coeffCO(7)
      end if

c 3.- Hot giants
      if(it.eq.3) then
        call read_coeffCO(3)
      end if
c 4.- Warm giants
      if(it.eq.4) then
        call read_coeffCO(4)
      end if
c 5.- Cool giants
      if(it.eq.5) then
        call read_coeffCO(5)
      end if
c 5.- Cold giants
      if(it.eq.6) then
        call read_coeffCO(6)
      end if
c 8.- AGB stars
      if(it.eq.8) then
        call read_coeffCO(8)
      end if

c Now it computes the index
      indt=0.D0
      j=0
      do i=1,ncomax
         if(flcoebox(i)) then 
            indt=indt+coefbox(i)*xf(i)
            j=j+1
            x(j)=xf(i)
         end if
      end do
      nter=j

      if(flpolybox) then
         fin=real(indt) ! polinomial fit
      else
         fin=real(cte+dexp(indt))
      end if            ! exponential fit

c and the error
      call comperrCO(nter,x,covar,derr)
      ein=real(sr*dsqrt(derr))
      if(.not.flpolybox) ein=ein*(fin-real(cte))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comperrCO(n,x,covar,v)
      implicit none

      integer ncomax
      parameter(ncomax=20)

      integer n
      double precision x(ncomax),v,covar(ncomax,ncomax)

      integer i,j
      double precision x1(ncomax,ncomax),x2(ncomax,ncomax)
      double precision temp(ncomax,ncomax),fin(ncomax,ncomax)

c First it reconstructs the whole v-c matrix
      do i=2,n
         do j=1,i-1
            covar(i,j)=covar(j,i)
         end do
      end do

      do j=1,n
         x1(1,j)=x(j)
         x2(j,1)=x(j)
      end do

      call multmatrixCO(x1,covar,temp,1,n,n,n)

      call multmatrixCO(temp,x2,fin,1,n,n,1)
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
 
      subroutine multmatrixCO(a,b,c,n1,n2,n3,n4)
      implicit none

      integer ncomax
      parameter(ncomax=20)

      integer n1,n2,n3,n4
      integer i,j,k
      double precision a(ncomax,ncomax),b(ncomax,ncomax)
      double precision c(ncomax,ncomax)

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
      

c******************************************************************************
      subroutine read_coeffCO(it)
      implicit none

      integer ncomax
      parameter(ncomax=20)
      
      integer it
      integer i,ii,j,jj
      double precision coefbox(ncomax),cte,sr
      double precision c(ncomax,ncomax),covar(ncomax,ncomax)
      logical flpolybox,flweigbox
      logical flcoebox(ncomax)
                                                     
      common/blkbox2/coefbox,cte,sr
      common/blkcovar/covar
      common/blkfitbox/flweigbox,flpolybox,flcoebox
      
      flweigbox=.false.
      flpolybox=.false.
      do i=1,ncomax
        flcoebox(i)=.false.
        coefbox(i)=0.D0
      end do
      do i=1,ncomax
        do ii=1,ncomax
          c(i,ii)=0.D0
        end do
      end do

C 1: Hot dwarfs   
      if(it.eq.1)then
c        print*,"Box 1: hot dwarfs"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.00557 
        
        flcoebox(1) = .true.
        coefbox(1) =  0.4990E-01 
        
        c( 1, 1) =  0.3613275463615035E-01
      end if

C 2: Cool dwarfs
      if(it.eq.2)then
c        print*,"Box 2: cool dwarfs"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.01254

        flcoebox(1) = .true.
        flcoebox(2) = .true.
        flcoebox(4) = .true.

        coefbox(1) = -0.2919E-01
        coefbox(2) =  0.1006E+00
        coefbox(4) =  0.1740E-01

        c( 1, 1) =  0.6933933560618123E+01
        c( 1, 2) = -0.6772316706182383E+01
        c( 1, 4) =  0.9876773279330903E-01
        c( 2, 2) =  0.6865432466423429E+01
        c( 2, 4) =  0.9029410053379856E-01
        c( 4, 4) =  0.2506507940307781E+00
      end if

C 7: Cold dwarfs
      if(it.eq.7)then
c        print*,"Box 7: cold dwarfs"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.01195 
        
        flcoebox(1) = .true.
        coefbox(1) =   0.1025E+00
        
        c( 1, 1) =  0.1467084569095854E+00
      end if

C 3: Hot giants     
      if(it.eq.3)then
c        print*,"Box 3: hot giants" 
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.00398 
        
        flcoebox(1) = .true.
        coefbox(1) =  0.4593E-01 
        
        c( 1, 1) =  0.6683713759125805E-01
      end if
      
C 4: Warm giants
      if(it.eq.4)then
c        print*,"Box 4: warm giants"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.00428 

        flcoebox( 1) = .true.
        flcoebox( 2) = .true.
        flcoebox( 4) = .true.
        flcoebox( 5) = .true.
        flcoebox( 9) = .true.

        coefbox( 1) = -0.3073E+00
        coefbox( 2) =  0.3876E+00
        coefbox( 4) = -0.1016E+00
        coefbox( 5) =  0.1072E+00
        coefbox( 9) = -0.2334E-02

        c( 1, 1) =  0.1155714142626352E+01
        c( 1, 2) = -0.1086914746336309E+01
        c( 1, 4) =  0.6681691917102971E+00
        c( 1, 5) = -0.6295587274116440E+00
        c( 1, 9) = -0.9242526688430821E-03
        c( 2, 2) =  0.1029177287355958E+01
        c( 2, 4) = -0.6303994760801731E+00
        c( 2, 5) =  0.5960077000253381E+00
        c( 2, 9) = -0.5021569136132354E-03
        c( 4, 4) =  0.9546218748120007E+00
        c( 4, 5) = -0.8726802323703131E+00
        c( 4, 9) =  0.1786108624822258E-01
        c( 5, 5) =  0.8263159725831999E+00
        c( 5, 9) = -0.4091026319139710E-03
        c( 9, 9) =  0.1148851448353762E-01
       end if  

C 5: Cool giants
      if(it.eq.5)then
c        print*,"Box 5: cool giants"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.00890

        flcoebox( 1) = .true.
        flcoebox( 2) = .true.
        flcoebox( 4) = .true.
        flcoebox( 5) = .true.
        flcoebox( 7) = .true.
        flcoebox( 9) = .true.
       
        coefbox( 1) = -0.5224E+00
        coefbox( 2) =  0.8257E+00
        coefbox( 4) =  0.6741E-01
        coefbox( 5) = -0.4441E-01
        coefbox( 7) = -0.2200E+00
        coefbox( 9) = -0.2329E-02

        c( 1, 1) =  0.1187933343425294E+03
        c( 1, 2) = -0.1733806743923548E+03
        c( 1, 4) =  0.5225759189817347E+00
        c( 1, 5) = -0.5468898975749468E+00
        c( 1, 7) =  0.6207323233644743E+02
        c( 1, 9) = -0.3047838926003409E+00
        c( 2, 2) =  0.2535468526860264E+03
        c( 2, 4) = -0.3388600000069259E+00
        c( 2, 5) =  0.5292778221401269E+00
        c( 2, 7) = -0.9092364682762397E+02
        c( 2, 9) =  0.4565099774687261E+00
        c( 4, 4) =  0.1296249767222967E+01
        c( 4, 5) = -0.8107019893274662E+00
        c( 4, 7) = -0.7915064713622872E-02
        c( 4, 9) =  0.3060567881074111E-01
        c( 5, 5) =  0.5341954813762932E+00
        c( 5, 7) = -0.1059988667569693E+00
        c( 5, 9) =  0.3198723098165772E-02
        c( 7, 7) =  0.3265359033826216E+02
        c( 7, 9) = -0.1679684057759219E+00
        c( 9, 9) =  0.2415698842904209E-01
      end if

C 6: Cold giants
      if(it.eq.6)then
c        print*,"Box 6: cold giants"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.02698 
        
        flcoebox(1) = .true.
        coefbox(1) =   0.2397E+00
        
        c( 1, 1) =  0.1582881020191905E+00
      end if

C 8: AGB stars  
      if(it.eq.8)then
c        print*,"Box 8: AGB cool giants"
        flweigbox=.true.
        flpolybox=.false.
        cte = 0.0
        sr =  0.00892 
        
        flcoebox(1) = .true.
        flcoebox(2) = .true.
        flcoebox(7) = .true.
        
        coefbox(1) =  -0.8893E+00
        coefbox(2) =   0.1495E+01
        coefbox(7) =  -0.4816E+00
        
        c( 1, 1) =  0.6077575384643615E+03
        c( 1, 2) = -0.9973549039048630E+03
        c( 1, 7) =  0.4013994054338606E+03
        c( 2, 2) =  0.1642444026682195E+04
        c( 2, 7) = -0.6631230718155189E+03
        c( 7, 7) =  0.2685519539487602E+03
      end if

c And it moves the coefficients different from zero to the upper side of the
C bidimensional matrix (necessary for the correct computation of the error
c index)
      j=0
      jj=0
      do i=1,ncomax
        if(flcoebox(i)) j=j+1
        jj=j-1
        do ii=i,ncomax
          if(c(i,ii).ne.0.0) then
            jj=jj+1
            covar(j,jj)=c(i,ii)
          end if
        end do
      end do
c
      return
      end
