      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (nstk=400) !IRTF
      DIMENSION aak(nstk,4)
      DIMENSION stark(nstk)
      DIMENSION nscb(8),nsc(8),nbande(8)
      DIMENSION sigbox(8,3)
      CHARACTER*20 stark
      CHARACTER*6 anom
      nstark=0
      OPEN(99,FILE=
c     &'/home/vazdekis/vazdekis/scratch/D/INPUT/PARAM_IRTF_ALL'
c     &,STATUS='OLD')
     &'/home/vazdekis/vazdekis/scratch/D/INPUT/PARAM_IRTF'
     &,STATUS='OLD')
      do k=1,99999
         read(99,*,end=24)anom,(aak(k,l),l=1,4)
         nstark=nstark+1
      enddo
24    CLOSE(99)
      iiak=0
      nbxstt=0 !numero total estrellas en las 8 subcajas
      re10=0.6827d0*0.5d0 !frac. subbox max. norm.Q (0.34=>1sigma)
      re05=0.3829d0*0.5d0 !frac. subbox max. norm.Q (0.19=>0.5sigma)
      re15=0.8664d0*0.5d0 !frac. subbox max. norm.Q (0.43=>1.5sigma)
      ultra0=1.0d-12
      t50=5040.0d0
      tz0_i=3250.0d0 !temp mas baja para la que distinguimos[M/H]
      tz0_f=7000.0d0 !temp mas alta para la que distinguimos [M/H]
      ddt0=0.008d0
      tlinmi=50.0d0
      ddt0M=0.15d0
      tlinma=525.0d0
      DO iy=1,nstark !bucle estrellas IRTF (ahi debe maximizar dens)
        tiss=aak(iy,1)
        giss=aak(iy,2)
        fiss=aak(iy,3)
	if((ddt0*tiss*tiss/t50).lt.tlinmi)then
	  ddt=abs(t50*tlinmi/(tiss*tiss))
	elseif((ddt0*tiss*tiss/t50).gt.tlinma)then
	  ddt=abs(t50*tlinma/(tiss*tiss))
	else
	  ddt=ddt0
	endif
	if(abs(ddt0M*tiss*tiss/t50).gt.tlinma)then
	  ddtmax=abs(t50*tlinma/(tiss*tiss))
	else
	  ddtmax=ddt0M
	endif
	tisst=t50/tiss
	ddz=0.09d0
	ddg=0.11d0
	ddgmax=0.59d0
	ddzmax=0.30d0
	do isbox=1,8
	 sigbox(isbox,1)=ddtmax
	 sigbox(isbox,2)=ddgmax
	 sigbox(isbox,3)=ddzmax
	enddo
c	DMA=300.0d0 !PARAM_IRTF
c	DMA=700.0d0 !PARAM_IRTF_ALL
	sigtim=3.0d0
	rtimes=0.0d0
	xs=1.5d0
        nttmax=10
	fxs=0.5d0
	snamax=100.0d0
	snamin=25.0d0
	xnumno=0.0d0 !numero pixels en banda normalizacion
      sigmat=0.0d0
      sigmag=0.0d0
      sigmaz=0.0d0
	DO nnnns=1,1000
	  rtimes=rtimes+1.0d0
	  starsi=0.0d0
	  dens=0.0d0
	  dddt=rtimes*sigtim*ddt
	  dddg=rtimes*sigtim*ddg
	  dddz=rtimes*sigtim*ddz
	  do l=1,nstark
           IF(tiss.le.tz0_i.or.tiss.ge.tz0_f)THEN
	if(abs(t50/aak(l,1)-tisst).le.dddt.and.
     &  abs(aak(l,2)-giss).le.dddg)then
     	      starsi=starsi+1.0d0
	endif
	   ELSE
	    if(abs(t50/aak(l,1)-tisst).le.dddt
     &       .and.abs(aak(l,2)-giss).le.dddg
     &       .and.abs(aak(l,3)-fiss).le.dddz)then
     	     starsi=starsi+1.0d0
	    endif
	   ENDIF
	  enddo
	  if(starsi.gt.0.0d0)then
	    dens=starsi/(2.0d0*dddt*2.0d0*dddg*2.0d0*dddz)
	    if(dens.gt.DMA)then
       	     sigmat=ddt
             sigmag=ddg
             sigmaz=ddz
            else
	     dense=((dens-DMA)/DMA)**2.      
             sigmat=ddt*exp(dense*dlog(ddtmax/ddt))
             sigmag=ddg*exp(dense*dlog(ddgmax/ddg))
             sigmaz=ddz*exp(dense*dlog(ddzmax/ddz))
           endif      
	    goto 222
	  endif
	ENDDO
222	continue
	do lp=1,8
          nsc(lp)=0 !numero de estrellas por cubo
          nscb(lp)=0 !numero de estrellas por cubo provisional
	  nbande(lp)=0 !bandera: 0=no estrellas 1=si estrellas
	enddo
	write(44,'(F5.0)')dens
      ENDDO
      END
