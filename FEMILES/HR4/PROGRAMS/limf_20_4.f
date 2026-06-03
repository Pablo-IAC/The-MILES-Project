c	---------Subrutina que calcula la IMF-------------------------
	SUBROUTINE limf(xxmasa,fimm)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	CHARACTER*1 shape
	DIMENSION sson(4)
c	double precision sson(4)
	COMMON/lmu0/bicl,bicp,bich
	COMMON/lmu1/ZMU
	COMMON/lshape/shape
	COMMON/lson/sson
	COMMON/lbarb/xlplow,xlpmed,ZMU1,ZMU2,ZMU3
      COMMON/turnp/xtpoint
c Turning point IMF:Ferreras,Weidner,Vazdekis,LaBarbera+15
c xtpoint=0.5d0
	alfx1=4.0d0 !Ferreras U peaked IMF (x option)
	alfx3=1.9d0 !Ferreras U peaked IMF (x option)
c      xnorm=4.d0/( 2.653238076d0+(1.d0-(.2d0)**(-ZMU))/(1.d0-ZMU))
      xnorm=1.0d0
c La Barbera
c      xlplow=0.4d0
c      xlpmed=0.7d0
c Kroupa universal
	xxmak0=0.08d0
	xxmak1=0.5d0
	alfu0=0.30d0
	alfu1=1.30d0
	alfu2=2.30d0
c Kroupa revisada
	xxmar0=0.08d0
	xxmar1=0.5d0
	xxmar2=1.0d0
	alfr0=0.30d0
	alfr1=1.80d0
	alfr2=2.70d0
	alfr3=2.30d0
c Chabrier
	c008=dlog10(0.08d0)
	cconst=2.0d0
	csigma=cconst*(0.69d0**cconst)
	cmasa=1.0d0 !masa cambio pendiente para normalizacion (fimmc)
      fimmc=(1.d0/cmasa)*exp(-((dlog10(cmasa)-c008)**2.d0)/csigma)
cccccccccccccccccc
	fimm=0.d0
      IF(shape.eq.'b')THEN
	 if(xxmasa.le.bicl)then
		fimm=bicp**(-ZMU+1.0d0)/xxmasa !nota:se pone +1.0 =>dividir
	 elseif((xxmasa.gt.bicl).and.(xxmasa.lt.bich))then
	 fimm=(sson(1)*xxmasa*xxmasa+sson(2)*xxmasa+sson(3))
     &  *1.0d0+sson(4)/xxmasa
	 else
		fimm=xxmasa**(-ZMU)
	 endif
	ELSEIF(shape.eq.'k')THEN
	 if(xxmasa.lt.xxmak0)then
	fimm=(xxmasa/xxmak0)**(-alfu0)
	 elseif(xxmasa.ge.xxmak0.and.xxmasa.lt.xxmak1)then
	fimm=(xxmasa/xxmak0)**(-alfu1)
	 else
	fimm=(xxmak1/xxmak0)**(-alfu1)*(xxmasa/xxmak1)**(-alfu2)
	 endif	 
	ELSEIF(shape.eq.'r')THEN
	 if(xxmasa.lt.xxmar0)then
	fimm=(xxmasa/xxmar0)**(-alfr0)
	 elseif(xxmasa.ge.xxmar0.and.xxmasa.lt.xxmar1)then
	fimm=(xxmasa/xxmar0)**(-alfr1)
	 elseif(xxmasa.ge.xxmar1.and.xxmasa.lt.xxmar2)then
	fimm=(xxmar1/xxmar0)**(-alfr1)*(xxmasa/xxmar1)**(-alfr2)
	 else
	fimm=(xxmar1/xxmar0)**(-alfr1)*(xxmar2/xxmar1)**(-alfr2)*
     & (xxmasa/xxmar2)**(-alfr3)
	 endif	 
	ELSEIF(shape.eq.'f')THEN
	 x=xxmasa/xtpoint
c	 if(xxmasa.lt.xtpoint)then
	 if(x.lt.1.0d0)then
c	  fimm=xxmasa**(-ZMU)
        fimm=xnorm*(x**(-ZMU))
	 else
c	  fimm=xxmasa**(-2.3d0)
        fimm=xnorm*(x**(-2.3d0))
	 endif
	ELSEIF(shape.eq.'x')THEN
	 x=xxmasa/xtpoint
	 if(x.lt.1.0d0)then
        fimm=x**(-alfx1)
	 else
        fimm=x**(-alfx3)
	 endif
	ELSEIF(shape.eq.'l')THEN
	 if(xxmasa.le.xlplow)then
		fimm=xxmasa**(-ZMU1)
	 elseif((xxmasa.gt.xlplow).and.(xxmasa.le.xlpmed))then
	      fimm=xxmasa**(-ZMU2)
	 else
		fimm=xxmasa**(-ZMU3)
	 endif
	ELSEIF(shape.eq.'c')THEN
ccccc	c008=dlog10(0.08d0)
ccccc	cconst=2.0d0
ccccc	csigma=cconst*(0.69d0**cconst)
ccccc	cmasa=1.0d0 !masa cambio pendiente para normalizacion (fimmc)
ccccc      fimmc=(1.d0/cmasa)*exp(-((dlog10(cmasa)-c008)**2.d0)/csigma)
	 if(xxmasa.le.cmasa)then
c      fimm=0.158d0*
c     &(1.d0/xxmasa)*exp(-((dlog10(xxmasa)-c008)**2.0d0)/csigma)
c        fimm=exp(-((dlog10(xxmasa)-c008)**2.0d0)/csigma)
ccccccccccc
c      fimm=exp(-(((dlog10(xxmasa)-c008)**2.d0)/csigma))
      fimm=(1.d0/xxmasa)*exp(-((dlog10(xxmasa)-c008)**2.d0)/csigma)
	fimm=fimm/fimmc
	 else
	  fimm=xxmasa**(-2.3d0)
	 endif
	ELSE
	 fimm=xxmasa**(-ZMU)
	ENDIF
	return
	end
