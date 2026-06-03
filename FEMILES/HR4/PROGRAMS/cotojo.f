c Subrutina que transforma V-I Cousin a V-I Johnson
c Se pasan giss,VR,VI Cousin y devuelve VI Johnson
c Llamadas desde hrsl.f,STU.f
      SUBROUTINE cotojo(giss,vrkc,vikc,vij)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
ccccccccccccccc44444ccccccccccccccccccccccccccc
c Besell 1983, PASP,95,480 combinado con transf. MAKEISO.f
c Fernie 1983, PASP,95,782 (V-I via R-I, V-R, transf. inversa)
      rikc=vikc-vrkc
      vij=0.0d0
      rij=0.0d0
      vrj=0.0d0
c      rikcg1=-(.1d0+.034d0/.845d0)/(1.d0-1.d0/.845d0)
      rikcg1=-(.095d0+.034d0/.845d0)/(1.d0-1.d0/.845d0)
      rikcg2=-(.15d0+1.045d0*.034d0/.845d0)/(1.d0-1.045d0/.845d0)
      rikcdd=-(.094d0+1.045d0*.034d0/.845d0)/(1.d0-1.045d0/.845d0)
      IF(giss.ge.3.5d0)THEN
       if(rikc.le.rikcdd)then
        rij=(rikc-0.034d0)/0.845d0
       else
        rij=(rikc+0.094d0)/1.045d0
       endif
c       if(rikc.gt.0.5d0.and.rikc.lt.1.1d0)then
c        rij=rikc+0.04d0
c       endif
      ELSE
       if(rikc.le.rikcg2)then
c       if(rikc.le.rikcg1)then
        rij=(rikc-0.034d0)/0.845d0
       else
        rij=(rikc+0.15d0)/1.045d0
c        rij=rikc+0.095d0
       endif
c       if(rikc.gt.rikcg1.and.rikc.lt.1.0d0)then
c        rij=rikc+0.10d0
c       endif
      ENDIF
      vrkcdd=(0.12d0+0.6d0*0.024d0/0.73d0)/(1.d0-0.6d0/0.73d0)
      vrkcgg=(0.24d0+0.5d0*0.024d0/0.73d0)/(1.d0-0.5d0/0.73d0)
      IF(giss.ge.3.5d0)THEN
       if(vrkc.le.vrkcdd)then
        vrj=(vrkc+0.024d0)/0.73d0
       else
        vrj=(vrkc-0.12d0)/0.6d0
       endif
      ELSE
       if(rikc.le.vrkcgg)then
        vrj=(vrkc+0.024d0)/0.73d0
       else
        vrj=(vrkc-0.24d0)/0.5d0
       endif
      ENDIF
      vij=vrj+rij
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc11111ccccccccccccccccccccccccccc
c transf. inversa cousin-johnson aplicada en MAKEISOC.f
cc transf cousin johnson R-I
c Fernie 1983, PASP,95,782 (V-I via R-I, V-R, transf. inversa)
c      vij=0.0d0
c      ricut=0.034d0+0.845d0*0.8d0
c      if(rikc.le.ricut)then
c	rij=(rikc-0.034d0)/0.845d0
c      elseif(rikc.gt.1.42d0)then
c	rij=0.205d0+0.733d0*rikc+0.171d0*rikc*rikc
c      else
c	ckc=(-rikc-0.239d0)
c	bkc=1.315d0
c	akc=-0.175d0
c	dkc=sqrt(bkc*bkc-4.d0*akc*ckc)
c	rij=(-bkc+dkc)/(2.d0*akc)
c      endif
c      vrcut=-0.024d0+.73d0*1.1d0
c      if(vrkc.le.vrcut)then
c	vrj=(vrkc+0.024d0)/0.73d0
c      else
c	vrj=(vrkc-0.218d0)/0.522d0
c      endif
c      vij=vrj+rij
ccccccccccccccc22222cccccccccccccccccccccccccccc
cc espectro o (record brillante)
cc Fernie 1983, PASP,95,782 (V-I directo)
c      vij2=0.0d0
c      if(vikc.le.1.5d0)then
c	vij2=-0.005d0+1.273d0*vikc
c      else
c	vij2=0.723d0+0.486d0*vikc+0.215d0*vikc*vikc
c      endif
cccccccccccccccc33333ccccccccccccccccccccccccccc
cc espectro medio
cc Fernie 1983, PASP,95,782 (V-I via R-I, V-R)
c      vij3=0.0d0
c      if(vrkc.le.0.8d0)then
c	vrj3=0.034d0+1.364d0*vrkc
c      else
c	vrj3=-0.311d0+1.803d0*vrkc
c      endif
c      if(rikc.le.0.7d0)then
c	rij3=-0.04d0+1.176d0*rikc
c      else
c	rij3=0.205d0+0.733d0*rikc+0.171d0*rikc*rikc
c      endif
c      vij3=vrj3+rij3
cccccccccccccccc44444ccccccccccccccccccccccccccc
cc Besell 1979, PASP,91,589
c      vij4=0.0d0
c      if(vikc.le.0.0d0)then
c	vij4=vikc/0.713d0
c      elseif(vikc.gt.0.0d0.and.vikc.le.(0.778d0*2.0d0))then 
c	vij4=vikc/0.778d0
c      else
c	vij4=(vikc+0.13d0)/0.835d0
c      endif
cccccccccccccccc44444ccccccccccccccccccccccccccc
cc Besell 1983, PASP,95,480
c	 grav=1.
c      vij5=0.0d0
c      if(vrkc.le.0.8523d0)then
c	vrj5=(vrkc+0.02d0)/0.715d0
c      else
c	if(grav.le.3.0d0)then
c	 vrj5=(vrkc-0.24d0)/0.5d0
c	else
c	 vrj5=(vrkc-0.12d0)/0.6d0
c	endif
c      endif
c      rikcd=(0.035d0+0.84d0*0.094d0/1.045d0)/(1.d0-0.84d0/1.045d0)
c      rikcg=(0.84d0*0.1d0+0.035d0)/(1.d0-0.84d0)
c      IF(grav.gt.3.0d0)then
c	if(rikc.le.rikcd)then
c	 rij5=(rikc-0.035d0)/0.84d0
c	else
c	 rij5=(rikc+0.094d0)/1.045d0
c	endif
c      ELSE
c	if(rikc.le.rikcg)then
c	 rij5=(rikc-0.035d0)/0.84d0
c	else
c	 rij5=rikc+0.1d0
c	endif
c      ENDIF
c      vij5=vrj5+rij5
