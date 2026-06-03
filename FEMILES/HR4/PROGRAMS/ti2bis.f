      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	 rikc=-0.2d0
	 vrkc=-0.6d0
	 vikc=0.d0
      DO i=1,50
	rikc=rikc+0.05
	vrkc=vrkc+0.05
	vikc=vrkc+rikc
ccccccccccccccc11111ccccccccccccccccccccccccccc
c transf. inversa cousin-johnson aplicada en MAKEISOC.f
cc transf cousin johnson R-I
c Fernie 1983, PASP,95,782 (V-I via R-I, V-R, transf. inversa)
      vij=0.0d0
      ricut=0.034d0+0.845d0*0.8d0
      if(rikc.le.ricut)then
       rij=(rikc-0.034d0)/0.845d0
      elseif(rikc.gt.1.42d0)then
       rij=0.205d0+0.733d0*rikc+0.171d0*rikc*rikc
      else
       ckc=(-rikc-0.239d0)
       bkc=1.315d0
       akc=-0.175d0
       dkc=sqrt(bkc*bkc-4.d0*akc*ckc)
       rij=(-bkc+dkc)/(2.d0*akc)
      endif
      vrcut=-0.024d0+.73d0*1.1d0
      if(vrkc.le.vrcut)then
       vrj=(vrkc+0.024d0)/0.73d0
      else
       vrj=(vrkc-0.218d0)/0.522d0
      endif
      vij=vrj+rij
cccccccccccccc22222cccccccccccccccccccccccccccc
c espectro o (record brillante)
c Fernie 1983, PASP,95,782 (V-I directo)
      vij2=0.0d0
      if(vikc.le.1.5d0)then
       vij2=-0.005d0+1.273d0*vikc
      else
       vij2=0.723d0+0.486d0*vikc+0.215d0*vikc*vikc
      endif
ccccccccccccccc33333ccccccccccccccccccccccccccc
c espectro medio
c Fernie 1983, PASP,95,782 (V-I via R-I, V-R)
      vij3=0.0d0
      if(vrkc.le.0.8d0)then
       vrj3=0.034d0+1.364d0*vrkc
      else
       vrj3=-0.311d0+1.803d0*vrkc
      endif
      if(rikc.le.0.7d0)then
       rij3=-0.04d0+1.176d0*rikc
      else
       rij3=0.205d0+0.733d0*rikc+0.171d0*rikc*rikc
      endif
      vij3=vrj3+rij3
ccccccccccccccc44444ccccccccccccccccccccccccccc
c Besell 1979, PASP,91,589
      vij4=0.0d0
      if(vikc.le.0.0d0)then
       vij4=vikc/0.713d0
      elseif(vikc.gt.0.0d0.and.vikc.le.(0.778d0*2.0d0))then 
       vij4=vikc/0.778d0
      else
       vij4=(vikc+0.13d0)/0.835d0
      endif
ccccccccccccccc44444ccccccccccccccccccccccccccc
c Besell 1983, PASP,95,480
	grav=1.
      vij5=0.0d0
      if(vrkc.le.0.8523d0)then
       vrj5=(vrkc+0.02d0)/0.715d0
      else
       if(grav.le.3.0d0)then
        vrj5=(vrkc-0.24d0)/0.5d0
       else
        vrj5=(vrkc-0.12d0)/0.6d0
       endif
      endif
      rikcd=(0.035d0+0.84d0*0.094d0/1.045d0)/(1.d0-0.84d0/1.045d0)
      rikcg=(0.84d0*0.1d0+0.035d0)/(1.d0-0.84d0)
      IF(grav.gt.3.0d0)then
       if(rikc.le.rikcd)then
        rij5=(rikc-0.035d0)/0.84d0
       else
        rij5=(rikc+0.094d0)/1.045d0
       endif
      ELSE
       if(rikc.le.rikcg)then
        rij5=(rikc-0.035d0)/0.84d0
       else
        rij5=rikc+0.1d0
       endif
      ENDIF
      vij5=vrj5+rij5
ccccccccccccccc44444ccccccccccccccccccccccccccc
c Besell 1983, PASP,95,480 combinado con
c Fernie 1983, PASP,95,782 (V-I via R-I, V-R, transf. inversa)
      vij0d=0.0d0
      vij0g1=0.0d0
      vij0g2=0.0d0
      rij0d=0.0d0
      rij0g1=0.0d0
      rij0g2=0.0d0
      vrj0d=0.0d0
      vrj0g=0.0d0
      rikcg1=-(0.1d0+0.034d0/0.845d0)/(1.d0-1.d0/0.845d0)
      rikcg2=-(0.15d0+1.045d0*0.034d0/0.845d0)/(1.d0-1.045d0/0.845d0)
      rikcdd=-(0.094d0+1.045d0*0.034d0/0.845d0)/(1.d0-1.045d0/0.845d0)
      DO k=1,2
       if(k.eq.1)giss=1.d0
       if(k.eq.2)giss=4.d0
      IF(giss.ge.3.5)THEN
       if(rikc.le.rikcdd)then
        rij0d=(rikc-0.034d0)/0.845d0
       else
        rij0d=(rikc+0.094d0)/1.045d0
       endif
      ELSE
       if(rikc.le.rikcg2)then
        rij0g2=(rikc-0.034d0)/0.845d0
       else
        rij0g2=(rikc+0.15d0)/1.045d0
       endif
      ENDIF
      vrkcdd=(0.12d0+0.6d0*0.024d0/0.73d0)/(1.d0-0.6d0/0.73d0)
      vrkcgg=(0.24d0+0.5d0*0.024d0/0.73d0)/(1.d0-0.5d0/0.73d0)
      IF(giss.ge.3.5)THEN
       if(vrkc.le.vrkcdd)then
        vrj0d=(vrkc+0.024d0)/0.73d0
       else
        vrj0d=(vrkc-0.12d0)/0.6d0
       endif
      ELSE
       if(rikc.le.vrkcgg)then
        vrj0g=(vrkc+0.024d0)/0.73d0
       else
        vrj0g=(vrkc-0.24d0)/0.5d0
       endif
      ENDIF
      ENDDO
      vij0d=vrj0d+rij0d
      vij0g2=vrj0g+rij0g2
c      vij0g1=vrj0g+rij0g1
      WRITE(37,*)vikc,vij,vij2,vij3,vij4,vij5,vij0d,vij0g2
      ENDDO
      end
