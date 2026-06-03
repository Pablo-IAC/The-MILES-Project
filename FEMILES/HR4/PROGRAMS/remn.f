c**************************************************************************
c	---------Subrutina que calcula los remanentes estelares------------
	SUBROUTINE remn(m,xm,r)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON/REMNAN/re(15,2,99),nr(15)
	dimension xmari1(nr(m)),xmari2(nr(m))
	do i=1,nr(m)
	 xmari1(i)=re(m,1,i)
	 xmari2(i)=abs(re(m,1,i)-re(m,2,i))
	enddo
	if(xm.ge.xmari1(1).and.xm.le.5.0d0)then
c Receta de Marigo (2001)
 	  call hunt(xmari1,nr(m),xm,k1)
	  call polint(xmari1(k1-1),xmari2(k1-1),3,xm,r,dy)
c	 write(61,'(4(f7.4,1x))')xmari1(k1-1),xmari1(k1),xm,r
c	elseif(xm.gt.5.0d0.and.xm.lt.8.5d0)then
c Receta de Renzini & Ciotti (1993)
c	     r=0.077d0*xm+0.48d0
c	elseif(xm.ge.8.5d0.and.xm.lt.40.0d0)then
	elseif(xm.gt.5.0d0.and.xm.lt.40.0d0)then
	     r=1.4d0
	elseif(xm.ge.40.0d0)then
	     r=0.5d0*xm
	else
	     r=0.0d0
	endif
	return
	end
