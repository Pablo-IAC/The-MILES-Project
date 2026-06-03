c*****SUBRUTINA QUE ESCRIBE LOS RESULTADOS******************************
      SUBROUTINE e_XSL(alpha,izo,realeo,izy,realey)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a90i(3)
      CHARACTER*1 shape
      CHARACTER*2 shapa
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
      COMMON/lshape/shape
      COMMON/fezsol/fez(15)
      COMMON/qua/quaqn,quaq10,quaq05,quaq15
      COMMON/escI/a90i !a,e
cc    COMMON/Uflujs/fxUtm,fxUti,fxUt,fxUtms,fxUtis,fxUts,fxUtis004 !a,e
      if(shape.eq.'u')then
       shapa='UN'
      elseif(shape.eq.'b')then
       shapa='BI'
      elseif(shape.eq.'t')then
       shapa='Ba'
      elseif(shape.eq.'k')then
       shapa='Ku'
      elseif(shape.eq.'r')then
       shapa='Kb'
      elseif(shape.eq.'c')then
       shapa='Ch'
      endif
      if(alpha.eq.1.0)then
       zetay=9.99d0
       reay=9.99d0
      else
       zetay=fez(izy)	
       reay=realey
      endif	
      zetao=fez(izo)
c    escribe en out_ef
      write(50,'(A2,1x,(F5.2,1x),5(F7.4,1x),1(F6.0,2(1x,F6.3)))')
     &shapa,ZMU-1.d0,zetao,realeo,zetay,reay,alpha,(a90i(nt),nt=1,3)
c	escribe en out_u
cc        write(49,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F14.8,1x))')
cc     &  shapa,ZMU-1.d0,zetao,realeo,
cc     &  fxUtm/fxUt,fxUti/fxUt,fxUt,
cc     &  fxUtms/fxUts,fxUtis/fxUts,fxUtis004/fxUts,fxUts
c    escribe en out_qua
      write(399,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F12.8))')shapa,ZMU-1.0d0
     &,zetao,realeo,quaqn/quaq10,quaqn/quaq05
     &,quaqn/quaq15,quaqn,quaq10,quaq05,quaq15
      RETURN
      END
