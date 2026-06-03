c*****SUBRUTINA QUE ESCRIBE LOS RESULTADOS******************************
      SUBROUTINE e_IRTF(alpha,izo,realeo,izy,realey)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION a90k(3)
      CHARACTER*1 shape
      CHARACTER*2 shapa
      COMMON/lmu1/ZMU !aaa,a_*,limf,e_*
      COMMON/lshape/shape
      COMMON/fezsol/fez(15)
      COMMON/qua/quaqn,quaq10,quaq05,quaq15
      COMMON/escK/a90k
      if(shape.eq.'u')then
       shapa='UN'
      elseif(shape.eq.'b')then
       shapa='BI'
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
      write(52,'(A2,1x,(F5.2,1x),5(F7.4,1x),1(F6.0,2(1x,F6.3)))')
     &shapa,ZMU-1.d0,zetao,realeo,zetay,reay,alpha,(a90k(nt),nt=1,3)
c    escribe en out_qua
      write(39999,'(A2,1x,(F5.2,1x),2(F7.4,1x),7(F12.8))')shapa,
     &ZMU-1.0d0,zetao,realeo,quaqn/quaq10,quaqn/quaq05
     &,quaqn/quaq15,quaqn,quaq10,quaq05,quaq15
      RETURN
      END
