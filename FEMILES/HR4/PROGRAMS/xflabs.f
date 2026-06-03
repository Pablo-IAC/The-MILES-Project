c Subrutina que calcula flujos absolutos filtros V, I_Johnson:
c flujv,flujcj
c Se pasan biss,viss,vij
c Llamadas desde hrsl.f,STU.f
c COMMON FNOR: a.f,hrsl.f,STU.f,xflabs.f
      SUBROUTINE xflabs(uvbus,biss,viss,vij,vk,flujb,flujv,
     &flujcj,fluju,flujk)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FNOR/fsol,vijvega,zerovv,zerovi,uvvega,zerovu,
     &vkvega,zerovk
      flujv=10.**(-0.4d0*(viss+zerovv-0.03d0))
      flujcj=flujv*10.**(0.4d0*(vij-vijvega+zerovv-zerovi))
      flujk=flujv*10.**(0.4d0*(vk-vkvega+zerovv-zerovk))
      fluju=flujv*10.**(0.4d0*(-uvbus+uvvega+zerovv-zerovu))
      flujv=flujv/fsol 
      flujcj=flujcj/fsol 
      fluju=fluju/fsol
      flujb=(flujv*(10.**((-0.4d0)*(biss-viss))))
      flujk=flujk/fsol
      RETURN
      END
