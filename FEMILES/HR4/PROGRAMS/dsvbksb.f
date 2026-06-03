      SUBROUTINE dsvbksb(uk,wk,vk,m,n,mp,np,bk,xk)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION bk(mp),uk(mp,np),vk(np,np),wk(np),xk(np)
      PARAMETER (NMAX=500)
      INTEGER i,j,jj
      DOUBLE PRECISION s,tmp(NMAX)
      do 12 j=1,n
        s=0.0d0
        if(wk(j).ne.0.0d0)then
          do 11 i=1,m
            s=s+uk(i,j)*bk(i)
11        continue
          s=s/wk(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.0d0
        do 13 jj=1,n
          s=s+vk(j,jj)*tmp(jj)
13      continue
        xk(j)=s
14    continue
      return
      END
