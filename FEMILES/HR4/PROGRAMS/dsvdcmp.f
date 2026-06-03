      SUBROUTINE dsvdcmp(ak,m,n,mp,np,wk,vk,ntons)
c añado variable "ntons" por si falla dsvdcmp
      INTEGER ntons
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION ak(mp,np),vk(np,np),wk(np)
      PARAMETER (NMAX=500)
CU    USES dpythag
      INTEGER i,its,j,jj,k,l,nm
      DOUBLE PRECISION anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),dpythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+abs(ak(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              ak(k,i)=ak(k,i)/scale
              s=s+ak(k,i)*ak(k,i)
12          continue
            f=ak(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            ak(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+ak(k,i)*ak(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                ak(k,j)=ak(k,j)+f*ak(k,i)
14            continue
15          continue
            do 16 k=i,m
              ak(k,i)=scale*ak(k,i)
16          continue
          endif
        endif
        wk(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(ak(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              ak(i,k)=ak(i,k)/scale
              s=s+ak(i,k)*ak(i,k)
18          continue
            f=ak(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            ak(i,l)=f-g
            do 19 k=l,n
              rv1(k)=ak(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+ak(j,k)*ak(i,k)
21            continue
              do 22 k=l,n
                ak(j,k)=ak(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              ak(i,k)=scale*ak(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(wk(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              vk(j,i)=(ak(i,j)/ak(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+ak(i,k)*vk(k,j)
27            continue
              do 28 k=l,n
                vk(k,j)=vk(k,j)+s*vk(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            vk(i,j)=0.0d0
            vk(j,i)=0.0d0
31        continue
        endif
        vk(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=wk(i)
        do 33 j=l,n
          ak(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+ak(k,i)*ak(k,j)
34          continue
            f=(s/ak(i,i))*g
            do 35 k=i,m
              ak(k,j)=ak(k,j)+f*ak(k,i)
35          continue
36        continue
          do 37 j=i,m
            ak(j,i)=ak(j,i)*g
37        continue
        else
          do 38 j= i,m
            ak(j,i)=0.0d0
38        continue
        endif
        ak(i,i)=ak(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(wk(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=wk(i)
            h=dpythag(f,g)
            wk(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=ak(j,nm)
              z=ak(j,i)
              ak(j,nm)=(y*c)+(z*s)
              ak(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=wk(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              wk(k)=-z
              do 44 j=1,n
                vk(j,k)=-vk(j,k)
44            continue
            endif
            goto 3
          endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c lo añado yo por si falla dsvdcmp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          if(its.eq.30) STOP 'no convergence in svdcmp'
          if(its.eq.30) then
	   write(*,*)'no convergence in svdcmp'
	   ntons=1
	  else
	   ntons=0
	  endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          x=wk(l)
          nm=k-1
          y=wk(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=dpythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=wk(i)
            h=s*g
            g=c*g
            z=dpythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=vk(jj,j)
              z=vk(jj,i)
              vk(jj,j)= (x*c)+(z*s)
              vk(jj,i)=-(x*s)+(z*c)
45          continue
            z=dpythag(f,h)
            wk(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=ak(jj,j)
              z=ak(jj,i)
              ak(jj,j)= (y*c)+(z*s)
              ak(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          wk(k)=x
48      continue
3       continue
49    continue
      return
      END
