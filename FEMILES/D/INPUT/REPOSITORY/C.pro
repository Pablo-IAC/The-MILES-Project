pro C
kp=8264
kt=122310
jk12=fltarr(1)
jk12(0)=1.2
zt=fltarr(kt,6) & zp=fltarr(kp,6)
;
kp04=8486
ct=fltarr(kt,6) & cp=fltarr(kp04,6)
;
kp01=8422
mt=fltarr(kt,6) & mp=fltarr(kp01,6)
;
xtit='J-K' & ytit='Mbol'
y2=-5.5 & y1=-1.5
x1=0.7 & x2=1.55
;
xv=fltarr(1) & yv=fltarr(1)
xv(0)=-.2 & yv(0)=-2.
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=x1+0.1 & yn(0)=y2+1.
xv=fltarr(1) & yv=fltarr(1)
xv(0)=xn(0) & yv(0)=yn(0)+1.
ya=fltarr(1)
ya(0)=yn(0)+1.
;
;tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;DEVICE,DECOMPOSED=0
loadct,2
erase
device,filename='C.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0198T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li,Lj,Lh,Lk,mi,fm,vk,mbol
  zt(i,0)=a
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lj)
  zt(i,4)=-2.5*ALOG10(Lk)
  zt(i,5)=mbol
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 9.99 and zt(*,0) lt 10.01)
zt1=where(zt(*,0) gt 0.999 and zt(*,0) lt 1.001)
zt01=where(zt(*,0) gt 1.49 and zt(*,0) lt 1.51)
;
plot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),POSITION=[0.08,0.68,0.38,0.98],xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(101670:101910,3)-zt(101670:101910,4)),zt(101670:101910,5),PSYM=4,SYMSIZE=0.1,color=90
;xyouts,xv+1.5,yv-2.,'BaSTI',color=35
xyouts,xn,yn,'Z=0.0198',color=0
xyouts,xv,ya,'10 Gyr',color=0
plot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),POSITION=[0.08,0.38,0.38,0.68],xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(47310:47550,3)-zt(47310:47550,4)),zt(47310:47550,5),PSYM=5,SYMSIZE=0.1,color=90
xyouts,xv,ya,'1 Gyr',color=0
plot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),POSITION=[0.08,0.08,0.38,0.38],xcharsize=1.,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(51840:52080,3)-zt(51840:52080,4)),zt(51840:52080,5),PSYM=4,SYMSIZE=0.1,color=90
xyouts,xv,ya,'1.5 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0100T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li,Lj,Lh,Lk,mi,fm,vk,mbol
  zt(i,0)=a
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lj)
  zt(i,4)=-2.5*ALOG10(Lk)
  zt(i,5)=mbol
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 9.99 and zt(*,0) lt 10.01)
zt1=where(zt(*,0) gt 0.999 and zt(*,0) lt 1.001)
zt01=where(zt(*,0) gt 1.49 and zt(*,0) lt 1.51)
;
plot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),POSITION=[0.38,0.68,0.68,0.98],ycharsize=0.0001,xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(101670:101910,3)-zt(101670:101910,4)),zt(101670:101910,5),PSYM=4,SYMSIZE=0.1,color=90
xyouts,xn,yn,'Z=0.0100',color=0
plot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),POSITION=[0.38,0.38,0.68,0.68],ycharsize=0.0001,xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(47310:47550,3)-zt(47310:47550,4)),zt(47310:47550,5),PSYM=5,SYMSIZE=0.1,color=90
plot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),POSITION=[0.38,0.08,0.68,0.38],ycharsize=0.0001,xcharsize=1.,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(51840:52080,3)-zt(51840:52080,4)),zt(51840:52080,5),PSYM=4,SYMSIZE=0.1,color=90
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0080T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li,Lj,Lh,Lk,mi,fm,vk,mbol
  zt(i,0)=a
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lj)
  zt(i,4)=-2.5*ALOG10(Lk)
  zt(i,5)=mbol
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 9.99 and zt(*,0) lt 10.01)
zt1=where(zt(*,0) gt 0.999 and zt(*,0) lt 1.001)
zt01=where(zt(*,0) gt 1.49 and zt(*,0) lt 1.51)
;
plot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),POSITION=[0.68,0.68,0.98,0.98],ycharsize=0.0001,xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(zt10,3)-zt(zt10,4)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(101670:101910,3)-zt(101670:101910,4)),zt(101670:101910,5),PSYM=4,SYMSIZE=0.1,color=90
xyouts,xn,yn,'Z=0.0080',color=0
plot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),POSITION=[0.68,0.38,0.98,0.68],ycharsize=0.0001,xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt1,3)-zt(zt1,4)),zt(zt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(47310:47550,3)-zt(47310:47550,4)),zt(47310:47550,5),PSYM=5,SYMSIZE=0.1,color=90
plot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),POSITION=[0.68,0.08,0.98,0.38],ycharsize=0.0001,xcharsize=1.,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt01,3)-zt(zt01,4)),zt(zt01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zt(51840:52080,3)-zt(51840:52080,4)),zt(51840:52080,5),PSYM=4,SYMSIZE=0.1,color=90
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
;spawn,'gv C.eps &'
end
