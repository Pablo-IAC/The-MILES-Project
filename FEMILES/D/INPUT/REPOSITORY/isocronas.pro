pro isocronas
kt=122310
z1=fltarr(kt,6) & z2=fltarr(kt,6) & z3=fltarr(kt,6) ; a/Fe=0.0
r1=fltarr(kt,6) & r2=fltarr(kt,6) & r3=fltarr(kt,6) ; a/Fe=0.4
kp=8817
p1=fltarr(kp,6) & p2=fltarr(kp,6) & p3=fltarr(kp,6) ; a/Fe=0.0 Padova00
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0020T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  z1(i,0)=a
  z1(i,1)=t
  z1(i,2)=g
  z1(i,3)=-2.5*ALOG10(Lu)
  z1(i,4)=-2.5*ALOG10(Lb)
  z1(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0040T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  z2(i,0)=a
  z2(i,1)=t
  z2(i,2)=g
  z2(i,3)=-2.5*ALOG10(Lu)
  z2(i,4)=-2.5*ALOG10(Lb)
  z2(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0080T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  z3(i,0)=a
  z3(i,1)=t
  z3(i,2)=g
  z3(i,3)=-2.5*ALOG10(Lu)
  z3(i,4)=-2.5*ALOG10(Lb)
  z3(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0020T_aa',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  r1(i,0)=a
  r1(i,1)=t
  r1(i,2)=g
  r1(i,3)=-2.5*ALOG10(Lu)
  r1(i,4)=-2.5*ALOG10(Lb)
  r1(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0040T_aa',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  r2(i,0)=a
  r2(i,1)=t
  r2(i,2)=g
  r2(i,3)=-2.5*ALOG10(Lu)
  r2(i,4)=-2.5*ALOG10(Lb)
  r2(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0080T_aa',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  r3(i,0)=a
  r3(i,1)=t
  r3(i,2)=g
  r3(i,3)=-2.5*ALOG10(Lu)
  r3(i,4)=-2.5*ALOG10(Lb)
  r3(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0010_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  p1(i,0)=a
  p1(i,1)=t
  p1(i,2)=g
  p1(i,3)=-2.5*ALOG10(Lu)
  p1(i,4)=-2.5*ALOG10(Lb)
  p1(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0040_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  p2(i,0)=a
  p2(i,1)=t
  p2(i,2)=g
  p2(i,3)=-2.5*ALOG10(Lu)
  p2(i,4)=-2.5*ALOG10(Lb)
  p2(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;
openr,lun,'Z0080_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  p3(i,0)=a
  p3(i,1)=t
  p3(i,2)=g
  p3(i,3)=-2.5*ALOG10(Lu)
  p3(i,4)=-2.5*ALOG10(Lb)
  p3(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ageol=0.099 & ageoh=0.101 ; edad old a plotear
ageyl=0.049 & ageyh=0.051 ; edad young a plotear
;
z1o=where(z1(*,0) gt ageol and z1(*,0) lt ageoh)
z1y=where(z1(*,0) gt ageyl and z1(*,0) lt ageyh)
z2o=where(z2(*,0) gt ageol and z2(*,0) lt ageoh)
z2y=where(z2(*,0) gt ageyl and z2(*,0) lt ageyh)
z3o=where(z3(*,0) gt ageol and z3(*,0) lt ageoh)
z3y=where(z3(*,0) gt ageyl and z3(*,0) lt ageyh)
;
r1o=where(r1(*,0) gt ageol and r1(*,0) lt ageoh)
r1y=where(r1(*,0) gt ageyl and r1(*,0) lt ageyh)
r2o=where(r2(*,0) gt ageol and r2(*,0) lt ageoh)
r2y=where(r2(*,0) gt ageyl and r2(*,0) lt ageyh)
r3o=where(r3(*,0) gt ageol and r3(*,0) lt ageoh)
r3y=where(r3(*,0) gt ageyl and r3(*,0) lt ageyh)
;
ageol=7.79 & ageoh=7.81 ; edad old a plotear
ageyl=7.99 & ageyh=8.01 ; edad young a plotear
;
p1o=where(p1(*,0) gt ageol and p1(*,0) lt ageoh)
p1y=where(p1(*,0) gt ageyl and p1(*,0) lt ageyh)
p2o=where(p2(*,0) gt ageol and p2(*,0) lt ageoh)
p2y=where(p2(*,0) gt ageyl and p2(*,0) lt ageyh)
p3o=where(p3(*,0) gt ageol and p3(*,0) lt ageoh)
p3y=where(p3(*,0) gt ageyl and p3(*,0) lt ageyh)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
set_plot,'ps'
;xloadct
erase
device,filename='isocronas.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;
xtit='log Teff' & ytit='log g'
y2=-1.50 & y1=5.4
y2=-1.5 & y1=5.5
x1=4.35 & x2=3.35
xv=fltarr(1) & yv=fltarr(1)
xv(0)=4.2 & yv(0)=-0.3
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=4.2 & yn(0)=-0.3
;
ya=fltarr(1) & ya(0)=-0.3
;;;;;;;;;;;;;;;;;;;
plot,z1(z1o,1),z1(z1o,2),POSITION=[0.08,0.53,0.53,0.98],xcharsize=0.0001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,z2(z2o,1),z2(z2o,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,z3(z3o,1),z3(z3o,2),PSYM=5,SYMSIZE=0.1,color=3
;
oplot,p1(p1o,1),p2(p2o,2),PSYM=5,SYMSIZE=0.3,color=4
oplot,p2(p2o,1),p2(p2o,2),PSYM=5,SYMSIZE=0.3,color=2
oplot,p3(p3o,1),p3(p3o,2),PSYM=5,SYMSIZE=0.3,color=3
;
;xyouts,3.62,3.8,'[Mg/Fe]=0.0',color=2
;xyouts,3.62,3.3,'[Mg/Fe]=0.4',color=3
dd=0.45
xyouts,xn,yn+dd,'Z=0.0020',color=0
xyouts,xn,yn+2.*dd,'Z=0.0040',color=2
xyouts,xn,yn+3.*dd,'Z=0.0080',color=3
xyouts,xv,ya,'[Mg/Fe]=0.0, 0.1Gyr',color=0
;;;;;;;;;;;;;;;;;;;
plot,z1(z1y,1),z1(z1y,2),POSITION=[0.08,0.08,0.53,0.53],xcharsize=1.,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,z2(z2y,1),z2(z2y,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,z3(z3y,1),z3(z3y,2),PSYM=4,SYMSIZE=0.1,color=3
xyouts,xv,ya,'[Mg/Fe]=0.0, 0.05Gyr',color=0
;;;;;;;;;;;;;;;;;;;
plot,r1(r1o,1),r1(r1o,2),POSITION=[0.53,0.53,0.98,0.98],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,r2(r2o,1),r2(r2o,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,r3(r3o,1),r3(r3o,2),PSYM=5,SYMSIZE=0.1,color=3
xyouts,xv,ya,'[Mg/Fe]=0.4, 0.1Gyr',color=0
;;;;;;;;;;;;;;;;;;;;
plot,r1(r1y,1),r1(r1y,2),POSITION=[0.53,0.08,0.98,0.53],xcharsize=1.,ycharsize=.000001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,xtitle=xtit,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,r2(r2y,1),r2(r2y,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,r3(r3y,1),r3(r3y,2),PSYM=5,SYMSIZE=0.1,color=3
xyouts,xv,ya,'[Mg/Fe]=0.4, 0.05Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
end
