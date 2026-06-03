pro ISOC_SS_AA
kt=101925
kt=120045
kp=kt
zt=fltarr(kt,6) & zp=fltarr(kp,6)
;
kp04=kt
ct=fltarr(kt,6) & cp=fltarr(kp04,6)
;
kp01=kt
mt=fltarr(kt,6) & mp=fltarr(kp01,6)
;
xtit='log Teff' & ytit='log g'
y2=-0.99 & y1=5.4
x1=3.96 & x2=3.39
;
xv=fltarr(1) & yv=fltarr(1)
xv(0)=3.92 & yv(0)=-0.5
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=3.92 & yn(0)=-0.5
;
ya=fltarr(1)
ya(0)=0.
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='ISOC_SS_AA.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0240T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=a
  zt(i,1)=t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lu)
  zt(i,4)=-2.5*ALOG10(Lb)
  zt(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 11.99 and zt(*,0) lt 12.01)
zt1=where(zt(*,0) gt 1.999 and zt(*,0) lt 2.001)
openr,lun,'Z0240T_aa',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=a
  zp(i,1)=t
  zp(i,2)=g
  zp(i,3)=-2.5*ALOG10(Lu)
  zp(i,4)=-2.5*ALOG10(Lb)
  zp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zp10=where(zp(*,0) gt 11.99 and zp(*,0) lt 12.01)
zp1=where(zp(*,0) gt 1.999 and zp(*,0) lt 2.001)
;
plot,zt(zt10,1),zt(zt10,2),POSITION=[0.08,0.53,0.53,0.98],xcharsize=0.0001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,zt(zt10,1),zt(zt10,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,zp(zp10,1),zp(zp10,2),PSYM=4,SYMSIZE=0.1,color=3
;oplot,xv,yv,PSYM=5,SYMSIZE=0.1,color=0
xyouts,3.62,3.8,'[Mg/Fe]=0.0',color=2
xyouts,3.62,3.3,'[Mg/Fe]=0.4',color=3
xyouts,xn,yn,'Z=0.0240',color=0
xyouts,xv,ya,'12 Gyr',color=0
plot,zt(zt1,1),zt(zt1,2),POSITION=[0.08,0.08,0.53,0.53],xcharsize=1.,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,zt(zt1,1),zt(zt1,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,zp(zp1,1),zp(zp1,2),PSYM=4,SYMSIZE=0.1,color=3
xyouts,xv,ya,'2 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0010T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  ct(i,0)=a
  ct(i,1)=t
  ct(i,2)=g
  ct(i,3)=-2.5*ALOG10(Lu)
  ct(i,4)=-2.5*ALOG10(Lb)
  ct(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
ct10=where(ct(*,0) gt 11.99 and ct(*,0) lt 12.01)
ct1=where(ct(*,0) gt 1.999 and ct(*,0) lt 2.001)
openr,lun,'Z0010T_aa',/get_lun
 for i=0,kp04-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  cp(i,0)=a
  cp(i,1)=t
  cp(i,2)=g
  cp(i,3)=-2.5*ALOG10(Lu)
  cp(i,4)=-2.5*ALOG10(Lb)
  cp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
cp10=where(cp(*,0) gt 11.99 and cp(*,0) lt 12.01)
cp1=where(cp(*,0) gt 1.999 and cp(*,0) lt 2.001)
;
plot,ct(ct10,1),ct(ct10,2),POSITION=[0.53,0.53,0.98,0.98],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,ct(ct10,1),ct(ct10,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,cp(cp10,1),cp(cp10,2),PSYM=4,SYMSIZE=0.1,color=3
xyouts,xv,ya,'12 Gyr',color=0
xyouts,xn,yn,'Z=0.001',color=0
plot,ct(ct1,1),ct(ct1,2),POSITION=[0.53,0.08,0.98,0.53],xcharsize=1.,ycharsize=.000001,PSYM=5,SYMSIZE=0.1,color=0,xrange=[x1,x2],xstyle=1,xtitle=xtit,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,ct(ct1,1),ct(ct1,2),PSYM=5,SYMSIZE=0.1,color=2
oplot,cp(cp1,1),cp(cp1,2),PSYM=4,SYMSIZE=0.1,color=3
xyouts,xv,ya,'2 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
end
