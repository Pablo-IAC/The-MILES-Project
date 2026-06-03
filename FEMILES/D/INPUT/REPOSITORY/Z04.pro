pro Z04
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
xtit='B-V' & ytit='V'
y2=-0.75 & y1=16.5
;y2=-0.75 & y1=10.
y2=8. & y1=16.5
x1=0.62 & x2=1.95
x1=1.2 & x2=1.95

xv=fltarr(1) & yv=fltarr(1)
xv(0)=.62 & yv(0)=0.51
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=.62 & yn(0)=16.
;
ya=fltarr(1)
ya(0)=15.
;
;tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;DEVICE,DECOMPOSED=0
loadct,2
erase
device,filename='Z04.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0198T_ss',/get_lun
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
zt10=where(zt(*,0) gt 13.99 and zt(*,0) lt 14.01)
zt1=where(zt(*,0) gt 11.999 and zt(*,0) lt 12.001)
openr,lun,'Z0240T_ss',/get_lun
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
zp10=where(zp(*,0) gt 13.99 and zp(*,0) lt 14.01)
zp1=where(zp(*,0) gt 11.999 and zp(*,0) lt 12.001)
;
plot,(zt(zt10,4)-zt(zt10,5)),zt(zt10,5),POSITION=[0.53,0.53,0.98,0.98],xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt10,4)-zt(zt10,5)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zp(zp10,4)-zp(zp10,5)),zp(zp10,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xn,yn,'14 Gyr, Z=0.0198',color=35
xyouts,xn,yn-2.,'14 Gyr Z=0.0240',color=135
plot,(zt(zt1,4)-zt(zt1,5)),zt(zt1,5),POSITION=[0.53,0.08,0.98,0.53],xcharsize=1,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt1,4)-zt(zt1,5)),zt(zt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zp(zp1,4)-zp(zp1,5)),zp(zp1,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xn,yn,'12 Gyr, Z=0.0198',color=35
xyouts,xn,yn-2.,'12 Gyr, Z=0.0240',color=135
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0300T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  mt(i,0)=a
  mt(i,1)=t
  mt(i,2)=g
  mt(i,3)=-2.5*ALOG10(Lu)
  mt(i,4)=-2.5*ALOG10(Lb)
  mt(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
mt10=where(mt(*,0) gt 13.99 and mt(*,0) lt 14.01)
mt1=where(mt(*,0) gt 11.999 and mt(*,0) lt 12.001)
openr,lun,'Z0400T_ss',/get_lun
 for i=0,kp01-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  mp(i,0)=a
  mp(i,1)=t
  mp(i,2)=g
  mp(i,3)=-2.5*ALOG10(Lu)
  mp(i,4)=-2.5*ALOG10(Lb)
  mp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
mp10=where(mp(*,0) gt 13.99 and mp(*,0) lt 14.01)
mp1=where(mp(*,0) gt 11.99 and mp(*,0) lt 12.001)
;
plot,(mt(mt10,4)-mt(mt10,5)),mt(mt10,5),POSITION=[0.08,0.53,0.53,0.98],xcharsize=.00001,ycharsize=1,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,(mt(mt10,4)-mt(mt10,5)),mt(mt10,5),PSYM=5,SYMSIZE=0.91,color=35
oplot,(mp(mp10,4)-mp(mp10,5)),mp(mp10,5),PSYM=4,SYMSIZE=0.91,color=135
oplot,(zt(zt10,4)-zt(zt10,5)),zt(zt10,5),PSYM=5,SYMSIZE=0.91,color=0 ; Z=0.0198
oplot,(zp(zp10,4)-zp(zp10,5)),zp(zp10,5),PSYM=4,SYMSIZE=0.91,color=235 ; Z=0.0240
xyouts,xn,yn,'14 Gyr, Z=0.0300',color=35
xyouts,xn,yn-2.,'14 Gyr, Z=0.0400',color=135
;
x111=3.6 & x222=3.42
plot,mt(mt10,1),mt(mt10,5),POSITION=[0.08,0.08,0.53,0.53],xcharsize=1,ycharsize=1,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x111,x222],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,mt(mt10,1),mt(mt10,5),PSYM=5,SYMSIZE=0.91,color=35
oplot,mp(mp10,1),mp(mp10,5),PSYM=4,SYMSIZE=0.91,color=135
oplot,zt(zt10,1),zt(zt10,5),PSYM=5,SYMSIZE=0.91,color=0 ; Z=0.0198
oplot,zp(zp10,1),zp(zp10,5),PSYM=4,SYMSIZE=0.91,color=235 ; Z=0.0240
;plot,(mt(mt1,4)-mt(mt1,5)),mt(mt1,5),POSITION=[0.08,0.08,0.53,0.53],xcharsize=1,ycharsize=1,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
;oplot,(mt(mt1,4)-mt(mt1,5)),mt(mt1,5),PSYM=5,SYMSIZE=0.1,color=35
;oplot,(mp(mp1,4)-mp(mp1,5)),mp(mp1,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xn,yn,'12 Gyr, Z=0.0300',color=35
xyouts,xn,yn-2.,'12 Gyr, Z=0.0400',color=135
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
;spawn,'gv Z04.eps &'
end
