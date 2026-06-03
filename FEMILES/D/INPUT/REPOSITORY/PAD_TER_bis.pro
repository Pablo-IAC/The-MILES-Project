pro PAD_TER_bis
kt=101925
kp=8264
zt=fltarr(kt,6) & zp=fltarr(kp,6)
;
kp04=8486
ct=fltarr(kt,6) & cp=fltarr(kp04,6)
;
kp01=8422
mt=fltarr(kt,6) & mp=fltarr(kp01,6)
;
xtit='B-V' & ytit='V'
y2=-7. & y1=17.
x1=-.4 & x2=2.4
;
xv=fltarr(1) & yv=fltarr(1)
xv(0)=-.2 & yv(0)=-2.
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=-.2 & yn(0)=-2.
;
ya=fltarr(1)
ya(0)=15.
;
;tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;DEVICE,DECOMPOSED=0
loadct,2
erase
device,filename='PAD_TER_bis.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0198T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=a
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lu)
  zt(i,4)=-2.5*ALOG10(Lb)
  zt(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 9.99 and zt(*,0) lt 10.01)
zt1=where(zt(*,0) gt 0.999 and zt(*,0) lt 1.001)
zt01=where(zt(*,0) gt .0999 and zt(*,0) lt .101)
openr,lun,'Z0190_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=a
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=-2.5*ALOG10(Lu)
  zp(i,4)=-2.5*ALOG10(Lb)
  zp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zp10=where(zp(*,0) gt 9.99 and zp(*,0) lt 10.01)
zp1=where(zp(*,0) gt 8.999 and zp(*,0) lt 9.001)
zp01=where(zp(*,0) gt 7.999 and zp(*,0) lt 8.001)
;
plot,(zt(zt10,4)-zt(zt10,5)),zt(zt10,5),POSITION=[0.08,0.68,0.38,0.98],xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt10,4)-zt(zt10,5)),zt(zt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zp(zp10,4)-zp(zp10,5)),zp(zp10,5),PSYM=4,SYMSIZE=0.1,color=135
;oplot,xv,yv,PSYM=5,SYMSIZE=0.1,color=0
xyouts,xv+1.5,yv-2.,'BaSTI',color=35
;oplot,xv+0.05,yv-2.,PSYM=4,SYMSIZE=0.1,color=2
xyouts,xv+1.5,yv,'Padova00',color=135
xyouts,xn,yn,'Z=0.0190',color=0
xyouts,xn,yn-2.,'Z=0.0198',color=0
xyouts,xv,ya,'10 Gyr',color=0
plot,(zt(zt1,4)-zt(zt1,5)),zt(zt1,5),POSITION=[0.08,0.38,0.38,0.68],xcharsize=0.0001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt1,4)-zt(zt1,5)),zt(zt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zp(zp1,4)-zp(zp1,5)),zp(zp1,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'1 Gyr',color=0
plot,(zt(zt01,4)-zt(zt01,5)),zt(zt01,5),POSITION=[0.08,0.08,0.38,0.38],xcharsize=1.,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,(zt(zt01,4)-zt(zt01,5)),zt(zt01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(zp(zp01,4)-zp(zp01,5)),zp(zp01,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'0.1 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0040T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  ct(i,0)=a
  ct(i,1)=10^t
  ct(i,2)=g
  ct(i,3)=-2.5*ALOG10(Lu)
  ct(i,4)=-2.5*ALOG10(Lb)
  ct(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
ct10=where(ct(*,0) gt 9.99 and ct(*,0) lt 10.01)
ct1=where(ct(*,0) gt 0.999 and ct(*,0) lt 1.001)
ct01=where(ct(*,0) gt .0999 and ct(*,0) lt .101)
openr,lun,'Z0040_G',/get_lun
 for i=0,kp04-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  cp(i,0)=a
  cp(i,1)=10^t
  cp(i,2)=g
  cp(i,3)=-2.5*ALOG10(Lu)
  cp(i,4)=-2.5*ALOG10(Lb)
  cp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
cp10=where(cp(*,0) gt 9.99 and cp(*,0) lt 10.01)
cp1=where(cp(*,0) gt 8.999 and cp(*,0) lt 9.001)
cp01=where(cp(*,0) gt 7.999 and cp(*,0) lt 8.001)
;
plot,(ct(ct10,4)-ct(ct10,5)),ct(ct10,5),POSITION=[0.38,0.68,0.68,0.98],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,(ct(ct10,4)-ct(ct10,5)),ct(ct10,5),PSYM=4,SYMSIZE=0.1,color=35
oplot,(cp(cp10,4)-cp(cp10,5)),cp(cp10,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'10 Gyr',color=0
xyouts,xn,yn-2.,'Z=0.004',color=0
plot,(ct(ct1,4)-ct(ct1,5)),ct(ct1,5),POSITION=[0.38,0.38,0.68,0.68],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,(ct(ct1,4)-ct(ct1,5)),ct(ct1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(cp(cp1,4)-cp(cp1,5)),cp(cp1,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'1 Gyr',color=0
plot,(ct(ct01,4)-ct(ct01,5)),ct(ct01,5),POSITION=[0.38,0.08,0.68,0.38],xcharsize=1.,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,thick=1,charsize=1.,/NOERASE
oplot,(ct(ct01,4)-ct(ct01,5)),ct(ct01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(cp(cp01,4)-cp(cp01,5)),cp(cp01,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'0.1 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0010T_ss',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  mt(i,0)=a
  mt(i,1)=10^t
  mt(i,2)=g
  mt(i,3)=-2.5*ALOG10(Lu)
  mt(i,4)=-2.5*ALOG10(Lb)
  mt(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
mt10=where(mt(*,0) gt 9.99 and mt(*,0) lt 10.01)
mt1=where(mt(*,0) gt 0.999 and mt(*,0) lt 1.001)
mt01=where(mt(*,0) gt .0999 and mt(*,0) lt .101)
openr,lun,'Z0010_G',/get_lun
 for i=0,kp01-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  mp(i,0)=a
  mp(i,1)=10^t
  mp(i,2)=g
  mp(i,3)=-2.5*ALOG10(Lu)
  mp(i,4)=-2.5*ALOG10(Lb)
  mp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
mp10=where(mp(*,0) gt 9.99 and mp(*,0) lt 10.01)
mp1=where(mp(*,0) gt 8.999 and mp(*,0) lt 9.001)
mp01=where(mp(*,0) gt 7.999 and mp(*,0) lt 8.001)
;
plot,(mt(mt10,4)-mt(mt10,5)),mt(mt10,5),POSITION=[0.68,0.68,0.98,0.98],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,(mt(mt10,4)-mt(mt10,5)),mt(mt10,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(mp(mp10,4)-mp(mp10,5)),mp(mp10,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xn,yn-2.,'Z=0.001',color=0
xyouts,xv,ya,'10 Gyr',color=0
plot,(mt(mt1,4)-mt(mt1,5)),mt(mt1,5),POSITION=[0.68,0.38,0.98,0.68],xcharsize=0.0001,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
oplot,(mt(mt1,4)-mt(mt1,5)),mt(mt1,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(mp(mp1,4)-mp(mp1,5)),mp(mp1,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'1 Gyr',color=0
plot,(mt(mt01,4)-mt(mt01,5)),mt(mt01,5),POSITION=[0.68,0.08,0.98,0.38],xcharsize=1.,ycharsize=.000001,PSYM=5,SYMSIZE=0.0001,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,thick=1,charsize=1.,/NOERASE
oplot,(mt(mt01,4)-mt(mt01,5)),mt(mt01,5),PSYM=5,SYMSIZE=0.1,color=35
oplot,(mp(mp01,4)-mp(mp01,5)),mp(mp01,5),PSYM=4,SYMSIZE=0.1,color=135
xyouts,xv,ya,'0.1 Gyr',color=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
;spawn,'gv PAD_TER_bis.eps &'
end
