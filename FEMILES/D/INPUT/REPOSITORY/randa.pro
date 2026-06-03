pro randa
kt=10609
kt=10503
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
y2=-0.699 & y1=5.499
x1=4.499 & x2=3.399
;
xv=fltarr(1) & yv=fltarr(1)
xv(0)=3.92 & yv(0)=-0.5
;
xn=fltarr(1) & yn=fltarr(1)
xn(0)=4.42 & yn(0)=-0.2
;
ya=fltarr(1)
ya(0)=0.4
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='randa.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z008NNN',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=10.^a/1.e9
  zt(i,1)=t
  zt(i,2)=g
  zt(i,3)=-2.5*ALOG10(Lu)
  zt(i,4)=-2.5*ALOG10(Lb)
  zt(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zt10=where(zt(*,0) gt 0.0099 and zt(*,0) lt 0.0101)
openr,lun,'Z02NNN',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=10.^a/1.e9
  zp(i,1)=t
  zp(i,2)=g
  zp(i,3)=-2.5*ALOG10(Lu)
  zp(i,4)=-2.5*ALOG10(Lb)
  zp(i,5)=-2.5*ALOG10(Lv)
 endfor
free_lun,lun
zp10=where(zp(*,0) gt 0.0099 and zp(*,0) lt 0.0101)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
rsgt=fltarr(3)
rsgg=fltarr(3)
rsgt(0)=ALOG10(3550.) & rsgg(0)=-.135 ; rsg(0,2)=0.03 ; m0199V HD039801
rsgt(1)=ALOG10(3614.) & rsgg(1)=0.00 ; rsg(1,2)=0.42 ; m0211V HD042543
rsgt(2)=ALOG10(4117.) & rsgg(2)=0.20 ; rsg(2,2)=-0.2 ; m0239V HD052005
hd=strarr(3)
hd(0)='HD39801'
hd(1)='HD42543'
hd(2)='HD52005'
;
plot,zt(zt10,1),zt(zt10,2),POSITION=[0.15,0.15,0.95,0.95],PSYM=6,SYMSIZE=1.0,color=0,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=2,charsize=1.8,/NOERASE
oplot,zp(zp10,1),zp(zp10,2),PSYM=5,SYMSIZE=1.0,color=2
oplot,zt(zt10,1),zt(zt10,2),PSYM=6,SYMSIZE=1.0,color=0
oplot,xn,yn,PSYM=6,SYMSIZE=1.0,color=0
oplot,xn,ya,PSYM=5,SYMSIZE=1.0,color=2
;xyouts,3.62,3.8,'[Mg/Fe]=0.0',color=2
;xyouts,3.62,3.3,'[Mg/Fe]=0.4',color=3
xyouts,xn,yn+0.1,'  [M/H]=-0.4',color=0,charsize=1.8
xyouts,xn,ya+0.1,'  [M/H]=0.0',color=2,charsize=1.8
for jj=0,2 do begin
 xyouts,rsgt(jj)-0.012,rsgg(jj)+0.045,hd(jj),color=3
endfor
oplot,rsgt,rsgg,PSYM=2,SYMSIZE=1,color=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
end
