pro VLMS
ml=0.655
kt=24
kp=11
ke=93
zt=fltarr(kt,6)
zp=fltarr(kp,6)
ze=fltarr(ke,6)
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='VLMS.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
xtit='M/Mo'
ytit='Teff'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x1=0.095
x2=0.55
y1=2801.
y2=4500.
openr,lun,'Z0198T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=m
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=Lu
  zt(i,4)=Lb
  zt(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z0190_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=m
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=Lu
  zp(i,4)=Lb
  zp(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z02NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  ze(i,0)=m
  ze(i,1)=10^t
  ze(i,2)=g
  ze(i,3)=Lu
  ze(i,4)=Lb
  ze(i,5)=Lv
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vlp=where(zp(*,0) lt ml)
vlt=where(zt(*,0) lt ml)
vle=where(ze(*,0) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot,zp(vlp,0),zp(vlp,1),xtitle=xtit,ytitle=ytit,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,linestyle=0,thick=1,charsize=1.,/NOERASE
oplot,ze(vle,0),ze(vle,1),linestyle=1
oplot,zt(vlt,0),zt(vlt,1),linestyle=2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0080T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=m
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=Lu
  zt(i,4)=Lb
  zt(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z0080_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=m
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=Lu
  zp(i,4)=Lb
  zp(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z008NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  ze(i,0)=m
  ze(i,1)=10^t
  ze(i,2)=g
  ze(i,3)=Lu
  ze(i,4)=Lb
  ze(i,5)=Lv
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vlp=where(zp(*,0) lt ml)
vlt=where(zt(*,0) lt ml)
vle=where(ze(*,0) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
oplot,zp(vlp,0),zp(vlp,1),color=2,linestyle=0
oplot,ze(vle,0),ze(vle,1),color=2,linestyle=1
oplot,zt(vlt,0),zt(vlt,1),color=2,linestyle=2
xyouts,x1+0.019,4400.,'_____ Girardi et al. (2000)'
xyouts,x1+0.019,4300.,'_ _ _ Cassisi et al. (2000)'
xyouts,x1+0.019,4200.,'......... Pols et al. (1995)'
xyouts,0.43,y1+100.,'_____ Z=0.020'
xyouts,0.43,y1+200.,'_____ Z=0.008',color=2
xyouts,0.43,y1+300.,'_____ Z=0.004',color=3
openr,lun,'Z0040T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zt(i,0)=m
  zt(i,1)=10^t
  zt(i,2)=g
  zt(i,3)=Lu
  zt(i,4)=Lb
  zt(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z0040_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=m
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=Lu
  zp(i,4)=Lb
  zp(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z004NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  ze(i,0)=m
  ze(i,1)=10^t
  ze(i,2)=g
  ze(i,3)=Lu
  ze(i,4)=Lb
  ze(i,5)=Lv
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vlp=where(zp(*,0) lt ml)
vlt=where(zt(*,0) lt ml)
vle=where(ze(*,0) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
oplot,zp(vlp,0),zp(vlp,1),color=3,linestyle=0
oplot,ze(vle,0),ze(vle,1),color=3,linestyle=1
oplot,zt(vlt,0),zt(vlt,1),color=3,linestyle=2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
spawn,'gv VLMS.eps'
end
