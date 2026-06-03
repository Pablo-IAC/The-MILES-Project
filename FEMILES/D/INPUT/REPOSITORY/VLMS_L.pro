pro VLMS_L
ml=0.655
kt=24
kp=11
ke=93
zt=fltarr(kt,8,3)
zp=fltarr(kp,8,3)
ze=fltarr(ke,8,3)
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
x1=0.1
x2=0.55
device,filename='VLMS_L.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,2]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'Z0198T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zt(i,0,0)=m
  zt(i,1,0)=10^t
  zt(i,2,0)=g
  zt(i,3,0)=-2.5*alog10(Lu)
  zt(i,4,0)=-2.5*alog10(Lb)
  zt(i,5,0)=-2.5*alog10(Lv)
  zt(i,6,0)=-2.5*alog10(Lr)
  zt(i,7,0)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z0190_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zp(i,0,0)=m
  zp(i,1,0)=10^t
  zp(i,2,0)=g
  zp(i,3,0)=-2.5*alog10(Lu)
  zp(i,4,0)=-2.5*alog10(Lb)
  zp(i,5,0)=-2.5*alog10(Lv)
  zp(i,6,0)=-2.5*alog10(Lr)
  zp(i,7,0)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z02NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  ze(i,0,0)=m
  ze(i,1,0)=10^t
  ze(i,2,0)=g
  ze(i,3,0)=-2.5*alog10(Lu)
  ze(i,4,0)=-2.5*alog10(Lb)
  ze(i,5,0)=-2.5*alog10(Lv)
  ze(i,6,0)=-2.5*alog10(Lr)
  ze(i,7,0)=-2.5*alog10(Li)
 endfor
free_lun,lun
;
openr,lun,'Z0080T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zt(i,0,1)=m
  zt(i,1,1)=10^t
  zt(i,2,1)=g
  zt(i,3,1)=-2.5*alog10(Lu)
  zt(i,4,1)=-2.5*alog10(Lb)
  zt(i,5,1)=-2.5*alog10(Lv)
  zt(i,6,1)=-2.5*alog10(Lr)
  zt(i,7,1)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z0080_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zp(i,0,1)=m
  zp(i,1,1)=10^t
  zp(i,2,1)=g
  zp(i,3,1)=-2.5*alog10(Lu)
  zp(i,4,1)=-2.5*alog10(Lb)
  zp(i,5,1)=-2.5*alog10(Lv)
  zp(i,6,1)=-2.5*alog10(Lr)
  zp(i,7,1)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z008NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  ze(i,0,1)=m
  ze(i,1,1)=10^t
  ze(i,2,1)=g
  ze(i,3,1)=-2.5*alog10(Lu)
  ze(i,4,1)=-2.5*alog10(Lb)
  ze(i,5,1)=-2.5*alog10(Lv)
  ze(i,6,1)=-2.5*alog10(Lr)
  ze(i,7,1)=-2.5*alog10(Li)
 endfor
free_lun,lun
;
openr,lun,'Z0040T',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zt(i,0,2)=m
  zt(i,1,2)=10^t
  zt(i,2,2)=g
  zt(i,3,2)=-2.5*alog10(Lu)
  zt(i,4,2)=-2.5*alog10(Lb)
  zt(i,5,2)=-2.5*alog10(Lv)
  zt(i,6,2)=-2.5*alog10(Lr)
  zt(i,7,2)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z0040_G',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  zp(i,0,2)=m
  zp(i,1,2)=10^t
  zp(i,2,2)=g
  zp(i,3,2)=-2.5*alog10(Lu)
  zp(i,4,2)=-2.5*alog10(Lb)
  zp(i,5,2)=-2.5*alog10(Lv)
  zp(i,6,2)=-2.5*alog10(Lr)
  zp(i,7,2)=-2.5*alog10(Li)
 endfor
free_lun,lun
openr,lun,'Z004NNN',/get_lun
 for i=0,ke-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv,Lr,Li
  ze(i,0,2)=m
  ze(i,1,2)=10^t
  ze(i,2,2)=g
  ze(i,3,2)=-2.5*alog10(Lu)
  ze(i,4,2)=-2.5*alog10(Lb)
  ze(i,5,2)=-2.5*alog10(Lv)
  ze(i,6,2)=-2.5*alog10(Lr)
  ze(i,7,2)=-2.5*alog10(Li)
 endfor
free_lun,lun
;
vlp0=where(zp(*,0,0) lt ml)
vlt0=where(zt(*,0,0) lt ml)
vle0=where(ze(*,0,0) lt ml)
vlp1=where(zp(*,0,1) lt ml)
vlt1=where(zt(*,0,1) lt ml)
vle1=where(ze(*,0,1) lt ml)
vlp2=where(zp(*,0,2) lt ml)
vlt2=where(zt(*,0,2) lt ml)
vle2=where(ze(*,0,2) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ytit='T!Deff!N!X'
y1=2801.
y2=4500.
;
plot,zp(vlp0,0,0),zp(vlp0,1,0),POSITION=[0.1,0.5,0.9,0.9],xcharsize=0.0001,ytitle=ytit,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,linestyle=0,thick=1,charsize=1.,/NOERASE
oplot,ze(vle0,0,0),ze(vle0,1,0),linestyle=1
oplot,zt(vlt0,0,0),zt(vlt0,1,0),linestyle=2
xyouts,x1+0.019,4350.,charsize=0.85,'___ Girardi et al. (2000)'
xyouts,x1+0.019,4250.,charsize=0.85,'_ _ Cassisi et al. (2000)'
xyouts,x1+0.019,4150.,charsize=0.85,'..... Pols et al. (1995)'
oplot,zp(vlp1,0,1),zp(vlp1,1,1),color=2,linestyle=0
oplot,ze(vle1,0,1),ze(vle1,1,1),color=2,linestyle=1
oplot,zt(vlt1,0,1),zt(vlt1,1,1),color=2,linestyle=2
oplot,zp(vlp2,0,2),zp(vlp2,1,2),color=3,linestyle=0
oplot,ze(vle2,0,2),ze(vle2,1,2),color=3,linestyle=1
oplot,zt(vlt2,0,2),zt(vlt2,1,2),color=3,linestyle=2
xyouts,0.46,y1+100.,charsize=0.85,'___ Z=0.020'
xyouts,0.46,y1+200.,charsize=0.85,'___ Z=0.008',color=2
xyouts,0.46,y1+300.,charsize=0.85,'___ Z=0.004',color=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
xtit='!20!S!DO!R!I ' + string(183b) + '!X!N'
xtit='M/M!D!9n!N!X'
;xtit='M/Mo'
ytit='M!DV!N!X'
y2=7.75
y1=16.
;
 plot,zp(vlp0,0,0),zp(vlp0,5,0),POSITION=[0.1,0.1,0.9,0.5],xtitle=xtit,ytitle=ytit,yrange=[y1,y2],ystyle=1,xrange=[x1,x2],xstyle=1,linestyle=0,thick=1,charsize=1.,/NOERASE
oplot,ze(vle0,0,0),ze(vle0,5,0),linestyle=1
oplot,zt(vlt0,0,0),zt(vlt0,5,0),linestyle=2
oplot,zp(vlp1,0,1),zp(vlp1,5,1),color=2,linestyle=0
oplot,ze(vle1,0,1),ze(vle1,5,1),color=2,linestyle=1
oplot,zt(vlt1,0,1),zt(vlt1,5,1),color=2,linestyle=2
oplot,zp(vlp2,0,2),zp(vlp2,5,2),color=3,linestyle=0
oplot,ze(vle2,0,2),ze(vle2,5,2),color=3,linestyle=1
oplot,zt(vlt2,0,2),zt(vlt2,5,2),color=3,linestyle=2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
spawn,'gv VLMS_L.eps'
end
;
;
;
;xtit='M/Mo'
;ytit='M!sB!n'
;y2=9.01
;y1=19.99
;
; plot,zp(vlp0,0,0),zp(vlp0,3,0),xtitle=xtit,ytitle=ytit,yrange=[y1,y2],ystyle=1,xrange=[x1,x2],xstyle=1,linestyle=0,thick=1,charsize=1.
;oplot,ze(vle0,0,0),ze(vle0,3,0),linestyle=1
;oplot,zt(vlt0,0,0),zt(vlt0,3,0),linestyle=2
;oplot,zp(vlp1,0,1),zp(vlp1,3,1),color=2,linestyle=0
;oplot,ze(vle1,0,1),ze(vle1,3,1),color=2,linestyle=1
;oplot,zt(vlt1,0,1),zt(vlt1,3,1),color=2,linestyle=2
;oplot,zp(vlp2,0,2),zp(vlp2,3,2),color=3,linestyle=0
;oplot,ze(vle2,0,2),ze(vle2,3,2),color=3,linestyle=1
;oplot,zt(vlt2,0,2),zt(vlt2,3,2),color=3,linestyle=2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;xtit='M/Mo'
;ytit='M!sR!n'
;y2=8.5
;y1=18.5
;;
; plot,zp(vlp0,0,0),zp(vlp0,5,0),xtitle=xtit,ytitle=ytit,yrange=[y1,y2],ystyle=1,xrange=[x1,x2],xstyle=1,linestyle=0,thick=1,charsize=1.
;oplot,ze(vle0,0,0),ze(vle0,5,0),linestyle=1
;oplot,zt(vlt0,0,0),zt(vlt0,5,0),linestyle=2
;oplot,zp(vlp1,0,1),zp(vlp1,5,1),color=2,linestyle=0
;oplot,ze(vle1,0,1),ze(vle1,5,1),color=2,linestyle=1
;oplot,zt(vlt1,0,1),zt(vlt1,5,1),color=2,linestyle=2
;oplot,zp(vlp2,0,2),zp(vlp2,5,2),color=3,linestyle=0
;oplot,ze(vle2,0,2),ze(vle2,5,2),color=3,linestyle=1
;oplot,zt(vlt2,0,2),zt(vlt2,5,2),color=3,linestyle=2

;print,ze(vle0,0,0),ze(vle0,1,0),vle0
