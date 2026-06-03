pro aa14
kt=2265
zt=fltarr(kt,6)
zp=fltarr(kt,6)
ze=fltarr(kt,6)
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='aa14.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
xtit='Teff'
ytit='logg'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;x1=0.095
;x2=0.55
;y1=2801.
;y2=4500.
openr,lun,'Z0100T_aa_14',/get_lun
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
openr,lun,'Z0198T_aa_14',/get_lun
 for i=0,kt-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=m
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=Lu
  zp(i,4)=Lb
  zp(i,5)=Lv
 endfor
free_lun,lun
openr,lun,'Z0080T_aa_14',/get_lun
 for i=0,kt-1 do begin
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
;vlp=where(zp(*,0) lt ml)
;vlt=where(zt(*,0) lt ml)
;vle=where(ze(*,0) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x1=6000.
x2=2800.
y1=5.5
y2=-0.5
 plot,zp(*,1),zp(*,2),xtitle=xtit,ytitle=ytit,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,linestyle=0,thick=1,charsize=1.,/NOERASE
oplot,ze(*,1),ze(*,2),linestyle=1,color=3 ; Z0080
oplot,zt(*,1),zt(*,2),linestyle=2,color=2 ; Z0198
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
y11=3.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;oplot,zp(vlp,0),zp(vlp,1),color=2,linestyle=0
;oplot,ze(vle,0),ze(vle,1),color=2,linestyle=1
;oplot,zt(vlt,0),zt(vlt,1),color=2,linestyle=2
;xyouts,x1+0.019,4400.,'_____ Girardi et al. (2000)'
;xyouts,x1+0.019,4300.,'_ _ _ Cassisi et al. (2000)'
;xyouts,x1+0.019,4200.,'......... Pols et al. (1995)'
xyouts,4000,y11+.100,'_____ Z=0.0198 ae'
xyouts,4000,y11+.400,'_ _ _ Z=0.0100 ae',color=2
xyouts,4000,y11+.700,'......... Z=0.0080 ae',color=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
;spawn,'gv 14.eps &'
end
