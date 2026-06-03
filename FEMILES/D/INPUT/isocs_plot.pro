pro isocs_plot
kt=2265
kp=206
zt=fltarr(kt,6)
zp=fltarr(kp,6)
;ze=fltarr(kt,6)
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='1.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,1]
xtit='Teff'
ytit='logg'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;x1=0.095
;x2=0.55
;y1=2801.
;y2=4500.
openr,lun,'Z0010T_ss_1Gyr',/get_lun
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
openr,lun,'Z0010_G_1Gyr',/get_lun
 for i=0,kp-1 do begin
  readf,lun,a,m,t,g,Lu,Lb,Lv
  zp(i,0)=m
  zp(i,1)=10^t
  zp(i,2)=g
  zp(i,3)=Lu
  zp(i,4)=Lb
  zp(i,5)=Lv
;  print,zp(i,1),zp(i,2)
 endfor
free_lun,lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;vlp=where(zp(*,0) lt ml)
;vlt=where(zt(*,0) lt ml)
;vle=where(ze(*,0) lt ml)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x1=11500
x2=2000
y1=5.7
y2=-1.2
;plot,zp(*,1),zp(*,2),xtitle=xtit,ytitle=ytit,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,linestyle=0,thick=1,charsize=1.,/NOERASE
plot,zt(*,1),zt(*,2),xtitle=xtit,ytitle=ytit,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,psym=4,symsize=0.5,charsize=1.6,/NOERASE
oplot,zp(*,1),zp(*,2),psym=5,symsize=0.5,color=3
xyouts,x1-700,0.5,'Padova00,1Gyr,Z=0.001',color=3
xyouts,x1-700,0.0,'BaSTI04,1Gyr,Z=0.001'

device,/close
loadct,0
;spawn,'gv 14.eps &'
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;oplot,zp(vlp,0),zp(vlp,1),color=2,linestyle=0
;oplot,ze(vle,0),ze(vle,1),color=2,linestyle=1
;oplot,zt(vlt,0),zt(vlt,1),color=2,linestyle=2
;xyouts,x1+0.019,4400.,'_____ Girardi et al. (2000)'
;xyouts,x1+0.019,4300.,'_ _ _ Cassisi et al. (2000)'
;xyouts,x1+0.019,4200.,'......... Pols et al. (1995)'
;xyouts,0.43,y1+100.,'_____ Z=0.020'
;xyouts,0.43,y1+200.,'_____ Z=0.008',color=2
;xyouts,0.43,y1+300.,'_____ Z=0.004',color=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
