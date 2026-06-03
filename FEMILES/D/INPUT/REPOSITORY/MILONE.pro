pro MILONE
kt=925
a=' '
z=fltarr(kt,4)
;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
;xloadct
erase
device,filename='MILONE.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
!p.multi=[1,1,0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'PARAM_MILES.idl',/get_lun
 for i=0,kt-1 do begin
  readf,lun,t,g,fe,sn,mgfe
  z(i,0)=t
  z(i,1)=g
  z(i,2)=fe
  z(i,3)=mgfe
 endfor
free_lun,lun
ssgiants=where(z(*,1) lt 3.01 and z(*,3) lt 0.201)
ssdwarfs=where(z(*,1) ge 3.0 and z(*,3) lt 0.201)
aagiants=where(z(*,1) lt 3.01 and z(*,3) ge 0.20)
aadwarfs=where(z(*,1) ge 3.0 and z(*,3) ge 0.20)
dwarfs=where(z(*,1) gt 3.)
giants=where(z(*,1) le 3.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
yt='T!Deff!N!X'
;yt='Teff!N!X'
xt='[Fe/H]'
y1=2501.
y2=7999.
;y2=39999.
x1=-2.75
x2=0.79
x11=fltarr(1)
y11=fltarr(1)
x11(0)=x1
y11(0)=y1
xx=-2.55
yy=7600.
;yy=7500.
xx0=-2.55
yy0=3099.
yy00=2750.

plot,x11,y11,POSITION=[0.5725,0.10,0.9985,0.996],ycharsize=0.0001,symsize=0.01,psym=4,xtitle=xt,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.515,/NOERASE
oplot,z(ssdwarfs,2),z(ssdwarfs,0),symsize=0.6,psym=4,color=2
oplot,z(aadwarfs,2),z(aadwarfs,0),symsize=0.6,psym=5,color=3
xyouts,xx,yy,'Dwarfs',color=0,charsize=1.8
plot,x11,y11,POSITION=[0.1465,0.10,0.5725,0.996],ytitle=yt,symsize=0.1,psym=4,xtitle=xt,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.515,/NOERASE
oplot,z(ssgiants,2),z(ssgiants,0),symsize=0.6,psym=4,color=2
oplot,z(aagiants,2),z(aagiants,0),symsize=0.6,psym=5,color=3
xyouts,xx,yy,'Giants',color=0,charsize=1.8
xyouts,xx0,yy0,'[Mg/Fe]<0.2',color=2,charsize=1.6
xyouts,xx0,yy00,'[Mg/Fe]>0.2',color=3,charsize=1.6
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
device,/close
loadct,0
;spawn,'gv MILONE.eps'
end
;plot,z(aadwarfs,2),z(aadwarfs,0),POSITION=[0.10,0.08,0.54,0.53],xtitle=xt,ytitle=yt,symsize=0.1,psym=4,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
;oplot,z(aadwarfs,2),z(aadwarfs,0),symsize=0.2,psym=4,color=3
;plot,z(aagiants,2),z(aagiants,0),POSITION=[0.54,0.08,0.98,0.53],xtitle=xt,ycharsize=0.0001,symsize=0.1,psym=4,ytitle=yt,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,thick=1,charsize=1.,/NOERASE
;oplot,z(aagiants,2),z(aagiants,0),symsize=0.2,psym=4,color=3
;xyouts,x1+0.019,4350.,charsize=0.85,'___ Girardi et al. (2000)'
