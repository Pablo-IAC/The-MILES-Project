pro NGSL
l=378
l=317
lp=66
lp=5
z=fltarr(l,3)
zp=fltarr(lp,3)
openr,lun,'PARAM_NGSL_2012.idl',/get_lun
for i=0,l-1 do begin
 readf,lun,t,g,f
 z(i,0)=t
 z(i,1)=g
 z(i,2)=f
endfor
free_lun,lun
zp=fltarr(l,3)
openr,lun,'PARAM_NGSL_2012_weight_10',/get_lun
for i=0,lp-1 do begin
 readf,lun,t,g,f
 zp(i,0)=t
 zp(i,1)=g
 zp(i,2)=f
endfor
free_lun,lun
;
giants=where(z(*,1) lt 3.0 and z(*,1) ge 1.5)
Pgiants=where(zp(*,1) lt 3.0 and zp(*,1) ge 1.5)
supergiants=where(z(*,1) lt 1.5)
Psupergiants=where(zp(*,1) lt 1.5)
dwarfs=where(z(*,1) ge 3.0)
Pdwarfs=where(zp(*,1) ge 3.0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xsc=[-1.,1.,1.,-1.,-1.] & ysc=[1.,1.,-1.,-1.,1.]
usersym,xsc,ysc,/fill
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!p.multi=[4,2,2]
set_plot,'ps'
device,filename='NGSL.eps',xsize=20,ysize=24,yoffset=1.,/color
erase
;
x1=-2.3
x2=1.0
y1=10000.
y2=33000.
xtit='[Fe/H](Dwarfs)'
ytit='[Teff](Dwarfs)'
plot,z(dwarfs,2),z(dwarfs,0),xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,$
ytitle=ytit,xtitle=xtit,psym=4,symsize=0.65
oplot,zp(Pdwarfs,2),zp(Pdwarfs,0),psym=4,color=2
;
y1=3000.
y2=10000.
xtit='[Fe/H](Dwarfs)'
ytit='[Teff](Dwarfs)'
plot,z(dwarfs,2),z(dwarfs,0),xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,$
ytitle=ytit,xtitle=xtit,psym=4,symsize=0.65
oplot,zp(Pdwarfs,2),zp(Pdwarfs,0),psym=4,color=2
;
;
y1=3000.
y2=9000.
xtit='[Fe/H](Giants)'
ytit='[Teff](Giants)'
plot,z(giants,2),z(giants,0),xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,$
ytitle=ytit,xtitle=xtit,symsize=0.65,psym=4
oplot,zp(Pgiants,2),zp(Pgiants,0),psym=4,color=2
;
;
y1=3000.
y2=7000.
xtit='[Fe/H](Supergiants)'
ytit='[Teff](Supergiants)'
plot,z(supergiants,2),z(supergiants,0),xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,$
ytitle=ytit,xtitle=xtit,symsize=0.65,psym=4
oplot,zp(Psupergiants,2),zp(Psupergiants,0),psym=4,color=2
;
!p.multi=0
device,/close
set_plot,'x'
end
