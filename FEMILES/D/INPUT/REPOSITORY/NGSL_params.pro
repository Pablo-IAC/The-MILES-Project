pro NGSL_params
km=135
k311=311
ko=176
z=fltarr(km,6)
r=fltarr(km,3)
zz=fltarr(ko,3)
rr=fltarr(ko,3)
ttf=[-2.25,-2.0,-1.75,-1.5,-1.25,-1.,-.75,-.5,-.25,.0,.25,.5]
ttt=[3000.,4000.,5000.,6000.,7000.,8000.]
xtit='Teff' & ytit=''
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
DEVICE,DECOMPOSED=0
;loadct,2
erase
device,filename='NGSL_params.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
;!p.multi=[1,2,2]
!p.multi=[1,1,1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
openr,lun,'PARAM_NGSL_23_08_2012.idl',/get_lun
 for i=0,km-1 do begin
  readf,lun,t,g,f,s1,s2,s3,s4,to,go,fo
  z(i,0)=t
  z(i,1)=g
  z(i,2)=f
  z(i,3)=to
  z(i,4)=go
  z(i,5)=fo
 endfor
free_lun,lun
giants=where(z(*,1) le 3.0 and z(*,3) lt 8000.)
dwarfs=where(z(*,1) gt 3.0 and z(*,3) lt 10000.)
;
openr,lun,'PARAM_NGSL_23_08_2012.idl',/get_lun
 for i=km,k311-1 do begin
  readf,lun,t,g,f,s1,s2,s3,s4,to,go,fo
  zz(i-km,0)=to
  zz(i-km,1)=go
  zz(i-km,2)=fo
 endfor
free_lun,lun
xtit='Teff' & ytit='Teff - Teff'
;
;print,z(dwarfs,3)-z(dwarfs,0)
;
y2=500. & y1=-400.
x1=-2.25 & x2=0.75
xtit='[Fe/H]' & ytit='Teff(NGSL) - Teff(MILES)'
plot,z(dwarfs,2),(z(dwarfs,3)-z(dwarfs,0)),position=[0.5,0.5,0.9,0.9],PSYM=5,SYMSIZE=0.5,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,thick=1,charsize=1.,/NOERASE
oplot,z(giants,2),(z(giants,3)-z(giants,0)),PSYM=5,SYMSIZE=0.5,color=3
coefd = poly_fit(z(dwarfs,2),(z(dwarfs,3)-z(dwarfs,0)),1)
oplot,ttf,poly(ttf,coefd),linestyle=1,color=0
coefg = poly_fit(z(giants,2),(z(giants,3)-z(giants,0)),1)
oplot,ttf,poly(ttf,coefg),linestyle=1,color=3
print,'coefd',coefd(0),coefd(1)
print,'-----------------------'
print,'coefg',coefg(0),coefg(1)
for i=0,ko-1 do begin
 rr(i,2)=zz(i,2)
 rr(i,1)=zz(i,1)
 rr(i,0)=zz(i,0)
endfor
for i=0,km-1 do begin
 r(i,2)=z(i,2)
 r(i,1)=z(i,1)
 r(i,0)=z(i,3)
endfor
giantsrr=where(rr(*,1) le 3.0)
dwarfsrr=where(rr(*,1) gt 3.0)
giantsm=where(r(*,1) le 3.0)
dwarfsm=where(r(*,1) gt 3.0)
rr(giantsrr,0)=rr(giantsrr,0)-(coefg(0)+rr(giantsrr,2)*coefg(1))
oplot,rr(giantsrr,2),(rr(giantsrr,0)-zz(giantsrr,0)),PSYM=4,SYMSIZE=0.5,color=2
rr(dwarfsrr,0)=rr(dwarfsrr,0)-(coefd(0)+rr(dwarfsrr,2)*coefd(1))
oplot,rr(dwarfsrr,2),(rr(dwarfsrr,0)-zz(dwarfsrr,0)),PSYM=3,SYMSIZE=0.5,color=4
r(giantsm,0)=r(giantsm,0)-(coefg(0)+r(giantsm,2)*coefg(1))
oplot,r(giantsm,2),(r(giantsm,0)-z(giantsm,0)),PSYM=4,SYMSIZE=0.5,color=2
r(dwarfsm,0)=r(dwarfsm,0)-(coefd(0)+r(dwarfsm,2)*coefd(1))
oplot,r(dwarfsm,2),(r(dwarfsm,0)-z(dwarfsm,0)),PSYM=3,SYMSIZE=0.5,color=4
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x1=2500. & x2=10000.
plot,z(dwarfs,3),(z(dwarfs,3)-z(dwarfs,0)),position=[0.1,0.5,0.5,0.9],PSYM=5,SYMSIZE=0.5,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,z(giants,3),(z(giants,3)-z(giants,0)),PSYM=5,SYMSIZE=0.5,color=3
coef = poly_fit(z(dwarfs,3),(z(dwarfs,3)-z(dwarfs,0)),1)
oplot,ttt,poly(ttt,coef),linestyle=1,color=0
coef = poly_fit(z(giants,3),(z(giants,3)-z(giants,0)),1)
oplot,ttt,poly(ttt,coef),linestyle=1,color=3
oplot,zz(giantsrr,0),(rr(giantsrr,0)-zz(giantsrr,0)),PSYM=5,SYMSIZE=0.5,color=2
oplot,zz(dwarfsrr,0),(rr(dwarfsrr,0)-zz(dwarfsrr,0)),PSYM=5,SYMSIZE=0.5,color=4
;;;;;;;;;;;;;;;;;;CORREGIDO Fe/H;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
y2=300. & y1=-300.
x1=-2.25 & x2=0.75
xtit='[Fe/H]' & ytit='Teff(NGSL) - Teff(MILES)'
plot,z(dwarfs,2),(r(dwarfs,0)-z(dwarfs,0)),position=[0.5,0.1,0.9,0.5],PSYM=5,SYMSIZE=0.5,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,thick=1,charsize=1.,/NOERASE
oplot,z(giants,2),(r(giants,0)-z(giants,0)),PSYM=5,SYMSIZE=0.5,color=3
coef = poly_fit(r(dwarfs,2),(r(dwarfs,0)-z(dwarfs,0)),1)
oplot,ttf,poly(ttf,coef),linestyle=1,color=0
coef = poly_fit(r(giants,2),(r(giants,0)-z(giants,0)),1)
oplot,ttf,poly(ttf,coef),linestyle=1,color=3
oplot,z(dwarfs,2),(r(dwarfs,0)-z(dwarfs,0)),position=[0.5,0.1,0.9,0.5],PSYM=5,SYMSIZE=0.5,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,thick=1,charsize=1.,/NOERASE
oplot,z(giants,2),(r(giants,0)-z(giants,0)),PSYM=5,SYMSIZE=0.5,color=3
coef = poly_fit(r(dwarfs,2),(r(dwarfs,0)-z(dwarfs,0)),1)
oplot,ttf,poly(ttf,coef),linestyle=1,color=0
coef = poly_fit(r(giants,2),(r(giants,0)-z(giants,0)),1)
oplot,ttf,poly(ttf,coef),linestyle=1,color=3
;;;;;;;;;;;;;;;;;;CORREGIDO Teff;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
y2=300. & y1=-300.
x1=2500. & x2=10000.
plot,z(dwarfs,0),(r(dwarfs,0)-z(dwarfs,0)),position=[0.1,0.1,0.5,0.5],PSYM=5,SYMSIZE=0.5,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,ytitle=ytit,thick=1,charsize=1.,/NOERASE
oplot,z(giants,0),(r(giants,0)-z(giants,0)),PSYM=5,SYMSIZE=0.5,color=3
coef = poly_fit(r(dwarfs,0),(r(dwarfs,0)-z(dwarfs,0)),1)
oplot,ttt,poly(ttt,coef),linestyle=1,color=0
coef = poly_fit(r(giants,0),(r(giants,0)-z(giants,0)),1)
oplot,ttt,poly(ttt,coef),linestyle=1,color=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
!p.multi=0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; print,z(giants,3)-z(giants,0)
; print,zz(giantsrr,0)-rr(giantsrr,0)
device,/close
loadct,0
;spawn,'gv NGSL_params.eps &'
end
