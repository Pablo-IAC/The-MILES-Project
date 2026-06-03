pro NGSL_partransf
km=135
k311=311
ko=176
z=fltarr(km,6)
r=fltarr(km,3)
zz=fltarr(ko,3)
rr=fltarr(ko,3)
ttf=[-2.2,-2.0,-1.75,-1.5,-1.25,-1.,-.75,-.5,-.25,.0,.25,.7]
ttt=[3000.,4000.,5000.,6000.,7000.,8000.]
xtit='Teff' & ytit=''
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
set_plot,'ps'
DEVICE,DECOMPOSED=0
;loadct,2
erase
device,filename='NGSL_partransf.eps',xsize=17,ysize=17,xoffset=2.,yoffset=2.,/color
;!p.multi=[1,2,2]
!p.multi=[1,1,1]
;
xsc=[-1.,1.,1.,-1.,-1.] & ysc=[1.,1.,-1.,-1.,1.]
usersym,xsc,ysc,/fill
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
x1=-2.2 & x2=0.7
;xtit='[Fe/H]' & ytit='Teff(NGSL) - Teff(MILES)'
;xtit='[Fe/H]' & ytit=(' + cgGreek('delta', /CAPITAL) + 'M)'
xtit='[Fe/H]' & ytit='!7D!3T!deff!n (KV12 - MILES)'
;
xyo=fltarr(1) & xyo(0)=-0.2
yyo=fltarr(1) & yyo(0)=430.
;
plot,z(dwarfs,2),(z(dwarfs,3)-z(dwarfs,0)),PSYM=5,SYMSIZE=1.2,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,thick=1,charsize=1.5,/NOERASE
oplot,z(giants,2),(z(giants,3)-z(giants,0)),PSYM=2,SYMSIZE=1.2,color=0
;
oplot,xyo,yyo+10.,PSYM=2,SYMSIZE=1.2,color=0
xyouts,xyo+0.17,yyo,'Giants',charsize=1.9,color=0
oplot,xyo,yyo-50.,PSYM=5,SYMSIZE=1.2
xyouts,xyo+0.17,yyo-60.,'Dwarfs',charsize=1.9
xx1=[-0.3,0.7]
yy1=[350,350]
oplot,xx1,yy1,linestyle=0
xx2=[-0.3,-0.3]
yy2=[350,500]
oplot,xx2,yy2,linestyle=0
;
coefd = poly_fit(z(dwarfs,2),(z(dwarfs,3)-z(dwarfs,0)),1)
oplot,ttf,poly(ttf,coefd),linestyle=0,color=0,thick=3
coefg = poly_fit(z(giants,2),(z(giants,3)-z(giants,0)),1)
oplot,ttf,poly(ttf,coefg),linestyle=2,color=0,thick=3
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
!p.multi=0
device,/close
loadct,0
end
