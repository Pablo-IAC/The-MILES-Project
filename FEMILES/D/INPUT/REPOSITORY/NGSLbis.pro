pro NGSLbis
l=378
l=317
l=303
lp=66
lp=5
z=fltarr(l,3)
openr,lun,'PARAM_NGSL.idl',/get_lun
for i=0,l-1 do begin
 readf,lun,t,g,f
 z(i,0)=t
 z(i,1)=g
 z(i,2)=f
endfor
free_lun,lun
;
giants=where(z(*,1) lt 3.0)
dwarfs=where(z(*,1) ge 3.0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xsc=[-1.,1.,1.,-1.,-1.] & ysc=[1.,1.,-1.,-1.,1.]
usersym,xsc,ysc,/fill
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;!p.multi=0
set_plot,'ps'
device,filename='NGSLbis.eps',xsize=20,ysize=24,yoffset=1.,/color
erase
;
x1=-2.2
x2=0.9999
y1=3000.
y2=40000.
xyo=-2.
yyo=29000.
xtit='[Fe/H]'
;ytit='Teff'
tval=strarr(17)
tval=[3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40]
tvals=intarr(17)
tvals=[3000,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,$
18000,20000,25000,30000,35000,40000]
;
plot,z(dwarfs,2),z(dwarfs,0),POSITION=[0.088,0.07,0.543,0.99],/ylog,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1$
,xtitle=xtit,psym=4,symsize=0.65,charsize=1.5,yticks=16,ytickv=tvals,ytickn=tval,/NOERASE
xyouts,xyo,yyo,'Dwarfs',charsize=1.8
xyouts,-2.6,8530.,'T!deff!n (x 10!e3!n K)',orientation=90,charsize=1.8
;
plot,z(giants,2),z(giants,0),POSITION=[0.543,0.07,0.998,0.99],/ylog,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1$
,xtitle=xtit,symsize=0.65,psym=4,charsize=1.5,yticks=16,ytickv=tvals,ycharsize=0.0001,/NOERASE
xyouts,xyo,yyo,'Giants',charsize=1.8
;
!p.multi=0
device,/close
set_plot,'x'
end
