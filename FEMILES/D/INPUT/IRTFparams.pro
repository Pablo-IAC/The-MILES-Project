pro IRTFparams
l=380
z=fltarr(l,3)
openr,lun,'PARAM_IRTF.idl',/get_lun
for i=0,l-1 do begin
 readf,lun,t,g,f
 z(i,0)=t
 z(i,1)=g
 z(i,2)=f
endfor
free_lun,lun
;
l=205
s=fltarr(l,3)
openr,lun,'PARAM_IRTF.idl',/get_lun
for i=0,l-1 do begin
 readf,lun,t,g,f
 s(i,0)=t
 s(i,1)=g
 s(i,2)=f
endfor
free_lun,lun
;
giants=where(z(*,1) lt 3.0)
dwarfs=where(z(*,1) ge 3.0)
giantss=where(s(*,1) lt 3.0)
dwarfss=where(s(*,1) ge 3.0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xsc=[-1.,1.,1.,-1.,-1.] & ysc=[1.,1.,-1.,-1.,1.]
usersym,xsc,ysc,/fill
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;!p.multi=0
set_plot,'ps'
device,filename='IRTFparams.eps',xsize=24,ysize=26,yoffset=1.,/color
erase
;
x1=-2.499
x2=0.9999
y1=2000.
y2=15000.
xyo=-2.2
yyo=13800.
xtit='[Fe/H]'
xtit='[Fe/H]'
ytit='T!deff!n (x 10!e3!n K)'
y_tickn=['2.0','2.5','3','3.5','4','4.5','5','6','7','8','9','10','12','14','16']
y_tickv=[2000,2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000,12000,14000,16000]
;
yyyticks=13
;
plot,z(dwarfs,2),z(dwarfs,0),POSITION=[0.088,0.07,0.543,0.99],/ylog,xrange=[x1,x2],xstyle=1,$
yrange=[y1,y2],ystyle=1,xtitle=xtit,ytitle=ytit,psym=4,symsize=0.65,charsize=1.6,/NOERASE,$
yticks=n_elements(y_tickv)-1, ytickv=y_tickv, ytickn=y_tickn
;
oplot,s(dwarfss,2),s(dwarfss,0),color=2,psym=4,symsize=0.65
xyouts,xyo,yyo,'Dwarfs',charsize=1.8
;
plot,z(giants,2),z(giants,0),POSITION=[0.543,0.07,0.998,0.99],/ylog,xrange=[x1,x2],xstyle=1,$
yrange=[y1,y2],ystyle=1,xtitle=xtit,symsize=0.65,psym=4,charsize=1.6,ycharsize=0.000001,/NOERASE,$
yticks=n_elements(y_tickv)-1,ytickv=y_tickv, ytickn=y_tickn
;
oplot,s(giantss,2),z(giantss,0),color=2,psym=4,symsize=0.65
xyouts,xyo,yyo,'Giants',charsize=1.8
;
!p.multi=0
device,/close
set_plot,'x'
end
