pro FUVparams_oneplot
l=476
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
l=173
s=fltarr(l,3)
openr,lun,'PARAM_NGSL.idl',/get_lun
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
device,filename='FUVparams_oneplot.ps',xsize=24,ysize=26,yoffset=1.,/color
erase
;
x1=-2.499
x2=0.9999
y1=2500.
y2=89999.
;xyo=-2.
;yyo=29000.
xyo=-2.2
yyo=72000.
xtit='[Fe/H]'
ytit='T!deff!n (x 10!e3!n K)'
y_tickn=['2.5','3','3.5','4','4.5','5','6','7','8','9','10','12','14','16','18',$
'20','25','30','40','50','60','70','80']
y_tickv=[2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,$
20000,25000,30000,40000,50000,60000,70000,80000]
;
plot,z(dwarfs,2),z(dwarfs,0),POSITION=[0.099,0.099,0.99,0.99],/ylog,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1$
,xtitle=xtit,ytitle=ytit,psym=4,symsize=0.9,charsize=1.6,/NOERASE,$
yticks=n_elements(y_tickv)-1, ytickv=y_tickv, ytickn=y_tickn
;
xyouts,xyo,yyo,'Dwarfs',charsize=1.8
oplot,z(giants,2),z(giants,0),color=3,psym=4,symsize=0.9
xyouts,xyo,yyo-10000,'Giants',charsize=2,color=3
;
!p.multi=0
device,/close
set_plot,'x'
end
