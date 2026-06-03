pro MILESparams_magH
  l=1881
  z=fltarr(l,4)
  a1 = intarr(l)
  a2 = intarr(l)
  a3 = fltarr(l)
  d1 = intarr(l)
  d2 = intarr(l)
  d3 = fltarr(l)
  name = strarr(l)

  openr,lun,'lista_GTC_Mh',/get_lun

  for i=0,l-1 do begin
    ; Temporary variables
    f_h = 0.0
    f_t = 0.0
    f_g = 0.0
    f_f = 0.0
    f_a1 = 0
    f_a2 = 0
    f_a3 = 0.0
    f_d1 = 0
    f_d2 = 0
    f_d3 = 0.0
    f_name = ''

    ; Read the line into temporary scalars
;    readf, lun, format='(F7.3, F8.2, F6.3, F6.3, I2, I2, F13.10, I3, I2, F12.9, A13)', $
;          f_h, f_t, f_g, f_f, f_a1, f_a2, f_a3, f_d1, f_d2, f_d3, f_name
    readf, lun, f_h, f_t, f_g, f_f, f_a1, f_a2, f_a3, f_d1, f_d2, f_d3, f_name

    ; Assign to arrays
    z[i,0] = f_t
    z[i,1] = f_g
    z[i,2] = f_f
    z[i,3] = f_h
    a1[i] = f_a1
    a2[i] = f_a2
    a3[i] = f_a3
    d1[i] = f_d1
    d2[i] = f_d2
    d3[i] = f_d3
    name[i] = f_name
  endfor

  free_lun,lun
;;;;
gi=where(z(*,1) lt 3.0)
dw=where(z(*,1) ge 3.0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xsc=[-1.,1.,1.,-1.,-1.] & ysc=[1.,1.,-1.,-1.,1.]
usersym,xsc,ysc,/fill
tvlct,[255,0,255,0],[255,255,0,0],[255,0,0,255],1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;!p.multi=0
set_plot,'ps'
device,filename='MILESGTC.eps',xsize=24,ysize=26,yoffset=1.,/color
erase
;
x1=-2.499
x2=0.999
y1=2500.
y2=40000.
xyo=-2.2
yyo=32500.
xtit='[Fe/H]'
; Log-spaced tick positions and names (in actual scale, not log10)
y_tickv=[2500,3000,3500,4000,4500,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000,25000,30000,40000]
y_tickn=['2.5','3','3.5','4','4.5','5','6','7','8','9','10','12','14','16','18','20','25','30','40']

; In the plot command:
yticks = n_elements(y_tickv)-1
;
ytit='T!deff!n (x 10!e3!n K)'
;
plot,z(dw,2),z(dw,0),POSITION=[0.09,0.1,0.54,0.99],/ylog,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1$
,xtitle=xtit,ytitle=ytit,psym=4,symsize=0.2,charsize=1.6,/NOERASE,$
yticks=n_elements(y_tickv)-1, ytickv=y_tickv, ytickn=y_tickn
xyouts,xyo,yyo,'Dwarfs',charsize=1.8
;
OPENW, lun, 'GTC_lista_selected', /GET_LUN
;;;metal poor dwarfs;;;
dp=where(z(*,1) ge 3.0 and z(*,0) ge 4850. and z(*,0) le 6600. and z(*,2) le -1.7 and z(*,3) ge 8.5)
;   PRINTF, lun, format='(F7.3, F8.2, F6.3, F6.3, I2, I2, F13.10, I3, I2, F12.9, A)', $
;   z[dp,3],z[dp,0],z[dp,1],z[dp,2],a1[dp],a2[dp],a3[dp],d1[dp],d2[dp],d3[dp],name[dp]
oplot,z(dp,2),z(dp,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dp)-1 DO BEGIN
     idx = dp[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;intermediate metallicity dwarfs;;;
dii=where(z(*,1) ge 3.0 and z(*,0) ge 4850. and z(*,0) le 6600. and z(*,2) gt -1.7 and  z(*,2) le -1.2 and z(*,3) ge 8.3)
oplot,z(dii,2),z(dii,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dii)-1 DO BEGIN
     idx = dii[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;subsolar dwarfs;;;
dis=where(z(*,1) ge 3.0 and z(*,0) ge 4850. and z(*,0) le 6600. and z(*,2) gt -1.2 and z(*,2) lt -0.3 and z(*,3) ge 7.0)
oplot,z(dis,2),z(dis,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dis)-1 DO BEGIN
     idx = dis[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;metal-rich dwarfs;;;
dir=where(z(*,1) ge 3.0 and z(*,0) ge 4850. and z(*,0) le 6600. and z(*,2) ge -0.3 and z(*,2) le 0.5 and z(*,3) ge 6.35)
oplot,z(dir,2),z(dir,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dir)-1 DO BEGIN
     idx = dir[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;M dwarfs;;;
dl=where(z(*,1) ge 3.0 and z(*,2) lt 0.6 and z(*,0) lt 3600. and z(*,3) ge 6.2)
oplot,z(dl,2),z(dl,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dl)-1 DO BEGIN
     idx = dl[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;Cool dwarfs;;;
dcl=where(z(*,1) ge 3.0 and z(*,2) lt 0.6 and z(*,0) lt 4850. and z(*,0) ge 3600. and z(*,3) ge 5.3)
oplot,z(dcl,2),z(dcl,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dcl)-1 DO BEGIN
     idx = dcl[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;hot dwarfs;;;
dh=where(z(*,1) ge 3.0 and z(*,2) lt 0.6 and z(*,0) gt 6600. and z(*,0) le 20000. and z(*,3) ge 5.8)
oplot,z(dh,2),z(dh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dh)-1 DO BEGIN
     idx = dh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;super hot dwarfs;;;
dsh=where(z(*,1) ge 3.0 and z(*,2) lt 0.6 and z(*,0) gt 20000. and z(*,3) ge 8.0)
oplot,z(dsh,2),z(dsh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(dsh)-1 DO BEGIN
     idx = dsh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot,z(gi,2),z(gi,0),POSITION=[0.54,0.1,0.99,0.99],/ylog,xrange=[x1,x2],xstyle=1,yrange=[y1,y2],ystyle=1$
,xtitle=xtit,symsize=0.2,psym=4,charsize=1.6,yticks=n_elements(y_tickv)-1, ytickv=y_tickv, ytickn=y_tickn,ycharsize=0.0001,/NOERASE
xyouts,xyo,yyo,'Giants',charsize=1.8
;;;metal-poor giants;;;
gph=where(z(*,1) lt 3.0 and z(*,2) le -1.5 and z(*,0) lt 9999 and z(*,3) ge 6.1)
oplot,z(gph,2),z(gph,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gph)-1 DO BEGIN
     idx = gph[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;poor giants;;;
gsh=where(z(*,1) lt 3.0 and z(*,2) gt -1.5 and z(*,2) le -0.5 and z(*,0) ge 4200. and z(*,3) ge 5.5) 
oplot,z(gsh,2),z(gsh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gsh)-1 DO BEGIN
     idx = gsh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;sub-solar giants;;;
grh=where(z(*,1) lt 3.0 and z(*,2) gt -0.5 and z(*,2) le -0.1 and z(*,0) ge 4200. and z(*,3) ge 4.95) 
oplot,z(grh,2),z(grh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(grh)-1 DO BEGIN
     idx = grh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;hot metal-rich giants;;;
ghmr=where(z(*,1) lt 3.0 and z(*,2) gt -0.1 and z(*,2) le 0.6 and z(*,0) lt 9999. and z(*,0) ge 5200. and z(*,3) ge 5.0) 
oplot,z(ghmr,2),z(ghmr,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(ghmr)-1 DO BEGIN
     idx = ghmr[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;metal-rich giants;;;
grrh=where(z(*,1) lt 3.0 and z(*,2) gt -0.1 and z(*,2) le 0.6 and z(*,0) lt 5200. and z(*,0) ge 4200. and z(*,3) ge 4.0) 
oplot,z(grrh,2),z(grrh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(grrh)-1 DO BEGIN
     idx = grrh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;cool giants;;;
gch=where(z(*,1) lt 3.0 and z(*,0) lt 4200. and z(*,0) ge 3800. and z(*,2) gt -1.2 and z(*,3) ge 2.5) 
oplot,z(gch,2),z(gch,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gch)-1 DO BEGIN
     idx = gch[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;cool cool giants;;;
gcch=where(z(*,1) lt 3.0 and z(*,0) lt 4200. and z(*,0) ge 3800. and z(*,2) gt -1.2 and z(*,3) ge 2.5) 
oplot,z(gcch,2),z(gcch,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gcch)-1 DO BEGIN
     idx = gcch[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;supergiants;;;
gsgh=where(z(*,1) lt 3.0 and z(*,0) lt 3800. and z(*,0) ge 3100. and z(*,2) gt -1.2 and z(*,3) ge 1.0) 
oplot,z(gsgh,2),z(gsgh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gsgh)-1 DO BEGIN
     idx = gsgh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
;;;cool supergiants;;;
gcsgh=where(z(*,1) lt 3.0 and z(*,0) lt 3100. and z(*,3) ge -1.5) 
oplot,z(gcsgh,2),z(gcsgh,0),color=2,psym=4,symsize=0.55
  FOR i = 0, N_ELEMENTS(gcsgh)-1 DO BEGIN
     idx = gcsgh[i]
     PRINTF, lun, format='(F10.3, F10.2, F8.3, F8.3, I4, I4, F15.10, I5, I4, F15.9, A15)', $
        z[idx,3], z[idx,0], z[idx,1], z[idx,2], a1[idx], a2[idx], a3[idx], d1[idx], d2[idx], d3[idx], name[idx]
  ENDFOR
;
FREE_LUN, lun
!p.multi=0
device,/close
set_plot,'x'
END

;;;;;;
;OPENW, lun, 'GTC_lista_selected', /GET_LUN
;FOR i = 0, N_ELEMENTS(v_sel)-1 DO BEGIN
;   PRINTF, lun, v_sel[i]
;ENDFOR
;FREE_LUN, lun
;hh=[8.5,8.0,7.5,7.0,6.5,6.0,5.5,5.0,4.5,4.0,3.5,3.0]
;s85=8.5 & s80=8.0 & s75=7.5 & s70=7.0 & s65=6.5 & s60=6.0 & s55=5.5 & S50=5.0 & s45=4.5 & s40=4.0 & s35=3.5 & s30=3.0
;d85=where(z(*,1) ge 3.0 and z(*,3) ge s85)
;d8085=where(z(*,1) ge 3.0 and z(*,3) lt s85 and z(*,3) ge s80 and z(*,2) lt 0.65)
;d7580=where(z(*,1) ge 3.0 and z(*,3) lt s80 and z(*,3) ge s75 and z(*,2) lt 0.65)
;d7075=where(z(*,1) ge 3.0 and z(*,3) lt s75 and z(*,3) ge s70 and z(*,2) lt 0.65)
;d6570=where(z(*,1) ge 3.0 and z(*,3) lt s70 and z(*,3) ge s65 and z(*,2) lt 0.65)
;d6065=where(z(*,1) ge 3.0 and z(*,3) lt s65 and z(*,3) ge s60 and z(*,2) lt 0.65)
;d5560=where(z(*,1) ge 3.0 and z(*,3) lt s60 and z(*,3) ge s55 and z(*,2) lt 0.65)
;d5055=where(z(*,1) ge 3.0 and z(*,3) lt s55 and z(*,3) ge s50 and z(*,2) lt 0.65)
;d4550=where(z(*,1) ge 3.0 and z(*,3) lt s50 and z(*,3) ge s45 and z(*,2) lt 0.65)
;oplot,z(d85,2),z(d85,0),color=3,psym=4,symsize=0.55
;oplot,z(d8085,2),z(d8085,0),color=3,psym=4,symsize=0.55
;oplot,z(d7580,2),z(d7580,0),color=3,psym=4,symsize=0.55
;oplot,z(d7075,2),z(d7075,0),color=3,psym=4,symsize=0.55
;oplot,z(d6570,2),z(d6570,0),color=3,psym=4,symsize=0.55
;oplot,z(d6065,2),z(d6065,0),color=3,psym=4,symsize=0.55
;oplot,z(d5560,2),z(d5560,0),color=3,psym=4,symsize=0.55
;oplot,z(d5055,2),z(d5055,0),color=3,psym=4,symsize=0.55
;oplot,z(d4550,2),z(d4550,0),color=3,psym=4,symsize=0.55
