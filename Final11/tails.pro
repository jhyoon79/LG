PRO DENnVr,x,y,z,vx,vy,vz,Etot,ra,dec,tail=tail,Npts=Npts,$
          m_phi,surf_den,surf_den_err,mean_Vr,sigma_Vr,phi_init,Vr,Vrcen

;=== estimate density along the orbit ===
IF (tail EQ 'leading') THEN BEGIN
  tailpart = where(ra-ra[0] lt 0)
ENDIF ELSE IF (tail EQ 'trailing') THEN BEGIN
  tailpart = where(ra-ra[0] gt 0)
ENDIF

Npts_inbin = Npts
Ntailpart = N_ELEMENTS(tailpart)
Nbin = (Ntailpart/Npts_inbin eq Ntailpart/double(Npts_inbin)) ? Ntailpart/Npts_inbin:Ntailpart/Npts_inbin+1

ra_tail = ra[tailpart]
dec_tail = dec[tailpart]
phi_init = angsep(ra[0],dec[0],ra_tail,dec_tail)

;=== estimate sigma_v along the orbit ===
x_sun = -8.d
vx_sun = 10.00	; from Dehnen&Binney 1998MNRAS.298.387
vy_sun = 5.25+220
vz_sun = 7.17
vx = vx-vx_sun
vy = vy-vy_sun
vz = vz-vz_sun
v = sqrt(vx^2.+vy^2.+vz^2.)
cos_theta = ((x-x_sun)*vx+y*vy+z*vz) / (sqrt((x-x_sun)^2.+y^2.+z^2.)*abs(v))
Vradial = v*cos_theta
Vtangen = v*sqrt(1-cos_theta^2)
Vrcen = Vradial[0]
Vr = Vradial[tailpart]
Vr_tail = Vr

m_phi = dblarr(Nbin)
surf_den = dblarr(Nbin)
surf_den_err = dblarr(Nbin)
mean_Vr = dblarr(Nbin)
sigma_Vr = dblarr(Nbin)
;plot,ra,dec,psym=3
for i=0,Nbin-1 do begin
  sort_phi = sort(phi_init)
  init = i*Npts+1
  final = (i lt Nbin-1)? (i+1)*Npts:Ntailpart-1
  Num = final-init+1
  if Num ge 2 then begin
    m_phi[i] = mean(phi_init[sort_phi[init:final]])
    mean_Vr[i] = mean(Vr_tail[sort_phi[init:final]])
    sigma_Vr[i] = stddev(Vr_tail[sort_phi[init:final]])
  endif else begin
    m_phi[i] = phi_init[sort_phi[init:final]]
    mean_Vr[i] = Vr_tail[sort_phi[init:final]]
    sigma_Vr[i] = 0;stddev(Vr_tail[sort_phi[init:final]]
  endelse
  area = (phi_init[sort_phi[final]]-phi_init[sort_phi[init]])^2.
  surf_den[i] = Num/area
  surf_den_err[i] = sqrt(Num)/area
  if area eq 0 then surf_den[i] = -999
;print,m_phi[i],surf_den[i],Num
;oplot,ra_tail[sort_phi],dec_tail[sort_phi],psym=1,color=255-i*15,symsize=.5
;wait,0.5
endfor
if (where(surf_den lt 0) ne -1) then begin
  m_phi = m_phi[0:Nbin-2]
  mean_Vr = mean_Vr[0:Nbin-2]
  sigma_Vr = sigma_Vr[0:Nbin-2]
  surf_den = surf_den[0:Nbin-2]
  surf_den_err = surf_den_err[0:Nbin-2]
endif
END

pro tails

dir = '/media/SEADISK/LG/Final11/'
spawn,'ls '+dir+'snap* > tmp'
chr_rdtbl,'tmp',0,fname
spawn,'rm -f tmp'
Nfile = N_elements(fname)

@plot_setting
loadct,0
device,file='tails_all.eps',/cmyk,/enc,xsize=24,ysize=18
!p.charsize=1.7
!p.multi=[0,4,3]
clr = 150
for i=0,11 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
  arr = double(arr)
    arr = double(arr)
    t_pts = (arr[0,0])[0]
    x_pts = reform(arr[1,*])
    y_pts = reform(arr[2,*])
    z_pts = reform(arr[3,*])
    vx_pts = reform(arr[4,*])
    vy_pts = reform(arr[5,*])
    vz_pts = reform(arr[6,*])
    Etot_pts = reform(arr[8,*])
    Nline = 232
    if tname eq '08770' or tname eq '09110' then Nline=290
    chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
    arr = double(arr)
    t_orbit = reform(arr[0,*])
    x_orbit = reform(arr[1,*])
    y_orbit = reform(arr[2,*])
    z_orbit = reform(arr[3,*])
  ;=== convert XYZ into RA & Dec ===
    xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
    xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  
    sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  
  xt = ''  &  yt = ''
  if (i/4 eq i/4.) and (i lt 8) then begin
    xt = '' 
    yt = textoidl('\Delta y [kpc]')
  endif else if i eq 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = textoidl('\Delta y [kpc]')
  endif else if i gt 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = ''
  endif
  xcen = x_pts[0]
  ycen = y_pts[0]

  if i eq 0 then begin
    xmin = -9  &  xmax = 4
    ymin = -9  &  ymax = 4
  endif else if i eq 1 then begin
    xmin = -7  &  xmax = 3
    ymin = -6  &  ymax = 4
  endif else if i eq 2 then begin
    xmin = -5  &  xmax = 5
    ymin = -4  &  ymax = 6
  endif else if i eq 3 then begin
    xmin = -6  &  xmax = 5
    ymin = -2  &  ymax = 9
  endif else if i eq 4 then begin
    xmin = -2  &  xmax = 18
    ymin = -10  &  ymax = 10
  endif else if i eq 5 then begin
    xmin = -5  &  xmax = 4
    ymin = -5  &  ymax = 4
  endif else if i eq 6 then begin
    xmin = -19 &  xmax = 2
    ymin = -15 &  ymax = 6
  endif else if i eq 7 then begin
    xmin = -10 &  xmax = 4
    ymin = -9  &  ymax = 5
  endif else if i eq 8 then begin
    xmin = -24 &  xmax = 1
    ymin = -18 &  ymax = 7
  endif else if i eq 9 then begin
    xmin = -6  &  xmax = 1
    ymin = -6  &  ymax = 1
  endif else if i eq 10 then begin
    xmin = -6  &  xmax = 3
    ymin = -6  &  ymax = 3
  endif else if i eq 11 then begin
    xmin = -6.5 &  xmax = 1.5
    ymin = -4. &  ymax = 4
  endif
  dx = (i mod 4)*0.23
  dy = i/4*0.31
  x1 = 0.08+dx   &  y1 = 0.70-dy
  x2 = 0.28+dx   &  y2 = 0.97-dy
print,i,dx,dy,x1,y1,x2,y2

;  plot,x_pts-xcen,y_pts-ycen,psym=3,xtitle=xt,ytitle=yt,/isotropic,pos=[x1,y1,x2,y2],title=t
  plot,x_pts-xcen,y_pts-ycen,psym=3,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=xt,ytitle=yt,/isotropic,pos=[x1,y1,x2,y2]
  oplot,[0],[0],psym=1,color=150,symsize=1.5,thick=5
  legend,t,box=0,charsize=.9,/bottom,/right
endfor
multiplot,/reset
erase
device,/close

;=====================
; pts contour map
;=====================
device,file='tails_cntr.eps',/cmyk,/enc,xsize=24,ysize=18
!p.charsize=1.7
!p.multi=[0,4,3]
clr = 150
for i=0,11 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
  arr = double(arr)
    arr = double(arr)
    t_pts = (arr[0,0])[0]
    x_pts = reform(arr[1,*])
    y_pts = reform(arr[2,*])
    z_pts = reform(arr[3,*])
    vx_pts = reform(arr[4,*])
    vy_pts = reform(arr[5,*])
    vz_pts = reform(arr[6,*])
    Etot_pts = reform(arr[8,*])
    Nline = 232
    if tname eq '08770' or tname eq '09110' then Nline=290
    chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
    arr = double(arr)
    t_orbit = reform(arr[0,*])
    x_orbit = reform(arr[1,*])
    y_orbit = reform(arr[2,*])
    z_orbit = reform(arr[3,*])
  ;=== convert XYZ into RA & Dec ===
    xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
    xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  
    sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  
  xt = ''  &  yt = ''
  if (i/4 eq i/4.) and (i lt 8) then begin
    xt = '' 
    yt = textoidl('\Delta y [kpc]')
  endif else if i eq 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = textoidl('\Delta y [kpc]')
  endif else if i gt 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = ''
  endif
  xcen = x_pts[0]
  ycen = y_pts[0]

  if i eq 0 then begin
    xmin = -9  &  xmax = 4
    ymin = -9  &  ymax = 4
  endif else if i eq 1 then begin
    xmin = -7  &  xmax = 3
    ymin = -6  &  ymax = 4
  endif else if i eq 2 then begin
    xmin = -5  &  xmax = 5
    ymin = -4  &  ymax = 6
  endif else if i eq 3 then begin
    xmin = -6  &  xmax = 5
    ymin = -2  &  ymax = 9
  endif else if i eq 4 then begin
    xmin = -2  &  xmax = 18
    ymin = -10  &  ymax = 10
  endif else if i eq 5 then begin
    xmin = -5  &  xmax = 4
    ymin = -5  &  ymax = 4
  endif else if i eq 6 then begin
    xmin = -19 &  xmax = 2
    ymin = -15 &  ymax = 6
  endif else if i eq 7 then begin
    xmin = -10 &  xmax = 4
    ymin = -9  &  ymax = 5
  endif else if i eq 8 then begin
    xmin = -24 &  xmax = 1
    ymin = -18 &  ymax = 7
  endif else if i eq 9 then begin
    xmin = -6  &  xmax = 1
    ymin = -6  &  ymax = 1
  endif else if i eq 10 then begin
    xmin = -6  &  xmax = 3
    ymin = -6  &  ymax = 3
  endif else if i eq 11 then begin
    xmin = -6.5 &  xmax = 1.5
    ymin = -4. &  ymax = 4
  endif
  dx = (i mod 4)*0.23
  dy = i/4*0.31
  x1 = 0.08+dx   &  y1 = 0.70-dy
  x2 = 0.28+dx   &  y2 = 0.97-dy
  xr = max(x_pts)-min(x_pts)
  yr = max(y_pts)-min(y_pts)
  xgn = round(xr)*10;15
  ygn = round(yr)*10;15
  lev = [2,5,10,20,30,40,50,60,70,80,90]
  pts_contour,x_pts-xcen,y_pts-ycen,xgn,ygn,gv,x,y
  contour,gv,x,y,levels=lev,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=xt,ytitle=yt,/isotropic,pos=[x1,y1,x2,y2],title=strtrim(xgn,2)+' '+strtrim(ygn,2),thick=1.
  oplot,[0],[0],psym=1,color=150,symsize=1.5,thick=5
  legend,t,box=0,charsize=.9,/bottom,/right
endfor
multiplot,/reset
erase

;=====================
; radec contour map
;=====================
device,file='tails_cntr_radec.eps',/cmyk,/enc,xsize=25,ysize=14
!p.charsize=1.2
!p.multi=[0,3,2]
ind = [0,2,5,10]
for ii=1,3 do begin
  i = ind[ii]
  if ii eq 0 then chr_rdtbl,'/media/SEADISK/LG/realistic/snapshot/snap08440',0,arr,/silent
  if ii gt 0 then chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
  arr = double(arr)
    t_pts = (arr[0,0])[0]
    x_pts = reform(arr[1,*])
    y_pts = reform(arr[2,*])
    z_pts = reform(arr[3,*])
    vx_pts = reform(arr[4,*])
    vy_pts = reform(arr[5,*])
    vz_pts = reform(arr[6,*])
    Etot_pts = reform(arr[8,*])
    Nline = 232
    if tname eq '08770' or tname eq '09110' then Nline=290
    chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
    arr = double(arr)
    t_orbit = reform(arr[0,*])
    x_orbit = reform(arr[1,*])
    y_orbit = reform(arr[2,*])
    z_orbit = reform(arr[3,*])
  ;=== convert XYZ into RA & Dec ===
    xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
    xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  
    sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  
  xmin = min(ra_pts)-1
  xmax = max(ra_pts)+1
  ymin = min(dec_pts)-1
  ymax = max(dec_pts)+1
  xt = 'RA [deg]'  &  yt = ''
  if ii eq 1 then yt = 'Dec [deg]'
  xr = xmax-xmin
  yr = ymax-ymin

  multiplot,xgap=0.02,/doyaxis,/rowmajor
plot,ra_pts,dec_pts,psym=3,xr=[xmax,xmin],yr=[ymin,ymax],ytit=yt
  if ii eq 2 then arrow,235,0,232.5,0,/data
;  if ii eq 3 then arrow,232,0.25,230,0.25,/data
  if ii eq 3 then arrow,244,2,244,3.5,/data
  legend,t,box=0,/bottom,charsize=1.2

  ss = 24/60.	; 24arcmin
  xgn = round(xr*1/ss)
  ygn = round(yr*1/ss)
  pts_contour,ra_pts,dec_pts,xgn,ygn,gv,x,y
  dy = (ii-1)*0.3
  x1 = 0.15
  x2 = 0.7
  y1 = 0.7 - dy
  y2 = 0.95 - dy
  mgv = max(gv)
  sig = 0.12*(ss*60)^2./5. *(3000.*10./25.)/1650.
  lev = [1.5*sig,2*sig,3*sig,5*sig,7*sig,9*sig,11*sig,13*sig,15*sig,17*sig]
  multiplot,xgap=0.02,/doyaxis,/rowmajor
  contour,gv,x,y,levels=lev,xtit=xt,ytit=yt,thick=1.,xr=[xmax,xmin],yr=[ymin,ymax]
  oplot,[ra_pts[0]],dec_pts[[0]],psym=1,color=150,symsize=2,thick=7
  if ii eq 2 then arrow,235,0,232.5,0,/data
;  if ii eq 3 then arrow,232,0.25,230,0.25,/data
  if ii eq 3 then arrow,244,2,244,3.5,/data
endfor
multiplot,/reset

;========================
; E&J plot
;========================
s = 0.0069943216
epsilon = 120.29439
device,file='tails_all_qj.eps',/cmyk,/enc,xsize=24,ysize=18
!p.multi=[0,4,3]
clr = 150
for i=0,11 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
    arr = double(arr)
    t_pts = (arr[0,0])[0]
    x_pts = reform(arr[1,*])
    y_pts = reform(arr[2,*])
    z_pts = reform(arr[3,*])
    vx_pts = reform(arr[4,*])
    vy_pts = reform(arr[5,*])
    vz_pts = reform(arr[6,*])
Etotpts = reform(arr[8,*])
Jpts = reform(arr[9,*])
delE = Etotpts-Etotpts[0]
qpts = delE/epsilon
dJpts = Jpts-Jpts[0]
dJsJ = dJpts / (s * Jpts)
    Nline = 232
    if tname eq '08770' or tname eq '09110' then Nline=290
    chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
    arr = double(arr)
    t_orbit = reform(arr[0,*])
    x_orbit = reform(arr[1,*])
    y_orbit = reform(arr[2,*])
    z_orbit = reform(arr[3,*])
  ;=== convert XYZ into RA & Dec ===
    xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
    xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  
    sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  
  xt = ''  &  yt = ''
  if (i/4 eq i/4.) and (i lt 8) then begin
    xt = '' 
    yt = textoidl('\Delta y [kpc]')
  endif else if i eq 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = textoidl('\Delta y [kpc]')
  endif else if i gt 8 then begin
    xt = textoidl('\Delta x [kpc]')
    yt = ''
  endif
  xcen = x_pts[0]
  ycen = y_pts[0]

  dx = (i mod 4)*0.23
  dy = i/4*0.31
  x1 = 0.08+dx   &  y1 = 0.70-dy
  x2 = 0.28+dx   &  y2 = 0.97-dy
print,i,dx,dy,x1,y1,x2,y2

  plot,qpts,dJsJ,psym=3,xr=[5,-5],yr=[-5,5]
  oplot,[0],[0],psym=1,color=150
  legend,t,box=0,charsize=.9,/bottom,/right
endfor
multiplot,/reset
erase
device,/close

;========================
; bifurcation
;========================
loadct,13
device,file='tails_bifur.eps',/cmyk,/enc,xsize=18,ysize=18,/color
!p.multi=[0,2,2]
!p.charsize=1.2
for i=6,11,5 do begin
  if i eq 6 then begin 
    chr_rdtbl,dir+'00000_08770',0,arr
    arr = double(arr)
    Etot_init = reform(arr[8,*])
    delE_init = Etot_init-Etot_init[0]
    q_init = delE_init/epsilon
    qmin = min(q_init);*0.9
    qmax = max(q_init);*0.9
    qcolor = round((q_init-qmin)/(qmax-qmin)*245.)+10
  endif else begin
    chr_rdtbl,dir+'00000_08770',0,arr
;    chr_rdtbl,dir+'00000_10453',0,arr
    arr = double(arr)
    Etot_init = reform(arr[8,*])
    delE_init = Etot_init-Etot_init[0]
    q_init = delE_init/epsilon
    qmin = min(q_init);*0.9
    qmax = max(q_init);*0.9
    qcolor = round((q_init-qmin)/(qmax-qmin)*245.)+10
  endelse
  Nc = N_elements(qcolor)
  chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])
  vx_pts = reform(arr[4,*])
  vy_pts = reform(arr[5,*])
  vz_pts = reform(arr[6,*])
  Etotpts = reform(arr[8,*])
  Jpts = reform(arr[9,*])
  delE = Etotpts-Etotpts[0]
  qpts = delE/epsilon
  dJpts = Jpts-Jpts[0]
  dJsJ = dJpts / (s * Jpts)
  Nline = 232
  if tname eq '08770' or tname eq '09110' then Nline=290
  chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
  arr = double(arr)
  t_orbit = reform(arr[0,*])
  x_orbit = reform(arr[1,*])
  y_orbit = reform(arr[2,*])
  z_orbit = reform(arr[3,*])
;=== convert XYZ into RA & Dec ===
  xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
  xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit

  sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
  ra_orbit_sub = ra_orbit[sub_orbit]
  dec_orbit_sub = dec_orbit[sub_orbit]
  
  xt = textoidl('\Delta x [kpc]')
  yt = textoidl('\Delta y [kpc]')
  if i gt 6 then yt = ' '
  xcen = x_pts[0]
  ycen = y_pts[0]
  if i eq 6 then begin
    xmin = -19  &  xmax = 2
    ymin = -15  &  ymax = 6
    qmin = -12  &  qmax = 5
    jmin = -14  &  jmax = 3
    x1 = 0  & x2 = -2
    y1 = -2 & y2 = -4.4
  endif else begin
    xmin = -6.7 &  xmax = 1.3
    ymin = -4   &  ymax = 4
    qmin = -6   &  qmax = 4
    jmin = -5   &  jmax = 5
    x1 = 0  & x2 = -1
    y1 = -0.5  & y2 = -1.0
  endelse
  xg = 0.03
;  if i gt 7 then xg = 0.1
  yg = 0.05
;  if i gt 6 then xg = 0.08
  multiplot,/rowmajor,xgap=xg,ygap=yg,/doxaxis,/doyaxis
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=xt,ytitle=yt;,/isotropic
  for k=0,Nc-1 do oplot,[x_pts[k]-xcen],[y_pts[k]-ycen],psym=3,color=qcolor[k]
  oplot,[0],[0],psym=1,symsize=1.5,thick=5
  arrow,x1,y1,x2,y2,/data
  legend,t,box=0,charsize=.9,/right,/bottom

  xt = 'q'
  yt = textoidl('\DeltaJ / sJ')
  if i gt 6 then yt = ' '
  multiplot,/rowmajor,xgap=xg,ygap=yg,/doxaxis,/doyaxis
  plot,[0],[0],/nodata,xr=[qmax,qmin],yr=[jmin,jmax],xtitle=xt,ytitle=yt;,/isotropic
  vline,0,linestyle=2
  hline,0,linestyle=2
  for k=0,Nc-1 do oplot,[qpts[k]],[dJsJ[k]],psym=3,color=qcolor[k]
endfor
multiplot,/reset
erase
device,/close

END
