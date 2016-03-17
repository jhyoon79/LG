PRO vxyz2vrvt,x,y,z,vx,vy,vz,vr=vr,vt=vt

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
Vr = Vradial
Vt = Vtangen
END

;===
FUNCTION pwr_pts,s,xu,xl,Npts

ss = double(1.-s)
b = xu^ss/abs(ss)
c = xl^ss/abs(ss)
a = abs(c-b)
xnew = a*randomu(626,Npts)+min([b,c])
ynew = (xnew*abs(ss))^(1d/ss)

RETURN,ynew
END

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
slope = dblarr(Nbin)
slope_err = dblarr(Nbin)
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


PRO den_close

mass = ['e6e7_30000','e9e10']
opt = 'd'
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,'/media/SEADISK/LG/realistic/part_peri',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
Msat = 10000.
s = (Msat/mr)^(1./3.)
epsilon = s*(G*mr/r_peri)


@plot_setting
device,file='den_close.eps',/color,/enc,/cmyk,xsize=10
loadct,0
!p.charsize=2.
!p.multi=[0,4,2]
clr = 150
incr = 0

for j=0,1 do begin
  list = mass[j]+['',['_20p_']+opt]
  Nlist = N_elements(list)
  
  leg = ['e6e7','e9e10']
  dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[0]+'/'
  chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  chr_rdtbl,dir+'snapshot/snap00000',0,arr0,/silent
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])
  Etot_pts = reform(arr[8,*])
  Jtot_pts = reform(arr[9,*])
  dE = Etot_pts - Etot_pts[0]
  dJ = Jtot_pts - Jtot_pts[0]
  q = dE/epsilon
  dJsJ = dJ/(s*Jtot_pts)

  Etot0 = double(reform(arr0[8,0]))
  Jtot0 = double(reform(arr0[9,0]))
  qinit = (Etot_pts[0] - Etot0)/epsilon
  dJsJinit = (Jtot_pts[0]-Jtot0)/(s*Jtot_pts[0])
print,'q,J-init',qinit,dJsJinit
  
  chr_rdtbl,dir+'part001',0,arr,N_LINES=950,/silent
  arr = double(arr)
  t_orbit = reform(arr[0,*])
  x_orbit = reform(arr[1,*])
  y_orbit = reform(arr[2,*])
  z_orbit = reform(arr[3,*])

  dir = '/media/SEADISK/LG/FinalRun_close/inVL_'+list[1]+'/'
  chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  chr_rdtbl,dir+'snapshot/snap00000',0,arr0,/silent
  arr = double(arr)
  t_pts20p = (arr[0,0])[0]
  x_pts20p = reform(arr[1,*])
  y_pts20p = reform(arr[2,*])
  z_pts20p = reform(arr[3,*])
  Etot_pts20p = reform(arr[8,*])
  Jtot_pts20p = reform(arr[9,*])
  dE20p = Etot_pts20p - Etot_pts20p[0]
  dJ20p = Jtot_pts20p - Jtot_pts20p[0]
  q20p = dE20p/epsilon
  dJsJ20p = dJ20p/(s*Jtot_pts20p)

  Etot0 = double(reform(arr0[8,0]))
  Jtot0 = double(reform(arr0[9,0]))
  qinit = (Etot_pts20p[0] - Etot0)/epsilon
  dJsJinit = (Jtot_pts20p[0]-Jtot0)/(s*Jtot_pts20p[0])
print,'q,J-init 20%',qinit,dJsJinit

  chr_rdtbl,dir+'part001',0,arr,N_LINES=950,/silent
  arr = double(arr)
  t_orbit20p = reform(arr[0,*])
  x_orbit20p = reform(arr[1,*])
  y_orbit20p = reform(arr[2,*])
  z_orbit20p = reform(arr[3,*])

  ;=== convert XYZ into RA & Dec ===
  xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
  xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  xyz2radec,x_pts20p,y_pts20p,z_pts20p,ra_pts20p,dec_pts20p
  xyz2radec,x_orbit20p,y_orbit20p,z_orbit20p,ra_orbit20p,dec_orbit20p

  sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
  ra_orbit_sub = ra_orbit[sub_orbit]
  dec_orbit_sub = dec_orbit[sub_orbit]
  sub_orbit = where(abs(t_orbit20p-t_pts20p) lt 120)	; orbit within 120Myr
  ra_orbit20p_sub = ra_orbit20p[sub_orbit]
  dec_orbit20p_sub = dec_orbit20p[sub_orbit]

  if j eq 0 then begin  
    xmin = 213 &  xmax = 243
    ymin = -23.9 &  ymax = 12
  endif else begin
    xmin = 230 &  xmax = 317
    ymin = -23.9 &  ymax = 12
  endelse
  di = incr*0.4
  x1 = 0.17+di &  x2 = 0.57+di
  y1 = 0.60   &  y2 = 0.95
  x3 = x1     &  x4 = x2
  y3 = 0.12   &  y4 = 0.47
  xt = 'RA [deg]'
  yt = ''
  if incr eq 0 then yt = 'Dec [deg]'
  ytick = replicate(' ',10)
  dy = 10
  if j eq 0 then begin
    xtick = [' ','235',' ','225',' ','215']
    plot,[0],[0],/nodata,xtitle=xt,ytitle=yt,xtickname=xtick,xr=[xmax,xmin],yr=[ymin-dy,ymax],psym=3,pos=[x1,y1,x2,y2]
  endif else begin
    plot,[0],[0],/nodata,xtitle=xt,ytitle=yt,ytickname=ytick,xr=[xmax,xmin],yr=[ymin-dy,ymax],psym=3,pos=[x1,y1,x2,y2]
  endelse
  oplot,ra_orbit_sub,dec_orbit_sub,color=clr
  oplot,ra_pts,dec_pts,psym=3
  oplot,[ra_pts[0]],[dec_pts[0]],psym=1
  legend,leg[incr],box=0,/bottom,charsize=1.2
  if j eq 0 then arrow,227,-17,223,-23,/data
  if j eq 1 then arrow,255,-5,243,-6,/data
  oplot,ra_orbit_sub,dec_orbit_sub-dy,color=clr
  oplot,ra_pts,dec_pts-dy,psym=3
  oplot,[ra_pts[0]],[dec_pts[0]]-dy,psym=1
  if j eq 0 then xyouts,240,-11,'20%',charsize=1.0
  if j eq 0 then xyouts,240,0,'All',charsize=1.0
  if j eq 1 then xyouts,300,-17,'20%',charsize=1.0
  if j eq 1 then xyouts,300,-8,'All',charsize=1.0

  dy = 4
  xt = 'q'
  yt = '' 
  if incr eq 0 then yt = textoidl('\DeltaJ / sJ')
  if incr eq 0 then begin
    plot,q,dJsJ,xr=[5,-5],yr=[-5-dy,5],xtitle=xt,ytitle=yt,psym=3,pos=[x3,y3,x4,y4]
  endif else begin
    plot,q,dJsJ,xr=[5,-5],yr=[-5-dy,5],xtitle=xt,ytitle=yt,ytickname=ytick,psym=3,pos=[x3,y3,x4,y4]
  endelse
  oplot,[0],[0],psym=1
  oplot,q20p,dJsJ20p-dy,psym=3
  oplot,[0],[0]-dy,psym=1
  if j eq 0 then xyouts,4,-5,'20%',charsize=1.0
  if j eq 0 then xyouts,4,1,'All',charsize=1.0
  if j eq 1 then xyouts,4,-5,'20%',charsize=1.0
  if j eq 1 then xyouts,4,1,'All',charsize=1.0
  incr += 1
endfor
ERASE
multiplot,/reset
device,/close

print,'===== BEGIN 2nd plot ====='

device,file='tail_varyingmass.eps',/color,/enc,/cmyk,xsize=15,ysize=25,yoffset=0.5,/portrait
!p.multi=[0,1,3]
!p.charsize=1.5
list = ['nosub','e5e6','e6e7','e7e8','e8e9','e9e10']
Nlist = N_elements(list)
xmin = -15 &  xmax = 15
ymin = -60 &  ymax = 10

; RA & Dec
multiplot
plot,[0],[0],/nodata,ytitle=textoidl('\DeltaDec [deg]'),xr=[xmax,xmin],yr=[ymin,ymax],pos=[0.22,0.71,0.94,0.98]
for j=0,Nlist-2 do begin
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
  xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
  dy = 10*j
  racen = ra_pts[0]
  deccen = dec_pts[0]
  oplot,ra_pts-racen,dec_pts-deccen-dy,psym=3
  oplot,[0],[0]-dy,psym=1,color=clr,thick=5
  xyouts,13,3-dy,list[j],charsize=1.2
endfor

; dv
ymin = -170
ymax = 90
multiplot,/doxaxis
plot,[0],[0],/nodata,xtitle=textoidl('\DeltaRA [deg]'),ytitle=textoidl('\DeltaV_r [km/s]'),xr=[xmax,xmin],yr=[ymin,ymax],pos=[0.22,0.44,0.94,0.71]
for j=0,Nlist-2 do begin
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])
  vx_pts = reform(arr[1,*])
  vy_pts = reform(arr[2,*])
  vz_pts = reform(arr[3,*])

  vxyz2vrvt,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,vr=vr,vt=vt
;=== convert XYZ into RA & Dec ===
  xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
  dy = 25*j
  racen = ra_pts[0]
  oplot,ra_pts-racen,vr-vr[0]-dy,psym=3
  oplot,[0],[0]-dy,psym=1,color=clr,thick=5
  xyouts,13,-45-dy,list[j],charsize=1.2
endfor


; q & dJ
qmin = 6 &  qmax = -7
jmin = -27 &  jmax = 5
multiplot
plot,[0],[0],/nodata,xtitle=textoidl('q'),ytitle=textoidl('\DeltaJ / sJ'),xr=[qmin,qmax],yr=[jmin,jmax],pos=[0.22,0.08,0.94,0.36]
for j=0,Nlist-2 do begin
  chr_rdtbl,'/media/SEADISK/LG/realistic/part_peri',0,arr,/silent
  arr = double(arr)
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  Msat = 10000.
  s = (Msat/mr)^(1./3.)
  epsilon = s*(G*mr/r_peri)
  
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  Etot_pts = reform(arr[8,*])
  Jtot_pts = reform(arr[9,*])
  dE = Etot_pts - Etot_pts[0]
  dJ = Jtot_pts - Jtot_pts[0]
  q = dE/epsilon
  dJsJ = dJ/(s*Jtot_pts)
  dy = 5*j
  oplot,q,dJsJ-dy,psym=3
  oplot,[0],[0]-dy,psym=1,thick=5,color=clr
  xyouts,5,2-dy,list[j],charsize=1.2
endfor
device,/close   

;=======================
; contour plot
;=======================
device,file='tail_varyingmass_cntr.eps',/color,/enc,/cmyk,ysize=25,/portrait
!p.multi=[0,1,2]
!p.charsize=1.5
list = ['nosub','e5e6','e6e7','e7e8','e8e9','e9e10']
Nlist = N_elements(list)
xmin = -15 &  xmax = 15
ymin = -60 &  ymax = 10

plot,[0],[0],/nodata,xtitle=textoidl('\DeltaRA'),ytitle=textoidl('\DeltaDec'),xr=[xmax,xmin],yr=[ymin,ymax]
for j=0,Nlist-2 do begin
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
  xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
  dy = 10*j
  racen = ra_pts[0]
  deccen = dec_pts[0]

  pts_contour,ra_pts,dec_pts,25,25,gv,x,y
  contour,gv,x-racen,y-deccen-dy,/overplot
  oplot,[0],[0]-dy,psym=1,color=clr,thick=5
  xyouts,13,3-dy,list[j],charsize=1.2
endfor

qmin = 6 &  qmax = -7
jmin = -27 &  jmax = 5
plot,[0],[0],/nodata,xtitle=textoidl('q'),ytitle=textoidl('\DeltaJ / sJ'),xr=[qmin,qmax],yr=[jmin,jmax],/noerase
for j=0,Nlist-2 do begin
  chr_rdtbl,'/media/SEADISK/LG/realistic/part_peri',0,arr,/silent
  arr = double(arr)
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  Msat = 10000.
  s = (Msat/mr)^(1./3.)
  epsilon = s*(G*mr/r_peri)
  
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  Etot_pts = reform(arr[8,*])
  Jtot_pts = reform(arr[9,*])
  dE = Etot_pts - Etot_pts[0]
  dJ = Jtot_pts - Jtot_pts[0]
  q = dE/epsilon
  dJsJ = dJ/(s*Jtot_pts)
  dy = 5*j
  oplot,q,dJsJ-dy,psym=3
  oplot,[0],[0]-dy,psym=1,thick=5,color=clr
  xyouts,5,2-dy,list[j],charsize=1.2
endfor




print,'===== BEGIN 3rd plot ====='
device,file='den_profile.eps',/color,/enc,/cmyk,xsize=24,ysize=12,/portrait
loadct,0
!p.multi=[0,5,3]
!p.charsize=1.2

G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
  
list = ['nosub','e5e6','e6e7','e7e8','e8e9','e9e10']
Nlist = N_elements(list)
for j=0,Nlist-2 do begin
  if j eq 0 then begin
    dir = '/media/SEADISK/LG/realistic/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else if j eq 1 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_inner_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08410',0,arr,/silent
  endif else if j eq 2 then begin 
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'_30000/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endelse
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])
  vx_pts = reform(arr[4,*])
  vy_pts = reform(arr[5,*])
  vz_pts = reform(arr[6,*])
  Etot_pts = reform(arr[8,*])
  Jtot_pts = reform(arr[9,*])
  dE = Etot_pts - Etot_pts[0]
  dJ = Jtot_pts - Jtot_pts[0]
  q = dE/epsilon
  dJsJ = dJ/(s*Jtot_pts)

  chr_rdtbl,dir+'part001',0,arr,N_LINES=920,/silent
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

  Npts = 100
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,tail='leading',Npts=Npts,$
        m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err,mean_Vr_pts_l,sigma_Vr_pts_l,phi_l,Vr_l,Vrcen
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,tail='trailing',Npts=Npts,$
        m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,mean_Vr_pts_t,sigma_Vr_pts_t,phi_t,Vr_t,Vrcen

m_phi_pts_t = -1.*m_phi_pts_t

  xmin = -14.9  &  xmax = 24.9
  ymin = 0   &  ymax = 1200
  yt = ' '
  if j eq 0 then yt = textoidl('\Sigma [#/deg^2]')
  multiplot,/rowmajor
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],ytitle=yt,title=list[j],charsize=1
  oplot,m_phi_pts_l,surf_den_pts_l
  oplot,m_phi_pts_t,surf_den_pts_t
  mean_err = sqrt((total(surf_den_pts_l_err^2.)+total(surf_den_pts_t_err^2.))/N_elements([surf_den_pts_l_err,surf_den_pts_t_err]))
  oploterror,[20],[ymax*0.8],[mean_err],psym=3
;  oploterror,m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err
;  oploterror,m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err
;  oploterror,m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err
;  oploterror,[22],[100],sqrt(total([surf_den_pts_t_err^2.,surf_den_pts_l_err^2.])/N_elements([m_phi_pts_l,m_phi_pts_t]))
; legend,list[j],box=0,/bottom,charsize=1.2

  ymin = 0  &  ymax = 13
  if j eq 0 then yt = textoidl('\sigma_{V_r} [km/s]')
  multiplot,/rowmajor
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],ytitle=yt
  oplot,m_phi_pts_l,sigma_Vr_pts_l
  oplot,m_phi_pts_t,sigma_Vr_pts_t

  ymin = -16  &  ymax = 35
  xt = textoidl('\lambda [deg]')
  if j eq 0 then yt = textoidl('\DeltaV_r [km/s]')
  multiplot,/rowmajor
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=xt,ytitle=yt
  oplot,m_phi_pts_l,mean_Vr_pts_l-Vrcen
  oplot,m_phi_pts_t,mean_Vr_pts_t-Vrcen

  Nphi = 360*2
  k = (findgen(Nphi)+1)/Nphi*(360d)
phi_all = phi_t
;phi_all = [phi_t,phi_l]
dphi = abs(phi_all - shift(phi_all,1))
nyq = min(dphi)/2.

  Bk = dblarr(Nphi)
;  Bk_t = dblarr(Nphi)
  max_phi = max(abs(phi_all))
  min_phi = min(abs(phi_all))
print,'min phi',min_phi
  for m=0,Nphi-1 do begin
    fact = k[m]*abs(phi_all)*!dtor
;   fact_t = k[m]*phi_t/max_t;*!dtor
    Bk[m] = abs(1./Nphi*total(complex(cos(fact),sin(fact))))
;    Bk_t[m] = abs(1./Nphi*total(complex(cos(fact_t),sin(fact_t))))
  endfor
  k2 = 360./k

print,nyq

endfor
multiplot,/reset
ERASE
device,/close

END
