;===
PRO convert_xyz_radec,x,y,z,ra=ra,dec=dec

l = atan(y/(x+8.))*!radeg 
b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
lcosb = l*cos(b*!dtor)
dec = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
ra = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(dec*!dtor) )*!radeg + 282.25

END


PRO DENnVr,x,y,z,vx,vy,vz,Etot,ra,dec,ra_proj,dec_proj,tail=tail,Npts=Npts,$
          m_phi,surf_den,surf_den_err,sigma_Vr,phi,Vr,slope,slope_err

;=== estimate density along the orbit ===
IF (tail EQ 'leading') THEN BEGIN
  tailpart = where(Etot-Etot[0] lt 0)
ENDIF ELSE IF (tail EQ 'trailing') THEN BEGIN
  tailpart = where(Etot-Etot[0] gt 0)
ENDIF

Npts_inbin = Npts
Ntailpart = N_ELEMENTS(tailpart)
Nbin = (Ntailpart/Npts_inbin eq Ntailpart/double(Npts_inbin)) ? Ntailpart/Npts_inbin:Ntailpart/Npts_inbin+1
phi = angsep(ra[0],dec[0],ra_proj[tailpart],dec_proj[tailpart])
sort_proj = sort(phi)
al = ra_proj[tailpart[sort_proj]]
dl = dec_proj[tailpart[sort_proj]]
ra_part = ra[tailpart]
dec_part = dec[tailpart]

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
Vr = Vradial[tailpart]

m_phi = dblarr(Nbin)
surf_den = dblarr(Nbin)
surf_den_err = dblarr(Nbin)
sigma_Vr = dblarr(Nbin)
slope = dblarr(Nbin)
slope_err = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntailpart-1) then final = Ntailpart-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin = final-initial+1
  if final-initial ge 1 then begin 
    m_phi[i] = mean(phi[sort_proj[initial:final]])
    sigma_Vr[i] = stddev(Vr[sort_proj[initial:final]])
  endif else begin
    m_phi[i] = (phi[sort_proj[initial:final]])
    sigma_Vr[i] = (Vr[sort_proj[initial:final]])
  endelse
  area = angsep(ai,di,af,df)^2.
  surf_den[i] = N_bin/area
  surf_den_err[i] = sqrt(N_bin)/area

  fit = poly_fit(ra_part[initial:final],dec_part[initial:final],1,sigma=sigma)
  slope[i] = fit[1]
  slope_err[i] = sigma[1] 
endfor

END


PRO den_profile_Mclass

mass = ['e5e6','e6e7_30000','e7e8','e8e9','e9e10']
opt = 'd'

@plot_setting
DEVICE,file='den_profile_Mclass_'+opt+'.ps',/color,ysize=24,yoffset=0.5

for j=0,N_elements(mass)-1 do begin

list = mass[j]+['',['_20p_','_15p_','_10p_','_5p_']+opt]
Nlist = N_elements(list)

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

;=== N-body ===

!p.multi=[0,4,Nlist]
for k=0,Nlist-1 do begin
  if k eq 0 then begin
    dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[k]+'/'
    chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
  endif else begin
    dir = '/media/SEADISK/LG/FinalRun_close/inVL_'+list[k]+'/'
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

  chr_rdtbl,dir+'part001',0,arr,N_LINES=986;,/silent
  arr = double(arr)
  t_orbit = reform(arr[0,*])
  x_orbit = reform(arr[1,*])
  y_orbit = reform(arr[2,*])
  z_orbit = reform(arr[3,*])
;=== convert XYZ into RA & Dec ===
  convert_xyz_radec,x_pts,y_pts,z_pts,ra=ra_pts,dec=dec_pts
  convert_xyz_radec,x_orbit,y_orbit,z_orbit,ra=ra_orbit,dec=dec_orbit

  sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
  ra_orbit_sub = ra_orbit[sub_orbit]
  dec_orbit_sub = dec_orbit[sub_orbit]

  xmin = 215 &  xmax = 241
  ymin = -25 &  ymax = 10
  plot,[0],[0],/nodata,xr=[xmax,xmin],yr=[ymin,ymax],psym=3,/isotropic,title=list[k]
  oplot,ra_orbit_sub,dec_orbit_sub,color=100
  oplot,ra_pts,dec_pts,psym=3
  min_d_curve,ra_orbit_sub,dec_orbit_sub,ra_pts,dec_pts,x_proj=ra_proj,y_proj=dec_proj

  plot,q,dJsJ,xr=[5,-5],yr=[-5,5],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),psym=3

  Npts = 150 
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='leading',Npts=Npts,$
        m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err,sigma_Vr_pts_l,phi_l,Vr_l,slope_pts_l,slope_pts_l_err
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='trailing',Npts=Npts,$
        m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,sigma_Vr_pts_t,phi_t,Vr_t,slope_pts_t,slope_pts_t_err

  xmin = 0  &  xmax = 26;max([m_phi_pts_l,m_phi_pts_t])
  ymin = .1  &  ymax = 1000
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
  oploterror,m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err
  oploterror,m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,linestyle=0,color=255,errcolor=255

  IF (k EQ 0) THEN  legend,['leading','trailing'],linestyle=[0,2],color=[0,255],box=0,charsize=0.8,/bottom
  
  ymin = 0  &  ymax = 11
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
  oplot,m_phi_pts_l,sigma_Vr_pts_l
  oplot,m_phi_pts_t,sigma_Vr_pts_t,linestyle=0,color=255
endfor
multiplot,/reset
erase
endfor

stop

;=== with 30,000 subhalo ===
!p.multi=[0,3,5]
list = ['all','e5e6','e5e6_300000','e6e7','e6e7_30000']
for k=0,4 do begin
  dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[k]+'/'
  if k eq 0 then dir = '/media/SEADISK/LG/inVL_12790/'
  if k eq 1 or k eq 2 then begin	; erase when 300000 is done
    chr_rdtbl,dir+'snapshot/snap01000',0,arr,/silent
  endif else begin
    chr_rdtbl,dir+'snapshot/snap12790',0,arr,/silent
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

  chr_rdtbl,dir+'part001',0,arr,N_lines=1045,/silent
  arr = double(arr)
  t_orbit = reform(arr[0,*])
  x_orbit = reform(arr[1,*])
  y_orbit = reform(arr[2,*])
  z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
  convert_xyz_radec,x_pts,y_pts,z_pts,ra=ra_pts,dec=dec_pts
  convert_xyz_radec,x_orbit,y_orbit,z_orbit,ra=ra_orbit,dec=dec_orbit

  sub_orbit = where(abs(t_orbit-t_pts) lt 120,count)	;orbit within 120Myr
  if (count gt 0) then begin
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  endif else begin
    ra_orbit_sub = ra_orbit;[sub_orbit]		; erase when 300000 is done
    dec_orbit_sub = dec_orbit;[sub_orbit]
  endelse

  plot,[0],[0],/nodata,xr=[max(ra_pts),min(ra_pts)],yr=[min(dec_pts),max(dec_pts)],psym=3,/isotropic,title=list[k]
  oplot,ra_orbit_sub,dec_orbit_sub,color=100
  oplot,ra_pts,dec_pts,psym=3
  min_d_curve,ra_orbit_sub,dec_orbit_sub,ra_pts,dec_pts,x_proj=ra_proj,y_proj=dec_proj

  Npts = 100 
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='leading',Npts=Npts,$
        m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err,sigma_Vr_pts_l,phi_l,Vr_l,slope_pts_l,slope_pts_l_err
  
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='trailing',Npts=Npts,$
        m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,sigma_Vr_pts_t,phi_t,Vr_t,slope_pts_t,slope_pts_t_err

  xmin = 0  &  xmax = 40;max([m_phi_pts_l,m_phi_pts_t])
  case k of
    0: begin
         ymin = 1  &  ymax = 100000
       end
    1: begin
         ymin = 1  &  ymax = 51000
       end
    2: begin
         ymin = 1  &  ymax = 51000
       end
    3: begin
         ymin = 1  &  ymax = 5100
       end
    4: begin
         ymin = 1  &  ymax = 5100
       end
  endcase
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
  oploterror,m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err
  oploterror,m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,linestyle=0,color=255,errcolor=255
  IF (k EQ 0) THEN  legend,['leading','trailing'],linestyle=[0,2],color=[0,255],box=0,charsize=0.8,/bottom
  
  case k of
    0: begin
         ymin = 0  &  ymax = 36
       end
    1: begin
         ymin = 0  &  ymax = 36
       end
    2: begin
         ymin = 0  &  ymax = 36
       end
    3: begin
         ymin = 0  &  ymax = 36
       end
    4: begin
         ymin = 0  &  ymax = 36
       end
  endcase
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
  oplot,m_phi_pts_l,sigma_Vr_pts_l
  oplot,m_phi_pts_t,sigma_Vr_pts_t,linestyle=0,color=255
endfor
multiplot,/reset
ERASE
device,/close

END
