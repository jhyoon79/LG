pro den_profile


mass_name = 'e5e6'
dir_Mbin = '/media/SEADISK/LG/FinalRun/inVL_'+mass_name
dir_close = '/media/SEADISK/LG/FinalRun_close/inVL_'+mass_name

dir_5p_part = dir_close+'_5p_d/snapshot/'
dir_10p_part = dir_close+'_10p_d/snapshot/'
dir_15p_part = dir_close+'_15p_d/snapshot/'
dir_20p_part = dir_close+'_20p_d/snapshot/'
dir_5p_orbit = dir_close+'_5p_d/'
dir_10p_orbit = dir_close+'_10p_d/'
dir_15p_orbit = dir_close+'_15p_d/'
dir_20p_orbit = dir_close+'_20p_d/'

dir_nosub = '/scratch/jhyoon/Research/LG/realistic/'
dir_nbody = '/media/SEADISK/LG/Nbody/pal5_final/'

x_sun = -8.d

openr,1,'../frogin_initial_Pal5'
readf,1,tmp
readf,1,t_final
close,1
t_final = t_final+3
fname_final = 'snap'+string(t_final,f='(i5.5)')
;=== no subhalo ===
chr_rdtbl,dir_nosub+'snapshot/'+fname_final,0,arr,/silent
arr = double(arr)
t_nosub = (arr[0,0])[0]
x_nosub = reform(arr[1,*])
y_nosub = reform(arr[2,*])
z_nosub = reform(arr[3,*])
vx_nosub = reform(arr[4,*])
vy_nosub = reform(arr[5,*])
vz_nosub = reform(arr[6,*])
Etot_nosub = reform(arr[8,*])
v_nosub = sqrt(Vx_nosub^2.+Vy_nosub^2.+Vz_nosub^2.)

chr_rdtbl,dir_nosub+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
t_orbit_nosub = t_orbit
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_nosub = atan(y_nosub/(x_nosub+8.))*!radeg 
b_nosub = atan(z_nosub/sqrt((x_nosub+8.)^2.+y_nosub^2.))*!radeg
lcosb_nosub = l_nosub*cos(b_nosub*!dtor)
delta_nosub = asin( cos(b_nosub*!dtor)*sin((l_nosub-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nosub*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_nosub = asin( (cos(b_nosub*!dtor)*sin((l_nosub-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nosub*!dtor)*sin(62.6*!dtor))/cos(delta_nosub*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_nosub) lt 120)	; orbit within 500Myr
alpha_orbit_ns = alpha_orbit
delta_orbit_ns = delta_orbit
alpha_orbit_nosub = alpha_orbit[sub_orbit]
delta_orbit_nosub = delta_orbit[sub_orbit]


set_plot,'ps'
@plot_setting
!p.charsize=1.5
device,file='den_profile.ps',/color,/landscape,xsize=25
;plot,alpha_nosub,delta_nosub,psym=3,xr=[200,260],yr=[-83,17],/isotropic

plot,alpha_nosub,delta_nosub,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
min_d_curve,alpha_orbit_nosub,delta_orbit_nosub,alpha_nosub,delta_nosub,x_proj=alpha_proj,y_proj=delta_proj
for i=0,N_elements(alpha_nosub)-1 do oplot,[alpha_proj[i],alpha_nosub[i]],[delta_proj[i],delta_nosub[i]],color=200

plot,alpha_nosub,delta_nosub,psym=3,xr=[219.5,223.5],yr=[-15,-8],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
min_d_curve,alpha_orbit_nosub,delta_orbit_nosub,alpha_nosub,delta_nosub,x_proj=alpha_proj,y_proj=delta_proj
for i=0,N_elements(alpha_nosub)-1 do oplot,[alpha_proj[i],alpha_nosub[i]],[delta_proj[i],delta_nosub[i]],color=200




;=== estimate density along the orbit ===
leading_nosub = where(Etot_nosub-Etot_nosub[0] lt 0)
trailing_nosub = where(Etot_nosub-Etot_nosub[0] gt 0)
help,leading_nosub,trailing_nosub
Npts_inbin = 150

Nlead = N_elements(leading_nosub)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_nosub_l = angsep(alpha_nosub[0],delta_nosub[0],alpha_proj[leading_nosub],delta_proj[leading_nosub])
sort_proj = sort(phi_nosub_l)
al = alpha_proj[leading_nosub[sort_proj]]
dl = delta_proj[leading_nosub[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_nosub_l = alpha_nosub[leading_nosub]
delta_nosub_l = delta_nosub[leading_nosub]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_nosub-x_sun)*vx_nosub+y_nosub*vy_nosub+z_nosub*vz_nosub) / (sqrt((x_nosub-x_sun)^2.+y_nosub^2.+z_nosub^2.)*abs(v_nosub))
Vradial = v_nosub*cos_theta
Vtangen = v_nosub*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_nosub]

m_phi_nosub_l = dblarr(Nbin)
surf_den_nosub_l = dblarr(Nbin)
surf_den_nosub_err_l = dblarr(Nbin)
sigma_Vr_nosub_l = dblarr(Nbin)
slope_nosub_l = dblarr(Nbin)
slope_err_nosub_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_nosub_l[i] = mean(phi_nosub_l[sort_proj[initial:final]])
  sigma_Vr_nosub_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_nosub_l[i] = N_bin_l/area
  surf_den_nosub_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_nosub_l[initial:final],delta_nosub_l[initial:final],1,sigma=sigma)
  slope_nosub_l[i] = fit[1]
  slope_err_nosub_l[i] = sigma[1] 
endfor


Ntrail = N_elements(trailing_nosub)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_nosub_t = angsep(alpha_nosub[0],delta_nosub[0],alpha_proj[trailing_nosub],delta_proj[trailing_nosub])
sort_proj = sort(phi_nosub_t)
at = alpha_proj[trailing_nosub[sort_proj]]
dt = delta_proj[trailing_nosub[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_nosub_t = alpha_nosub[trailing_nosub]
delta_nosub_t = delta_nosub[trailing_nosub]

Vr_t = Vradial[trailing_nosub]
m_phi_nosub_t = dblarr(Nbin)
surf_den_nosub_t = dblarr(Nbin)
surf_den_nosub_err_t = dblarr(Nbin)
sigma_Vr_nosub_t = dblarr(Nbin)
slope_nosub_t = dblarr(Nbin)
slope_err_nosub_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_nosub_t[i] = mean(phi_nosub_t[sort_proj[initial:final]])
  sigma_Vr_nosub_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_nosub_t[i] = N_bin_t/area
  surf_den_nosub_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_nosub_t[initial:final],delta_nosub_t[initial:final],1,sigma=sigma)
  slope_nosub_t[i] = fit[1]
  slope_err_nosub_t[i] = sigma[1] 
endfor

;===========
;=== 5p ===
;===========
chr_rdtbl,dir_5p_part+fname_final,0,arr,/silent
arr = double(arr)
t_5p = (arr[0,0])[0]
x_5p = reform(arr[1,*])
y_5p = reform(arr[2,*])
z_5p = reform(arr[3,*])
vx_5p = reform(arr[4,*])
vy_5p = reform(arr[5,*])
vz_5p = reform(arr[6,*])
Etot_5p = reform(arr[8,*])
v_5p = sqrt(Vx_5p^2.+Vy_5p^2.+Vz_5p^2.)

chr_rdtbl,dir_5p_orbit+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_5p = atan(y_5p/(x_5p+8.))*!radeg 
b_5p = atan(z_5p/sqrt((x_5p+8.)^2.+y_5p^2.))*!radeg
lcosb_5p = l_5p*cos(b_5p*!dtor)
delta_5p = asin( cos(b_5p*!dtor)*sin((l_5p-33.)*!dtor)*sin(62.6*!dtor)+sin(b_5p*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_5p = asin( (cos(b_5p*!dtor)*sin((l_5p-33.)*!dtor)*cos(62.6*!dtor)-sin(b_5p*!dtor)*sin(62.6*!dtor))/cos(delta_5p*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_5p) lt 120)	; orbit within 500Myr
alpha_orbit_5p = alpha_orbit[sub_orbit]
delta_orbit_5p = delta_orbit[sub_orbit]

plot,alpha_5p,delta_5p,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_5p,delta_orbit_5p,color=100
min_d_curve,alpha_orbit_5p,delta_orbit_5p,alpha_5p,delta_5p,x_proj=alpha_proj,y_proj=delta_proj


;=== estimate density along the orbit ===
leading_5p = where(alpha_5p-alpha_5p[0] lt 0)
trailing_5p = where(alpha_5p-alpha_5p[0] gt 0)
;leading_5p = where(Etot_5p-Etot_5p[0] lt 0)
;trailing_5p = where(Etot_5p-Etot_5p[0] gt 0)
help,leading_5p,trailing_5p
Npts_inbin = 150

Nlead = N_elements(leading_5p)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_5p_l = angsep(alpha_5p[0],delta_5p[0],alpha_proj[leading_5p],delta_proj[leading_5p])
sort_proj = sort(phi_5p_l)
al = alpha_proj[leading_5p[sort_proj]]
dl = delta_proj[leading_5p[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_5p_l = alpha_5p[leading_5p]
delta_5p_l = delta_5p[leading_5p]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_5p-x_sun)*vx_5p+y_5p*vy_5p+z_5p*vz_5p) / (sqrt((x_5p-x_sun)^2.+y_5p^2.+z_5p^2.)*abs(v_5p))
Vradial = v_5p*cos_theta
Vtangen = v_5p*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_5p]

m_phi_5p_l = dblarr(Nbin)
surf_den_5p_l = dblarr(Nbin)
surf_den_5p_err_l = dblarr(Nbin)
sigma_Vr_5p_l = dblarr(Nbin)
slope_5p_l = dblarr(Nbin)
slope_err_5p_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_5p_l[i] = mean(phi_5p_l[sort_proj[initial:final]])
  sigma_Vr_5p_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_5p_l[i] = N_bin_l/area
  surf_den_5p_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_5p_l[initial:final],delta_5p_l[initial:final],1,sigma=sigma)
  slope_5p_l[i] = fit[1]
  slope_err_5p_l[i] = sigma[1] 
endfor

Ntrail = N_elements(trailing_5p)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_5p_t = angsep(alpha_5p[0],delta_5p[0],alpha_proj[trailing_5p],delta_proj[trailing_5p])
sort_proj = sort(phi_5p_t)
at = alpha_proj[trailing_5p[sort_proj]]
dt = delta_proj[trailing_5p[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_5p_t = alpha_5p[trailing_5p]
delta_5p_t = delta_5p[trailing_5p]

Vr_t = Vradial[trailing_5p]
m_phi_5p_t = dblarr(Nbin)
surf_den_5p_t = dblarr(Nbin)
surf_den_5p_err_t = dblarr(Nbin)
sigma_Vr_5p_t = dblarr(Nbin)
slope_5p_t = dblarr(Nbin)
slope_err_5p_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_5p_t[i] = mean(phi_5p_t[sort_proj[initial:final]])
  sigma_Vr_5p_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_5p_t[i] = N_bin_t/area
  surf_den_5p_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_5p_t[initial:final],delta_5p_t[initial:final],1,sigma=sigma)
  slope_5p_t[i] = fit[1]
  slope_err_5p_t[i] = sigma[1] 
endfor


;===========
;=== 10p ===
;===========
chr_rdtbl,dir_10p_part+fname_final,0,arr,/silent
arr = double(arr)
t_10p = (arr[0,0])[0]
x_10p = reform(arr[1,*])
y_10p = reform(arr[2,*])
z_10p = reform(arr[3,*])
vx_10p = reform(arr[4,*])
vy_10p = reform(arr[5,*])
vz_10p = reform(arr[6,*])
Etot_10p = reform(arr[8,*])
v_10p = sqrt(Vx_10p^2.+Vy_10p^2.+Vz_10p^2.)

chr_rdtbl,dir_10p_orbit+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_10p = atan(y_10p/(x_10p+8.))*!radeg 
b_10p = atan(z_10p/sqrt((x_10p+8.)^2.+y_10p^2.))*!radeg
lcosb_10p = l_10p*cos(b_10p*!dtor)
delta_10p = asin( cos(b_10p*!dtor)*sin((l_10p-33.)*!dtor)*sin(62.6*!dtor)+sin(b_10p*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_10p = asin( (cos(b_10p*!dtor)*sin((l_10p-33.)*!dtor)*cos(62.6*!dtor)-sin(b_10p*!dtor)*sin(62.6*!dtor))/cos(delta_10p*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_10p) lt 120)	; orbit within 500Myr
alpha_orbit_10p = alpha_orbit[sub_orbit]
delta_orbit_10p = delta_orbit[sub_orbit]

plot,alpha_10p,delta_10p,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_10p,delta_orbit_10p,color=100
min_d_curve,alpha_orbit_10p,delta_orbit_10p,alpha_10p,delta_10p,x_proj=alpha_proj,y_proj=delta_proj


;=== estimate density along the orbit ===
leading_10p = where(alpha_10p-alpha_10p[0] lt 0)
trailing_10p = where(alpha_10p-alpha_10p[0] gt 0)
;leading_10p = where(Etot_10p-Etot_10p[0] lt 0)
;trailing_10p = where(Etot_10p-Etot_10p[0] gt 0)
help,leading_10p,trailing_10p
Npts_inbin = 150

Nlead = N_elements(leading_10p)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_10p_l = angsep(alpha_10p[0],delta_10p[0],alpha_proj[leading_10p],delta_proj[leading_10p])
sort_proj = sort(phi_10p_l)
al = alpha_proj[leading_10p[sort_proj]]
dl = delta_proj[leading_10p[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_10p_l = alpha_10p[leading_10p]
delta_10p_l = delta_10p[leading_10p]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_10p-x_sun)*vx_10p+y_10p*vy_10p+z_10p*vz_10p) / (sqrt((x_10p-x_sun)^2.+y_10p^2.+z_10p^2.)*abs(v_10p))
Vradial = v_10p*cos_theta
Vtangen = v_10p*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_10p]

m_phi_10p_l = dblarr(Nbin)
surf_den_10p_l = dblarr(Nbin)
surf_den_10p_err_l = dblarr(Nbin)
sigma_Vr_10p_l = dblarr(Nbin)
slope_10p_l = dblarr(Nbin)
slope_err_10p_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_10p_l[i] = mean(phi_10p_l[sort_proj[initial:final]])
  sigma_Vr_10p_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_10p_l[i] = N_bin_l/area
  surf_den_10p_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_10p_l[initial:final],delta_10p_l[initial:final],1,sigma=sigma)
  slope_10p_l[i] = fit[1]
  slope_err_10p_l[i] = sigma[1] 
endfor

Ntrail = N_elements(trailing_10p)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_10p_t = angsep(alpha_10p[0],delta_10p[0],alpha_proj[trailing_10p],delta_proj[trailing_10p])
sort_proj = sort(phi_10p_t)
at = alpha_proj[trailing_10p[sort_proj]]
dt = delta_proj[trailing_10p[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_10p_t = alpha_10p[trailing_10p]
delta_10p_t = delta_10p[trailing_10p]

Vr_t = Vradial[trailing_10p]
m_phi_10p_t = dblarr(Nbin)
surf_den_10p_t = dblarr(Nbin)
surf_den_10p_err_t = dblarr(Nbin)
sigma_Vr_10p_t = dblarr(Nbin)
slope_10p_t = dblarr(Nbin)
slope_err_10p_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_10p_t[i] = mean(phi_10p_t[sort_proj[initial:final]])
  sigma_Vr_10p_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_10p_t[i] = N_bin_t/area
  surf_den_10p_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_10p_t[initial:final],delta_10p_t[initial:final],1,sigma=sigma)
  slope_10p_t[i] = fit[1]
  slope_err_10p_t[i] = sigma[1] 
endfor

;=== ALL subhalo ===
chr_rdtbl,dir_BI5_part+fname_final,0,arr,/silent
arr = double(arr)
t_allsub = (arr[0,0])[0]
x_allsub = reform(arr[1,*])
y_allsub = reform(arr[2,*])
z_allsub = reform(arr[3,*])
vx_allsub = reform(arr[4,*])
vy_allsub = reform(arr[5,*])
vz_allsub = reform(arr[6,*])
Etot_allsub = reform(arr[8,*])
v_allsub = sqrt(Vx_allsub^2.+Vy_allsub^2.+Vz_allsub^2.)

chr_rdtbl,dir_BI5_orbit+'part001_e7_1000sub',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_allsub = atan(y_allsub/(x_allsub+8.))*!radeg 
b_allsub = atan(z_allsub/sqrt((x_allsub+8.)^2.+y_allsub^2.))*!radeg
lcosb_allsub = l_allsub*cos(b_allsub*!dtor)
delta_allsub = asin( cos(b_allsub*!dtor)*sin((l_allsub-33.)*!dtor)*sin(62.6*!dtor)+sin(b_allsub*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_allsub = asin( (cos(b_allsub*!dtor)*sin((l_allsub-33.)*!dtor)*cos(62.6*!dtor)-sin(b_allsub*!dtor)*sin(62.6*!dtor))/cos(delta_allsub*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_allsub) lt 120)	; orbit within 500Myr
alpha_orbit_allsub = alpha_orbit[sub_orbit]
delta_orbit_allsub = delta_orbit[sub_orbit]

plot,alpha_allsub,delta_allsub,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_allsub,delta_orbit_allsub,color=100
min_d_curve,alpha_orbit_allsub,delta_orbit_allsub,alpha_allsub,delta_allsub,x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
leading_allsub = where(alpha_allsub-alpha_allsub[0] lt 0)
trailing_allsub = where(alpha_allsub-alpha_allsub[0] gt 0)
;leading_allsub = where(Etot_allsub-Etot_allsub[0] lt 0)
;trailing_allsub = where(Etot_allsub-Etot_allsub[0] gt 0)
help,leading_allsub,trailing_allsub
Npts_inbin = 150

Nlead = N_elements(leading_allsub)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_allsub_l = angsep(alpha_allsub[0],delta_allsub[0],alpha_proj[leading_allsub],delta_proj[leading_allsub])
sort_proj = sort(phi_allsub_l)
al = alpha_proj[leading_allsub[sort_proj]]
dl = delta_proj[leading_allsub[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_allsub_l = alpha_allsub[leading_allsub]
delta_allsub_l = delta_allsub[leading_allsub]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_allsub-x_sun)*vx_allsub+y_allsub*vy_allsub+z_allsub*vz_allsub) / (sqrt((x_allsub-x_sun)^2.+y_allsub^2.+z_allsub^2.)*abs(v_allsub))
Vradial = v_allsub*cos_theta
Vtangen = v_allsub*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_allsub]

m_phi_allsub_l = dblarr(Nbin)
surf_den_allsub_l = dblarr(Nbin)
surf_den_allsub_err_l = dblarr(Nbin)
sigma_Vr_allsub_l = dblarr(Nbin)
slope_allsub_l = dblarr(Nbin)
slope_err_allsub_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_allsub_l[i] = mean(phi_allsub_l[sort_proj[initial:final]])
  sigma_Vr_allsub_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_allsub_l[i] = N_bin_l/area
  surf_den_allsub_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_allsub_l[initial:final],delta_allsub_l[initial:final],1,sigma=sigma)
  slope_allsub_l[i] = fit[1]
  slope_err_allsub_l[i] = sigma[1] 
endfor

Ntrail = N_elements(trailing_allsub)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_allsub_t = angsep(alpha_allsub[0],delta_allsub[0],alpha_proj[trailing_allsub],delta_proj[trailing_allsub])
sort_proj = sort(phi_allsub_t)
at = alpha_proj[trailing_allsub[sort_proj]]
dt = delta_proj[trailing_allsub[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_allsub_t = alpha_allsub[trailing_allsub]
delta_allsub_t = delta_allsub[trailing_allsub]

Vr_t = Vradial[trailing_allsub]
m_phi_allsub_t = dblarr(Nbin)
surf_den_allsub_t = dblarr(Nbin)
surf_den_allsub_err_t = dblarr(Nbin)
sigma_Vr_allsub_t = dblarr(Nbin)
slope_allsub_t = dblarr(Nbin)
slope_err_allsub_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_allsub_t[i] = mean(phi_allsub_t[sort_proj[initial:final]])
  sigma_Vr_allsub_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_allsub_t[i] = N_bin_t/area
  surf_den_allsub_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_allsub_t[initial:final],delta_allsub_t[initial:final],1,sigma=sigma)
  slope_allsub_t[i] = fit[1]
  slope_err_allsub_t[i] = sigma[1] 
endfor

;=== in VL 13000 ===
chr_rdtbl,dir_Mbin+'snap13000',0,arr,/silent
arr = double(arr)
t_13000 = (arr[0,0])[0]
x_13000 = reform(arr[1,*])
y_13000 = reform(arr[2,*])
z_13000 = reform(arr[3,*])
vx_13000 = reform(arr[4,*])
vy_13000 = reform(arr[5,*])
vz_13000 = reform(arr[6,*])
Etot_13000 = reform(arr[8,*])
v_13000 = sqrt(Vx_13000^2.+Vy_13000^2.+Vz_13000^2.)

chr_rdtbl,dir_Mbin+'part001_13000',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_13000 = atan(y_13000/(x_13000+8.))*!radeg 
b_13000 = atan(z_13000/sqrt((x_13000+8.)^2.+y_13000^2.))*!radeg
lcosb_13000 = l_13000*cos(b_13000*!dtor)
delta_13000 = asin( cos(b_13000*!dtor)*sin((l_13000-33.)*!dtor)*sin(62.6*!dtor)+sin(b_13000*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_13000 = asin( (cos(b_13000*!dtor)*sin((l_13000-33.)*!dtor)*cos(62.6*!dtor)-sin(b_13000*!dtor)*sin(62.6*!dtor))/cos(delta_13000*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_13000) lt 180)	; orbit within 500Myr
alpha_orbit_13000 = alpha_orbit[sub_orbit]
delta_orbit_13000 = delta_orbit[sub_orbit]

print,'start!!!'
;plot,alpha_13000,delta_13000,psym=3,xr=[200,280],yr=[-16,10],/isotropic
plot,alpha_13000,delta_13000,psym=3,xr=[237.0,237.1],yr=[5.5,6.5]
oplot,alpha_orbit_13000,delta_orbit_13000,color=100
oplot,alpha_orbit_13000,delta_orbit_13000,color=100,psym=8
min_d_curve,alpha_orbit_13000,delta_orbit_13000,alpha_13000,delta_13000,x_proj=alpha_proj,y_proj=delta_proj
for i=0,N_elements(alpha_13000)-1 do oplot,[alpha_proj[i],alpha_13000[i]],[delta_proj[i],delta_13000[i]],color=200

;=== estimate density along the orbit ===
leading_13000 = where(alpha_13000-alpha_13000[0] lt 0)
trailing_13000 = where(alpha_13000-alpha_13000[0] gt 0)
help,leading_13000,trailing_13000
Npts_inbin = 150

Nlead = N_elements(leading_13000)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_13000_l = angsep(alpha_13000[0],delta_13000[0],alpha_proj[leading_13000],delta_proj[leading_13000])
sort_proj = sort(phi_13000_l)
al = alpha_proj[leading_13000[sort_proj]]
dl = delta_proj[leading_13000[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_13000_l = alpha_13000[leading_13000]
delta_13000_l = delta_13000[leading_13000]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_13000-x_sun)*vx_13000+y_13000*vy_13000+z_13000*vz_13000) / (sqrt((x_13000-x_sun)^2.+y_13000^2.+z_13000^2.)*abs(v_13000))
Vradial = v_13000*cos_theta
Vtangen = v_13000*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_13000]

m_phi_13000_l = dblarr(Nbin)
surf_den_13000_l = dblarr(Nbin)
surf_den_13000_err_l = dblarr(Nbin)
sigma_Vr_13000_l = dblarr(Nbin)
slope_13000_l = dblarr(Nbin)
slope_err_13000_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_13000_l[i] = mean(phi_13000_l[sort_proj[initial:final]])
  sigma_Vr_13000_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_13000_l[i] = N_bin_l/area
  surf_den_13000_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_13000_l[initial:final],delta_13000_l[initial:final],1,sigma=sigma)
  slope_13000_l[i] = fit[1]
  slope_err_13000_l[i] = sigma[1] 
print,'f13000l',fit[1],sigma[1]
endfor

Ntrail = N_elements(trailing_13000)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
Nbin = Nbin-1
phi_13000_t = angsep(alpha_13000[0],delta_13000[0],alpha_proj[trailing_13000],delta_proj[trailing_13000])
sort_proj = sort(phi_13000_t)
at = alpha_proj[trailing_13000[sort_proj]]
dt = delta_proj[trailing_13000[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_13000_t = alpha_13000[trailing_13000]
delta_13000_t = delta_13000[trailing_13000]

Vr_t = Vradial[trailing_13000]
m_phi_13000_t = dblarr(Nbin)
surf_den_13000_t = dblarr(Nbin)
surf_den_13000_err_t = dblarr(Nbin)
sigma_Vr_13000_t = dblarr(Nbin)
slope_13000_t = dblarr(Nbin)
slope_err_13000_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  if initial ge 1500 then goto,out
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_13000_t[i] = mean(phi_13000_t[sort_proj[initial:final]])
  sigma_Vr_13000_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_13000_t[i] = N_bin_t/area
  surf_den_13000_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  xyouts,ai,di-3,strtrim(i,2),alignment=0.5
  fit = poly_fit(alpha_13000_t[initial:final],delta_13000_t[initial:final],1,sigma=sigma)
  slope_13000_t[i] = fit[1]
  slope_err_13000_t[i] = sigma[1] 
print,'f13000t',fit[1],sigma[1]
endfor
out:

;=== in VL 13010 ===
chr_rdtbl,dir_Mbin+'snap13010',0,arr,/silent
arr = double(arr)
t_13010 = (arr[0,0])[0]
x_13010 = reform(arr[1,*])
y_13010 = reform(arr[2,*])
z_13010 = reform(arr[3,*])
vx_13010 = reform(arr[4,*])
vy_13010 = reform(arr[5,*])
vz_13010 = reform(arr[6,*])
Etot_13010 = reform(arr[8,*])
v_13010 = sqrt(Vx_13010^2.+Vy_13010^2.+Vz_13010^2.)

chr_rdtbl,dir_Mbin+'part001_13010',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_13010 = atan(y_13010/(x_13010+8.))*!radeg 
b_13010 = atan(z_13010/sqrt((x_13010+8.)^2.+y_13010^2.))*!radeg
lcosb_13010 = l_13010*cos(b_13010*!dtor)
delta_13010 = asin( cos(b_13010*!dtor)*sin((l_13010-33.)*!dtor)*sin(62.6*!dtor)+sin(b_13010*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_13010 = asin( (cos(b_13010*!dtor)*sin((l_13010-33.)*!dtor)*cos(62.6*!dtor)-sin(b_13010*!dtor)*sin(62.6*!dtor))/cos(delta_13010*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_13010) lt 210)	; orbit within 500Myr
alpha_orbit_13010 = alpha_orbit[sub_orbit]
delta_orbit_13010 = delta_orbit[sub_orbit]

plot,alpha_13010,delta_13010,psym=3,xr=[215,390],yr=[-60,20],/isotropic
oplot,alpha_orbit_13010,delta_orbit_13010,color=100
min_d_curve,alpha_orbit_13010,delta_orbit_13010,alpha_13010,delta_13010,x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
leading_13010 = where(alpha_13010-alpha_13010[0] lt 0)
trailing_13010 = where(alpha_13010-alpha_13010[0] gt 0)
help,leading_13010,trailing_13010
Npts_inbin = 150

Nlead = N_elements(leading_13010)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_13010_l = angsep(alpha_13010[0],delta_13010[0],alpha_proj[leading_13010],delta_proj[leading_13010])
sort_proj = sort(phi_13010_l)
al = alpha_proj[leading_13010[sort_proj]]
dl = delta_proj[leading_13010[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_13010_l = alpha_13010[leading_13010]
delta_13010_l = delta_13010[leading_13010]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_13010-x_sun)*vx_13010+y_13010*vy_13010+z_13010*vz_13010) / (sqrt((x_13010-x_sun)^2.+y_13010^2.+z_13010^2.)*abs(v_13010))
Vradial = v_13010*cos_theta
Vtangen = v_13010*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_13010]

m_phi_13010_l = dblarr(Nbin)
surf_den_13010_l = dblarr(Nbin)
surf_den_13010_err_l = dblarr(Nbin)
sigma_Vr_13010_l = dblarr(Nbin)
slope_13010_l = dblarr(Nbin)
slope_err_13010_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_13010_l[i] = mean(phi_13010_l[sort_proj[initial:final]])
  sigma_Vr_13010_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_13010_l[i] = N_bin_l/area
  surf_den_13010_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_13010_l[initial:final],delta_13010_l[initial:final],1,sigma=sigma)
  slope_13010_l[i] = fit[1]
  slope_err_13010_l[i] = sigma[1] 
print,'f13010l',fit[1],sigma[1]
endfor

Ntrail = N_elements(trailing_13010)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_13010_t = angsep(alpha_13010[0],delta_13010[0],alpha_proj[trailing_13010],delta_proj[trailing_13010])
sort_proj = sort(phi_13010_t)
at = alpha_proj[trailing_13010[sort_proj]]
dt = delta_proj[trailing_13010[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_13010_t = alpha_13010[trailing_13010]
delta_13010_t = delta_13010[trailing_13010]

Vr_t = Vradial[trailing_13010]
m_phi_13010_t = dblarr(Nbin)
surf_den_13010_t = dblarr(Nbin)
surf_den_13010_err_t = dblarr(Nbin)
sigma_Vr_13010_t = dblarr(Nbin)
slope_13010_t = dblarr(Nbin)
slope_err_13010_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_13010_t[i] = mean(phi_13010_t[sort_proj[initial:final]])
  sigma_Vr_13010_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_13010_t[i] = N_bin_t/area
  surf_den_13010_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_13010_t[initial:final],delta_13010_t[initial:final],1,sigma=sigma)
  slope_13010_t[i] = fit[1]
  slope_err_13010_t[i] = sigma[1] 
print,'f13010t',fit[1],sigma[1]
endfor


;=== N-body ===
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
mu = 1.e4
ru = 0.0075
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec
;
;openr,1,dir_nbody+'SNAP081'
;readf,1,Npart,time
;close,1
;t_nbody = round(time*tu*1000.)/1000.
;
;chr_rdtbl,dir_nbody+'SNAP081',1,arr,/silent
;arr = double(arr)
;x_nbody = reform(arr[1,*])*ru
;y_nbody = reform(arr[2,*])*ru
;z_nbody = reform(arr[3,*])*ru
;vx_nbody = reform(arr[4,*])*vu
;vy_nbody = reform(arr[5,*])*vu
;vz_nbody = reform(arr[6,*])*vu
;Etot_nbody = reform(arr[8,*])*vu^2.
;v_nbody = sqrt(Vx_nbody^2.+Vy_nbody^2.+Vz_nbody^2.)*vu^2.
;tunbound_nbody = reform(arr[9,*])
;bounded = where(tunbound_nbody eq 0)
;escaped = where(tunbound_nbody ne 0)

chr_rdtbl,dir_nbody+'backward/snap12760',1,arr,/silent
arr = double(arr)
t_nbody = 12796+(arr[0,0])[0]
x_nbody = reform(arr[1,*])
y_nbody = reform(arr[2,*])
z_nbody = reform(arr[3,*])
vx_nbody = reform(arr[4,*])
vy_nbody = reform(arr[5,*])
vz_nbody = reform(arr[6,*])
Etot_nbody = reform(arr[8,*])
v_nbody = sqrt(Vx_nbody^2.+Vy_nbody^2.+Vz_nbody^2.)
tunbound_nbody = reform(arr[12,*])
bounded = where(tunbound_nbody eq 0)
escaped = where(tunbound_nbody ne 0)


iseed = 1
N_pre = 4050
N_sub = 0
while (N_sub ne 3000) do begin
  sub_selected = round(randomu(iseed,N_pre)*(N_elements(escaped)-1))
  sub_selected = redundant(sub_selected,/non) ; 3000 orbits selected.
  N_sub = N_elements(sub_selected)
  N_pre = (N_sub gt 3000)? N_pre-1:N_pre+1
endwhile
help,sub_selected,escaped

x_nbody = x_nbody[escaped[sub_selected]]
y_nbody = y_nbody[escaped[sub_selected]]
z_nbody = z_nbody[escaped[sub_selected]]
vx_nbody = vx_nbody[escaped[sub_selected]]
vy_nbody = vy_nbody[escaped[sub_selected]]
vz_nbody = vz_nbody[escaped[sub_selected]]
v_nbody = v_nbody[escaped[sub_selected]]

chr_rdtbl,dir_nbody+'SCFCEN',0,arr,/silent
arr = double(arr)
t_orbit = round(reform(arr[0,*])*tu*1000)
x_orbit = reform(arr[2,*])*ru
y_orbit = reform(arr[3,*])*ru
z_orbit = reform(arr[4,*])*ru
;nbody_cen = where(fix(t_orbit*1000) eq fix(t_nbody))
;x_nbody_cen = (x_orbit[nbody_cen])[0]
;y_nbody_cen = (y_orbit[nbody_cen])[0]
;z_nbody_cen = (z_orbit[nbody_cen])[0]
chr_rdtbl,dir_nbody+'backward/SCFCEN_snap12760',0,arr,/silent
x_nbody_cen = (double(arr[1,*]))[0]
y_nbody_cen = (double(arr[2,*]))[0]
z_nbody_cen = (double(arr[3,*]))[0]

;=== convert XYZ into RA & Dec ===
l_nbody = atan(y_nbody/(x_nbody+8.))*!radeg 
b_nbody = atan(z_nbody/sqrt((x_nbody+8.)^2.+y_nbody^2.))*!radeg
lcosb_nbody = l_nbody*cos(b_nbody*!dtor)
delta_nbody = asin( cos(b_nbody*!dtor)*sin((l_nbody-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_nbody = asin( (cos(b_nbody*!dtor)*sin((l_nbody-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody*!dtor)*sin(62.6*!dtor))/cos(delta_nbody*!dtor) )*!radeg + 282.25

l_nbody_cen = atan(y_nbody_cen/(x_nbody_cen+8.))*!radeg 
b_nbody_cen = atan(z_nbody_cen/sqrt((x_nbody_cen+8.)^2.+y_nbody_cen^2.))*!radeg
lcosb_nbody_cen = l_nbody_cen*cos(b_nbody_cen*!dtor)
delta_nbody_cen = asin( cos(b_nbody_cen*!dtor)*sin((l_nbody_cen-33.)*!dtor)*sin(62.6*!dtor)+sin(b_nbody_cen*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_nbody_cen = asin( (cos(b_nbody_cen*!dtor)*sin((l_nbody_cen-33.)*!dtor)*cos(62.6*!dtor)-sin(b_nbody_cen*!dtor)*sin(62.6*!dtor))/cos(delta_nbody_cen*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

t_orbit = t_orbit_nosub
sub_orbit = where(abs(t_orbit-t_nbody) lt 00250)	; orbit within 120Myr
help,t_orbit,t_nbody,sub_orbit
;alpha_orbit_nbody = alpha_orbit[sub_orbit]
;delta_orbit_nbody = delta_orbit[sub_orbit]
alpha_orbit_nbody = alpha_orbit_ns[sub_orbit]
delta_orbit_nbody = delta_orbit_ns[sub_orbit]

plot,alpha_nbody,delta_nbody,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_nbody,delta_orbit_nbody,color=100
min_d_curve,alpha_orbit_nbody,delta_orbit_nbody,alpha_nbody,delta_nbody,x_proj=alpha_proj,y_proj=delta_proj


;=== estimate density along the orbit ===
leading_nbody = where(delta_nbody-delta_nbody_cen lt 0)
trailing_nbody = where(delta_nbody-delta_nbody_cen gt 0)
;leading_nbody = where(Etot_nbody-Etot_nbody[0] lt 0)
;trailing_nbody = where(Etot_nbody-Etot_nbody[0] gt 0)

help,leading_nbody,trailing_nbody
Npts_inbin = 150

Nlead = N_elements(leading_nbody)
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
phi_nbody_l = angsep(alpha_nbody_cen,delta_nbody_cen,alpha_proj[leading_nbody],delta_proj[leading_nbody])
sort_proj = sort(phi_nbody_l)
al = alpha_proj[leading_nbody[sort_proj]]
dl = delta_proj[leading_nbody[sort_proj]]
oplot,al,dl,psym=7,color=255
alpha_nbody_l = alpha_nbody[leading_nbody]
delta_nbody_l = delta_nbody[leading_nbody]

;=== estimate sigma_v along the orbit ===
cos_theta = ((x_nbody-x_sun)*vx_nbody+y_nbody*vy_nbody+z_nbody*vz_nbody) / (sqrt((x_nbody-x_sun)^2.+y_nbody^2.+z_nbody^2.)*abs(v_nbody))
Vradial = v_nbody*cos_theta
Vtangen = v_nbody*sqrt(1-cos_theta^2)
Vr_l = Vradial[leading_nbody]

m_phi_nbody_l = dblarr(Nbin)
surf_den_nbody_l = dblarr(Nbin)
surf_den_nbody_err_l = dblarr(Nbin)
sigma_Vr_nbody_l = dblarr(Nbin)
slope_nbody_l = dblarr(Nbin)
slope_err_nbody_l = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_nbody_l[i] = mean(phi_nbody_l[sort_proj[initial:final]])
  sigma_Vr_nbody_l[i] = stddev(Vr_l[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_nbody_l[i] = N_bin_l/area
  surf_den_nbody_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_nbody_l[initial:final],delta_nbody_l[initial:final],1,sigma=sigma)
  slope_nbody_l[i] = fit[1]
  slope_err_nbody_l[i] = sigma[1] 
endfor

Ntrail = N_elements(trailing_nbody)
Nbin = (Ntrail/Npts_inbin eq Ntrail/double(Npts_inbin)) ? Ntrail/Npts_inbin:Ntrail/Npts_inbin+1
phi_nbody_t = angsep(alpha_nbody_cen,delta_nbody_cen,alpha_proj[trailing_nbody],delta_proj[trailing_nbody])
sort_proj = sort(phi_nbody_t)
at = alpha_proj[trailing_nbody[sort_proj]]
dt = delta_proj[trailing_nbody[sort_proj]]
oplot,at,dt,psym=7,color=215
alpha_nbody_t = alpha_nbody[trailing_nbody]
delta_nbody_t = delta_nbody[trailing_nbody]

Vr_t = Vradial[trailing_nbody]
m_phi_nbody_t = dblarr(Nbin)
surf_den_nbody_t = dblarr(Nbin)
surf_den_nbody_err_t = dblarr(Nbin)
sigma_Vr_nbody_t = dblarr(Nbin)
slope_nbody_t = dblarr(Nbin)
slope_err_nbody_t = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntrail-1) then final = Ntrail-1
  ai = at[initial]
  di = dt[initial]
  af = at[final]
  df = dt[final]
  N_bin_t = final-initial+1
  m_phi_nbody_t[i] = mean(phi_nbody_t[sort_proj[initial:final]])
  sigma_Vr_nbody_t[i] = stddev(Vr_t[sort_proj[initial:final]])
  area = angsep(ai,di,af,df)^2.
  surf_den_nbody_t[i] = N_bin_t/area
  surf_den_nbody_err_t[i] = sqrt(N_bin_t)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
  fit = poly_fit(alpha_nbody_t[initial:final],delta_nbody_t[initial:final],1,sigma=sigma)
  slope_nbody_t[i] = fit[1]
  slope_err_nbody_t[i] = sigma[1] 
endfor
oplot,[alpha_nbody_cen],[delta_nbody_cen],psym=1,color=70,thick=8

plot,alpha_nbody,delta_nbody,psym=3,xr=[205,265],yr=[-73,17],/isotropic
oplot,alpha_orbit_nbody,delta_orbit_nbody,color=100
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=200,linestyle=2
oplot,[alpha_nbody_cen],[delta_nbody_cen],psym=1,color=70,thick=8


;============
;=== PLOT ===
;============
!p.multi=[0,7,4]
!p.charsize=1.3
RAmin = 255
RAmax = 205
Decmin = -23
Decmax = 11
erase
multiplot
plot,alpha_nosub,delta_nosub,psym=3,ytitle='Dec',xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
legend,['nosub'],box=0,/bottom,charsize=0.8
for i=0,14 do begin
  xl = ['240','230','220']
  xyouts,240-i*10,-27,xl[round((i/3.-i/3)*3)],alignment=0.5
endfor
for i=0,4 do xyouts,230-i*30,-31,'RA',alignment=0.5

multiplot
plot,alpha_5p,delta_5p,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_5p,delta_orbit_5p,color=100
legend,['5p: 27'],box=0,/bottom,charsize=0.8
multiplot
plot,alpha_10p,delta_10p,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_10p,delta_orbit_10p,color=100
legend,['10p: 64'],box=0,/bottom,charsize=0.8
multiplot
plot,alpha_allsub,delta_allsub,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_allsub,delta_orbit_allsub,color=100
legend,['allsub: 1000'],box=0,/bottom,charsize=0.8
multiplot
plot,alpha_13000,delta_13000,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_13000,delta_orbit_13000,color=100
legend,['in VL: t=13.000Gyr'],box=0,/bottom,charsize=0.8
multiplot
plot,alpha_13010,delta_13010,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_13010,delta_orbit_13010,color=100
legend,['in VL: t=13.010Gyr'],box=0,/bottom,charsize=0.8
multiplot
plot,alpha_nbody,delta_nbody,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_nbody,delta_orbit_nbody,color=100
legend,['nbody nosub'],box=0,/bottom,charsize=0.8

phimin = 0
phimax = 34
denmin = 0.4
denmax = 9000
Vr_min = 0
Vr_max = 36
;=== nosub density profile ===
multiplot
ploterror,m_phi_nosub_l,surf_den_nosub_l,surf_den_nosub_err_l,xr=[phimin,phimax],yr=[denmin,denmax],ytitle=textoidl('\Sigma [n/deg^2]'),/ylog;,pos=[0.115,0.34,0.285,0.58]
oploterror,m_phi_nosub_t,surf_den_nosub_t,surf_den_nosub_err_t,linestyle=2,color=255,errcolor=255
legend,['leading','trailing'],linestyle=[0,2],color=[0,255],box=0,charsize=0.8

multiplot
ploterror,m_phi_5p_l,surf_den_5p_l,surf_den_5p_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.285,0.34,0.455,0.58]
oploterror,m_phi_5p_t,surf_den_5p_t,surf_den_5p_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_10p_l,surf_den_10p_l,surf_den_10p_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.455,0.34,0.625,0.58]
oploterror,m_phi_10p_t,surf_den_10p_t,surf_den_10p_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_allsub_l,surf_den_allsub_l,surf_den_allsub_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_allsub_t,surf_den_allsub_t,surf_den_allsub_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_13000_l,surf_den_13000_l,surf_den_13000_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_13000_t,surf_den_13000_t,surf_den_13000_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_13010_l,surf_den_13010_l,surf_den_13010_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_13010_t,surf_den_13010_t,surf_den_13010_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_nbody_l,surf_den_nbody_l,surf_den_nbody_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog;,pos=[0.795,0.34,0.965,0.58]
oploterror,m_phi_nbody_t,surf_den_nbody_t,surf_den_nbody_err_t,linestyle=2,color=255,errcolor=255


;=== nosub sigma_v profile ===
multiplot
plot,m_phi_nosub_l,sigma_Vr_nosub_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]');,pos=[0.115,0.10,0.285,0.34]
oplot,m_phi_nosub_t,sigma_Vr_nosub_t,linestyle=2,color=255

multiplot
plot,m_phi_5p_l,sigma_Vr_5p_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.285,0.10,0.455,0.34]
oplot,m_phi_5p_t,sigma_Vr_5p_t,linestyle=2,color=255

multiplot
plot,m_phi_10p_l,sigma_Vr_10p_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.455,0.10,0.625,0.34]
oplot,m_phi_10p_t,sigma_Vr_10p_t,linestyle=2,color=255

multiplot
plot,m_phi_allsub_l,sigma_Vr_allsub_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.625,0.10,0.795,0.34]
oplot,m_phi_allsub_t,sigma_Vr_allsub_t,linestyle=2,color=255

multiplot
plot,m_phi_13000_l,sigma_Vr_13000_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.625,0.10,0.795,0.34]
oplot,m_phi_13000_t,sigma_Vr_13000_t,linestyle=0,color=255

multiplot
plot,m_phi_13010_l,sigma_Vr_13010_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.625,0.10,0.795,0.34]
oplot,m_phi_13010_t,sigma_Vr_13010_t,linestyle=0,color=255

multiplot
plot,m_phi_nbody_l,sigma_Vr_nbody_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]');,pos=[0.795,0.10,0.965,0.34]
oplot,m_phi_nbody_t,sigma_Vr_nbody_t,linestyle=2,color=255


;=== phi direction profile ===
dphimin = -1
dphimax = 4
multiplot
ploterror,m_phi_nosub_l,slope_nosub_l,slope_err_nosub_l,xr=[phimin,phimax],yr=[dphimin,dphimax],ytitle=textoidl('tan \theta');,pos=[0.115,0.34,0.285,0.58]
oploterror,m_phi_nosub_t,slope_nosub_t,slope_err_nosub_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_5p_l,slope_5p_l,slope_err_5p_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.285,0.34,0.455,0.58]
oploterror,m_phi_5p_t,slope_5p_t,slope_err_5p_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_10p_l,slope_10p_l,slope_err_10p_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.455,0.34,0.625,0.58]
oploterror,m_phi_10p_t,slope_10p_t,slope_err_10p_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_allsub_l,slope_allsub_l,slope_err_allsub_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_allsub_t,slope_allsub_t,slope_err_allsub_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_13000_l,slope_13000_l,slope_err_13000_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_13000_t,slope_13000_t,slope_err_13000_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_13010_l,slope_13010_l,slope_err_13010_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.625,0.34,0.795,0.58]
oploterror,m_phi_13010_t,slope_13010_t,slope_err_13010_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_nbody_l,slope_nbody_l,slope_err_nbody_l,xr=[phimin,phimax],yr=[dphimin,dphimax];,pos=[0.795,0.34,0.965,0.58]
oploterror,m_phi_nbody_t,slope_nbody_t,slope_err_nbody_t,linestyle=2,color=255,errcolor=255



multiplot,/reset
erase

!p.multi=[0,2,2]
plot,alpha_nosub,v_nosub,psym=3,xr=[221,234],yr=[65,160]
plot,alpha_5p,v_5p,psym=3,xr=[221,234],yr=[65,160]
plot,alpha_10p,v_10p,psym=3,xr=[221,234],yr=[65,160]
plot,alpha_allsub,v_allsub,psym=3,xr=[221,234],yr=[65,160]
device,/close


device,file='den_profile_kathryn.ps',/color,/portrait,ysize=18,yoffset=3.0
!p.charsize=1.2
phimin = 0
phimax = 34
denmin = 0.4
denmax = 10000
Vr_min = 0
Vr_max = 19
!p.multi=[0,3,3]
multiplot
plot,alpha_nosub,delta_nosub,psym=3,ytitle='Dec',xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
legend,['nosub'],box=0,/bottom,charsize=0.8
for i=0,14 do begin
  xl = ['250','240','230','220','210']
  xyouts,250-i*10,-29,xl[round((i/5.-i/5)*5)],alignment=0.5
endfor
for i=0,2 do xyouts,230-i*50,-35,'RA',alignment=0.5

multiplot
plot,alpha_13000,delta_13000,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_13000,delta_orbit_13000,color=100
legend,['in VL: t=13.000Gyr'],box=0,/bottom,charsize=0.8

multiplot
plot,alpha_nbody,delta_nbody,psym=3,xr=[RAmin,RAmax],yr=[Decmin,Decmax],/isotropic
oplot,alpha_orbit_nbody,delta_orbit_nbody,color=100
legend,['nbody nosub'],box=0,/bottom,charsize=0.8

multiplot
ploterror,m_phi_nosub_l,surf_den_nosub_l,surf_den_nosub_err_l,xr=[phimin,phimax],yr=[denmin,denmax],ytitle=textoidl('\Sigma [n/deg^2]'),/ylog,pos=[0.1505,0.335,0.4190,0.58]
oploterror,m_phi_nosub_t,surf_den_nosub_t,surf_den_nosub_err_t,linestyle=2,color=255,errcolor=255
legend,['leading','trailing'],linestyle=[0,2],color=[0,255],box=0,charsize=0.8

multiplot
ploterror,m_phi_13000_l,surf_den_13000_l,surf_den_13000_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog,pos=[0.4190,0.335,0.6875,0.58]
oploterror,m_phi_13000_t,surf_den_13000_t,surf_den_13000_err_t,linestyle=2,color=255,errcolor=255

multiplot
ploterror,m_phi_nbody_l,surf_den_nbody_l,surf_den_nbody_err_l,xr=[phimin,phimax],yr=[denmin,denmax],/ylog,pos=[0.6875,0.335,0.9555,0.58]
oploterror,m_phi_nbody_t,surf_den_nbody_t,surf_den_nbody_err_t,linestyle=2,color=255,errcolor=255

multiplot
plot,m_phi_nosub_l,sigma_Vr_nosub_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]'),pos=[0.1505,0.09,0.4190,0.335]
oplot,m_phi_nosub_t,sigma_Vr_nosub_t,linestyle=2,color=255

multiplot
plot,m_phi_13000_l,sigma_Vr_13000_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),pos=[0.4190,0.09,0.6875,0.335]
oplot,m_phi_13000_t,sigma_Vr_13000_t,linestyle=2,color=255


multiplot
plot,m_phi_nbody_l,sigma_Vr_nbody_l,xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),pos=[0.6875,0.09,0.9555,0.335]
oplot,m_phi_nbody_t,sigma_Vr_nbody_t,linestyle=2,color=255
multiplot,/reset
device,/close
END
