pro den_profile


mass_name = 'e7'
dir_close = '/scratch/jhyoon/Research/LG/close_encounter/'

dir_1Rs_part = dir_close+'encounter_1Rs_'+mass_name+'/snapshot/'
dir_5Rs_part = dir_close+'encounter_5Rs_'+mass_name+'/snapshot/'
dir_10Rs_part = dir_close+'encounter_10Rs_'+mass_name+'/snapshot/'
dir_20Rs_part = dir_close+'encounter_20Rs_'+mass_name+'/snapshot/'
dir_1Rs_orbit = dir_close+'encounter_1Rs_'+mass_name+'/'
dir_5Rs_orbit = dir_close+'encounter_5Rs_'+mass_name+'/'
dir_10Rs_orbit = dir_close+'encounter_10Rs_'+mass_name+'/'
dir_20Rs_orbit = dir_close+'encounter_20Rs_'+mass_name+'/'
dir_BI5_part = '/scratch/jhyoon/Research/LG/BI5_test/snapshot_e7_1000sub/'
dir_BI5_orbit = '/scratch/jhyoon/Research/LG/BI5_test/'

dir_nosub = '/scratch/jhyoon/Research/LG/realistic/snapshot_nosub/'
dir_orbit = '/scratch/jhyoon/Research/LG/realistic/'

dir_nbody = '/scratch/jhyoon/Research/LG/Nbody_new/pal5/ten4_new/'

x_sun = -8.d

;=== no subhalo ===
chr_rdtbl,dir_nosub+'snap12760',0,arr,/silent
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

chr_rdtbl,dir_orbit+'part001_nosub',0,arr,/silent
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
device,file='den_profile.ps',/color,/landscape
;plot,alpha_nosub,delta_nosub,psym=3,xr=[200,260],yr=[-83,17],/isotropic

plot,alpha_nosub,delta_nosub,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
min_d_curve,alpha_orbit_nosub,delta_orbit_nosub,alpha_nosub,delta_nosub,x_proj=alpha_proj,y_proj=delta_proj


;=== estimate density along the orbit ===
leading_nosub = where(Etot_nosub-Etot_nosub[0] lt 0)
trailing_nosub = where(Etot_nosub-Etot_nosub[0] gt 0)
help,leading_nosub,trailing_nosub

;Npts = N_elements(x_nosub)
Nlead = N_elements(leading_nosub)
Npts_inbin = 150
Nbin = (Nlead/Npts_inbin eq Nlead/double(Npts_inbin)) ? Nlead/Npts_inbin:Nlead/Npts_inbin+1
;Nbin = Npts/Npts_inbin+1

;sort_proj = sort(alpha_proj)
;eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line
;oplot,m_alpha,m_delta,psym=1
;vline,bin_line
;phi_nosub = angsep(alpha_nosub[0],delta_nosub[0],m_alpha,m_delta)
phi_nosub_l = angsep(alpha_nosub[0],delta_nosub[0],alpha_proj[leading_nosub],delta_proj[leading_nosub])
sort_proj = sort(phi_nosub_l)
al = alpha_proj[leading_nosub[sort_proj]]
dl = delta_proj[leading_nosub[sort_proj]]
oplot,al,dl,psym=7,color=255
m_phi_nosub_l = dblarr(Nbin)
surf_den_nosub_l = dblarr(Nbin)
surf_den_nosub_err_l = dblarr(Nbin)

for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Nlead-1) then final = Nlead-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin_l = final-initial+1
  m_phi_nosub_l[i] = mean(phi_nosub_l[sort_proj[leading_nosub[initial:final]]])
  area = angsep(ai,di,af,df)^2.
  surf_den_nosub_l[i] = N_bin_l/area
  surf_den_nosub_err_l[i] = sqrt(N_bin_l)/area
  oplot,[ai],[di],psym=6
  oplot,[af],[df],psym=5
;  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.
;  oplot,[alpha_proj[sort_proj[initial]]],[delta_proj[sort_proj[initial]]],psym=6
;  oplot,[alpha_proj[sort_proj[final]]],[delta_proj[sort_proj[final]]],psym=5
endfor


;surf_den_nosub = Npts_inbin/area
;surf_den_nosub_err = sqrt(Npts_inbin)/area
;leading_nosub = where(m_alpha-alpha_nosub[0] lt 0,complement=trailing_nosub)


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_nosub-x_sun)*vx_nosub+y_nosub*vy_nosub+z_nosub*vz_nosub) / (sqrt((x_nosub-x_sun)^2.+y_nosub^2.+z_nosub^2.)*abs(v_nosub))
Vradial = v_nosub*cos_theta
Vtangen = v_nosub*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_nosub,num=Nbin,bin_line=bin_line


;===========
;=== 1Rs ===
;===========
chr_rdtbl,dir_1Rs_part+'snap12760',0,arr,/silent
arr = double(arr)
t_1Rs = (arr[0,0])[0]
x_1Rs = reform(arr[1,*])
y_1Rs = reform(arr[2,*])
z_1Rs = reform(arr[3,*])
vx_1Rs = reform(arr[4,*])
vy_1Rs = reform(arr[5,*])
vz_1Rs = reform(arr[6,*])
Etot_1Rs = reform(arr[8,*])
v_1Rs = sqrt(Vx_1Rs^2.+Vy_1Rs^2.+Vz_1Rs^2.)

chr_rdtbl,dir_1Rs_orbit+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_1Rs = atan(y_1Rs/(x_1Rs+8.))*!radeg 
b_1Rs = atan(z_1Rs/sqrt((x_1Rs+8.)^2.+y_1Rs^2.))*!radeg
lcosb_1Rs = l_1Rs*cos(b_1Rs*!dtor)
delta_1Rs = asin( cos(b_1Rs*!dtor)*sin((l_1Rs-33.)*!dtor)*sin(62.6*!dtor)+sin(b_1Rs*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_1Rs = asin( (cos(b_1Rs*!dtor)*sin((l_1Rs-33.)*!dtor)*cos(62.6*!dtor)-sin(b_1Rs*!dtor)*sin(62.6*!dtor))/cos(delta_1Rs*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_1Rs) lt 120)	; orbit within 500Myr
alpha_orbit_1Rs = alpha_orbit[sub_orbit]
delta_orbit_1Rs = delta_orbit[sub_orbit]

min_d_curve,alpha_orbit_1Rs,delta_orbit_1Rs,alpha_1Rs[1:*],delta_1Rs[1:*],x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
Npts = N_elements(x_1Rs)
Npts_inbin = 150
Nbin = Npts/Npts_inbin
sort_proj = sort(alpha_proj)
eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line

phi_1Rs = angsep(alpha_1Rs[0],delta_1Rs[0],m_alpha,m_delta)
area = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (i eq Nbin-1) then final = initial + Npts_inbin-2
  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.
endfor

surf_den_1Rs = Npts_inbin/area
surf_den_1Rs_err = sqrt(Npts_inbin)/area
leading_1Rs = where(m_alpha-alpha_1Rs[0] lt 0,complement=trailing_1Rs)


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_1Rs-x_sun)*vx_1Rs+y_1Rs*vy_1Rs+z_1Rs*vz_1Rs) / (sqrt((x_1Rs-x_sun)^2.+y_1Rs^2.+z_1Rs^2.)*abs(v_1Rs))
vradial = v_1Rs*cos_theta
vtangen = v_1Rs*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_1Rs,num=Nbin,bin_line=bin_line

plot,alpha_1Rs,delta_1Rs,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_1Rs,delta_orbit_1Rs,color=100
oplot,m_alpha,m_delta,psym=1
vline,bin_line


;===========
;=== 5Rs ===
;===========
chr_rdtbl,dir_5Rs_part+'snap12760',0,arr,/silent
arr = double(arr)
t_5Rs = (arr[0,0])[0]
x_5Rs = reform(arr[1,*])
y_5Rs = reform(arr[2,*])
z_5Rs = reform(arr[3,*])
vx_5Rs = reform(arr[4,*])
vy_5Rs = reform(arr[5,*])
vz_5Rs = reform(arr[6,*])
Etot_5Rs = reform(arr[8,*])
v_5Rs = sqrt(Vx_5Rs^2.+Vy_5Rs^2.+Vz_5Rs^2.)

chr_rdtbl,dir_5Rs_orbit+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_5Rs = atan(y_5Rs/(x_5Rs+8.))*!radeg 
b_5Rs = atan(z_5Rs/sqrt((x_5Rs+8.)^2.+y_5Rs^2.))*!radeg
lcosb_5Rs = l_5Rs*cos(b_5Rs*!dtor)
delta_5Rs = asin( cos(b_5Rs*!dtor)*sin((l_5Rs-33.)*!dtor)*sin(62.6*!dtor)+sin(b_5Rs*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_5Rs = asin( (cos(b_5Rs*!dtor)*sin((l_5Rs-33.)*!dtor)*cos(62.6*!dtor)-sin(b_5Rs*!dtor)*sin(62.6*!dtor))/cos(delta_5Rs*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_5Rs) lt 120)	; orbit within 500Myr
alpha_orbit_5Rs = alpha_orbit[sub_orbit]
delta_orbit_5Rs = delta_orbit[sub_orbit]

min_d_curve,alpha_orbit_5Rs,delta_orbit_5Rs,alpha_5Rs[1:*],delta_5Rs[1:*],x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
Npts = N_elements(x_5Rs)
Npts_inbin = 150
Nbin = Npts/Npts_inbin
sort_proj = sort(alpha_proj)
eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line

phi_5Rs = angsep(alpha_5Rs[0],delta_5Rs[0],m_alpha,m_delta)
area = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (i eq Nbin-1) then final = initial + Npts_inbin-2
  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.
endfor

surf_den_5Rs = Npts_inbin/area
surf_den_5Rs_err = sqrt(Npts_inbin)/area
leading_5Rs = where(m_alpha-alpha_5Rs[0] lt 0,complement=trailing_5Rs)


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_5Rs-x_sun)*vx_5Rs+y_5Rs*vy_5Rs+z_5Rs*vz_5Rs) / (sqrt((x_5Rs-x_sun)^2.+y_5Rs^2.+z_5Rs^2.)*abs(v_5Rs))
vradial = v_5Rs*cos_theta
vtangen = v_5Rs*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_5Rs,num=Nbin,bin_line=bin_line

plot,alpha_5Rs,delta_5Rs,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_5Rs,delta_orbit_5Rs,color=100
oplot,m_alpha,m_delta,psym=1
vline,bin_line


;===========
;=== 10Rs ===
;===========
chr_rdtbl,dir_10Rs_part+'snap12760',0,arr,/silent
arr = double(arr)
t_10Rs = (arr[0,0])[0]
x_10Rs = reform(arr[1,*])
y_10Rs = reform(arr[2,*])
z_10Rs = reform(arr[3,*])
vx_10Rs = reform(arr[4,*])
vy_10Rs = reform(arr[5,*])
vz_10Rs = reform(arr[6,*])
Etot_10Rs = reform(arr[8,*])
v_10Rs = sqrt(Vx_10Rs^2.+Vy_10Rs^2.+Vz_10Rs^2.)

chr_rdtbl,dir_10Rs_orbit+'part001',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])
x_orbit = reform(arr[1,*])
y_orbit = reform(arr[2,*])
z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
l_10Rs = atan(y_10Rs/(x_10Rs+8.))*!radeg 
b_10Rs = atan(z_10Rs/sqrt((x_10Rs+8.)^2.+y_10Rs^2.))*!radeg
lcosb_10Rs = l_10Rs*cos(b_10Rs*!dtor)
delta_10Rs = asin( cos(b_10Rs*!dtor)*sin((l_10Rs-33.)*!dtor)*sin(62.6*!dtor)+sin(b_10Rs*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_10Rs = asin( (cos(b_10Rs*!dtor)*sin((l_10Rs-33.)*!dtor)*cos(62.6*!dtor)-sin(b_10Rs*!dtor)*sin(62.6*!dtor))/cos(delta_10Rs*!dtor) )*!radeg + 282.25

l_orbit = atan(y_orbit/(x_orbit+8.))*!radeg 
b_orbit = atan(z_orbit/sqrt((x_orbit+8.)^2.+y_orbit^2.))*!radeg
lcosb_orbit = l_orbit*cos(b_orbit*!dtor)
delta_orbit = asin( cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*sin(62.6*!dtor)+sin(b_orbit*!dtor)*cos(62.6*!dtor) )*!radeg
alpha_orbit = asin( (cos(b_orbit*!dtor)*sin((l_orbit-33.)*!dtor)*cos(62.6*!dtor)-sin(b_orbit*!dtor)*sin(62.6*!dtor))/cos(delta_orbit*!dtor) )*!radeg + 282.25

sub_orbit = where(abs(t_orbit-t_10Rs) lt 120)	; orbit within 500Myr
alpha_orbit_10Rs = alpha_orbit[sub_orbit]
delta_orbit_10Rs = delta_orbit[sub_orbit]

min_d_curve,alpha_orbit_10Rs,delta_orbit_10Rs,alpha_10Rs[1:*],delta_10Rs[1:*],x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
Npts = N_elements(x_10Rs)
Npts_inbin = 150
Nbin = Npts/Npts_inbin
sort_proj = sort(alpha_proj)
eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line

phi_10Rs = angsep(alpha_10Rs[0],delta_10Rs[0],m_alpha,m_delta)
area = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (i eq Nbin-1) then final = initial + Npts_inbin-2
  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.
endfor

surf_den_10Rs = Npts_inbin/area
surf_den_10Rs_err = sqrt(Npts_inbin)/area
leading_10Rs = where(m_alpha-alpha_10Rs[0] lt 0,complement=trailing_10Rs)


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_10Rs-x_sun)*vx_10Rs+y_10Rs*vy_10Rs+z_10Rs*vz_10Rs) / (sqrt((x_10Rs-x_sun)^2.+y_10Rs^2.+z_10Rs^2.)*abs(v_10Rs))
vradial = v_10Rs*cos_theta
vtangen = v_10Rs*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_10Rs,num=Nbin,bin_line=bin_line

;plot,alpha_10Rs,delta_10Rs,psym=3,xr=[215,245],yr=[-23,7],/isotropic
;oplot,alpha_orbit_10Rs,delta_orbit_10Rs,color=100
;oplot,m_alpha,m_delta,psym=1
;vline,bin_line


;=== ALL subhalo ===
chr_rdtbl,dir_BI5_part+'snap12760',0,arr,/silent
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

min_d_curve,alpha_orbit_allsub,delta_orbit_allsub,alpha_allsub[1:*],delta_allsub[1:*],x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
Npts = N_elements(x_allsub)
Npts_inbin = 150
Nbin = Npts/Npts_inbin
sort_proj = sort(alpha_proj)
eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line

phi_allsub = angsep(alpha_allsub[0],delta_allsub[0],m_alpha,m_delta)
area = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (i eq Nbin-1) then final = initial + Npts_inbin-2
  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.
endfor

surf_den_allsub = Npts_inbin/area
surf_den_allsub_err = sqrt(Npts_inbin)/area
leading_allsub = where(m_alpha-alpha_allsub[0] lt 0,complement=trailing_allsub)


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_allsub-x_sun)*vx_allsub+y_allsub*vy_allsub+z_allsub*vz_allsub) / (sqrt((x_allsub-x_sun)^2.+y_allsub^2.+z_allsub^2.)*abs(v_allsub))
Vradial = v_allsub*cos_theta
Vtangen = v_allsub*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_allsub,num=Nbin,bin_line=bin_line

plot,alpha_allsub,delta_allsub,psym=3,xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_allsub,delta_orbit_allsub,color=100
oplot,m_alpha,m_delta,psym=1
vline,bin_line


;=== N-body ===
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
mu = 1.e4
ru = 0.0075
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec

openr,1,dir_nbody+'SNAP081'
readf,1,Npart,time
close,1
t_nbody = round(time*tu*1000.)/1000.

chr_rdtbl,dir_nbody+'SNAP081',1,arr,/silent
arr = double(arr)
x_nbody = reform(arr[1,*])*ru
y_nbody = reform(arr[2,*])*ru
z_nbody = reform(arr[3,*])*ru
vx_nbody = reform(arr[4,*])*vu
vy_nbody = reform(arr[5,*])*vu
vz_nbody = reform(arr[6,*])*vu
Etot_nbody = reform(arr[8,*])*vu^2.
v_nbody = sqrt(Vx_nbody^2.+Vy_nbody^2.+Vz_nbody^2.)*vu^2.
tunbound_nbody = reform(arr[9,*])
bounded = where(tunbound_nbody eq 0)
escaped = where(tunbound_nbody ne 0)
help,escaped
x_nbody = x_nbody[escaped]
y_nbody = y_nbody[escaped]
z_nbody = z_nbody[escaped]
vx_nbody = vx_nbody[escaped]
vy_nbody = vy_nbody[escaped]
vz_nbody = vz_nbody[escaped]
v_nbody = v_nbody[escaped]

chr_rdtbl,dir_nbody+'SCFCEN',0,arr,/silent
arr = double(arr)
t_orbit = reform(arr[0,*])*tu
x_orbit = reform(arr[2,*])*ru
y_orbit = reform(arr[3,*])*ru
z_orbit = reform(arr[4,*])*ru
nbody_cen = where(fix(t_orbit*1000) eq fix(t_nbody*1000))
x_nbody_cen = x_orbit[nbody_cen]
y_nbody_cen = y_orbit[nbody_cen]
z_nbody_cen = z_orbit[nbody_cen]

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

t_orbit = t_orbit_nosub/1000
sub_orbit = where(abs(t_orbit-t_nbody) lt 0.250)	; orbit within 120Myr
;alpha_orbit_nbody = alpha_orbit[sub_orbit]
;delta_orbit_nbody = delta_orbit[sub_orbit]
alpha_orbit_nbody = alpha_orbit_ns[sub_orbit]
delta_orbit_nbody = delta_orbit_ns[sub_orbit]

min_d_curve,alpha_orbit_nbody,delta_orbit_nbody,alpha_nbody,delta_nbody,x_proj=alpha_proj,y_proj=delta_proj

;=== estimate density along the orbit ===
Npts = N_elements(x_nbody)
Npts_inbin = 300;150
Nbin = Npts/Npts_inbin
;eqnum_bin,alpha_proj[sort_proj],delta_proj[sort_proj],m_alpha,m_delta,err_alpha,err_delta,num=Nbin,bin_line=bin_line
;eqnum_bin,delta_proj[sort_proj],alpha_proj[sort_proj],m_delta,m_alpha,err_delta,err_alpha,num=Nbin,bin_line=bin_line

phi_pts = angsep(alpha_nbody_cen,delta_nbody_cen,alpha_nbody,delta_nbody)
sort_proj = sort(phi_nbody)

phi_nbody = angsep(alpha_nbody_cen,delta_nbody_cen,m_alpha,m_delta)
leading_nbody = where(m_delta-delta_nbody_cen lt 0,complement=trailing_nbody)

area = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (i eq Nbin-1) then final = initial + Npts_inbin-2
  area[i] = angsep(alpha_proj[sort_proj[initial]],delta_proj[sort_proj[initial]],alpha_proj[sort_proj[final]],delta_proj[sort_proj[final]])^2.

endfor

surf_den_nbody = Npts_inbin/area
surf_den_nbody_err = sqrt(Npts_inbin)/area


;=== estimate sigma_v along the orbit ===
cos_theta = ((x_nbody-x_sun)*vx_nbody+y_nbody*vy_nbody+z_nbody*vz_nbody) / (sqrt((x_nbody-x_sun)^2.+y_nbody^2.+z_nbody^2.)*abs(v_nbody))
Vradial = v_nbody*cos_theta
Vtangen = v_nbody*sqrt(1-cos_theta^2)
eqnum_bin,alpha_proj[sort_proj],Vradial[sort_proj],m_alpha,m_V,err_alpha,sigma_Vr_nbody,num=Nbin,bin_line=bin_line

;plot,alpha_nbody,delta_nbody,psym=3,xr=[215,245],yr=[-23,7],/isotropic
plot,alpha_nbody,delta_nbody,psym=3,xr=[205,265],yr=[-73,17],/isotropic
oplot,alpha_nbody_cen,delta_nbody_cen,psym=1,color=255,thick=8
oplot,alpha_orbit_nbody,delta_orbit_nbody,color=100
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=200,linestyle=2
oplot,m_alpha,m_delta,psym=1


;============
;=== PLOT ===
;============
!p.multi=[0,4,3]
plot,alpha_nosub,delta_nosub,psym=3,xtitle='RA',ytitle='Dec',xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_nosub,delta_orbit_nosub,color=100
legend,['nosub'],box=0,/bottom,/right,charsize=0.8
plot,alpha_1Rs,delta_1Rs,psym=3,xtitle='RA',ytitle='Dec',xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_1Rs,delta_orbit_1Rs,color=100
legend,['1Rs: 27'],box=0,/bottom,/right,charsize=0.8
plot,alpha_5Rs,delta_5Rs,psym=3,xtitle='RA',ytitle='Dec',xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_5Rs,delta_orbit_5Rs,color=100
legend,['5Rs: 64'],box=0,/bottom,/right,charsize=0.8
;plot,alpha_10Rs,delta_10Rs,psym=3,xtitle='RA',ytitle='Dec',xr=[215,245],yr=[-23,7],/isotropic
;oplot,alpha_orbit_10Rs,delta_orbit_10Rs,color=100
;legend,['10Rs: 102'],box=0,/bottom,/right
plot,alpha_allsub,delta_allsub,psym=3,xtitle='RA',ytitle='Dec',xr=[215,245],yr=[-23,7],/isotropic
oplot,alpha_orbit_allsub,delta_orbit_allsub,color=100
legend,['allsub: 1000'],box=0,/bottom,/right,charsize=0.8


phimin = 0
phimax = 22
;=== nosub density profile ===
denmin = 1
denmax = 5000
ploterror,phi_nosub[leading_nosub],surf_den_nosub[leading_nosub],surf_den_nosub_err[leading_nosub],xr=[phimin,phimax],yr=[denmin,denmax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
oploterror,phi_nosub[trailing_nosub],surf_den_nosub[trailing_nosub],surf_den_nosub_err[trailing_nosub],linestyle=2
legend,['leading','trailing'],linestyle=[0,2],box=0,charsize=0.7

;=== 1Rs density profile ===
ploterror,phi_1Rs[leading_1Rs],surf_den_1Rs[leading_1Rs],surf_den_1Rs_err[leading_1Rs],xr=[phimin,phimax],yr=[denmin,denmax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
help,trailing_1Rs,phi_1Rs
oploterror,phi_1Rs[trailing_1Rs],surf_den_1Rs[trailing_1Rs],surf_den_1Rs_err[trailing_1Rs],linestyle=2

;=== 5Rs density profile ===
ploterror,phi_5Rs[leading_5Rs],surf_den_5Rs[leading_5Rs],surf_den_5Rs_err[leading_5Rs],xr=[phimin,phimax],yr=[denmin,denmax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
oploterror,phi_5Rs[trailing_5Rs],surf_den_5Rs[trailing_5Rs],surf_den_5Rs_err[trailing_5Rs],linestyle=2

;;=== 10Rs density profile ===
;ploterror,phi_10Rs[leading_10Rs],surf_den_10Rs[leading_10Rs],surf_den_10Rs_err[leading_10Rs],xr=[phimin,phimax],yr=[0,denmax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]')
;oploterror,phi_10Rs[trailing_10Rs],surf_den_10Rs[trailing_10Rs],surf_den_10Rs_err[trailing_10Rs],linestyle=2

;=== ALL density profile ===
ploterror,phi_allsub[leading_allsub],surf_den_allsub[leading_allsub],surf_den_allsub_err[leading_allsub],xr=[phimin,phimax],yr=[denmin,denmax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
oploterror,phi_allsub[trailing_allsub],surf_den_allsub[trailing_allsub],surf_den_allsub_err[trailing_allsub],linestyle=2


;=== nosub sigma_v profile ===
Vr_min = 10
Vr_max = 60
plot,phi_nosub[leading_nosub],sigma_Vr_nosub[leading_nosub],xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
oplot,phi_nosub[leading_nosub],sigma_Vr_nosub[leading_nosub],psym=7
oplot,phi_nosub[trailing_nosub],sigma_Vr_nosub[trailing_nosub],linestyle=2
oplot,phi_nosub[trailing_nosub],sigma_Vr_nosub[trailing_nosub],psym=7

;=== 1Rs sigma_v profile ===
plot,phi_1Rs[leading_1Rs],sigma_Vr_1Rs[leading_1Rs],xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
oplot,phi_1Rs[leading_1Rs],sigma_Vr_1Rs[leading_1Rs],psym=7
oplot,phi_1Rs[trailing_1Rs],sigma_Vr_1Rs[trailing_1Rs],linestyle=2
oplot,phi_1Rs[trailing_1Rs],sigma_Vr_1Rs[trailing_1Rs],psym=7

;=== 5Rs sigma_v profile ===
plot,phi_5Rs[leading_5Rs],sigma_Vr_5Rs[leading_5Rs],xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
oplot,phi_5Rs[leading_5Rs],sigma_Vr_5Rs[leading_5Rs],psym=7
oplot,phi_5Rs[trailing_5Rs],sigma_Vr_5Rs[trailing_5Rs],linestyle=2
oplot,phi_5Rs[trailing_5Rs],sigma_Vr_5Rs[trailing_5Rs],psym=7

;;=== 10Rs sigma_v profile ===
;plot,phi_10Rs[leading_10Rs],sigma_Vr_10Rs[leading_10Rs],xr=[phimin,phimax],yr=[0,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
;oplot,phi_10Rs[trailing_10Rs],sigma_Vr_10Rs[trailing_10Rs],linestyle=2

;=== allsub sigma_v profile ===
plot,phi_allsub[leading_allsub],sigma_Vr_allsub[leading_allsub],xr=[phimin,phimax],yr=[Vr_min,Vr_max],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
oplot,phi_allsub[leading_allsub],sigma_Vr_allsub[leading_allsub],psym=7
oplot,phi_allsub[trailing_allsub],sigma_Vr_allsub[trailing_allsub],linestyle=2
oplot,phi_allsub[trailing_allsub],sigma_Vr_allsub[trailing_allsub],psym=7

device,/close

END
