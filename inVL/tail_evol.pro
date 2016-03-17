pro tail_evol

dir_nosub12916 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_12916/snapshot_nosub/'
dir_nosub13000 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13000/snapshot_nosub/'
dir_nosub13011 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13011/snapshot_nosub/'
dir_inVLallsub = '/scratch/jhyoon/Research/LG/inVL_allsub/'
dir_13000 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13000/'
dir_part13000 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13000/snapshot/'
dir_subhalo13000 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13000/snapshot_subhalo/'
dir_12916 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_12916/'
dir_part12916 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_12916/snapshot/'
dir_subhalo12916 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_12916/snapshot_subhalo/'
dir_13011 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13011/'
dir_part13011 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13011/snapshot/'
dir_subhalo13011 = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_13011/snapshot_subhalo/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
r_tide = 0.108653
Msat = 10000.

;=== perigee ===
chr_rdtbl, dir_12916+'part_peri',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
r_peri12916 = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
p = r_peri12916/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
epsilon_12916 = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri12916)

chr_rdtbl, dir_13000+'part_peri',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
r_peri13000 = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
p = r_peri13000/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
epsilon_13000 = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri13000)

chr_rdtbl, dir_13011+'part_peri',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
r_peri13011 = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
p = r_peri13011/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
epsilon_13011 = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri13011)


chr_rdtbl,dir_part13000+'snap00000',0,arr,/silent
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon_13000
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (r_tide/r_peri13000 * J_initial)
qmin = min(q_initial)
qmax = max(q_initial)
dJmin = min(scale_dJ_initial)
dJmax = max(scale_dJ_initial)
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)


spawn, 'ls '+dir_part12916+' > temp9'
chr_rdtbl,'temp9',0,fname_part12916,/silent
fname_part12916 = dir_part12916+reform(fname_part12916)
spawn, 'ls '+dir_subhalo12916+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo12916,/silent
fname_subhalo12916 = dir_subhalo12916+reform(fname_subhalo12916)

spawn, 'ls '+dir_part13000+' > temp9'
chr_rdtbl,'temp9',0,fname_part13000,/silent
fname_part13000 = dir_part13000+reform(fname_part13000)
spawn, 'ls '+dir_subhalo13000+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo13000,/silent
fname_subhalo13000 = dir_subhalo13000+reform(fname_subhalo13000)

spawn, 'ls '+dir_part13011+' > temp9'
chr_rdtbl,'temp9',0,fname_part13011,/silent
fname_part13011 = dir_part13011+reform(fname_part13011)
spawn, 'ls '+dir_subhalo13011+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo13011,/silent
fname_subhalo13011 = dir_subhalo13011+reform(fname_subhalo13011)

spawn, 'ls '+dir_nosub12916+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub12916,/silent
fname_nosub12916 = dir_nosub12916+reform(fname_nosub12916)
spawn, 'ls '+dir_nosub13000+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub13000,/silent
fname_nosub13000 = dir_nosub13000+reform(fname_nosub13000)
spawn, 'ls '+dir_nosub13011+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub13011,/silent
fname_nosub13011 = dir_nosub13011+reform(fname_nosub13011)
spawn,'rm -f temp9'

chr_rdtbl,'../frog_VLsubhalo_all.dat',0,arr,/silent
arr = double(arr)
Mtidal = reform(arr[6,*])
Rs_subhalo = reform(arr[7,*])
Rtidal = reform(arr[8,*])

N_time = N_elements(fname_part13011)


@plot_setting
!p.charsize=1.1
!p.multi=[0,3,3]
out_name = dir_inVLallsub+'tail_evol.ps'
device,file=out_name,/color,xoffset=0.2,yoffset=0.2,ysize=18

for i=0L,N_time-1,2 do begin
;=== read out data of the tail particles ===
  chr_rdtbl,fname_part12916[i],0,arr,/silent
  arr = double(arr)
  t_part12916 = reform(arr[0,*])
  x_part12916 = reform(arr[1,*])
  y_part12916 = reform(arr[2,*])
  z_part12916 = reform(arr[3,*])
  Vx_part12916 = reform(arr[4,*])
  Vy_part12916 = reform(arr[5,*])
  Vz_part12916 = reform(arr[6,*])
  V_part12916 = sqrt(Vx_part12916^2.+Vy_part12916^2.+Vz_part12916^2.)
  dE_part12916 = reform(arr[7,*])
  Etot_part12916 = reform(arr[8,*])
  J_part12916 = reform(arr[9,*])
  delE = Etot_part12916-Etot_part12916[0]
  q_part12916 = delE/epsilon_12916
  dJ_part12916 = J_part12916-J_part12916[0]
  scale_dJ_part12916 = dJ_part12916 / (r_tide/r_peri12916 * J_part12916)
  t_step12916 = string(t_part12916[0]/1000.,f='(f6.3)')
  if (i/50. eq i/50) then print,t_step12916+'Gyr'

  chr_rdtbl,fname_part13000[i],0,arr,/silent
  arr = double(arr)
  t_part13000 = reform(arr[0,*])
  x_part13000 = reform(arr[1,*])
  y_part13000 = reform(arr[2,*])
  z_part13000 = reform(arr[3,*])
  Vx_part13000 = reform(arr[4,*])
  Vy_part13000 = reform(arr[5,*])
  Vz_part13000 = reform(arr[6,*])
  V_part13000 = sqrt(Vx_part13000^2.+Vy_part13000^2.+Vz_part13000^2.)
  dE_part13000 = reform(arr[7,*])
  Etot_part13000 = reform(arr[8,*])
  J_part13000 = reform(arr[9,*])
  delE = Etot_part13000-Etot_part13000[0]
  q_part13000 = delE/epsilon_13000
  dJ_part13000 = J_part13000-J_part13000[0]
  scale_dJ_part13000 = dJ_part13000 / (r_tide/r_peri13000 * J_part13000)
  t_step13000 = string(t_part13000[0]/1000.,f='(f6.3)')


  chr_rdtbl,fname_part13011[i],0,arr,/silent
  arr = double(arr)
  t_part13011 = reform(arr[0,*])
  x_part13011 = reform(arr[1,*])
  y_part13011 = reform(arr[2,*])
  z_part13011 = reform(arr[3,*])
  Vx_part13011 = reform(arr[4,*])
  Vy_part13011 = reform(arr[5,*])
  Vz_part13011 = reform(arr[6,*])
  V_part13011 = sqrt(Vx_part13011^2.+Vy_part13011^2.+Vz_part13011^2.)
  dE_part13011 = reform(arr[7,*])
  Etot_part13011 = reform(arr[8,*])
  J_part13011 = reform(arr[9,*])
  delE = Etot_part13011-Etot_part13011[0]
  q_part13011 = delE/epsilon_13011
  dJ_part13011 = J_part13011-J_part13011[0]
  scale_dJ_part13011 = dJ_part13011 / (r_tide/r_peri13011 * J_part13011)
  t_step13011 = string(t_part13011[0]/1000.,f='(f6.3)')

;=== read out data of the tail without subhalo ===
  chr_rdtbl,fname_nosub12916[i],0,arr,/silent
  arr = double(arr)
;  x_nosub12916 = reform(arr[1,*])
;  y_nosub12916 = reform(arr[2,*])
;  z_nosub12916 = reform(arr[3,*])
;  Vx_nosub12916 = reform(arr[4,*])
;  Vy_nosub12916 = reform(arr[5,*])
;  Vz_nosub12916 = reform(arr[6,*])
  Etot_nosub12916 = reform(arr[8,*])

  chr_rdtbl,fname_nosub13000[i],0,arr,/silent
  arr = double(arr)
  Etot_nosub13000 = reform(arr[8,*])

  chr_rdtbl,fname_nosub13011[i],0,arr,/silent
  arr = double(arr)
  Etot_nosub13011 = reform(arr[8,*])


;=== read out subhalo data ===
  chr_rdtbl,fname_subhalo12916[i],0,arr,/silent
  arr = double(arr)
  t_subhalo12916 = reform(arr[0,*])
  x_subhalo12916 = reform(arr[1,*])
  y_subhalo12916 = reform(arr[2,*])
  z_subhalo12916 = reform(arr[3,*])

  chr_rdtbl,fname_subhalo13000[i],0,arr,/silent
  arr = double(arr)
  t_subhalo13000 = reform(arr[0,*])
  x_subhalo13000 = reform(arr[1,*])
  y_subhalo13000 = reform(arr[2,*])
  z_subhalo13000 = reform(arr[3,*])

  chr_rdtbl,fname_subhalo13011[i],0,arr,/silent
  arr = double(arr)
  t_subhalo13011 = reform(arr[0,*])
  x_subhalo13011 = reform(arr[1,*])
  y_subhalo13011 = reform(arr[2,*])
  z_subhalo13011 = reform(arr[3,*])

;=== read out data of the tail without subhalo ===
;=== estimate b,delv,delE ===
  b_impact12916 = 999999.d
  b_impact13000 = 999999.d
  b_impact13011 = 999999.d
  for j=0L,N_particle-1 do begin
    d2subhalo12916 = sqrt((x_part12916[j]-x_subhalo12916)^2.+(y_part12916[j]-y_subhalo12916)^2.+(z_part12916[j]-z_subhalo12916)^2.)
    tmp_min12916= min(d2subhalo12916,tmp_sub)
    if (tmp_min12916 lt b_impact12916) then begin
      d2subhalo_close_pts12916 = d2subhalo12916
      b_impact12916 = tmp_min12916
    endif

    d2subhalo13000 = sqrt((x_part13000[j]-x_subhalo13000)^2.+(y_part13000[j]-y_subhalo13000)^2.+(z_part13000[j]-z_subhalo13000)^2.)
    tmp_min13000 = min(d2subhalo13000,tmp_sub)
    if (tmp_min13000 lt b_impact13000) then begin
      d2subhalo_close_pts13000 = d2subhalo13000
      b_impact13000 = tmp_min13000
    endif

    d2subhalo13011 = sqrt((x_part13011[j]-x_subhalo13011)^2.+(y_part13011[j]-y_subhalo13011)^2.+(z_part13011[j]-z_subhalo13011)^2.)
    tmp_min13011 = min(d2subhalo13011,tmp_sub)
    if (tmp_min13011 lt b_impact13011) then begin
      d2subhalo_close_pts13011 = d2subhalo13011
      b_impact13011 = tmp_min13011
    endif
  endfor

  xmin = -20   &  xmax = 20
  ymin = -20   &  ymax = 20
  qmin = -16.4 &  qmax = 16.4
  jmin = -5.9  &  jmax = 7.4
  dEmin = -22  &  dEmax = 16
  Rs_cri = Rs_subhalo*10
;=============
;=== 12916 ===
;=============
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x',ytitle='y',title=t_step12916+'Gyr',/isotropic
  for k=0,N_particle-1 do oplot,[x_part12916[k]],[y_part12916[k]],psym=3,color=qcolor[k]
  sub_Rs = where(d2subhalo_close_pts12916 lt Rs_cri) 
  if (sub_Rs[0] ne -1) then begin
    for l=0L,N_elements(sub_Rs)-1 do begin
      Rs_xy_circle = circle(x_subhalo12916[sub_Rs[l]],y_subhalo12916[sub_Rs[l]],double(Rs_subhalo[sub_Rs[l]]))
      oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
    endfor
  endif

;=== energy gained by subhalos ===
  Ediff = (Etot_part12916-Etot_nosub12916)/epsilon_12916
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],xtitle='q',ytitle=textoidl('\Deltaq')
  for k=0,N_particle-1 do oplot,[q_part12916[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy & angular momentum ===
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part12916[k]],[scale_dJ_part12916[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1


;=============
;=== 13000 ===
;=============
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x',ytitle='y',/isotropic
  for k=0,N_particle-1 do oplot,[x_part13000[k]],[y_part13000[k]],psym=3,color=qcolor[k]
  sub_Rs = where(d2subhalo_close_pts13000 lt Rs_cri) 
  if (sub_Rs[0] ne -1) then begin
    for l=0L,N_elements(sub_Rs)-1 do begin
      Rs_xy_circle = circle(x_subhalo13000[sub_Rs[l]],y_subhalo13000[sub_Rs[l]],double(Rs_subhalo[sub_Rs[l]]))
      oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
    endfor
  endif

;=== energy gained by subhalos ===
  Ediff = (Etot_part13000-Etot_nosub13000)/epsilon_13000
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],xtitle='q',ytitle=textoidl('\Deltaq')
  for k=0,N_particle-1 do oplot,[q_part13000[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy & angular momentum ===
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part13000[k]],[scale_dJ_part13000[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

;=============
;=== 13011 ===
;=============
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle='x',ytitle='y',/isotropic
  for k=0,N_particle-1 do oplot,[x_part13011[k]],[y_part13011[k]],psym=3,color=qcolor[k]
  sub_Rs = where(d2subhalo_close_pts13011 lt Rs_cri) 
  if (sub_Rs[0] ne -1) then begin
    for l=0L,N_elements(sub_Rs)-1 do begin
      Rs_xy_circle = circle(x_subhalo13011[sub_Rs[l]],y_subhalo13011[sub_Rs[l]],double(Rs_subhalo[sub_Rs[l]]))
      oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
    endfor
  endif

;=== energy gained by subhalos ===
  Ediff = (Etot_part13011-Etot_nosub13011)/epsilon_13011
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],xtitle='q',ytitle=textoidl('\Deltaq')
  for k=0,N_particle-1 do oplot,[q_part13011[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy & angular momentum ===
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part13011[k]],[scale_dJ_part13011[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

endfor
erase
device,/close

END
