pro one_encounter_snap_movie

run = 7
dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/'
dir_real = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_nosub/'
if run eq 1 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy/'
endif else if run eq 2 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy2/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy2/'
endif else if run eq 3 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_0.5b/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_0.5b/'
endif else if run eq 4 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_1b/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_1b/'
endif else if run eq 5 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_kink/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_kink/'
endif else if run eq 6 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_kink_e8/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_kink_e8/'
endif else if run eq 7 then begin
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_0.5b_v100/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_0.5b_v100/'
endif



G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
r_tide = 0.108653
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_real+'part_peri_nosub',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
V_peri = sqrt(Vx_peri^2.+Vy_peri^2.+Vz_peri^2.)
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)

;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part+'snap00000',0,arr,/silent
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (r_tide/r_peri * J_initial)
qmin = min(q_initial)
qmax = max(q_initial)
dJmin = min(scale_dJ_initial)
dJmax = max(scale_dJ_initial)
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)
qmin = qmin-0.2
qmax = qmax+0.2
dJmin = dJmin-0.05
dJmax = dJmax+0.05

;=== read out simulation info ===
chr_rdtbl,'orbit_subhalo_e7.dat',0,arr,/silent
;chr_rdtbl,'orbit_subhalo_e8.dat',0,arr,/silent
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fname_part)


spawn,'ls '+dir_part+' > temp'
chr_rdtbl,'temp',0,fname_part,/silent
spawn,'ls '+dir_subhalo+' > temp'
chr_rdtbl,'temp',0,fname_subhalo,/silent
spawn,'rm -f temp'


set_plot,'ps'
@plot_setting
;p.charsize=2.0
!p.multi=[0,2,3]
if run eq 1 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif else if run eq 2 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie2.ps'
  p1 = 2989  &  p3 = 2936
  p2 = 2935  &  p4 = 2937
endif else if run eq 3 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie_0.5b.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif else if run eq 4 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie_1b.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif else if run eq 5 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie_kink.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif else if run eq 6 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie_kink_e8.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif else if run eq 7 then begin
  file_out = dir_out+'one_encounter_xy_snap_movie_0.5b_v100.ps'
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif

device,file=file_out,/color,ysize=22,xoffset=0.2,yoffset=0.2
t_p2 = dblarr(9999)    & t_p3 = dblarr(9999)    & t_p4 = dblarr(9999)
r_p2 = dblarr(9999)    & r_p3 = dblarr(9999)    & r_p4 = dblarr(9999)
delv_p2 = dblarr(9999) & delv_p3 = dblarr(9999) & delv_p4 = dblarr(9999)
j = 0

for i=950,N_elements(fname_part)-1,10 do begin
;for i=950,1050-1,5 do begin
  if i/100 eq i/100. then print, i
  chr_rdtbl,dir_part+fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part = reform(arr[0,*])/1000.
  t_title = strtrim(string(reform(arr[0,0])/1000.-10.,f='(f8.3)'),2)+'Gyr'
  x_part = reform(arr[1,*])
  y_part = reform(arr[2,*])
  z_part = reform(arr[3,*])
  r_part = sqrt(x_part^2.+y_part^2.+z_part^2.)
  Vx_part = reform(arr[4,*])
  Vy_part = reform(arr[5,*])
  Vz_part = reform(arr[6,*])
  V_part = sqrt(Vx_part^2.+Vy_part^2.+Vz_part^2.)
  dE_part = reform(arr[7,*])
  Etot_part = reform(arr[8,*])
  J_part = reform(arr[9,*])
  delE = Etot_part-Etot_part[0]
  q_part = delE/epsilon
  dJ_part = J_part-J_part[0]
  scale_dJ_part = dJ_part / (r_tide/r_peri * J_part)

;=== nosub ===  
  chr_rdtbl,dir_nosub+fname_part[i],0,arr,/silent
  arr = double(arr)
  x_nosub = reform(arr[1,*])
  y_nosub = reform(arr[2,*])
  z_nosub = reform(arr[3,*])
  Vx_nosub = reform(arr[4,*])
  Vy_nosub = reform(arr[5,*])
  Vz_nosub = reform(arr[6,*])
  V_nosub = sqrt(Vx_nosub^2.+Vy_nosub^2.+Vz_nosub^2.)

  ; read out data of the subhalo 
  chr_rdtbl,dir_subhalo+fname_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_subhalo = reform(arr[1,*])
  y_subhalo = reform(arr[2,*])
  z_subhalo = reform(arr[3,*])
  Vx_subhalo = reform(arr[4,*])
  Vy_subhalo = reform(arr[5,*])
  Vz_subhalo = reform(arr[6,*])
  
  ; make a circle of a subhalo 
  Rs_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rs[0]))
  Rtidal_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rtidal[0]))
  Rs_yz_circle = circle(y_subhalo[0],z_subhalo[0],double(Rs[0]))
  Rtidal_yz_circle = circle(y_subhalo[0],z_subhalo[0],double(Rtidal[0]))
  Rs_xz_circle = circle(x_subhalo[0],z_subhalo[0],double(Rs[0]))
  Rtidal_xz_circle = circle(x_subhalo[0],z_subhalo[0],double(Rtidal[0]))
 
;============
;=== plot ===
;============
;=== 20Myr earlier than the closest approach ===
  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]],[y_part[k]],psym=3,color=qcolor[k]
  oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
  oplot,Rtidal_xy_circle[0,*],Rtidal_xy_circle[1,*],linestyle=2
  legend,['dt='+t_title],box=0,charsize=0.8,/bottom,/right

  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='z',/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]],[z_part[k]],psym=3,color=qcolor[k]
  oplot,Rs_xz_circle[0,*],Rs_xz_circle[1,*]
  oplot,Rtidal_xz_circle[0,*],Rtidal_xz_circle[1,*],linestyle=2

  x_cen = x_nosub[p1]
  y_cen = y_nosub[p1]
  z_cen = z_nosub[p1]
  xmin = -5.5  &  xmax = 5.5
  ymin = -5.5  &  ymax = 5.5
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]-x_cen],[y_part[k]-y_cen],psym=3,color=qcolor[k]
  oplot,Rs_xy_circle[0,*]-x_cen,Rs_xy_circle[1,*]-y_cen
  oplot,Rtidal_xy_circle[0,*]-x_cen,Rtidal_xy_circle[1,*]-y_cen,linestyle=2
  oplot,[x_part[p2]-x_cen],[y_part[p2]-y_cen],psym=4
  oplot,[x_part[p3]-x_cen],[y_part[p3]-y_cen],psym=5
  oplot,[x_part[p4]-x_cen],[y_part[p4]-y_cen],psym=6
  vline,0,linestyle=1
  hline,0,linestyle=1
  xyouts,[x_part[p2]-x_cen],[y_part[p2]-y_cen]+0.1,strtrim(string(V_part[p2],f='(f8.2)'),2),charsize=0.8,alignment=.5
  xyouts,[x_part[p3]-x_cen],[y_part[p3]-y_cen]+0.1,strtrim(string(V_part[p3],f='(f8.2)'),2),charsize=0.8,alignment=1
  xyouts,[x_part[p4]-x_cen],[y_part[p4]-y_cen]-0.1,strtrim(string(V_part[p4],f='(f8.2)'),2),charsize=0.8
  legend,[textoidl('\Deltav_{p2}=')+strtrim(string(V_part[p2]-V_nosub[p2],f='(f8.2)'),2), $
          textoidl('\Deltav_{p3}=')+strtrim(string(V_part[p3]-V_nosub[p3],f='(f8.2)'),2), $
          textoidl('\Deltav_{p4}=')+strtrim(string(V_part[p4]-V_nosub[p4],f='(f8.2)'),2)],box=0,charsize=0.8,psym=[4,5,6]

  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_nosub[k]-x_cen],[y_nosub[k]-y_cen],psym=3,color=qcolor[k]
  legend,['nosub'],box=0
  oplot,[x_nosub[p2]-x_cen],[y_nosub[p2]-y_cen],psym=4
  oplot,[x_nosub[p3]-x_cen],[y_nosub[p3]-y_cen],psym=5
  oplot,[x_nosub[p4]-x_cen],[y_nosub[p4]-y_cen],psym=6
  vline,0,linestyle=1
  hline,0,linestyle=1
  xyouts,[x_nosub[p2]-x_cen],[y_nosub[p2]-y_cen]+0.1,strtrim(string(V_nosub[p2],f='(f8.2)'),2),charsize=0.8,alignment=.5
  xyouts,[x_nosub[p3]-x_cen],[y_nosub[p3]-y_cen]+0.1,strtrim(string(V_nosub[p3],f='(f8.2)'),2),charsize=0.8,alignment=1
  xyouts,[x_nosub[p4]-x_cen],[y_nosub[p4]-y_cen]-0.1,strtrim(string(V_nosub[p4],f='(f8.2)'),2),charsize=0.8


  ; energy & angular momentum 
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part[k]],[scale_dJ_part[k]],psym=3,color=qcolor[k]
  oplot,[q_part[p2]],[scale_dJ_part[p2]],psym=4
  oplot,[q_part[p3]],[scale_dJ_part[p3]],psym=5
  oplot,[q_part[p4]],[scale_dJ_part[p4]],psym=6
  
  ; q histogram
  plothist,q_part,xr=[qmin,qmax],bin=0.1,yr=[0,210],xtitle='q'

  t_p2[j] = t_part[p2]
  r_p2[j] = r_part[p2]
  delv_p2[j] = (V_part-V_nosub)[p2]
  t_p3[j] = t_part[p3]
  r_p3[j] = r_part[p3]
  delv_p3[j] = (V_part-V_nosub)[p3]
  t_p4[j] = t_part[p4]
  r_p4[j] = r_part[p4]
  delv_p4[j] = (V_part-V_nosub)[p4]
  j += 1
endfor
erase
t_p2 = t_p2[0:j-1]
r_p2 = r_p2[0:j-1]
delv_p2 = delv_p2[0:j-1]
t_p3 = t_p3[0:j-1]
r_p3 = r_p3[0:j-1]
delv_p3 = delv_p3[0:j-1]
t_p4 = t_p4[0:j-1]
r_p4 = r_p4[0:j-1]
delv_p4 = delv_p4[0:j-1]

!p.multi=[0,1,3]
multiplot
plot,t_p2,r_p2/max(r_p2),xr=[9.8,18],ytitle='r [kpc]'
oplot,t_p2,delv_p2/max(delv_p2)/2.+0.5,color=255
hline,0.5,linestyle=1
multiplot
plot,t_p3,r_p3/max(r_p3),xr=[9.8,18],ytitle='r [kpc]'
oplot,t_p3,delv_p3/max(delv_p3)/2.+0.5,color=255
hline,0.5,linestyle=1
multiplot
plot,t_p4,r_p4/max(r_p4),xr=[9.8,18],xtitle='t [Gyr]',ytitle='r [kpc]'
oplot,t_p4,delv_p4/max(delv_p4)/2.+0.5,color=255
hline,0.5,linestyle=1
erase

!p.multi=[0,1,3]
multiplot
plot,t_p2,r_p2/max(r_p2),xr=[9.8,10.5],yr=[0.4,0.6],ytitle='r [kpc]'
oplot,t_p2,delv_p2/max(delv_p2)/2.+0.5,color=255
hline,0.5,linestyle=1
vline,10,linestyle=2
multiplot
plot,t_p3,r_p3/max(r_p3),xr=[9.8,10.5],yr=[0.4,0.6],ytitle='r [kpc]'
oplot,t_p3,delv_p3/max(delv_p3)/2.+0.5,color=255
hline,0.5,linestyle=1
vline,10,linestyle=2
multiplot
plot,t_p4,r_p4/max(r_p4),xr=[9.8,10.5],yr=[0.4,0.6],xtitle='t [Gyr]',ytitle='r [kpc]'
oplot,t_p4,delv_p4/max(delv_p4)/2.+0.5,color=255
hline,0.5,linestyle=1
vline,10,linestyle=2
multiplot,/reset


device,/close
;spawn,'convert '+dir_out+'one_encounter_xy_snap_movie.ps '+dir_out+'one_encounter_xy_snap_movie.gif'

END
