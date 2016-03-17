pro oscil

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
  dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy_1b_v50/'
  dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_subhalo_xy_1b_v50/'
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
;chr_rdtbl,'orbit_subhalo_e7.dat',0,arr,/silent
;chr_rdtbl,'orbit_subhalo_e8.dat',0,arr,/silent
chr_rdtbl,'orbit_subhalo_xy.dat',0,arr,/silent
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fname_part)

chr_rdtbl,dir_out+'part001_xy_0.5b',0,arr
t_orbit = double(arr[0,*])/1000.
x_orbit = double(arr[1,*])
y_orbit = double(arr[2,*])
z_orbit = double(arr[3,*])

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
  p1 = 99   &  p3 = 87
  p2 = 101  &  p4 = 111
endif

;device,file=file_out,/color,ysize=22,xoffset=0.2,yoffset=0.2
device,file='oscil.ps',/color,ysize=22,xoffset=0.2,yoffset=0.2
t_p2 = dblarr(9999)    & t_p3 = dblarr(9999)    & t_p4 = dblarr(9999)
r_p2 = dblarr(9999)    & r_p3 = dblarr(9999)    & r_p4 = dblarr(9999)
delv_p2 = dblarr(9999) & delv_p3 = dblarr(9999) & delv_p4 = dblarr(9999)
j = 0
tt = [0]
d1 = [0]
d2 = [0]
d3 = [0]
d4 = [0]
d5 = [0]
d6 = [0]
d7 = [0]
d_nosub2 = [0]

for i=950,N_elements(fname_part)-1,5 do begin
;for i=950,1050-1,5 do begin
  if i/100 eq i/100. then print, i
  chr_rdtbl,dir_part+fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part = reform(arr[0,0])/1000.
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

p5 = 2982
p6 = 2989
p7 = 2990
;nosub p2
sub_orbit = where(abs(t_part[0]-reform(t_orbit)) lt 0.5)
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p1],y_part[p1],x_proj=x_proj1,y_proj=y_proj1
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p2],y_part[p2],x_proj=x_proj2,y_proj=y_proj2
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p3],y_part[p3],x_proj=x_proj3,y_proj=y_proj3
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p4],y_part[p4],x_proj=x_proj4,y_proj=y_proj4
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p5],y_part[p5],x_proj=x_proj5,y_proj=y_proj5
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p6],y_part[p6],x_proj=x_proj6,y_proj=y_proj6
min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_part[p7],y_part[p7],x_proj=x_proj7,y_proj=y_proj7
tt = [tt,t_part[0]]
d1 = [d1,sqrt((x_proj1-x_part[p1])^2.+(y_proj1-y_part[p1])^2.)]
d2 = [d2,sqrt((x_proj2-x_part[p2])^2.+(y_proj2-y_part[p2])^2.)]
d3 = [d3,sqrt((x_proj3-x_part[p3])^2.+(y_proj3-y_part[p3])^2.)]
d4 = [d4,sqrt((x_proj4-x_part[p4])^2.+(y_proj4-y_part[p4])^2.)]
d5 = [d5,sqrt((x_proj5-x_part[p5])^2.+(y_proj5-y_part[p5])^2.)]
d6 = [d6,sqrt((x_proj6-x_part[p6])^2.+(y_proj6-y_part[p6])^2.)]
d7 = [d7,sqrt((x_proj7-x_part[p7])^2.+(y_proj7-y_part[p7])^2.)]

min_d_curve,x_orbit[sub_orbit],y_orbit[sub_orbit],x_nosub[p2],y_nosub[p2],x_proj=x_nosub2,y_proj=y_nosub2
d_nosub2 = [d_nosub2,sqrt((x_nosub2-x_nosub[p2])^2.+(y_nosub2-y_nosub[p2])^2.)]

endfor
tt = tt[1:*]-9
d1 = d1[1:*]
d2 = d2[1:*]
d3 = d3[1:*]
d4 = d4[1:*]
d5 = d5[1:*]
d6 = d6[1:*]
d7 = d7[1:*]
d_nosub2 = d_nosub2[1:*]

erase
!p.multi=[0,1,7]
ymax = 1;max([d1,d2,d3,d4,d6,d7])
multiplot
plot,tt,d1,xr=[0,9],yr=[0,ymax]
multiplot
plot,tt,d2,color=70,xr=[0,9],yr=[0,ymax]
oplot,tt,d_nosub2,color=100,linestyle=2
multiplot
plot,tt,d3,color=100,xr=[0,9],yr=[0,ymax]
multiplot
plot,tt,d4,ytitle=textoidl('\Deltad [kpc]'),color=150,xr=[0,9],yr=[0,ymax]
multiplot
plot,tt,d5,color=200,xr=[0,9],yr=[0,ymax]
multiplot
plot,tt,d6,color=220,xr=[0,9],yr=[0,ymax]
multiplot
plot,tt,d7,xtitle='t [Gyr]',color=255,xr=[0,9],yr=[0,ymax]
multiplot,/reset
erase

multiplot
plot,tt,d1/max(d1),xtitle='t [Gyr]',ytitle=textoidl('\Deltad [kpc]'),xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d2/max(d2),color=70,xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d3/max(d3),color=100,xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d4/max(d4),color=150,xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d5/max(d5),color=200,xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d6/max(d6),color=220,xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d7/max(d7),color=255,xr=[0,9],yr=[0,1.1]
multiplot,/reset
erase

nor = d7/max(d7)
multiplot
plot,tt,d1/max(d1)/nor,xtitle='t [Gyr]',ytitle=textoidl('\Deltad [kpc]'),xr=[0,9],yr=[0,1.1]
multiplot
plot,tt,d2/max(d2)/nor,color=70,xr=[0,9],yr=[0,2.1]
multiplot
plot,tt,d3/max(d3)/nor,color=100,xr=[0,9],yr=[0,2.1]
multiplot
plot,tt,d4/max(d4)/nor,color=150,xr=[0,9],yr=[0,2.1]
multiplot
plot,tt,d5/max(d5)/nor,color=200,xr=[0,9],yr=[0,2.1]
multiplot
plot,tt,d6/max(d6)/nor,color=220,xr=[0,9],yr=[0,2.1]
multiplot
plot,tt,d7/max(d7)/nor,color=255,xr=[0,9],yr=[0,2.1]
multiplot,/reset

device,/close

END
