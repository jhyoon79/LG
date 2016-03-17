pro compare_closesubs

mass_name = 'e7'
dir_close = '/scratch/jhyoon/Research/LG/close_encounter/'
dir_1Rs_part = dir_close+'encounter_1Rs_'+mass_name+'/snapshot/'
dir_1Rs_subhalo = dir_close+'encounter_1Rs_'+mass_name+'/snapshot_subhalo/'
dir_5Rs_part = dir_close+'encounter_5Rs_'+mass_name+'/snapshot/'
dir_5Rs_subhalo = dir_close+'encounter_5Rs_'+mass_name+'/snapshot_subhalo/'
dir_10Rs_part = dir_close+'encounter_10Rs_'+mass_name+'/snapshot/'
dir_10Rs_subhalo = dir_close+'encounter_10Rs_'+mass_name+'/snapshot_subhalo/'

dir_1000sub = '/scratch/jhyoon/Research/LG/BI5_test/'
dir_1000sub_part = dir_1000sub+'snapshot_'+mass_name+'_1000sub/'
dir_1000sub_subhalo = dir_1000sub+'snapshot_subhalo_'+mass_name+'_1000sub/'
dir_nosub = '/scratch/jhyoon/Research/LG/realistic/snapshot_nosub/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
chr_rdtbl, dir_1000sub+'part_peri_'+mass_name+'_1000sub',0,arr,/silent
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
r_tide = 0.108653
Msat = 20000.
epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)

chr_rdtbl, dir_1Rs_part+'snap00000',0, arr,/silent
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

spawn, 'ls '+dir_1Rs_part+' > temp9'
chr_rdtbl,'temp9',0,fname_1Rs_part,/silent
fname_1Rs_part = dir_1Rs_part+reform(fname_1Rs_part)

spawn, 'ls '+dir_1Rs_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_1Rs_subhalo,/silent
fname_1Rs_subhalo = dir_1Rs_subhalo+reform(fname_1Rs_subhalo)

spawn, 'ls '+dir_5Rs_part+' > temp9'
chr_rdtbl,'temp9',0,fname_5Rs_part,/silent
fname_5Rs_part = dir_5Rs_part+reform(fname_5Rs_part)

spawn, 'ls '+dir_5Rs_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_5Rs_subhalo,/silent
fname_5Rs_subhalo = dir_5Rs_subhalo+reform(fname_5Rs_subhalo)

spawn, 'ls '+dir_10Rs_part+' > temp9'
chr_rdtbl,'temp9',0,fname_10Rs_part,/silent
fname_10Rs_part = dir_10Rs_part+reform(fname_10Rs_part)

spawn, 'ls '+dir_10Rs_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_10Rs_subhalo,/silent
fname_10Rs_subhalo = dir_10Rs_subhalo+reform(fname_10Rs_subhalo)

spawn, 'ls '+dir_1000sub_part+' > temp9'
chr_rdtbl,'temp9',0,fname_1000sub_part,/silent
fname_1000sub_part = dir_1000sub_part+reform(fname_1000sub_part)

spawn, 'ls '+dir_1000sub_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_1000sub_subhalo,/silent
fname_1000sub_subhalo = dir_1000sub_subhalo+reform(fname_1000sub_subhalo)
spawn, 'rm -f temp9'


set_plot,'ps'
@plot_setting
;!p.charsize=1
!p.multi=[0,4,3]
out_name = dir_close+'compare_closesubs_'+mass_name+'_last'
device,file=out_name+'.ps',/color,/landscape,ysize=20;,xoffset=0.5,yoffset=0.5

N_time = 1277;N_elements(fname_1Rs_part)
for i=N_time-1,N_time-1 do begin
;for i=0L,N_time-1,1 do begin

;=== read out data of the tail particles ===
  chr_rdtbl,fname_1Rs_part[i],0,arr,/silent
  arr = double(arr)
  t_1Rs_part = reform(arr[0,*])
  x_1Rs_part = reform(arr[1,*])
  y_1Rs_part = reform(arr[2,*])
  z_1Rs_part = reform(arr[3,*])
  Vx_1Rs_part = reform(arr[4,*])
  Vy_1Rs_part = reform(arr[5,*])
  Vz_1Rs_part = reform(arr[6,*])
  V_1Rs_part = sqrt(Vx_1Rs_part^2.+Vy_1Rs_part^2.+Vz_1Rs_part^2.)
  dE_1Rs_part = reform(arr[7,*])
  Etot_1Rs_part = reform(arr[8,*])
  J_1Rs_part = reform(arr[9,*])
  delE = Etot_1Rs_part-Etot_1Rs_part[0]
  q_1Rs_part = delE/epsilon
  dJ_1Rs_part = J_1Rs_part-J_1Rs_part[0]
  scale_dJ_1Rs_part = dJ_1Rs_part / (r_tide/r_peri * J_1Rs_part)

;=== read out data of the tail particles ===
  chr_rdtbl,fname_5Rs_part[i],0,arr,/silent
  arr = double(arr)
  t_5Rs_part = reform(arr[0,*])
  x_5Rs_part = reform(arr[1,*])
  y_5Rs_part = reform(arr[2,*])
  z_5Rs_part = reform(arr[3,*])
  Vx_5Rs_part = reform(arr[4,*])
  Vy_5Rs_part = reform(arr[5,*])
  Vz_5Rs_part = reform(arr[6,*])
  V_5Rs_part = sqrt(Vx_5Rs_part^2.+Vy_5Rs_part^2.+Vz_5Rs_part^2.)
  dE_5Rs_part = reform(arr[7,*])
  Etot_5Rs_part = reform(arr[8,*])
  J_5Rs_part = reform(arr[9,*])
  delE = Etot_5Rs_part-Etot_5Rs_part[0]
  q_5Rs_part = delE/epsilon
  dJ_5Rs_part = J_5Rs_part-J_5Rs_part[0]
  scale_dJ_5Rs_part = dJ_5Rs_part / (r_tide/r_peri * J_5Rs_part)

;=== read out data of the tail particles ===
  chr_rdtbl,fname_10Rs_part[i],0,arr,/silent
  arr = double(arr)
  t_10Rs_part = reform(arr[0,*])
  x_10Rs_part = reform(arr[1,*])
  y_10Rs_part = reform(arr[2,*])
  z_10Rs_part = reform(arr[3,*])
  Vx_10Rs_part = reform(arr[4,*])
  Vy_10Rs_part = reform(arr[5,*])
  Vz_10Rs_part = reform(arr[6,*])
  V_10Rs_part = sqrt(Vx_10Rs_part^2.+Vy_10Rs_part^2.+Vz_10Rs_part^2.)
  dE_10Rs_part = reform(arr[7,*])
  Etot_10Rs_part = reform(arr[8,*])
  J_10Rs_part = reform(arr[9,*])
  delE = Etot_10Rs_part-Etot_10Rs_part[0]
  q_10Rs_part = delE/epsilon
  dJ_10Rs_part = J_10Rs_part-J_10Rs_part[0]
  scale_dJ_10Rs_part = dJ_10Rs_part / (r_tide/r_peri * J_10Rs_part)

;=== read out data of the tail particles ===
  chr_rdtbl,fname_1000sub_part[i],0,arr,/silent
  arr = double(arr)
  t_1000sub_part = reform(arr[0,*])
  x_1000sub_part = reform(arr[1,*])
  y_1000sub_part = reform(arr[2,*])
  z_1000sub_part = reform(arr[3,*])
  Vx_1000sub_part = reform(arr[4,*])
  Vy_1000sub_part = reform(arr[5,*])
  Vz_1000sub_part = reform(arr[6,*])
  dE_1000sub_part = reform(arr[7,*])
  Etot_1000sub_part = reform(arr[8,*])
  J_1000sub_part = reform(arr[9,*])
  delE_1000sub = Etot_1000sub_part-Etot_1000sub_part[0]
  q_1000sub_part = delE_1000sub/epsilon
  dJ_1000sub_part = J_1000sub_part-J_1000sub_part[0]
  scale_dJ_1000sub_part = dJ_1000sub_part / (r_tide/r_peri * J_1000sub_part)

  t_step = string(t_1Rs_part[0]/1000.,f='(f6.3)')
  if (i/10. eq i/10) then print,t_step+'Gyr'

;=== read out data of the subhalo ===
  chr_rdtbl,fname_1Rs_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_1Rs_subhalo = reform(arr[1,*])
  y_1Rs_subhalo = reform(arr[2,*])
  z_1Rs_subhalo = reform(arr[3,*])
  Vx_1Rs_subhalo = reform(arr[4,*])
  Vy_1Rs_subhalo = reform(arr[5,*])
  Vz_1Rs_subhalo = reform(arr[6,*])
  N_1Rs_subhalo = strtrim(N_elements(x_1Rs_subhalo),2)

;=== read out data of the subhalo ===
  chr_rdtbl,fname_5Rs_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_5Rs_subhalo = reform(arr[1,*])
  y_5Rs_subhalo = reform(arr[2,*])
  z_5Rs_subhalo = reform(arr[3,*])
  Vx_5Rs_subhalo = reform(arr[4,*])
  Vy_5Rs_subhalo = reform(arr[5,*])
  Vz_5Rs_subhalo = reform(arr[6,*])
  N_5Rs_subhalo = strtrim(N_elements(x_5Rs_subhalo),2)

;=== read out data of the subhalo ===
  chr_rdtbl,fname_10Rs_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_10Rs_subhalo = reform(arr[1,*])
  y_10Rs_subhalo = reform(arr[2,*])
  z_10Rs_subhalo = reform(arr[3,*])
  Vx_10Rs_subhalo = reform(arr[4,*])
  Vy_10Rs_subhalo = reform(arr[5,*])
  Vz_10Rs_subhalo = reform(arr[6,*])
  N_10Rs_subhalo = strtrim(N_elements(x_10Rs_subhalo),2)

;=== read out data of the subhalo ===
  chr_rdtbl,fname_1000sub_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_1000sub_subhalo = reform(arr[1,*])
  y_1000sub_subhalo = reform(arr[2,*])
  z_1000sub_subhalo = reform(arr[3,*])
  Vx_1000sub_subhalo = reform(arr[4,*])
  Vy_1000sub_subhalo = reform(arr[5,*])
  Vz_1000sub_subhalo = reform(arr[6,*])

;=== read out data of the tail without subhalo ===
print,dir_nosub+'snap'+strtrim(fix(t_1Rs_part[0]),2)
  chr_rdtbl,dir_nosub+'snap'+strtrim(fix(t_1Rs_part[0]),2),0,arr,/silent
  arr = double(arr)
  x_nosub = reform(arr[1,*])
  y_nosub = reform(arr[2,*])
  z_nosub = reform(arr[3,*])
  Vx_nosub = reform(arr[4,*])
  Vy_nosub = reform(arr[5,*])
  Vz_nosub = reform(arr[6,*])
  Etot_nosub = reform(arr[8,*])
  i_name = '_'+string(i,f='(i3.3)')

;=== estimate b,delv,delE ===
  b_impact = 999999.d
  for j=0L,N_particle-1 do begin
    d2subhalo = sqrt((x_1Rs_part[j]-x_1Rs_subhalo)^2.+(y_1Rs_part[j]-y_1Rs_subhalo)^2.+(z_1Rs_part[j]-z_1Rs_subhalo)^2.)
    tmp_min = min(d2subhalo,tmp_sub)
    if (tmp_min lt b_impact) then begin
      d2subhalo_close_pts = d2subhalo
      b_impact = float(tmp_min)
      min_sub = tmp_sub
      V_viewed_at_tail = sqrt((Vx_1Rs_subhalo[min_sub]-Vx_1Rs_part[j])^2.+$
             (Vy_1Rs_subhalo[min_sub]-Vy_1Rs_part[j])^2.+(Vz_1Rs_subhalo[min_sub]-Vz_1Rs_part[j])^2.)
      Vx_viewed_at_tail = Vx_1Rs_subhalo[min_sub] - Vx_1Rs_part[j]
      Vy_viewed_at_tail = Vy_1Rs_subhalo[min_sub] - Vy_1Rs_part[j]
      Vz_viewed_at_tail = Vz_1Rs_subhalo[min_sub] - Vz_1Rs_part[j]
    endif
  endfor


;=== plot begins ===
;=== in the tail's rest frame (xy plane ) 4 subhalos ===
  xmin = -5.99  &  xmax = 5.99
  ymin = -5.99  &  ymax = 5.99

;=== in the tail's rest frame (xy plane ) 1000 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_1000sub_part[k]-x_1000sub_part[0]],[y_1000sub_part[k]-y_1000sub_part[0]],psym=3,color=qcolor[k]
  legend,['ALL: 1000'],box=0,/bottom,/right,charsize=1.2
  legend,[t_step+'Gyr'],box=0,charsize=1.2
for j=0,23 do begin
  xl = ['','-4','-2','0','2','4','']
  xyouts,-6+j*2,-7.5,xl[round((j/6.-j/6)*6)],alignment=0.5
endfor
for j=0,3 do xyouts,j*12,-9.,textoidl('\Deltax'),alignment=0.5

  multiplot
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  for k=0,N_particle-1 do oplot,[x_1Rs_part[k]-x_1Rs_part[0]],[y_1Rs_part[k]-y_1Rs_part[0]],psym=3,color=qcolor[k]
  legend,['1Rs: '+N_1Rs_subhalo],box=0,/bottom,/right,charsize=1.2

;  for l=0,N_elements(x_1Rs_subhalo)-1 do begin
;    Rs_xy_circle = circle(x_1Rs_subhalo[l],y_1Rs_subhalo[l],double(Rs[l]))
;    Rtidal_xy_circle = circle(x_1Rs_subhalo[l],y_1Rs_subhalo[l],double(Rtidal[l]))
;    oplot,Rs_xy_circle[0,*]-x_1Rs_part[0],Rs_xy_circle[1,*]-y_1Rs_part[0]
;  endfor
;=== in the tail's rest frame (xy plane ) 10 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  for k=0,N_particle-1 do oplot,[x_5Rs_part[k]-x_5Rs_part[0]],[y_5Rs_part[k]-y_5Rs_part[0]],psym=3,color=qcolor[k]
  legend,['5Rs: '+N_5Rs_subhalo],box=0,/bottom,/right,charsize=1.2

;=== in the tail's rest frame (xy plane ) 10 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  for k=0,N_particle-1 do oplot,[x_10Rs_part[k]-x_10Rs_part[0]],[y_10Rs_part[k]-y_10Rs_part[0]],psym=3,color=qcolor[k]
  legend,['10Rs: '+N_10Rs_subhalo],box=0,/bottom,/right,charsize=1.2


  if (mass_name eq 'e6') then begin
    qmin = 7  &  qmax = -7
    jmin = -2   &  jmax = 2
    dEmin = -30  &  dEmax = 70
  endif else if (mass_name eq 'e7') then begin
    qmin = -2.9  &  qmax = 3.4
    jmin = -0.5   &  jmax = 2.4
    dEmin = -3  &  dEmax = 3
;    dEmin = -380  &  dEmax = 380
  endif else if (mass_name eq 'e8') then begin
    qmin = 22   &  qmax = -16
    jmin = -7.5 &  jmax = 3
    dEmin = -800  &  dEmax = 2600
  endif

;=== energy gained 4 subhalos ===
  Ediff = (Etot_1000sub_part-Etot_nosub)/epsilon
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],ytitle=textoidl('\Deltaq'),/isotropic,pos=[0.138,0.35,0.343,0.70]
  for k=0,N_particle-1 do oplot,[q_1000sub_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy gained 10 subhalos ===
  Ediff = (Etot_1Rs_part-Etot_nosub)/epsilon
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],/isotropic,pos=[0.343,0.35,0.548,0.70]
  for k=0,N_particle-1 do oplot,[q_1Rs_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy gained 10 subhalos ===
  Ediff = (Etot_5Rs_part-Etot_nosub)/epsilon
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],/isotropic,pos=[0.548,0.35,0.753,0.70]
  for k=0,N_particle-1 do oplot,[q_5Rs_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1

;=== energy gained 1000 subhalos ===
  Ediff = (Etot_10Rs_part-Etot_nosub)/epsilon
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],/isotropic,pos=[0.753,0.35,0.958,0.70]
  for k=0,N_particle-1 do oplot,[q_10Rs_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1


;=== energy & angular momentum 4 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ'),pos=[0.138,0.105,0.343,0.35]
  for k=0,N_particle-1 do oplot,[q_1000sub_part[k]],[scale_dJ_1000sub_part[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

;=== energy & angular momentum 10 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',pos=[0.343,0.105,0.548,0.35]
  for k=0,N_particle-1 do oplot,[q_1Rs_part[k]],[scale_dJ_1Rs_part[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

;=== energy & angular momentum 10 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',pos=[0.548,0.105,0.753,0.35]
  for k=0,N_particle-1 do oplot,[q_5Rs_part[k]],[scale_dJ_5Rs_part[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

;=== energy & angular momentum : 1000 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',pos=[0.753,0.105,0.958,0.35]
  for k=0,N_particle-1 do oplot,[q_10Rs_part[k]],[scale_dJ_10Rs_part[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1

  multiplot,/reset
  erase

endfor

device,/close
close,5

END
