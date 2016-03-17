pro tail_VLallsub

dir_inVL = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_0.5Myr/'
dir_part = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_0.5Myr/snapshot/'
dir_subhalo = '/scratch/jhyoon/Research/LG/inVL_allsub/inVL_0.5Myr/snapshot_subhalo/'
dir_nosub = '/scratch/jhyoon/Research/LG/realistic/snapshot_nosub/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl, dir_inVL+'part_peri',0,arr,/silent
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
Msat = 10000.
epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)

chr_rdtbl,dir_part+'snap00000',0,arr,/silent
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

spawn, 'ls '+dir_part+' > temp9'
chr_rdtbl,'temp9',0,fname_part,/silent
fname_part = dir_part+reform(fname_part)

spawn, 'ls '+dir_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo,/silent
fname_subhalo = dir_subhalo+reform(fname_subhalo)

spawn, 'ls '+dir_nosub+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub,/silent
fname_nosub = dir_nosub+reform(fname_nosub)

chr_rdtbl,'../frog_VLsubhalo_all.dat',0,arr,/silent
arr = double(arr)
Mtidal = reform(arr[6,*])
Rs = reform(arr[7,*])
Rtidal = reform(arr[8,*])

N_time = N_elements(fname_part)

;=== find snapshot of close encounter ===
select_snap = 14
snap_step = 5
t_snap = [0.484,1.490,2.050,3.500,4.100,5.270,5.455,5.620,5.785,7.832,8.620,10.395,11.470,11.640,12.760]
help,t_snap
print,t_snap

;for kk=0,n_elements(t_snap)-1 do begin
for kk=14,14 do begin

select_snap = kk
snap_file = long(t_snap[select_snap]*1000)
snap_begin = snap_file-80
snap_end = snap_file+80
if (snap_file eq 3500) then begin
  snap_begin = snap_file-300
  snap_end = snap_file+300
  snap_step = 1
endif else if (snap_file eq 4100) then begin
  snap_begin = snap_file-200
  snap_end = snap_file+100
  snap_step = 1
endif else if (snap_file eq 7832) then begin
  snap_begin = snap_file-100
  snap_end = snap_file+100
  snap_step = 1
endif else if (snap_file eq 12760) then begin
  snap_begin = snap_file
  snap_end = snap_file
  snap_step = 1
endif

print,snap_file,snap_begin,snap_end

@plot_setting
!p.charsize=1.1
!p.multi=[0,2,2]
;out_name = dir_inVL+'tail_VLallsub'
out_name = dir_inVL+'tail_VLallsub_'+string(snap_file,f='(i5.5)')
device,file=out_name+'.ps',/color;,ysize=25,yoffset=.5,xoffset=.5

;for i=0L,N_time-1,10 do begin
for i=snap_begin,snap_end,snap_step do begin
;=== read out data of the tail particles ===
  chr_rdtbl,fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part = reform(arr[0,*])
  x_part = reform(arr[1,*])
  y_part = reform(arr[2,*])
  z_part = reform(arr[3,*])
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

  if t_part[0] gt 12760 then goto,out
  t_step = string(t_part[0]/1000.,f='(f6.3)')
  if (i/50. eq i/50) then print,t_step+'Gyr'

;=== read out data of the subhalo ===
  chr_rdtbl,fname_subhalo[i],0,arr,/silent
  arr = double(arr)
  x_subhalo = reform(arr[1,*])
  y_subhalo = reform(arr[2,*])
  z_subhalo = reform(arr[3,*])
  Vx_subhalo = reform(arr[4,*])
  Vy_subhalo = reform(arr[5,*])
  Vz_subhalo = reform(arr[6,*])
  N_subhalo = strtrim(N_elements(x_subhalo),2)

;=== read out data of the tail without subhalo ===
  chr_rdtbl,fname_nosub[i],0,arr,/silent
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
    d2subhalo = sqrt((x_part[j]-x_subhalo)^2.+(y_part[j]-y_subhalo)^2.+(z_part[j]-z_subhalo)^2.)
    tmp_min = min(d2subhalo,tmp_sub)
    if (tmp_min lt b_impact) then begin
      d2subhalo_close_pts = d2subhalo
      b_impact = float(tmp_min)
      min_sub = tmp_sub
      V_viewed_at_tail = sqrt((Vx_subhalo[min_sub]-Vx_part[j])^2.+$
             (Vy_subhalo[min_sub]-Vy_part[j])^2.+(Vz_subhalo[min_sub]-Vz_part[j])^2.)
      Vx_viewed_at_tail = Vx_subhalo[min_sub] - Vx_part[j]
      Vy_viewed_at_tail = Vy_subhalo[min_sub] - Vy_part[j]
      Vz_viewed_at_tail = Vz_subhalo[min_sub] - Vz_part[j]
    endif
  endfor


;=== plot begins ===
;=== x,y coord. of the tail ===
;  Mtidal_name = strtrim(fix(alog10(Mtidal[min_sub])),2)
;  Rs_name = strtrim(string(Rs[min_sub],f='(f8.4)'),2)
;  ;multiplot
;  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',/isotropic
;  for k=0,N_particle-1 do oplot,[x_part[k]],[y_part[k]],psym=3,color=qcolor[k]
;
;  sub_Rs = where(d2subhalo_close_pts lt 10*Rs) 
;  if (sub_Rs[0] ne -1) then begin
;    for l=0L,N_elements(x_subhalo[sub_Rs])-1 do begin
;      Rs_xy_circle = circle(x_subhalo[sub_Rs[l]],y_subhalo[sub_Rs[l]],double(Rs[sub_Rs[l]]))
;      oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
;    endfor
;  endif
;  legend,[t_step+'Gyr',textoidl('M_{tidal}=e'+Mtidal_name),$
;         'R_s='+Rs_name+'kpc',$
;         'b='+strtrim(string(b_impact,f='(f9.4)'),2)+'kpc'],box=0
;
;  ;multiplot
;  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',/isotropic
;  for k=0,N_particle-1 do oplot,[x_nosub[k]],[y_nosub[k]],psym=3,color=qcolor[k]


;=== in the tail's rest frame (xy plane) ===

  plot,[0],[0],/nodata,xr=[-7,7],yr=[-7,7],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_nosub[k]-x_nosub[0]],[y_nosub[k]-y_nosub[0]],psym=3,color=qcolor[k]
  legend,['w/o subhalos'],box=0,charsize=1,pos=[-4.5,-4.5]
  legend,[t_step+'Gyr'],box=0,charsize=1

  ;multiplot
  plot,[0],[0],/nodata,xr=[-7,7],yr=[-7,7],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]-x_part[0]],[y_part[k]-y_part[0]],psym=3,color=qcolor[k]
  legend,[t_step+'Gyr'],box=0,charsize=1
  legend,['w/ VL subhalos'],box=0,charsize=1,pos=[-5.5,-4.5]
;  if (sub_Rs[0] ne -1) then begin
;    for l=0L,N_elements(x_subhalo[sub_Rs])-1 do begin
;      Rs_xy_circle = circle(x_subhalo[sub_Rs[l]],y_subhalo[sub_Rs[l]],double(Rs[sub_Rs[l]]))
;      oplot,Rs_xy_circle[0,*]-x_part[0],Rs_xy_circle[1,*]-y_part[0]
;    endfor
;  endif
;  legend,[textoidl('V_{rel}=')+strtrim(fix(V_viewed_at_tail),2)+'km/s',$
;          textoidl('V_{x,rel}=')+strtrim(fix(Vx_viewed_at_tail),2)+'km/s',$
;          textoidl('V_{y,rel}=')+strtrim(fix(Vy_viewed_at_tail),2)+'km/s',$
;          textoidl('V_{z,rel}=')+strtrim(fix(Vz_viewed_at_tail),2)+'km/s'],box=0

  ;multiplot

;=== energy gained by subhalos ===
  qmin = -2.2  &  qmax = 2.2
  jmin = -1.3  &  jmax = 1.3
  dEmin = -9  &  dEmax = 1
;  dEmin = -1200  &  dEmax = 490

  Ediff = (Etot_part-Etot_nosub)/epsilon
  ;multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],xtitle='q',ytitle=textoidl('\Deltaq')
  for k=0,N_particle-1 do oplot,[q_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
  hline,0,linestyle=1


;=== energy & angular momentum ===
  ;multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part[k]],[scale_dJ_part[k]],psym=3,color=qcolor[k]
  vline,0,linestyle=1
  hline,0,linestyle=1
  ;multiplot,/reset
  erase
endfor
out:device,/close
close,5

;spawn, 'convert '+out_name+'.ps '+out_name+'.gif'

endfor

END
