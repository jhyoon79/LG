pro tail_VL

;dir_inVL = '/scratch/jhyoon/Research/LG/inVL_0.25Myr/'
;dir_part = '/scratch/jhyoon/Research/LG/inVL_0.25Myr/snapshot/'
;dir_subhalo = '/scratch/jhyoon/Research/LG/inVL_0.25Myr/snapshot_subhalo/'
dir_inVL = '/scratch/jhyoon/Research/LG/inVL/'
dir_part = '/scratch/jhyoon/Research/LG/inVL/snapshot/'
dir_subhalo = '/scratch/jhyoon/Research/LG/inVL/snapshot_subhalo/'
;dir_part = '/scratch/jhyoon/Research/LG/inVL/snapshot_10Myr/'
;dir_subhalo = '/scratch/jhyoon/Research/LG/inVL/snapshot_subhalo_10Myr/'

;dir_nosub = '/scratch/jhyoon/Research/LG/snapshot_nosubhalo/'
;dir_nosub = '/scratch/jhyoon/Research/LG/snapshot_nosubhalo_1Myr/'
dir_nosub = '/scratch/jhyoon/Research/LG/realistic/snapshot/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl, dir_inVL+'part_peri',0,arr
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

chr_rdtbl, dir_part+'snap00000',0, arr
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
chr_rdtbl,'temp9',0,fname_part
fname_part = dir_part+reform(fname_part)

spawn, 'ls '+dir_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo
fname_subhalo = dir_subhalo+reform(fname_subhalo)

spawn, 'ls '+dir_nosub+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub
fname_nosub = dir_nosub+reform(fname_nosub)

;readcol,'close_10Rs_VL',sub_subhalo

chr_rdtbl,'../frog_VLsubhalo.dat',0,arr
;Mtidal = double(arr[6,sub_subhalo])
;Rs = double(arr[7,sub_subhalo])
;Rtidal = double(arr[8,sub_subhalo])
Mtidal = double(arr[6,*])
Rs = double(arr[7,*])
Rtidal = double(arr[8,*])

N_time = N_elements(fname_part)
select_snap = 7

t_snap = [6.220,7.110,7.790,7.970,8.290,10.000,12.419]
if select_snap le 6 then begin
  out_name = dir_inVL+'tail_VL13_'+strtrim(string(t_snap[select_snap],f='(f7.3)'),2)+'_w5000sub'
  start = fix(t_snap[select_snap]*1000)-30
  over = fix(t_snap[select_snap]*1000)+30
  skip = 1
endif else if select_snap eq 7 then begin
  out_name = dir_inVL+'tail_VL13_w5000sub_last'
  start = 0
  over = 8447;N_time-1
  skip = 8447;N_time-1
endif else begin
  out_name = dir_inVL+'tail_VL13_w5000sub_all'
  start = 0
  over = N_time-1
  skip = 20
endelse


set_plot,'ps'
@plot_setting
!p.charsize=1
!p.multi=[0,2,2]
;out_name = dir_inVL+'tail_VL_2.81Gyr'
;out_name = dir_inVL+'tail_VL_3.19_3.33Gyr'
device,file=out_name+'.ps',/color,/landscape,xoffset=0.4,yoffset=24

;for i=2800,2820 do begin
;for i=3180,3360 do begin
;for i=00,00 do begin
for i=start,over,skip do begin
;for i=0L,N_time-1 do begin

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

  t_step = string(t_part[0]/1000.,f='(f6.3)')
  if (i/50. eq i/50) then print,t_step+'Gyr'

;=== read out data of the subhalo ===
  chr_rdtbl,fname_subhalo[i],0,arr,/silent
  arr = double(arr)
;  arr = double(arr[*,sub_subhalo])
  x_subhalo = reform(arr[1,*])
  y_subhalo = reform(arr[2,*])
  z_subhalo = reform(arr[3,*])
  Vx_subhalo = reform(arr[4,*])
  Vy_subhalo = reform(arr[5,*])
  Vz_subhalo = reform(arr[6,*])
  N_subhalo = strtrim(N_elements(x_subhalo),2)

;=== read out data of the tail without subhalo ===
  ii = i;*10
if i ge 12758 then ii = 12758
  chr_rdtbl,fname_nosub[ii],0,arr,/silent
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


;;=== expected delV & delE ===
;  c = Rtidal[min_sub]/Rs[min_sub]
;  Mhalo = Mtidal[0]/(alog(c+1.)-c/(c+1.))
;  p = b_impact/Rs[min_sub]
;  mr_subhalo = Mhalo*(alog(p+1.)-p/(p+1.))
;  delV_exp = 2.*G*mr_subhalo / b_impact / V_viewed_at_tail
;  delE_exp = delV_exp*V_viewed_at_tail
;;print,V_viewed_at_tail,mr_subhalo/b_impact,delV_exp


;=== plot begins ===
;=== x,y coord. of the tail from 4 subhalos ===
  Mtidal_name = strtrim(fix(alog10(Mtidal[min_sub])),2)
  Rs_name = strtrim(string(Rs[min_sub],f='(f8.4)'),2)
  multiplot
  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]],[y_part[k]],psym=3,color=qcolor[k]
  for l=0L,N_elements(x_subhalo)-1 do begin
    Rs_xy_circle = circle(x_subhalo[l],y_subhalo[l],double(Rs[l]))
    Rtidal_xy_circle = circle(x_subhalo[l],y_subhalo[l],double(Rtidal[l]))
    oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
;    oplot,Rtidal_xy_circle[0,*],Rtidal_xy_circle[1,*],linestyle=2
  endfor
  legend,[t_step+'Gyr',textoidl('M_{tidal}=e'+Mtidal_name),$
         'R_s='+Rs_name+'kpc',$
         'b='+strtrim(string(b_impact,f='(f9.4)'),2)+'kpc'],box=0


;=== in the tail's rest frame (xy plane ) 4 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[-10,10],yr=[-10,10],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]-x_part[0]],[y_part[k]-y_part[0]],psym=3,color=qcolor[k]
  for l=0,N_elements(x_subhalo)-1 do begin
    Rs_xy_circle = circle(x_subhalo[l],y_subhalo[l],double(Rs[l]))
    Rtidal_xy_circle = circle(x_subhalo[l],y_subhalo[l],double(Rtidal[l]))
    oplot,Rs_xy_circle[0,*]-x_part[0],Rs_xy_circle[1,*]-y_part[0]
  endfor
  legend,[textoidl('V_{rel}=')+strtrim(fix(V_viewed_at_tail),2)+'km/s',$
          textoidl('V_{x,rel}=')+strtrim(fix(Vx_viewed_at_tail),2)+'km/s',$
          textoidl('V_{y,rel}=')+strtrim(fix(Vy_viewed_at_tail),2)+'km/s',$
          textoidl('V_{z,rel}=')+strtrim(fix(Vz_viewed_at_tail),2)+'km/s'],box=0

  qmin = 5  &  qmax = -5
  jmin = -2.  &  jmax = 2.
  dEmin = -990  &  dEmax = 990

;=== energy gained 4 subhalos ===
  Ediff = Etot_part-Etot_nosub
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dEmin,dEmax],xtitle='q',ytitle=textoidl('\DeltaE')
  for k=0,N_particle-1 do oplot,[q_part[k]],[Ediff[k]],psym=3,color=qcolor[k]
;  hline,delE_exp,linestyle=1

;=== energy & angular momentum 4 subhalos ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[jmin,jmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part[k]],[scale_dJ_part[k]],psym=3,color=qcolor[k]

  multiplot,/reset
  erase

endfor



device,/close
close,5

spawn, 'convert '+out_name+'.ps '+out_name+'.gif'

END
