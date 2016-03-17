pro movie_1encounter_xy

dir_out = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy/'
dir_subhalo = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_subhalo_xy/'
dir_nosub = '/scratch/jhyoon/Research/LG/one_encounter_test_real/snapshot_xy_nosub/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
chr_rdtbl,dir_out+'part_peri_xy',0, arr
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
V_peri = sqrt(Vx_peri^2.+Vy_peri^2.+Vz_peri^2.)
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
Msat = 10000.
s = (Msat/mr)^(1./3.)
epsilon = s*(4.3e-6*mr/r_peri)
;=== read out simulation info ===
chr_rdtbl,'orbit_subhalo_xy.dat',0,arr
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
eff_enc = arr[10,*]
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)


;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part+'snap00000',0, arr
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (s * J_initial)
qmin = min(q_initial)-0.3
qmax = max(q_initial)+0.3
dJmin = min(scale_dJ_initial)-0.2
dJmax = max(scale_dJ_initial)+0.2
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)


;=== define & read out the tail data ===
spawn, 'ls '+dir_part+' > temp9'
chr_rdtbl,'temp9',0,fname_part
fname_part = dir_part+reform(fname_part)

spawn, 'ls '+dir_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo
fname_subhalo = dir_subhalo+reform(fname_subhalo)

spawn, 'ls '+dir_nosub+' > temp9'
chr_rdtbl,'temp9',0,fname_nosub
fname_nosub = dir_nosub+reform(fname_nosub)
spawn, 'rm -f temp9'

set_plot,'ps'
@plot_setting
!p.charsize=1
!p.multi=[0,2,3]
file_out = dir_out+'movie_1encounter_xy'
device,file=file_out+'.ps',/color,xsize=20,ysize=20,xoffset=0.5,yoffset=0.5

N_time = N_elements(fname_part)
for i=0,N_time-1,10 do begin

;=== read out data of the tail particles ===
  chr_rdtbl,fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part = 5500+reform(arr[0,*])
  x_part = reform(arr[1,*])
  y_part = reform(arr[2,*])
  z_part = reform(arr[3,*])
  Vx_part = reform(arr[4,*])
  Vy_part = reform(arr[5,*])
  Vz_part = reform(arr[6,*])
  dE_part = reform(arr[7,*])
  Etot_part = reform(arr[8,*])
  J_part = reform(arr[9,*])
  delE = Etot_part-Etot_part[0]
  q_part = delE/epsilon
  dJ_part = J_part-J_part[0]
  scale_dJ_part = dJ_part / (s * J_part)
  t_step = string(t_part[0]/1000.,f='(f7.4)')
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
  max_delE_subhalo = reform(arr[8,*])


;=== draw a circle of a subhalo ===
  Rs_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rs[0]))
  Rtidal_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rtidal[0]))
  Rs_yz_circle = circle(y_subhalo[0],z_subhalo[0],double(Rs[0]))
  Rtidal_yz_circle = circle(y_subhalo[0],z_subhalo[0],double(Rtidal[0]))
  Rs_xz_circle = circle(x_subhalo[0],z_subhalo[0],double(Rs[0]))
  Rtidal_xz_circle = circle(x_subhalo[0],z_subhalo[0],double(Rtidal[0]))

;=== read out data of the tail without subhalo ===
  chr_rdtbl,fname_nosub[i],0,arr,/silent
  arr = double(arr)
  x_nosub = reform(arr[1,*])
  y_nosub = reform(arr[2,*])
  z_nosub = reform(arr[3,*])
  Vx_nosub = reform(arr[4,*])
  Vy_nosub = reform(arr[5,*])
  Vz_nosub = reform(arr[6,*])
  i_name = '_'+string(i,f='(i3.3)')

;=== estimate b,delv,delE ===
  d2subhalo = sqrt((x_part-x_subhalo[0])^2.+(y_part-y_subhalo[0])^2.+(z_part-z_subhalo[0])^2.)
  dxy2subhalo = sqrt((x_part-x_subhalo[0])^2.+(y_part-y_subhalo[0])^2.)
  b_impact = min(d2subhalo,min_sub)
  bxy_impact = min(dxy2subhalo,min_sub)
  delVx = (Vx_part[min_sub]-Vx_nosub[min_sub])
  delVy = (Vy_part[min_sub]-Vy_nosub[min_sub])
  delVz = (Vz_part[min_sub]-Vz_nosub[min_sub])
  delV = sqrt(delVx^2.+delVy^2.+delVz^2.)
  delEx = Vx_part[min_sub]*delVx
  delEy = Vy_part[min_sub]*delVy
  delEz = Vz_part[min_sub]*delVz
  delE = sqrt(delEx^2.+delEy^2.+delEz^2.)
  V_viewed_at_tail = sqrt(Vx_subhalo^2.+Vy_subhalo^2.+Vz_subhalo^2.)
;  V_viewed_at_tail = sqrt((Vx_subhalo-Vx_part[min_sub])^2.+(Vy_subhalo-Vy_part[min_sub])^2.+(Vz_subhalo-Vz_part[min_sub])^2.)
  delV_expected = G*double(Mtidal) / (b_impact*V_viewed_at_tail)
  if (t_step eq 1.60) then printf,5,delV,delV_expected,delE,f='(33317.5),1x)'

;  Vx_rel = Vx_subhalo-Vx_part[min_sub]
;  Vy_rel = Vy_subhalo-Vy_part[min_sub]
;  Vz_rel = Vz_subhalo-Vz_part[min_sub]
  Vx_rel = Vx_subhalo
  Vy_rel = Vy_subhalo
  Vz_rel = Vz_subhalo

  delE_expect = G*Mtidal/(b_impact*V_viewed_at_tail)


;=== plot begins ===
;=== x,y coord. of the tail ===
  multiplot
  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='y',title=t_step+'Gyr '+textoidl('V_{rel}=')+strtrim(fix(V_viewed_at_tail),2)+'km/s',/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]],[y_part[k]],psym=3,color=qcolor[k]
  oplot,Rs_xy_circle[0,*],Rs_xy_circle[1,*]
  oplot,Rtidal_xy_circle[0,*],Rtidal_xy_circle[1,*],linestyle=2
  legend,[textoidl('V_{x,rel}=')+strtrim(fix(Vx_rel),2)+'km/s',$
          textoidl('V_{y,rel}=')+strtrim(fix(Vy_rel),2)+'km/s',$
          textoidl('V_{z,rel}=')+strtrim(fix(Vz_rel),2)+'km/s'],box=0

;=== y,z coord. of the tail ===
  multiplot
  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='y',ytitle='z',/isotropic
  for k=0,N_particle-1 do oplot,[y_part[k]],[z_part[k]],psym=3,color=qcolor[k]
  oplot,Rs_yz_circle[0,*],Rs_yz_circle[1,*]
  oplot,Rtidal_yz_circle[0,*],Rtidal_yz_circle[1,*],linestyle=2
  legend,['b='+strtrim(b_impact,2),'bxy='+strtrim(bxy_impact,2),textoidl('\DeltaE_{max}=')+strtrim(max_delE_subhalo,2)],box=0

;=== x,z coord. of the tail ===
  multiplot
  plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x',ytitle='z',/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]],[z_part[k]],psym=3,color=qcolor[k]
  oplot,Rs_xz_circle[0,*],Rs_xz_circle[1,*]
  oplot,Rtidal_xz_circle[0,*],Rtidal_xz_circle[1,*],linestyle=2


;=== in the tail's rest frame ===
  multiplot
  plot,[0],[0],/nodata,xr=[-6,6],yr=[-6,6],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
  for k=0,N_particle-1 do oplot,[x_part[k]-x_part[0]],[y_part[k]-y_part[0]],psym=3,color=qcolor[k]
  oplot,Rs_xy_circle[0,*]-x_part[0],Rs_xy_circle[1,*]-y_part[0]
  oplot,Rtidal_xy_circle[0,*]-x_part[0],Rtidal_xy_circle[1,*]-y_part[0],linestyle=2

;;=== in the subhalos' rest frame ===
;  multiplot
;  plot,[0],[0],/nodata,xr=[-Rs[0]*1.1,Rs[0]*1.1],yr=[-Rs[0]*1.1,Rs[0]*1.1],xtitle=textoidl('\Deltax'),ytitle=textoidl('\Deltay'),/isotropic
;  for k=0,N_particle-1 do oplot,[x_part[k]-x_subhalo[0]],[y_part[k]-y_subhalo[0]],psym=3,color=qcolor[k]
;  oplot,Rs_xy_circle[0,*]-x_subhalo[0],Rs_xy_circle[1,*]-y_subhalo[0]
;  oplot,Rtidal_xy_circle[0,*]-x_subhalo[0],Rtidal_xy_circle[1,*]-y_subhalo[0],linestyle=2

;=== energy & angular momentum ===
  multiplot
  plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
  for k=0,N_particle-1 do oplot,[q_part[k]],[scale_dJ_part[k]],psym=3,$
      color=qcolor[k]

  multiplot
  plothist,q_part,bin=0.1,xtitle='q',/peak

  multiplot,/reset
  erase
endfor


device,/close
close,5

;spawn,'convert '+file_out+'.ps '+file_out+'.gif'
;readcol,'orbit_subhalo_xy.dat',x,y,z,vx,vy,vz
;name = '_v'+strtrim(fix(abs(vz)),2)
name = '_VxVy'
print,'name=',name
spawn,'cp '+file_out+'.ps '+file_out+name+'.ps'
spawn,'cp '+dir_out+'part001_xy '+dir_out+'part001_xy'+name
spawn,'cp '+dir_out+'part_peri_xy '+dir_out+'part_peri_xy'+name
spawn,'cp -r '+dir_out+'snapshot_xy '+dir_out+'snapshot_xy'+name
spawn,'cp -r '+dir_out+'snapshot_subhalo_xy '+dir_out+'snapshot_subhalo_xy'+name


END
