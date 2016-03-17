pro proceeding

dir_part = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/snapshot_xy/'
dir_real = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
dir_evol = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/'

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
r_tide = 0.108653
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
qcolor = (230-round((q_initial-qmin)/(qmax-qmin)*230.))
Npts = N_elements(qcolor)
qmin = qmin-0.2
qmax = qmax+0.2
dJmin = dJmin-0.05
dJmax = dJmax+0.05

chr_rdtbl,'../one_encounter_test_real/frog_rotate_0Gyr',0,arr
arr = double(arr)
t0 = arr[0,*]
x0 = arr[1,*]
y0 = arr[2,*]
z0 = arr[3,*]

chr_rdtbl,'../one_encounter_test_real/frog_rotate_10Gyr',0,arr
arr = double(arr)
t10 = arr[0,*]
x10 = arr[1,*]
y10 = arr[2,*]
z10 = arr[3,*]

chr_rdtbl,dir_evol+'part001_t0',0,arr
Nend = 389
arr = double(arr[*,0:Nend])
t_orbit = arr[0,*]
x_orbit = arr[1,*]
y_orbit = arr[2,*]
z_orbit = arr[3,*]

set_plot,'ps'
@plot_setting
loadct,0
!p.charsize=2.
device,file='proceeding_f1.eps',/enc,/color,/cmyk,/landscape,ysize=7.3,yoffset=24
!p.multi=[0,3,1]
plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x [kpc]', ytitle='y [kpc]',/isotropic
oplot,[0],[0],psym=1
for k=0,Npts-1 do oplot,[x0[k]],[y0[k]],psym=3,color=qcolor[k]
legend,['t=-10Gyr'],box=0,charsize=1.2
oplot,x_orbit,y_orbit,color=150,linestyle=2
arrow,x_orbit[Nend-1],y_orbit[Nend-1],x_orbit[Nend],y_orbit[Nend],/data,color=150
xmin = x0[0]-0.25
xmax = x0[0]+0.25
ymin = y0[0]-0.25
ymax = y0[0]+0.25
tname = replicate(' ',10)
oplot,[xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin]
oplot,[xmax,2.1],[ymax,-2.]
oplot,[xmax,2.1],[ymin,-19]
plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic,/noerase,position=[0.21,0.21,0.48,0.51],charsize=1,xtickname=tname,ytickname=tname
for k=0,Npts-1 do oplot,[x0[k]],[y0[k]],psym=3,color=qcolor[k]
oplot,x_orbit,y_orbit,color=150,linestyle=2


plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x [kpc]', ytitle='y [kpc]',/isotropic
for k=0,Npts-1 do oplot,[x10[k]],[y10[k]],psym=3,color=qcolor[k]
oplot,[x10[99]],[y10[99]],psym=7
legend,['t=0Gyr'],box=0,charsize=1.2


plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,Npts-1 do oplot,[q_initial[k]],[scale_dJ_initial[k]],psym=3,color=qcolor[k]
oplot,[q_initial[99]],[scale_dJ_initial[99]],psym=7

device,/close

chr_rdtbl,dir_part+'snap06000',0,arr,/silent
arr = double(arr)
t_final = arr[0,*]
x_final = arr[1,*]
y_final = arr[2,*]
z_final = arr[3,*]

Etot_final = arr[8,*]
J_final = arr[9,*]
delE = Etot_final-Etot_final[0]
q_final = delE/epsilon
dJ_final = J_final-J_final[0]
scale_dJ_final = dJ_final / (r_tide/r_peri * J_final)
qmin = min(q_final)
qmax = max(q_final)
dJmin = min(scale_dJ_final)
dJmax = max(scale_dJ_final)
qcolor = (230-round((q_final-qmin)/(qmax-qmin)*230.))
Npts = N_elements(qcolor)
qmin = qmin-0.2
qmax = qmax+0.2
dJmin = dJmin-0.05
dJmax = dJmax+0.05

device,file='proceeding_f2.eps',/enc,/color,/cmyk,/landscape,ysize=11,yoffset=24
!p.charsize=1.4
!p.multi=[0,2,1]
plot,[0],[0],/nodata,xr=[-20,20],yr=[-20,20],xtitle='x [kpc]', ytitle='y [kpc]',/isotropic
for k=0,Npts-1 do oplot,[x_final[k]],[y_final[k]],psym=3,color=qcolor[k]
oplot,[x_final[99]],[y_final[99]],psym=7
legend,['t=6Gyr'],box=0;charsize=1.2
xmin = x_final[99] - 12
ymin = y_final[99] - 5
xmax = x_final[99] + 12
ymax = y_final[99] + 5
print,xmin,xmax
;plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],/isotropic,/noerase,position=[0.21,0.21,0.48,0.51],xtickname=tname,ytickname=tname
;for k=0,Npts-1 do oplot,[x_final[k]],[y_final[k]],psym=3,color=qcolor[k]
plothist,x_final,xr=[xmin,xmax],bin=0.5,/noerase,position=[0.189,0.18,0.388,0.38],xtickname=tname,ytickname=tname,xsty=1,ysty=1
xyouts,-14.,27,'N',orientation=90,charsize=1

plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,Npts-1 do oplot,[q_final[k]],[scale_dJ_final[k]],psym=3,color=qcolor[k]
oplot,[q_final[99]],[scale_dJ_final[99]],psym=7
device,/close

END
