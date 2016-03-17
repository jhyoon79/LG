pro read_out,xpts,ypts,zpts,vr,qcolor,Rs_xy_circle,Rtidal_xy_circle,qpts,scale_dJpts,rssub=rssub,time=time

common share_read,dirpts,dir_subhalo,epsilon,s
common plot_read,Npts

;=== readout data for the center of Pal5 ===
chr_rdtbl, dirpts+'snap00000',0,arr,/silent
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (s * J_initial)
qmin = min(q_initial)*0.9
qmax = max(q_initial)*0.9
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
Npts = N_elements(qcolor)

;=== read out simulation info ===
if not keyword_set(rssub) then rssub='7'
chr_rdtbl,'orbit_subhalo_e'+rssub+'.dat',0,arr,/silent
Mtidal = arr[6,*]
Rs = strmid(arr[7,*],0,5)+'kpc'
Rtidal = arr[8,*]
b_impact_init = arr[9,*]+'kpc'
Mass_name = 'e'+strtrim(fix(alog10(double(Mtidal[0]))),2)
N_time = N_elements(fnamepts)


;=== vari ===
chr_rdtbl,dirpts+'snap'+time,0,arr
arr = double(arr)
tpts = reform(arr[0,*])
xpts = reform(arr[1,*])
ypts = reform(arr[2,*])
zpts = reform(arr[3,*])
vxpts = reform(arr[4,*])
vypts = reform(arr[5,*])
vzpts = reform(arr[6,*])
rpts = sqrt(xpts^2.+ypts^2.+zpts^2.)
vpts = sqrt(vxpts^2.+vypts^2.+vzpts^2.)
vr = (xpts*vxpts+ypts*vypts+zpts*vzpts)/(rpts*vpts)*vpts
dEpts = reform(arr[7,*])
Etotpts = reform(arr[8,*])
Jpts = reform(arr[9,*])
delE = Etotpts-Etotpts[0]
qpts = delE/epsilon
dJpts = Jpts-Jpts[0]
scale_dJpts = dJpts / (s * Jpts)

print,time,xpts[0:9]

; read out data of the subhalo 
chr_rdtbl,dir_subhalo+'subhalo'+time,0,arr
arr = double(arr)
x_subhalo = reform(arr[1,*])
y_subhalo = reform(arr[2,*])
z_subhalo = reform(arr[3,*])

; make a circle of a subhalo 
Rs_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rs[0]))
Rtidal_xy_circle = circle(x_subhalo[0],y_subhalo[0],double(Rtidal[0]))

RETURN
END

pro plot_row,leg

common plot_read,Npts
common share_plot,qmin,qmax,dJmin,dJmax
common share_plot1,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1

;=== at the moment of closest approach, 10Gyr ===
index=1673
xcen = xpts1[index]
ycen = ypts1[index]
plot,[0],[0],/nodata,xr=[-5,5],yr=[-5,5],xtitle=textoidl('\Deltax [kpc]'),ytitle=textoidl('\Deltay [kpc]'),/isotropic
for k=0,Npts-1 do oplot,[xpts1[k]-xcen],[ypts1[k]-ycen],psym=3,color=qcolor1[k]
oplot,[xpts1[index]-xcen],[ypts1[index]-ycen],psym=1

lcharsize=0.7

; energy & angular momentum 
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtitle='q',ytitle=textoidl('\DeltaJ / sJ')
for k=0,Npts-1 do oplot,[qpts1[k]],[scale_dJpts1[k]],psym=3,color=qcolor1[k]
oplot,[qpts1[index]],[scale_dJpts1[index]],psym=1
for i=0,N_elements(leg)-1 do xyouts,0,-0.4-i*3,leg[i],alignment=0.5,charsize=lcharsize

plot,xpts1-xcen,vr1-vr1[index],psym=3,xr=[-2,2],xtitle='x',ytitle='Vr'
oplot,[xpts1[index]-xcen],[vr1[index]-vr1[index]],psym=1

RETURN
END


pro evol

common share_read,dirpts,dir_subhalo,epsilon,s
common share_plot,qmin,qmax,dJmin,dJmax
common share_plot1,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1

dir_out = '/media/SEADISK/LG/various_encounter/'
dir_nosub = '/media/SEADISK/LG/one_encounter_test_real/'
G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_nosub+'part_peri_xy_nosub',0, arr
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s = (Msat/mr)^(1./3.)
epsilon = s*(4.3e-6*mr/r_peri)

;============
;=== plot ===
;============
set_plot,'ps'
@plot_setting
!p.charsize=1.8
!p.multi=[0,3,4]
file_out = dir_out+'evol.eps'
device,file=file_out,/color,/landscape,/enc,/cmyk,yoffset=24

dirpts = '/media/SEADISK/LG/various_encounter/evol/snapshot_xy/'
dir_subhalo = '/media/SEADISK/LG/various_encounter/evol/snapshot_subhalo_xy/'
read_out,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1,time='00400'
qmin = -4  &  qmax = 4
dJmin = -3 &  dJmax = 4
leg1 = ['100km/s','200km/s','400km/s']
plot_row,leg1

read_out,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1,time='02400'
plot_row,leg1

read_out,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1,time='04400'
plot_row,leg1
read_out,xpts1,ypts1,zpts1,vr1,qcolor1,Rs_xy_circle1,Rtidal_xy_circle1,qpts1,scale_dJpts1,time='06000'
plot_row,leg1



device,/close

END
