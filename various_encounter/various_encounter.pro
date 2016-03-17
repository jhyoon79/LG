pro read_out,xpts,ypts,xfinal,yfinal,qcolor,Rs_xy_circle_vari,Rtidal_xy_circle_vari,qpts,scale_dJpts,rssub=rssub

common share_read,dirpts,dir_subhalo_vari,epsilon,s
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
chr_rdtbl,dirpts+'snap00400',0,arr,/silent
arr = double(arr)
tpts = reform(arr[0,*])
xpts = reform(arr[1,*])
ypts = reform(arr[2,*])
zpts = reform(arr[3,*])
;chr_rdtbl,dirpts+'snap04560',0,arr,/silent
chr_rdtbl,dirpts+'snap04740',0,arr,/silent
arr = double(arr)
tfinal = reform(arr[0,*])
xfinal = reform(arr[1,*])
yfinal = reform(arr[2,*])
zfinal = reform(arr[3,*])
dEpts = reform(arr[7,*])
Etotpts = reform(arr[8,*])
Jpts = reform(arr[9,*])
delE = Etotpts-Etotpts[0]
qpts = delE/epsilon
dJpts = Jpts-Jpts[0]
scale_dJpts = dJpts / (s * Jpts)

; read out data of the subhalo 
chr_rdtbl,dir_subhalo_vari+'subhalo00400',0,arr,/silent
arr = double(arr)
x_subhalo_vari = reform(arr[1,*])
y_subhalo_vari = reform(arr[2,*])
z_subhalo_vari = reform(arr[3,*])

; make a circle of a subhalo 
Rs_xy_circle_vari = circle(x_subhalo_vari[0],y_subhalo_vari[0],double(Rs[0]))
Rtidal_xy_circle_vari = circle(x_subhalo_vari[0],y_subhalo_vari[0],double(Rtidal[0]))

RETURN
END

pro plot_row,leg,ndy

common plot_read,Npts
common share_plot,qmin,qmax,dJmin,dJmax,ii
common share_plot1,xpts1,ypts1,xfinal1,yfinal1,qcolor1,Rs_xy_circle_vari1,Rtidal_xy_circle_vari1,qpts1,scale_dJpts1
common share_plot2,xpts2,ypts2,xfinal2,yfinal2,qcolor2,Rs_xy_circle_vari2,Rtidal_xy_circle_vari2,qpts2,scale_dJpts2
common share_plot3,xpts3,ypts3,xfinal3,yfinal3,qcolor3,Rs_xy_circle_vari3,Rtidal_xy_circle_vari3,qpts3,scale_dJpts3

enc_index = 2025
;=== at the moment of closest approach, 10Gyr ===
if ii eq 1 then xt = textoidl('\Deltax [kpc]')
yt = textoidl('\Deltay [kpc]')
multiplot,xgap=0
xmin = -13
xmax = 8
plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[xmin,xmax],xtit=xt,ytit=yt,/isotropic
;for k=0,Npts-1 do oplot,[xpts1[k]-xpts1[0]],[ypts1[k]-ypts1[0]],psym=3,color=qcolor1[k]
oplot,xpts1-xpts1[0],ypts1-ypts1[0],psym=3
oplot,[0],[0],psym=1
;oplot,Rs_xy_circle_vari1[0,*]-xpts1[0],Rs_xy_circle_vari1[1,*]-ypts1[0],thick=4
oplot,Rtidal_xy_circle_vari1[0,*]-xpts1[0],Rtidal_xy_circle_vari1[1,*]-ypts1[0],linestyle=0,color=100,thick=2
;oplot,Rs_xy_circle_vari2[0,*]-xpts2[0],Rs_xy_circle_vari2[1,*]-ypts2[0],color=70,thick=4
oplot,Rtidal_xy_circle_vari2[0,*]-xpts2[0],Rtidal_xy_circle_vari2[1,*]-ypts2[0],linestyle=0,color=160,thick=2
;oplot,Rs_xy_circle_vari3[0,*]-xpts3[0],Rs_xy_circle_vari3[1,*]-ypts3[0],color=150,thick=4
oplot,Rtidal_xy_circle_vari3[0,*]-xpts3[0],Rtidal_xy_circle_vari3[1,*]-ypts3[0],linestyle=0,color=220,thick=2
lcharsize=1.1
legend,leg,box=0,color=[70,150,220],charsize=lcharsize,psym=8,/bottom
;print,'sub',Rs_xy_circle_vari1[0,*]-xpts1[0],Rs_xy_circle_vari1[1,*]-ypts1[0]

xmin = -11
xmax = 5
ymin = -11
ymax = 5
dy = 3
multiplot,/doyaxis,xgap=-0.02
plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtit=xt,/isotropic
;for k=0,Npts-1 do oplot,[xfinal1[k]-xfinal1[enc_index]],[yfinal1[k]-yfinal1[enc_index]],psym=3,color=qcolor1[k]
;for k=0,Npts-1 do oplot,[xfinal2[k]-xfinal2[enc_index]],[yfinal2[k]-yfinal2[enc_index]-dy],psym=3,color=qcolor1[k]
;for k=0,Npts-1 do oplot,[xfinal3[k]-xfinal3[enc_index]],[yfinal3[k]-yfinal3[enc_index]-dy*2],psym=3,color=qcolor1[k]
oplot,xfinal1-xfinal1[enc_index],yfinal1-yfinal1[enc_index],psym=3
oplot,xfinal2-xfinal2[enc_index],yfinal2-yfinal2[enc_index]-dy,psym=3
oplot,xfinal3-xfinal3[enc_index],yfinal3-yfinal3[enc_index]-dy*2,psym=3
xyouts,4.0,-4.1-ndy*dy,'(f)',alignment=0.5,charsize=lcharsize
legend,'+4.34Gyr',box=0,/bottom,charsize=lcharsize

; energy & angular momentum 
if ii eq 1 then xt = 'q'
yt = textoidl('\Delta J / sJ')
multiplot,/doyaxis,xgap=0.02
plot,[0],[0],/nodata,xr=[qmin,qmax],yr=[dJmin,dJmax],xtit=xt,ytit=yt
;for k=0,Npts-1 do oplot,[qpts1[k]],[scale_dJpts1[k]],psym=3,color=qcolor1[k]
;for k=0,Npts-1 do oplot,[qpts2[k]],[scale_dJpts2[k]]-3,psym=3,color=qcolor2[k]
;for k=0,Npts-1 do oplot,[qpts3[k]],[scale_dJpts3[k]]-6,psym=3,color=qcolor3[k]
oplot,qpts1,scale_dJpts1,psym=3
oplot,qpts2,scale_dJpts2-3,psym=3
oplot,qpts3,scale_dJpts3-6,psym=3
xyouts,-4.1,-3.5-ndy*dy,'(f)',alignment=0.5,charsize=lcharsize
legend,'+4.34Gyr',box=0,/bottom,charsize=lcharsize

RETURN
END


pro various_encounter

common share_read,dirpts,dir_subhalo_vari,epsilon,s
common share_plot,qmin,qmax,dJmin,dJmax,ii
common share_plot1,xpts1,ypts1,xfinal1,yfinal1,qcolor1,Rs_xy_circle_vari1,Rtidal_xy_circle_vari1,qpts1,scale_dJpts1
common share_plot2,xpts2,ypts2,xfinal2,yfinal2,qcolor2,Rs_xy_circle_vari2,Rtidal_xy_circle_vari2,qpts2,scale_dJpts2
common share_plot3,xpts3,ypts3,xfinal3,yfinal3,qcolor3,Rs_xy_circle_vari3,Rtidal_xy_circle_vari3,qpts3,scale_dJpts3

dir_out = '/media/SEADISK/LG/various_encounter/'
;dir_nosub = '/media/SEADISK/LG/one_encounter_test_real/'
G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

chr_rdtbl,dir_out+'part_peri_xy_nosub',0,arr,/silent
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
;print,s,epsilon,r_peri

;============
;=== plot ===
;============
set_plot,'ps'
@plot_setting
loadct,0
!p.charsize=1.2
!p.multi=[0,3,3]
file_out = dir_out+'various_encounter.eps'
device,file=file_out,/color,/landscape,/enc,/cmyk,yoffset=24

dirpts = '/media/SEADISK/LG/various_encounter/v100/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/v100/snapshot_subhalo_xy/'
;dirpts = '/media/SEADISK/LG/various_encounter/snapshot_xy/'
;dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/snapshot_subhalo_xy/'
read_out,xpts1,ypts1,xfinal1,yfinal1,qcolor1,Rs_xy_circle_vari1,Rtidal_xy_circle_vari1,qpts1,scale_dJpts1
dirpts = '/media/SEADISK/LG/various_encounter/v200/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/v200/snapshot_subhalo_xy/'
read_out,xpts2,ypts2,xfinal2,yfinal2,qcolor2,Rs_xy_circle_vari2,Rtidal_xy_circle_vari2,qpts2,scale_dJpts2
dirpts = '/media/SEADISK/LG/various_encounter/v400/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/v400/snapshot_subhalo_xy/'
read_out,xpts3,ypts3,xfinal3,yfinal3,qcolor3,Rs_xy_circle_vari3,Rtidal_xy_circle_vari3,qpts3,scale_dJpts3
qmin = 5  &  qmax = -5
dJmin = -10.9 &  dJmax = 4.9
leg1 = [textoidl('v_{enc}=100km/s'),textoidl('v_{enc}=200km/s (f)'),textoidl('v_{enc}=400km/s')]
leg2 = [textoidl('b=0R_s (f)'),textoidl('b=4R_s'),textoidl('b=8R_s')]
leg3 = [textoidl('M_{sub}=1e5.5'),textoidl('M_{sub}=1e6.5 (f)'),textoidl('M_{sub}=1e7.5')]
ii = 0
plot_row,leg1,1

dirpts = '/media/SEADISK/LG/various_encounter/v200/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/v200/snapshot_subhalo_xy/'
read_out,xpts1,ypts1,xfinal1,yfinal1,qcolor1,Rs_xy_circle_vari1,Rtidal_xy_circle_vari1,qpts1,scale_dJpts1
dirpts = '/media/SEADISK/LG/various_encounter/4Rs/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/4Rs/snapshot_subhalo_xy/'
read_out,xpts2,ypts2,xfinal2,yfinal2,qcolor2,Rs_xy_circle_vari2,Rtidal_xy_circle_vari2,qpts2,scale_dJpts2
dirpts = '/media/SEADISK/LG/various_encounter/8Rs/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/8Rs/snapshot_subhalo_xy/'
read_out,xpts3,ypts3,xfinal3,yfinal3,qcolor3,Rs_xy_circle_vari3,Rtidal_xy_circle_vari3,qpts3,scale_dJpts3
plot_row,leg2,0

dirpts = '/media/SEADISK/LG/various_encounter/e5.5/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/e5.5/snapshot_subhalo_xy/'
read_out,xpts1,ypts1,xfinal1,yfinal1,qcolor1,Rs_xy_circle_vari1,Rtidal_xy_circle_vari1,qpts1,scale_dJpts1,rssub='6'
dirpts = '/media/SEADISK/LG/various_encounter/v200/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/v200/snapshot_subhalo_xy/'
read_out,xpts2,ypts2,xfinal2,yfinal2,qcolor2,Rs_xy_circle_vari2,Rtidal_xy_circle_vari2,qpts2,scale_dJpts2
dirpts = '/media/SEADISK/LG/various_encounter/e7.5/snapshot_xy/'
dir_subhalo_vari = '/media/SEADISK/LG/various_encounter/e7.5/snapshot_subhalo_xy/'
read_out,xpts3,ypts3,xfinal3,yfinal3,qcolor3,Rs_xy_circle_vari3,Rtidal_xy_circle_vari3,qpts3,scale_dJpts3,rssub='8'
ii = 1
plot_row,leg3,1

device,/close
;dir = '/media/SEADISK/LG/various_encounter/'
;spawn,'convert '+dir+'various_encounter.eps -resize 350x450 '+dir+'various_encounter.small.eps'

END
