pro make_encounter_xy

;out_dir = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
;chr_rdtbl,out_dir+'part001_xy',0,arr
;;chr_rdtbl,out_dir+'part001_nosub_xy',0,arr
;arr = double(arr)
;t = arr[0,*]
;x = arr[1,*]
;y = arr[2,*]
;z = arr[3,*]
;
;t_meet = 10000;1411.	; meet at 0.1Gyr
;sub_meet = where(t eq t_meet)
;x_meet = x[sub_meet]
;y_meet = y[sub_meet]
;z_meet = z[sub_meet]
;print,sub_meet,x_meet,y_meet,z_meet
t_begin = 5500.
t_meet = 6000.
dir = '/scratch/jhyoon/Research/LG/one_encounter_test_real/'
fname = dir+'snapshot_xy_nosub/snap'+string(t_meet-t_begin,f='(i5.5)')
chr_rdtbl,fname,0,arr
arr = double(arr)
index = 1673
;t_meet = arr[0,index]
x_meet = arr[1,index]
y_meet = arr[2,index]
z_meet = arr[3,index]


;=== compute the properties of subhalo ===
logMtidal = 7.d
logRtidal = -3.8768103+0.63632018*logMtidal
logRs = -3.5173388+0.42683925*logMtidal
;Mtidal = 1e7
;Rtidal = 2.6942992
;Rs = 0.29628733
Mtidal = 10^logMtidal
Rtidal = 10^logRtidal
Rs = 10^logRs


kpc2km = 3.0827d16
Myr2sec = 3.1536d13
print,Mtidal,Rtidal,Rs


;=== parameters to vary ===
G = 4.3d-6
;r_peri = 5.53221
;r_tide = 0.108653
;v_circ = 200
;s = (1e4/1.8498359e10)^(1./3.);r_tide/r_peri
b_impact = 0*Rs;-0.05; impact parameter b=1pc
;b_cri = G*Mtidal/(s*v_circ^2.)
;A = 10.
;b_impact = A*2./3.*b_cri
print, 'b=',b_impact

t_meet2 = t_meet - t_begin
Vxsubhalo = 0
Vysubhalo = 0
Vzsubhalo = -200	; b test:v=400km/s
Vsubhalo = sqrt(Vxsubhalo^2.+Vysubhalo^2.+Vzsubhalo^2.)
;Vsubhalo = 200.
;Vxsubhalo = sqrt(Vsubhalo^2.*192./360.)
;Vysubhalo = -sqrt(Vsubhalo^2.*168./360.)
;Vzsubhalo = -sqrt(Vsubhalo^2.*0./360.)	; b test:v=400km/s
;print,sqrt(Vxsubhalo^2.+Vysubhalo^2.+Vzsubhalo^2.),Vsubhalo
x_subhalo = x_meet + b_impact/sqrt(2.) - Vxsubhalo/kpc2km*Myr2sec*t_meet2
y_subhalo = y_meet + b_impact/sqrt(2.) - Vysubhalo/kpc2km*Myr2sec*t_meet2
z_subhalo = z_meet - Vzsubhalo/kpc2km*Myr2sec*t_meet2

eff_enc = abs(b_impact / (G*Mtidal/Vsubhalo^2))
print,'eff_enc',eff_enc
;readcol,'frog_subrotate_1411Myr',t,x_subhalo,y_subhalo,z_subhalo,vxsubhalo,vysubhalo,vzsubhalo
;x_subhalo = x_subhalo - Vxsubhalo/kpc2km*Myr2sec*t_meet2
;y_subhalo = y_subhalo - Vysubhalo/kpc2km*Myr2sec*t_meet2
;z_subhalo = z_subhalo - Vzsubhalo/kpc2km*Myr2sec*t_meet2

openw,1,'orbit_subhalo_xy.dat'
printf,1,x_subhalo,y_subhalo,z_subhalo,Vxsubhalo,Vysubhalo,Vzsubhalo,Mtidal,Rs,Rtidal,b_impact,eff_enc,$
    f='(6(f14.8,1x),a16,2(1x,f13.8),f7.3,f8.2)'
close,1


END
