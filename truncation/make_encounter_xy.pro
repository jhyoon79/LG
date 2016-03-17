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
t_begin = 9000.

t_meet = 10.000*1000.
dir = '/scratch/jhyoon/Research/LG/one_encounter_test_evol/'
fname = dir+'snapshot_xy_nosub/snap'+string(t_meet-t_begin,f='(i5.5)')
chr_rdtbl,fname,0,arr
arr = double(arr)
;index = 2989
index = 99
t_meet = arr[0,index]
x_meet = arr[1,index]
y_meet = arr[2,index]
z_meet = arr[3,index]
print,fname,t_meet
print,x_meet,y_meet,z_meet


;=== compute the properties of subhalo ===
;logMtidal = 7.d
;logRtidal = -3.8768103+0.63632018*logMtidal
;logRs = -3.5173388+0.42683925*logMtidal
Mtidal = 1e7
Rtidal = 2.6942992
Rs = 0.29628733
;Mtidal = 10^logMtidal
;Rtidal = 10^logRtidal
;Rs = 10^logRs


kpc2km = 3.0827d16
Myr2sec = 3.1536d13
print,Mtidal,Rtidal,Rs

;=== parameters to vary ===
b_impact = -1*Rs;Rs*2.;0.;0.45;0.3;0.15; impact parameter b=1pc

t_meet2 = t_meet - t_begin
;=== for test with gap ===
Vxsubhalo = 0
Vysubhalo = 0
Vzsubhalo = -50	; b test:v=400km/s
;=== for test with kink(e7) ===
;Vxsubhalo = 50
;Vysubhalo = -10
;Vzsubhalo = 0	; b test:v=400km/s
;=== for test with kink(e8) ===
;Vxsubhalo = 50
;Vysubhalo = -20
;Vzsubhalo = 0	; b test:v=400km/s
x_subhalo = x_meet + b_impact/sqrt(2.) - Vxsubhalo/kpc2km*Myr2sec*t_meet2
y_subhalo = y_meet + b_impact/sqrt(2.) - Vysubhalo/kpc2km*Myr2sec*t_meet2
z_subhalo = z_meet - Vzsubhalo/kpc2km*Myr2sec*t_meet2

openw,1,'orbit_subhalo_xy.dat'
printf,1,x_subhalo,y_subhalo,z_subhalo,Vxsubhalo,Vysubhalo,Vzsubhalo,Mtidal,Rs,Rtidal,b_impact,$
    f='(6(f14.8,1x),a16,2(1x,f13.8),f7.3)'
close,1


END
