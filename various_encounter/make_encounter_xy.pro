pro make_encounter_xy

t_begin = 5500.
t_meet = 5900.
;dir = '/media/SEADISK/LG/one_encounter_test_real/'
dir = '/media/SEADISK/LG/various_encounter/'
fname = dir+'snapshot_xy_nosub6500/snap'+string(t_meet-t_begin,f='(i5.5)')
chr_rdtbl,fname,0,arr
arr = double(arr)
index = 2025;1673
;t_meet = arr[0,index]
x_meet = arr[1,index]
y_meet = arr[2,index]
z_meet = arr[3,index]
vx_meet = arr[4,index]
vy_meet = arr[5,index]
vz_meet = arr[6,index]


;=== compute the properties of subhalo ===
logMtidal = 5.5;[5.5,6.5,7.5,8.5,9.5]
;logRtidal = -3.8768103+0.63632018*logMtidal
;logRs = -3.5173388+0.42683925*logMtidal
;Mtidal = 1e7
;Rtidal = 2.6942992
;Rs = 0.29628733

;logRtidal = -3.1570183+0.46490731*logMtidal ; within Rgc < 50kpc
;logRs = -3.4927735+0.39870709*logMtidal
logRtidal = -3.1253187+0.44552343*logMtidal ; within Rgc < 25kpc
logRs = -3.5713973+0.40638194*logMtidal
Mtidal = 10d^logMtidal
Rtidal = 10d^logRtidal
Rs = 10d^logRs

kpc2km = 3.0827d16
Myr2sec = 3.1536d13
print,Mtidal,Rs,Rtidal

;=== parameters to vary ===
G = 4.3d-6
b_impact = 0*Rs;-0.05; impact parameter b=1pc
print, 'b=',b_impact

t_meet2 = t_meet - t_begin
Vxsubhalo = 0 + vx_meet
Vysubhalo = 0; + vy_meet
Vzsubhalo = -50	; b test:v=400km/s
Vsubhalo = sqrt(Vxsubhalo^2.+Vysubhalo^2.+Vzsubhalo^2.)
;Vsubhalo = 200.
;Vxsubhalo = sqrt(Vsubhalo^2.*192./360.)
;Vysubhalo = -sqrt(Vsubhalo^2.*168./360.)
;Vzsubhalo = -sqrt(Vsubhalo^2.*0./360.)	; b test:v=400km/s
;print,sqrt(Vxsubhalo^2.+Vysubhalo^2.+Vzsubhalo^2.),Vsubhalo
x_subhalo = x_meet + b_impact*0.8064 - Vxsubhalo/kpc2km*Myr2sec*t_meet2
y_subhalo = y_meet + b_impact*0.5914 - Vysubhalo/kpc2km*Myr2sec*t_meet2
z_subhalo = z_meet - Vzsubhalo/kpc2km*Myr2sec*t_meet2

eff_enc = abs(b_impact / (G*Mtidal/Vsubhalo^2))
print,'eff_enc',eff_enc

openw,1,'orbit_subhalo_xy.dat'
printf,1,x_subhalo,y_subhalo,z_subhalo,Vxsubhalo,Vysubhalo,Vzsubhalo,Mtidal,Rs,Rtidal,b_impact,eff_enc,$
    f='(6(f14.8,1x),a16,2(1x,f13.8),f7.3,f8.2)'
close,1


END
