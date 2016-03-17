pro make_evol_encounter_xy

t_begin = 5600.
t_meet = 6000.
dir = '/media/SEADISK/LG/various_encounter/evol_nosub/'
fname = dir+'snapshot_xy/snap'+string(t_meet-t_begin,f='(i5.5)')
chr_rdtbl,fname,0,arr
arr = double(arr)
index = 1673
;t_meet = arr[0,index]
x_meet = arr[1,index]
y_meet = arr[2,index]
z_meet = arr[3,index]
Vx_meet = arr[4,index]
Vy_meet = arr[5,index]
Vz_meet = arr[6,index]


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
b_impact = 0*Rs;-0.05; impact parameter b=1pc
print, 'b=',b_impact

t_meet2 = t_meet - t_begin
Vxsubhalo = Vx_meet
Vysubhalo = Vy_meet
Vzsubhalo = -200	; b test:v=400km/s
Vsubhalo = sqrt(Vxsubhalo^2.+Vysubhalo^2.+Vzsubhalo^2.)
x_subhalo = x_meet - b_impact/sqrt(10.) - Vxsubhalo/kpc2km*Myr2sec*t_meet2
y_subhalo = y_meet + b_impact/sqrt(10./9.) - Vysubhalo/kpc2km*Myr2sec*t_meet2
z_subhalo = z_meet - Vzsubhalo/kpc2km*Myr2sec*t_meet2

eff_enc = abs(b_impact / (G*Mtidal/Vsubhalo^2))
print,'eff_enc',eff_enc

openw,1,'orbit_subhalo_xy_evol.dat'
printf,1,x_subhalo,y_subhalo,z_subhalo,Vxsubhalo,Vysubhalo,Vzsubhalo,Mtidal,Rs,Rtidal,b_impact,eff_enc,$
    f='(6(f14.8,1x),a16,2(1x,f13.8),f7.3,f8.2)'
close,1


END
