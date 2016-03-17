pro make_Pal5_xy

;=== print out starting time ===
openw,31,'time_start'
printf,31,systime(1)
close,31

readcol, 'Pal5_xy_initial', t, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init

openr, 99, 'frogin'
readf, 99, N_stars
close, 99
N_stars = N_stars - 1
G = double(4.3e-6)

chr_rdtbl,'/scratch/jhyoon/Research/LG/part_peri_nosubhalo',0,arr,/silent
arr = double(arr)
R_peri = sqrt(arr[1,0]^2.+arr[2,0]^2.+arr[3,0]^2.)
openw,8,'R_peri'
print , 'R_peri=',R_peri
printf,8,R_peri

;;=== make ball like Pal5 ===
;A = 1.7d;0.88
;B = 0.28d;0.30
;;chr_rdtbl,'part_peri',0,arr,/silent
;close,8
;readcol,'M_Rperi',M_peri
;M_peri = M_peri[0]
;print, 'M_peri',M_peri
;;M_peri = 2.9856467e+10
;;R_peri = 7.036  &  M_peri = 3.12e+10
;
;Msat = 20000.d
;sigma_r = A*R_peri*(Msat/M_peri)^(1./3.)
;
;;sigma_r = 0.108653	; r_c = 0.024kpc, r_tide = 0.108653(16.1arcmin*23.2kpc / 60.*!dtor) from Oden+02
;;Msat = (sigma_r/R_peri)^3.*M_peri
;sigma_v = B*sqrt(G*Msat/sigma_r); + 1e-10
;;sigma_v = 0.22		; sigma_v = 0.22km/s(Dehnen), 1.1km/s(Oden)
;print, 'A=',A,'   B=',B,'   Msat=',Msat,'   sigma_r=',sigma_r,'   sigma_v=',sigma_v
;
;sigma_r = sigma_r/1.
;Vx_stars = randomn(91,N_stars)*sigma_v + Vx_init[0]
;Vy_stars = randomn(82,N_stars)*sigma_v + Vy_init[0]
;Vz_stars = randomn(73,N_stars)*sigma_v + Vz_init[0]
;x_stars = randomn(64,N_stars)*sigma_r + x_init[0]
;y_stars = randomn(55,N_stars)*sigma_r + y_init[0]
;z_stars = randomn(46,N_stars)*sigma_r + z_init[0]
;help, Vy_stars
;tmp1 = moment(Vy_stars)
;print, tmp1[0], sqrt(tmp1[1])

;;===========================


;=== make bar like Pal5 ===
r_tide = 0.108653d
sigma_r = r_tide/0.2
sph_coord = cv_coord(from_rect=[x_init,y_init,z_init],/to_sphere)
longitude = randomn(39,N_stars)*0.000 + sph_coord[0]
latitude = randomn(29,N_stars)*0.000 + sph_coord[1]
radial_out = randomu(19,N_stars/2)*sigma_r+r_tide*0.5+sph_coord[2]
radial_in = -1.*(randomu(9,N_stars/2+1))*sigma_r-r_tide*0.5+sph_coord[2]
radial = [radial_in,radial_out]

input_arr = transpose([[longitude],[latitude],[radial]])
rect_coord = cv_coord(from_sphere=input_arr,/to_rect)
x_stars = rect_coord[0,*]
y_stars = rect_coord[1,*]
z_stars = rect_coord[2,*]
sigma_v = 0.15
Vx_stars = randomn(91,N_stars)*sigma_v + Vx_init[0]
Vy_stars = randomn(82,N_stars)*sigma_v + Vy_init[0]
Vz_stars = randomn(73,N_stars)*sigma_v + Vz_init[0]

openw, 1, 'frog_Pal5_xy.dat'
printf, 1, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init, f='(6(f15.9,1x))'
for i=0L,N_stars-1 do printf, 1, x_stars[i], y_stars[i], z_stars[i], Vx_stars[i], Vy_stars[i], Vz_stars[i], f='(6(f15.9,1x))'
close, 1

END
