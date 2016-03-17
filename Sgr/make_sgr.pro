pro make_sgr

;=== print out starting time ===
openw,31,'time_start'
printf,31,systime(1)
close,31

readcol, 'sgr_initial', t, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init

openr, 99, 'frogin'
readf, 99, N_stars
close, 99
N_stars = N_stars - 1
G = double(4.3e-6)

;chr_rdtbl,'/scratch/jhyoon/Research/LG/part_peri_nosubhalo',0,arr,/silent
chr_rdtbl,'part_peri',0,arr,/silent
arr = double(arr)
R_peri = sqrt(arr[1,0]^2.+arr[2,0]^2.+arr[3,0]^2.)
openw,8,'R_peri'
print , 'R_peri=',R_peri
printf,8,R_peri


;=== make bar like Pal5 ===
r_tide = 3d
print, 'r_tide=',r_tide
sph_coord = cv_coord(from_rect=[x_init,y_init,z_init],/to_sphere)
longitude = randomn(39,N_stars)*0.000 + sph_coord[0]
latitude = randomn(29,N_stars)*0.000 + sph_coord[1]
radius = randomn(19,N_stars)*r_tide+sph_coord[2]
input_arr = transpose([[longitude],[latitude],[radius]])
rect_coord = cv_coord(from_sphere=input_arr,/to_rect)
x_stars = rect_coord[0,*]
y_stars = rect_coord[1,*]
z_stars = rect_coord[2,*]
help,rect_coord,x_stars,r_spread,longitude,input_arr
;Vx_stars = replicate(Vx_init,N_stars)
;Vy_stars = replicate(Vy_init,N_stars)
;Vz_stars = replicate(Vz_init,N_stars)
sigma_v = 11.4d
Vx_stars = randomn(91,N_stars)*sigma_v + Vx_init[0]
Vy_stars = randomn(82,N_stars)*sigma_v + Vy_init[0]
Vz_stars = randomn(73,N_stars)*sigma_v + Vz_init[0]

openw, 1, 'frog_sgr.dat'
printf, 1, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init, f='(6(f15.9,1x))'
for i=0L,N_stars-1 do printf, 1, x_stars[i], y_stars[i], z_stars[i], Vx_stars[i], Vy_stars[i], Vz_stars[i], f='(6(f15.9,1x))'
close, 1

END
