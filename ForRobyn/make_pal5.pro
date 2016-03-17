;==================================
pro make_pal5

;=== print out starting time ===
readcol, 'Pal5_initial', t, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init
openr, 99, 'frogin'
readf, 99, N_stars
close, 99
N_stars = N_stars; - 1
G = 4.3d-6

;dir = '/scratch/jhyoon/Research/LG/realistic/'
;chr_rdtbl,dir+'part001',0,arr,/silent
chr_rdtbl,'part0001',0,arr,/silent
arr = double(arr)
r_orbit = sqrt(arr[1,*]^2.+arr[2,*]^2.+arr[3,*]^2.)
v_orbit = sqrt(arr[4,*]^2.+arr[5,*]^2.+arr[6,*]^2.)
;=== make bar like Pal5 ===
sph_coord = cv_coord(from_rect=[x_init,y_init,z_init],/to_sphere)
longitude = randomn(39,N_stars)*0.000 + sph_coord[0]
latitude = randomn(29,N_stars)*0.000 + sph_coord[1]

rperi = min(r_orbit)
rapo =max(r_orbit)
mvir = 1.77d12
rvir = 389d
rhalo = 24.6d
Msat = 10000.d
cc = rvir/rhalo
mhalo = mvir/(alog(cc+1.)-cc/(cc+1.))
p = rperi/rhalo
mperi = mhalo*(alog(p+1.)-p/(p+1))
p = rapo/rhalo
mapo = mhalo*(alog(p+1)-p/(p+1))
s = (Msat/mperi)^(1./3.)
vperi = max(v_orbit)
vapo = min(v_orbit)

alpha = (findgen(1500)+1)/600+0.5d
beta = 2.5*alpha/3.d
alpha = [alpha,reverse(-alpha)]
beta = [beta,reverse(-beta)]
jay = s*rperi*vperi
dr = rapo*s*(alpha*G*mperi/rperi-beta*vapo^2.)/(G*mapo/rapo-vapo^2.)
dv = (beta*jay-vapo*dr)/rapo

;=== for the encounter test ===
radial = sph_coord[2]+dr
input_arr = transpose([[longitude],[latitude],[radial]])
rect_coord = cv_coord(from_sphere=input_arr,/to_rect)
x_stars = rect_coord[0,*]
y_stars = rect_coord[1,*]
z_stars = rect_coord[2,*]
help,input_arr,rect_coord,longitude

const = 0.6d
;sigma_v = const*s*sqrt(G*mperi/rperi)
J = vapo*rapo
sigma_v = const*s*J/rapo
V = sqrt(Vx_init[0]^2.+Vy_init[0]^2.+Vz_init[0]^2.)
sigma_v2 = randomn(13,N_stars)*sigma_v
vv = V+dv;+sigma_v2
vvx = vv*Vx_init[0]/V
vvy = vv*Vy_init[0]/V
vvz = vv*Vz_init[0]/V

Vx_stars = vvx+sigma_v2/sqrt(3.);(randomu(13,N_stars)*2.-1.)*sigma_v
Vy_stars = vvy+sigma_v2/sqrt(3.);(randomu(26,N_stars)*2.-1.)*sigma_v
Vz_stars = vvz+sigma_v2/sqrt(3.);(randomu(39,N_stars)*2.-1.)*sigma_v
Vtot = sqrt(Vx_stars^2.+Vy_stars^2.+Vz_stars^2.)

;openw, 2, 'Pal5_sigma.dat'
;printf, 2, sigma_v,stddev(Vtot)
;close, 2

openw, 1, 'frog_Pal5.dat'
printf, 1, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init, f='(6(f15.9,1x))'
for i=0L,N_stars-1 do printf, 1, x_stars[i], y_stars[i], z_stars[i], Vx_stars[i], Vy_stars[i], Vz_stars[i], f='(6(f15.9,1x))'
close, 1

!p.multi=[0,2,2]
plot,x_stars,y_stars,psym=1,/isotropic
plot,y_stars,z_stars,psym=1,/isotropic
plot,x_stars,z_stars,psym=1,/isotropic

END
