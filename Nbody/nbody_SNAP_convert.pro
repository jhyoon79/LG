pro nbody_SNAP_convert

dir_nbody = '/scratch/jhyoon/Research/LG/Nbody_new/pal5/ten4_new/'
;=== N-body ===
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
mu = 1.e4
ru = 0.0075
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec

openr,1,dir_nbody+'SNAP081'
readf,1,Npart,time
close,1
t_nbody = round(time*tu*1000.)

chr_rdtbl,dir_nbody+'SNAP081',1,arr,/silent
arr = double(arr)
m_nbody = reform(arr[0,*])*mu
x_nbody = reform(arr[1,*])*ru
y_nbody = reform(arr[2,*])*ru
z_nbody = reform(arr[3,*])*ru
vx_nbody = reform(arr[4,*])*vu
vy_nbody = reform(arr[5,*])*vu
vz_nbody = reform(arr[6,*])*vu
pot_int_nbody = reform(arr[7,*])*vu^2.
pot_ext_nbody = reform(arr[8,*])*vu^2.
tunbound_nbody = reform(arr[9,*])

chr_rdtbl,dir_nbody+'SCFCEN',1,arr
arr = double(arr)
orbit_t = round(arr[0,*]*tu*1000)
orbit_dt = arr[1,*] * tu
orbit_x = arr[2,*] * ru
orbit_y = arr[3,*] * ru
orbit_z = arr[4,*] * ru
orbit_vx = arr[5,*] * vu
orbit_vy = arr[6,*] * vu
orbit_vz = arr[7,*] * vu
sub12796 = where(orbit_t eq 12796)

openw,1,'nbody_snap'+strtrim(t_nbody,2)
for i=0,N_elements(x_nbody)-1 do printf,1,m_nbody[i],x_nbody[i],y_nbody[i],z_nbody[i],vx_nbody[i],vy_nbody[i],vz_nbody[i],pot_int_nbody[i],pot_ext_nbody[i],tunbound_nbody[i],f='(10(e14.6))'
close,1
openw,1,'nbody_SCFCEN'+strtrim(t_nbody,2)
printf,1,orbit_t[sub12796],orbit_x[sub12796],orbit_y[sub12796],orbit_z[sub12796],orbit_vx[sub12796],orbit_vy[sub12796],orbit_vz[sub12796],f='(7(e14.6))'
close,1

END
