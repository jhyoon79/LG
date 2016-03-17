pro nbody_orbit,dir

mu = 1.e4
ru = 0.0075
kmpkpc = 3.240756d-17
Gyr2sec = 3.1536d16
kpc2km = 3.0857d16
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
tu = sqrt((ru*kpc2km)^3./(mu*G*kpc2km))/Gyr2sec
vu = (ru*kpc2km)/tu/Gyr2sec

dir = '/media/SEADISK/LG/Nbody/'+dir

chr_rdtbl,dir+'SCFCEN',1,arr
arr = double(arr)
orbit_t = arr[0,*]
orbit_dt = arr[1,*]
orbit_x = arr[2,*] * ru
orbit_y = arr[3,*] * ru
orbit_z = arr[4,*] * ru
orbit_vx = arr[5,*] * vu
orbit_vy = arr[6,*] * vu
orbit_vz = arr[7,*] * vu
orbit_r = sqrt(orbit_x^2.+orbit_y^2.+orbit_z^2.)

print,'r_peri=',min(orbit_r)
print,'r_apo=',max(orbit_r)


spawn,'ls '+dir+'SNAP* > tmp5'
chr_rdtbl,'tmp5',0,arr,/silent
snap_fname = reform(arr[0,*])
spawn,'rm -f tmp5'

for i=180,N_elements(snap_fname)-1 do begin
;or i=196,197 do begin
  readcol,snap_fname[i],n,t,numline=1,/silent
  chr_rdtbl,snap_fname[i],1,arr,/silent
  arr = double(arr)
  x = reform(arr[1,*])*ru
  y = reform(arr[2,*])*ru
  z = reform(arr[3,*])*ru
  Vx = reform(arr[4,*])*vu
  Vy = reform(arr[5,*])*vu
  Vz = reform(arr[6,*])*vu
  dE = reform(arr[7,*])*vu^2
  Etot = reform(arr[8,*])*vu^2

  l = atan(y/(x+8.))*!radeg 
  b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
  lcosb = l*cos(b*!dtor)
  Dec = asin(cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor))*!radeg
  RA = asin((cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(Dec*!dtor))*!radeg + 282.25
  !p.multi=[0,2,2]
  sub_orbit = where(orbit_t le t[0]+orbit_dt[0]/2.,count)

  xcen = median(x)
  ycen = median(y)
  RAcen = median(RA)
  Deccen = median(Dec)
  plot,orbit_x[sub_orbit],orbit_y[sub_orbit],xr=[-20,20],yr=[-20,20],/isotropic
  oplot,[orbit_x[count-1]],[orbit_y[count-1]],psym=4
  plot,x-xcen,y-ycen,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
  plot,RA-RAcen,Dec-Deccen,psym=3,xr=[-15,15],yr=[-5,5],title=strmid(snap_fname[i],29,7),pos=[0.1,0.1,0.8,0.5]
wait,.5
endfor

stop
END
