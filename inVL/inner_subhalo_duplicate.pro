function v_vector, x_t, y_t, z_t, Vr, Vt, phi=phi
; PURPOSE: In order to estimate the velocity vector in Galactic rest xyz frame from the observed radial velocity and tangential velocity.
; PARAMETERS:
;	Vt - tangential velocity in the Galactic restframe viewed at the Sun
;	phi - orientation of the proper motion to the Northern Galactic Pole
;	Vr - radial velocity in the Galactic restframe viewed at the Sun, since the position of the Sun is (-8,0,0).

;=== Sun position ini xyz Galactic plane ===
x_sun = 0.d
y_sun = 0.d
z_sun = 0.d

d_3D = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.+(z_t-z_sun)^2.)
d_xy = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.)
Vr_x = Vr * (x_t-x_sun) / d_3D
Vr_y = Vr * y_t / d_3D
Vr_z = Vr * z_t / d_3D
Vt_x = -Vt*sin(phi) * y_t / d_xy - Vt*cos(phi)*z_t/d_3D * (x_t-x_sun)/d_xy
Vt_y = Vt*sin(phi) * (x_t-x_sun) / d_xy - Vt*cos(phi)*z_t/d_3D * y_t/d_xy
Vt_z = Vt*cos(phi)*d_xy/d_3D
V_x = Vr_x + Vt_x
V_y = Vr_y + Vt_y
V_z = Vr_z + Vt_z
Nx = N_elements(V_x)

array = dblarr(3,Nx)
array[0,*] = V_x
array[1,*] = V_y
array[2,*] = V_z
return,array
END

;=====================================
pro inner_subhalo_duplicate

chr_rdtbl,'../frog_VLsubhalo_all.dat',0,arr
Mtidal = double(arr[6,*])

sub1 = where(Mtidal ge 1e5 and Mtidal lt 1e6)
sub2 = where(Mtidal ge 1e6 and Mtidal lt 1e7)
sub3 = where(Mtidal ge 1e7 and Mtidal lt 1e8)
sub4 = where(Mtidal ge 1e8 and Mtidal lt 1e9)
sub5 = where(Mtidal ge 1e9 and Mtidal lt 1e10)
sub6 = where(Mtidal ge 1e10)

chr_rdtbl,'rmin.prt',1,arr1
inner1 = where(double(arr1[sub1]) lt 25)
inner2 = where(double(arr1[sub2]) lt 25)
inner3 = where(double(arr1[sub3]) lt 25)
inner4 = where(double(arr1[sub4]) lt 25)
inner5 = where(double(arr1[sub5]) lt 25)
inner6 = where(double(arr1[sub6]) lt 25)
help,inner6
forprint,arr1[sub1],textout='tmp'

;goto,jump

x = reform(double(arr[0,sub2[inner2]]))
y = reform(double(arr[1,sub2[inner2]]))
z = reform(double(arr[2,sub2[inner2]]))
vx = reform(double(arr[3,sub2[inner2]]))
vy = reform(double(arr[4,sub2[inner2]]))
vz = reform(double(arr[5,sub2[inner2]]))
Mtidal = reform(double(arr[6,sub2[inner2]]))
rs = reform(double(arr[7,sub2[inner2]]))
rtidal = reform(double(arr[8,sub2[inner2]]))
Nsub = N_elements(sub2[inner2])

r = reform(sqrt(x^2.+y^2.+z^2.))
V = reform(sqrt(vx^2.+vy^2.+vz^2.))
Vr = reform((x*vx+y*vy+z*vz)/r)
Vt = sqrt(V^2.-Vr^2.)
; sph[0]: longitude, theta - from x-axis, xy-plane,counterclock
; sph[1]: latitude, from xy-plane
; sph[2]: r
x30000 = x  &  vx30000 = vx
y30000 = y  &  vy30000 = vy
z30000 = z  &  vz30000 = vz
Mtidal30000 = Mtidal
rtidal30000 = rtidal
rs30000 = rs
for k=0,2 do begin
  theta = randomu(11+k,Nsub)*2*!dpi
  phi = (randomu(22+k,Nsub)-0.5)*!dpi
  Vt_phi = randomu(33+k,Nsub)*2*!dpi
  input_arr = dblarr(3,Nsub)
  input_arr[0,*] = theta
  input_arr[1,*] = phi
  input_arr[2,*] = r
  coord_new = cv_coord(from_sphere=input_arr,/to_rect)
  
  xnew = reform(coord_new[0,*])
  ynew = reform(coord_new[1,*])
  znew = reform(coord_new[2,*])
  x30000 = [x30000,xnew]
  y30000 = [y30000,ynew]
  z30000 = [z30000,znew]
  vnew = v_vector(xnew,ynew,znew,Vr,Vt,phi=Vt_phi)
  vx30000 = [vx30000,reform(vnew[0,*])]
  vy30000 = [vy30000,reform(vnew[1,*])]
  vz30000 = [vz30000,reform(vnew[2,*])]
  Mtidal30000 = [Mtidal30000,Mtidal]
  rtidal30000 = [rtidal30000,rtidal]
  rs30000 = [rs30000,rs]
  print,'r',sqrt(x[0]^2+y[0]^2+z[0]^2),sqrt(xnew[0]^2+ynew[0]^2+znew[0]^2)
  print,'v',sqrt(vx[0]^2+vy[0]^2+vz[0]^2),sqrt(vnew[0,0]^2+vnew[1,0]^2+vnew[2,0]^2)
endfor
help,x30000,Mtidal30000

x = reform(double(arr[0,sub1[inner1]]))
y = reform(double(arr[1,sub1[inner1]]))
z = reform(double(arr[2,sub1[inner1]]))
vx = reform(double(arr[3,sub1[inner1]]))
vy = reform(double(arr[4,sub1[inner1]]))
vz = reform(double(arr[5,sub1[inner1]]))
Mtidal = reform(double(arr[6,sub1[inner1]]))
rs = reform(double(arr[7,sub1[inner1]]))
rtidal = reform(double(arr[8,sub1[inner1]]))
Nsub = N_elements(sub1[inner1])
r = reform(sqrt(x^2.+y^2.+z^2.))
V = reform(sqrt(vx^2.+vy^2.+vz^2.))
Vr = reform((x*vx+y*vy+z*vz)/r)
Vt = sqrt(V^2.-Vr^2.)
; sph[0]: longitude, theta - from x-axis, xy-plane,counterclock
; sph[1]: latitude, from xy-plane
; sph[2]: r
x300000 = x  &  vx300000 = vx
y300000 = y  &  vy300000 = vy
z300000 = z  &  vz300000 = vz
Mtidal300000 = Mtidal
rtidal300000 = rtidal
rs300000 = rs
for k=0,47 do begin
  theta = randomu(100+k,Nsub)*2*!dpi
  phi = (randomu(200+k,Nsub)-0.5)*!dpi
  Vt_phi = randomu(300+k,Nsub)*2*!dpi
  input_arr = dblarr(3,Nsub)
  input_arr[0,*] = theta
  input_arr[1,*] = phi
  input_arr[2,*] = r
  coord_new = cv_coord(from_sphere=input_arr,/to_rect)
  
  xnew = reform(coord_new[0,*])
  ynew = reform(coord_new[1,*])
  znew = reform(coord_new[2,*])
  x300000 = [x300000,xnew]
  y300000 = [y300000,ynew]
  z300000 = [z300000,znew]
  vnew = v_vector(xnew,ynew,znew,Vr,Vt,phi=Vt_phi)
  vx300000 = [vx300000,reform(vnew[0,*])]
  vy300000 = [vy300000,reform(vnew[1,*])]
  vz300000 = [vz300000,reform(vnew[2,*])]
  Mtidal300000 = [Mtidal300000,Mtidal]
  rtidal300000 = [rtidal300000,rtidal]
  rs300000 = [rs300000,rs]
endfor
help,x300000,Mtidal300000
help,inner1,inner2,inner3,inner4,inner5
help,sub1,sub2,sub3,sub4,sub5

openw,1,'frog_VLsubhalo_inner.dat'
for i=0L,N_elements(x300000)-1 do printf,1,x300000[i],y300000[i],z300000[i],$
    vx300000[i],vy300000[i],vz300000[i],Mtidal300000[i],$
    rs300000[i],rtidal300000[i],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(x30000)-1 do printf,1,x30000[i],y30000[i],z30000[i],$
    vx30000[i],vy30000[i],vz30000[i],Mtidal30000[i],$
    rs30000[i],rtidal30000[i],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub3)-1 do printf,1,arr[*,sub3[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub4)-1 do printf,1,arr[*,sub4[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner.dat'

openw,1,'frog_VLsubhalo_inner_e5e6all.dat'
for i=0L,N_elements(x300000)-1 do printf,1,x300000[i],y300000[i],z300000[i],$
    vx300000[i],vy300000[i],vz300000[i],Mtidal300000[i],$
    rs300000[i],rtidal300000[i],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner_e5e6all.dat'


openw,1,'frog_VLsubhalo_inner_e6e10.dat'
for i=0L,N_elements(x30000)-1 do printf,1,x30000[i],y30000[i],z30000[i],$
    vx30000[i],vy30000[i],vz30000[i],Mtidal30000[i],$
    rs30000[i],rtidal30000[i],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub3)-1 do printf,1,arr[*,sub3[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub4)-1 do printf,1,arr[*,sub4[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner_e6e10.dat'

jump:

openw,1,'frog_VLsubhalo_inner_e7e10.dat'
for i=0L,N_elements(sub3)-1 do printf,1,arr[*,sub3[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub4)-1 do printf,1,arr[*,sub4[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner_e7e10.dat'

openw,1,'frog_VLsubhalo_inner_e8e10.dat'
for i=0L,N_elements(sub4)-1 do printf,1,arr[*,sub4[i]],f='(6(f14.5,1x),a16,2(f13.5))'
for i=0L,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner_e8e10.dat'

openw,1,'frog_VLsubhalo_inner_e9e10.dat'
for i=0L,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
spawn,'wc -l frog_VLsubhalo_inner_e9e10.dat'



END
