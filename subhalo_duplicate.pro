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


pro subhalo_duplicate

chr_rdtbl,'frog_VLsubhalo_all.dat',0,arr
Mtidal = double(arr[6,*])

sub1 = where(Mtidal ge 1e5 and Mtidal lt 1e6)
sub2 = where(Mtidal ge 1e6 and Mtidal lt 1e7)
sub3 = where(Mtidal ge 1e7 and Mtidal lt 1e8)
sub4 = where(Mtidal ge 1e8 and Mtidal lt 1e9)
sub5 = where(Mtidal ge 1e9 and Mtidal lt 1e10)
sub6 = where(Mtidal ge 1e10)

openw,1,'frog_VLsubhalo_e5e6.dat'
for i=0,N_elements(sub1)-1 do printf,1,arr[*,sub1[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
openw,1,'frog_VLsubhalo_e6e7.dat'
for i=0,N_elements(sub2)-1 do printf,1,arr[*,sub2[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
openw,1,'frog_VLsubhalo_e7e8.dat'
for i=0,N_elements(sub3)-1 do printf,1,arr[*,sub3[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
openw,1,'frog_VLsubhalo_e8e9.dat'
for i=0,N_elements(sub4)-1 do printf,1,arr[*,sub4[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
openw,1,'frog_VLsubhalo_e9e10.dat'
for i=0,N_elements(sub5)-1 do printf,1,arr[*,sub5[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1
openw,1,'frog_VLsubhalo_e10.dat'
for i=0,N_elements(sub6)-1 do printf,1,arr[*,sub6[i]],f='(6(f14.5,1x),a16,2(f13.5))'
close,1


test = 'on'
x = reform(double(arr[0,sub2]))
y = reform(double(arr[1,sub2]))
z = reform(double(arr[2,sub2]))
vx = reform(double(arr[3,sub2]))
vy = reform(double(arr[4,sub2]))
vz = reform(double(arr[5,sub2]))
Mtidal = reform(double(arr[6,sub2]))
rs = reform(double(arr[7,sub2]))
rtidal = reform(double(arr[8,sub2]))
Nsub = N_elements(sub2)

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
  print,'r',sqrt(x[0]^2+y[0]^2+z[0]^2),sqrt(xnew[0]^2+ynew[0]^2+znew[0]^2)
  print,'v',sqrt(vx[0]^2+vy[0]^2+vz[0]^2),sqrt(vnew[0,0]^2+vnew[1,0]^2+vnew[2,0]^2)
endfor

Mtidal30000 = [Mtidal,Mtidal,Mtidal,Mtidal]
rs30000 = [rs,rs,rs,rs]
rtidal30000 = [rtidal,rtidal,rtidal,rtidal]
help,x,x30000

help,x30000,vx30000,rtidal30000
openw,1,'frog_VLsubhalo_e6e7_30000.dat'
for i=0L,30000-1 do printf,1,x30000[i],y30000[i],z30000[i],$
    vx30000[i],vy30000[i],vz30000[i],Mtidal30000[i],$
    rs30000[i],rtidal30000[i],f='(6(f14.5,1x),a16,2(f13.5))'
close,1

; check duplicate coord.
testr = sqrt(x30000^2.+y30000^2.+z30000^2.)
testv = sqrt(vx30000^2.+vy30000^2.+vz30000^2.)
print,x30000[0],x30000[Nsub],x30000[2*Nsub],x30000[3*Nsub]
print,testr[0],testr[Nsub],testr[2*Nsub],testr[3*Nsub]
print,testv[0],testv[Nsub],testv[2*Nsub],testv[3*Nsub]


@plot_setting
device,file='subhalo_duplicate.ps',/color
!p.multi=[0,2,2]
plothist,r,bin=100,/peak,xtitle='r'
plothist,testr,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,x,bin=100,/peak,xtitle='x'
plothist,x30000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,y,bin=100,/peak,xtitle='y'
plothist,y30000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,z,bin=100,/peak,xtitle='z'
plothist,z30000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,V,bin=10,/peak,xtitle='v'
plothist,testv,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vx,bin=10,/peak,xtitle='vx'
plothist,vx30000,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vy,bin=10,/peak,xtitle='vy'
plothist,vy30000,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vz,bin=10,/peak,xtitle='vz'
plothist,vz30000,bin=10,/peak,/overplot,color=255,linestyle=1


x = reform(double(arr[0,sub1]))
y = reform(double(arr[1,sub1]))
z = reform(double(arr[2,sub1]))
vx = reform(double(arr[3,sub1]))
vy = reform(double(arr[4,sub1]))
vz = reform(double(arr[5,sub1]))
Mtidal = reform(double(arr[6,sub1]))
rs = reform(double(arr[7,sub1]))
rtidal = reform(double(arr[8,sub1]))
Nsub = N_elements(sub1)

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
endfor

Mtidal300000 = [Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal,Mtidal]
rs300000 = [rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs,rs]
rtidal300000 = [rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal,rtidal]
help,x,Mtidal300000

; check duplicate coord.
testr = sqrt(x300000^2.+y300000^2.+z300000^2.)
testv = sqrt(vx300000^2.+vy300000^2.+vz300000^2.)
print,x300000[0],x300000[Nsub],x300000[18*Nsub],x300000[36*Nsub]
print,testr[0],testr[Nsub],testr[18*Nsub],testr[36*Nsub]
print,testv[0],testv[Nsub],testv[18*Nsub],testv[36*Nsub]

openw,1,'frog_VLsubhalo_e5e6_300000.dat'
for i=0L,300000-1 do printf,1,x300000[i],y300000[i],z300000[i],$
    vx300000[i],vy300000[i],vz300000[i],Mtidal300000[i],$
    rs300000[i],rtidal300000[i],f='(6(f14.5,1x),a16,2(f13.5))'
close,1

plothist,r,bin=100,/peak,xtitle='r'
plothist,testr,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,x,bin=100,/peak,xtitle='x'
plothist,x300000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,y,bin=100,/peak,xtitle='y'
plothist,y300000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,z,bin=100,/peak,xtitle='z'
plothist,z300000,bin=100,/peak,/overplot,color=255,linestyle=1
plothist,V,bin=10,/peak,xtitle='v'
plothist,testv,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vx,bin=10,/peak,xtitle='vx'
plothist,vx300000,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vy,bin=10,/peak,xtitle='vy'
plothist,vy300000,bin=10,/peak,/overplot,color=255,linestyle=1
plothist,vz,bin=10,/peak,xtitle='vz'
plothist,vz300000,bin=10,/peak,/overplot,color=255,linestyle=1

END
