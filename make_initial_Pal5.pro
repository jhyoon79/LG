function v_vector, x_t, y_t, z_t, Vr, Vt, phi=phi
; PURPOSE: In order to estimate the velocity vector in Galactic rest xyz frame from the observed radial velocity and tangential velocity.
; PARAMETERS:
;	Vt - tangential velocity in the Galactic restframe viewed at the Sun
;	phi - orientation of the proper motion to the Northern Galactic Pole
;	Vr - radial velocity in the Galactic restframe viewed at the Sun, since the position of the Sun is (-8,0,0).

;=== Sun position ini xyz Galactic plane ===
x_sun = -8.d
y_sun = 0.d
z_sun = 0.d
;phi = (360.d - phi)*!dtor
phi = phi*!dtor

d_3D = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.+(z_t-z_sun)^2.)
d_xy = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.)
Vr_x = Vr * (x_t-x_sun) / d_3D
Vr_y = Vr * y_t / d_3D
Vr_z = Vr * z_t / d_3D
print,'Vr_x',Vr_x
print,'Vr_y',Vr_y
print,'Vr_z',Vr_z
print,'Vr=',Vr,sqrt(Vr_x^2.+Vr_y^2.+Vr_z^2.)	; correct

;daewon = asin(sin(d_3D*(!dpi/2.-acos(d_xy/d_3D)))*y_t/d_xy)
;a = asin(y_t/d_xy*d_3D*!dpi/2.)
;b = a*d_xy/d_3D
;print,daewon,b
;Vt_x = -Vt*sin(phi) * y_t / d_xy
;Vt_y = Vt*sin(phi) * (x_t-x_sun) / d_xy
;Vt_z = Vt*cos(phi)
Vt_x = -Vt*sin(phi) * y_t / d_xy - Vt*cos(phi)*z_t/d_3D * (x_t-x_sun)/d_xy
Vt_y = Vt*sin(phi) * (x_t-x_sun) / d_xy - Vt*cos(phi)*z_t/d_3D * y_t/d_xy
Vt_z = Vt*cos(phi)*d_xy/d_3D

print,'Vt1', -Vt*sin(phi) * y_t / d_xy
print,'Vt2', - Vt*cos(phi)*z_t/d_3D * (x_t-x_sun)/d_xy
print,'Vt_x',Vt_x
print,'Vt_y',Vt_y
print,'Vt_z',Vt_z
print,'Vt=',Vt,sqrt(Vt_x^2.+Vt_y^2.+Vt_z^2.)	; correct

V_x = Vr_x + Vt_x
V_y = Vr_y + Vt_y
V_z = Vr_z + Vt_z

return, [V_x, V_y, V_z]
END


pro make_initial_Pal5
; PURPOSE: convert the coordinate of the position to cartesian.
; d = 23.2kpc

;=== target position in xyz Galactic plane ===
x_init = 8.2d
y_init = 0.2d  
z_init = 16.6d 
Vr = -44.3d 	; radial velocity in Galactic rest frame viewed at Sun
Vt =  89.d		; old=69km/s  tangential velocity in Galactic rest frame viewed at Sun
V = v_vector(x_init, y_init, z_init, Vr, Vt, phi=280.d)
print, 'Vx, Vy, Vz of Pal5 =', V
Vx_init = V[0]
Vy_init = V[1]
Vz_init = V[2]

print,'VrVt=',sqrt(Vr^2.+Vt^2.),'    Vxyz=',sqrt(total(V^2.))

openw, 1, 'frog_Pal5_initial.dat'
printf, 1, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init, f='(6(f14.9))'
close, 1


END
