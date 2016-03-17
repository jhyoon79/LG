function v_vector, x_t, y_t, z_t, Vr, Vt, phi=phi
; PURPOSE: In order to estimate the velocity vector in Galactic rest xyz frame from the observed radial velocity and tangential velocity.
; PARAMETERS:
;	Vt - tangential velocity in the Galactic restframe viewed at the Sun
;	phi - orientation of the proper motion to the Northern Galactic Pole
;	Vr - radial velocity in the Galactic restframe viewed at the Sun, since the position of the Sun is (-8,0,0).

;=== Sun position ini xyz Galactic plane ===
x_sun = -8.
y_sun = 0.
z_sun = 0.
phi = (360.- phi)*!dtor

d_3D = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.+(z_t-z_sun)^2.)
d_xy = sqrt((x_t-x_sun)^2.+(y_t-y_sun)^2.)
Vr_x = Vr * (x_t-x_sun) / d_3D
Vr_y = Vr * y_t / d_3D
Vr_z = Vr * z_t / d_3D
Vt_x = Vt * sin(phi) * y_t / d_xy
Vt_y = -1. * Vt * sin(phi) * (x_t-x_sun) / d_xy
Vt_z = Vt * cos(phi)
print, 'V3',Vt_x,Vt_y,Vt_z

V_x = Vr_x + Vt_x
V_y = Vr_y + Vt_y
V_z = Vr_z + Vt_z

return, [V_x, V_y, V_z]
END

function v_vector2, x_t, y_t, z_t, Vr, Vt, phi=phi
; PURPOSE: In order to estimate the velocity vector in Galactic rest xyz frame from the observed radial velocity and tangential velocity.
; PARAMETERS:
;	Vt - tangential velocity in the Galactic restframe viewed at the Galactic center
;	Vr - radial velocity in the Galactic restframe viewed at the Galactic center since the position of the Sun is (-8,0,0).

;=== Sun position ini xyz Galactic plane ===
x_sun = -8.
y_sun = 0.
z_sun = 0.

d_3D = sqrt(x_t^2.+y_t^2.+z_t^2.)
d_xy = sqrt(x_t^2.+y_t^2.)

Vr_x = Vr * x_t / d_3D
Vr_y = Vr * y_t / d_3D
Vr_z = Vr * z_t / d_3D

Vt_x = Vt * sin(phi) * y_t / d_xy              ;232
Vt_y = -1. * Vt * sin(phi) * x_t / d_xy;0
Vt_z = Vt * cos(phi)                           ;194

V_x = Vr_x + Vt_x
V_y = Vr_y + Vt_y
V_z = Vr_z + Vt_z

return, [V_x, V_y, V_z]
END

function GC_coord,l,b,d,x_sun=x_sun
; PURPOSE: convert l,b coord. into Galactic cartesian coord.
; PARAMETERS:
;	l - Galactic longitude
;	b - Galactic latitude
;	d - heliocentric distance

if not keyword_set(x_sub) then x_sun = -8.
Xgc = d*cos(b)*cos(l) + x_sun
Ygc = d*cos(b)*sin(l)
Zgc = d*sin(b)

return, [Xgc,Ygc,Zgc]
END


pro make_initial_sgr
; PURPOSE: convert the coordinate of the position to cartesian.
; Sagittarius case
;	 v_tangential = km/s
;	 d = 26.3kpc	: Monaco et al. 2004MNRAS.353.874

l_sgr = 5.6d
b_sgr = -14.2d
d_sgr = 24.d
xyz = GC_coord(l_sgr,b_sgr,d_sgr)
print,xyz
;=== target position in xyz Galactic plane ===
x_init = 16.2;xyz[0]
y_init = 2.3;xyz[1]
z_init = -5.9;xyz[2]
Vr = 140.3	; using Ibata97 eq1, vgal = 171km/s
;171;140.3	; radial velocity in Galactic rest frame viewed at Sun
Vt = 250	; from Ibata97
;70		; tangential velocity in Galactic rest frame viewed at Sun
V = v_vector2(x_init, y_init, z_init, Vr, Vt, phi=0.)
print, 'Vx, Vy, Vz of Pal5 =', V
Vx_init = 215;V[0]
Vy_init = -5;V[1]
Vz_init = 205;V[2]

openw, 1, 'frog_sgr_initial.dat'
printf, 1, x_init, y_init, z_init, Vx_init, Vy_init, Vz_init, f='(6(f14.9))'
close, 1

END


