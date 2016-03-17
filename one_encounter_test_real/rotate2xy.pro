pro rotate2xy

dir_out = '/scratch/jhyoon/Research/LG/realistic/'
dir_part = '/scratch/jhyoon/Research/LG/realistic/snapshot/'

G = 4.3d-6
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
chr_rdtbl,dir_out+'part_peri',0, arr
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
V_peri = sqrt(Vx_peri^2.+Vy_peri^2.+Vz_peri^2.)
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
Msat = 10000.
s = (Msat/mr)^(1./3.)
epsilon = s*(4.3e-6*mr/r_peri)

;=== readout data for the center of Pal5 ===
chr_rdtbl, dir_part+'snap00000',0, arr
arr = double(arr)
Etot_initial = arr[8,*]
J_initial = arr[9,*]
delE = Etot_initial-Etot_initial[0]
q_initial = delE/epsilon
dJ_initial = J_initial-J_initial[0]
scale_dJ_initial = dJ_initial / (s * J_initial)
qmin = min(q_initial)
qmax = max(q_initial)
dJmin = min(scale_dJ_initial)
dJmax = max(scale_dJ_initial)
qcolor = round((q_initial-qmin)/(qmax-qmin)*245.)+10
N_particle = N_elements(qcolor)


;=== equation of the orbital plane ===
;chr_rdtbl,dir_part+'snap00000',0,arr,/silent
;arr = double(arr)
;x1 = reform(arr[1,0])
;y1 = reform(arr[2,0])
;z1 = reform(arr[3,0])
;chr_rdtbl,dir_part+'snap00100',0,arr,/silent
;arr = double(arr)
;x2 = reform(arr[1,0])
;y2 = reform(arr[2,0])
;z2 = reform(arr[3,0])
;chr_rdtbl,dir_part+'snap00200',0,arr,/silent
;arr = double(arr)
;x3 = reform(arr[1,0])
;y3 = reform(arr[2,0])
;z3 = reform(arr[3,0])
chr_rdtbl,dir_out+'part001',0,arr,/silent
arr = double(arr)
x1 = reform(arr[1,0])
y1 = reform(arr[2,0])
z1 = reform(arr[3,0])
x2 = reform(arr[1,100])
y2 = reform(arr[2,100])
z2 = reform(arr[3,100])
x3 = reform(arr[1,200])
y3 = reform(arr[2,200])
z3 = reform(arr[3,200])
A = determ(transpose([[1,1,1],[y1,y2,y3],[z1,z2,z3]]))
B = determ(transpose([[x1,x2,x3],[1,1,1],[z1,z2,z3]]))
C = determ(transpose([[x1,x2,x3],[y1,y2,y3],[1,1,1]]))
D = -determ(transpose([[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]))
D=0
theta = acos(A/sqrt(A^2.+B^2.))
phi = acos(C/sqrt(A^2.+B^2.+C^2.))
cos_theta = cos(theta)
sin_theta = sin(theta)
cos_phi = cos(phi)
sin_phi = sin(phi)

;=== define & read out the tail data ===
spawn, 'ls '+dir_part+' > temp9'
chr_rdtbl,'temp9',0,fname_part
fname_part = dir_part+reform(fname_part)
spawn, 'rm -f temp9'

openw,1,'frog_rotate2xy_5.5Gyr'
!p.multi=[0,2,2]
N_time = N_elements(fname_part)
for i=5400,6050,10 do begin
;for i=0,200,10 do begin
;for i=0,N_time-1,10 do begin

;=== read out data of the tail particles ===
  chr_rdtbl,fname_part[i],0,arr,/silent
  arr = double(arr)
  t_part = reform(arr[0,*])
  x_part = reform(arr[1,*])
  y_part = reform(arr[2,*])
  z_part = reform(arr[3,*])
  Vx_part = reform(arr[4,*])
  Vy_part = reform(arr[5,*])
  Vz_part = reform(arr[6,*])
  dE_part = reform(arr[7,*])
  Etot_part = reform(arr[8,*])
  J_part = reform(arr[9,*])
  delE = Etot_part-Etot_part[0]
  q_part = delE/epsilon
  dJ_part = J_part-J_part[0]
  scale_dJ_part = dJ_part / (s * J_part)
  t_step = string(t_part[0]/1000.,f='(f6.4)')
  if (i/50. eq i/50) then print,t_step+'Gyr'

plot,x_part,y_part,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
plot,y_part,z_part,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
r1=sqrt(x_part[0]^2+y_part[0]^2+z_part[0]^2)
;=== rotate the coord. ===
; | cos(-theta)  sin(-theta) | 
; | sin(-theta) -cos(-theta) | 
  xx = cos_theta*x_part-sin_theta*y_part
  yy = -sin_theta*x_part-cos_theta*y_part
  zz = z_part
  xxx = cos_phi*xx-sin_phi*zz
  yyy = yy
  zzz = -sin_phi*xx-cos_phi*zz
  x_part = xxx
  y_part = yyy
  z_part = zzz

  vxx = cos_theta*vx_part-sin_theta*vy_part
  vyy = -sin_theta*vx_part-cos_theta*vy_part
  vzz = vz_part
  vxxx = cos_phi*vxx-sin_phi*vz_part
  vyyy = vyy
  vzzz = -sin_phi*vxx-cos_phi*vzz
  vx_part = vxxx
  vy_part = vyyy
  vz_part = vzzz

if t_part[0] eq 5500. then for j=0,N_particle-1 do printf,1,t_part[j],x_part[j],y_part[j],z_part[j],vx_part[j],vy_part[j],vz_part[j],f='(i5.1,6(f14.8,1x))'
plot,x_part,y_part,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
plot,y_part,z_part,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
r2=sqrt(x_part[0]^2+y_part[0]^2+z_part[0]^2)
print,r1,r2
endfor
close,1

END
