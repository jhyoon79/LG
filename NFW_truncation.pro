pro NFW_truncation
; check the truncated NFW halo profiles
;  McMillan & Dehnen 2007MNRAS.378.541
;  Springel & White 1999MNRAS.307.162
n = 1000
r = findgen(n)+0.1
Mvir = 1.77d12
rvir = 389.d
rs = 24.6d
c = rvir/rs
delta_c = 200./3.*c^3./(alog(1.+c)-c/(1.+c))
H = 73. /1000. ; km/s/kpc
G = 4.3e-6 ; (m/s)^2 M^-1 kpc
rho_cri = 3.*H^2./(8.*!pi*G)
rho = delta_c / ((r/rs)*(1.+r/rs)^2.) * rho_cri
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
print,Mhalo,Mvir
phi0 = G*Mhalo/rs


;=== check taylor expansion of sech ===
x = r/rvir
y_anal = 1./cosh(x)
EulerNumber = dblarr(19)
EulerNumber[0] = 1.
EulerNumber[2] = -1.
EulerNumber[4] = 5.
EulerNumber[6] = -61
EulerNumber[8] = 1385
EulerNumber[10] = -50521
EulerNumber[12] = 2702765
EulerNumber[14] = -199360981
EulerNumber[16] = 19391512145
EulerNumber[18] = -2404879675441
;`print,EulerNumber
y_nu = dblarr(N_elements(x))
for i=0,N_elements(EulerNumber)-1 do begin
  term = EulerNumber[i]/factorial(i)*x^i
  y_nu += term
endfor
diff = (y_anal-y_nu)/y_anal
;plothist,diff,bin=stddev(diff)/5.
;========================================

r_t = rvir

;x = r/r_t
;sub_trun1 = where(x le 1./2.,complement=sub_trun2)
;truncation_function = dblarr(N_elements(x))
;truncation_function[sub_trun1] = 1.-x[sub_trun1]^2.
;truncation_function[sub_trun2] = 1./(16./9.*x[sub_trun2]+4./9.)
;rho_mine = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*truncation_function

rho_Mc = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*1./cosh(r/r_t)
rho_Mc2 = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*y_nu
;rho_mine = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*(atan((rvir-r)/45.)/!dpi+0.5)

a = c - (1.+3.*c)/(1.+c)
rho_Spr = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*(r^(-r))  *exp(-1.*(r-rvir)/rs)
rho_exp = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*exp(-1.*r/r_t)
sub = where(r gt rvir)

Mbnd = Mvir
rte = 1.02+1.38*alog(Mbnd)+0.37*(alog(Mbnd))^2.
ft = -0.007+0.35*alog(Mbnd)+0.39*(alog(Mbnd))^2.+0.23*(alog(Mbnd))^3.
rho_f = delta_c*rho_cri/((r/rs)*(1.+r/rs)^2.)*ft/(1.+(r/rte)^3.)/ (ft[0]/(1.+(r[0]/rte[0])^3.))
;print, rte,ft,rvir

mass_dr = dblarr(n)
mr_numerical = dblarr(n)
mass_dr_Mc = dblarr(n)
mr_numerical_Mc = dblarr(n)
for i=0L,n-1 do begin
  mass_dr[i] = (i eq 0) ? rho[i]*4./3.*!dpi*r[i]^3. : rho[i]*4./3.*!dpi*(r[i]^3.-r[i-1]^3.)
  mr_numerical[i] = total(mass_dr[0:i])
  mass_dr_Mc[i] = (i eq 0) ? rho_Mc[i]*4./3.*!dpi*r[i]^3. : rho_Mc[i]*4./3.*!dpi*(r[i]^3.-r[i-1]^3.)
  mr_numerical_Mc[i] = total(mass_dr_Mc[0:i])
endfor

set_plot,'ps'
@plot_setting
device,file='NFW_truncation.ps',/color,ysize=22,yoffset=1
;plot,r,truncation_function,/xlog

!p.multi=[0,1,2]
plot,r,rho,xtitle='r [kpc]',ytitle=textoidl('\rho / M'+sunsymbol()+'kpc^{-3}'),/xlog,/ylog
vline,rvir,1e-3,1e10,linestyle=2
oplot,r,rho_Mc,color=70
;oplot,r,rho_Mc2,color=255
;oplot,r[sub],rho_Spr[sub],color=255
;oplot,r,rho_exp,color=220
;oplot,r,rho_f,color=150
;oplot,r,rho_mine,color=30

p = r/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
mr22 = Mhalo*(alog(c+1.)-c/(c+1.))
forprint, r,mr,textout=2
print,mr22

mr_t = Mhalo*(alog(rvir/rs+1.)-rvir/rs/(rvir/rs+1.))
mr_Mc = Mhalo*(alog(p+1.)-p/(p+1.))/cosh(r/r_t)
mr_tanh = Mhalo*(alog(p+1.)-p/(p+1.)) *ft/(1.+(r/rte)^3.)/ (ft[0]/(1.+(r[0]/rte[0])^3.))

dr_smooth = 10.
sigma = dr_smooth*5.
Gaussian = Mvir*0.8*exp(-1.*(r-rvir)^2./(2.*sigma^2.))
mr_new = Mvir^(atan((r-rvir)/dr_smooth)/!dpi+0.5) + mr^(atan((rvir-r)/dr_smooth)/!dpi+0.5)
mr_new = mr_new + Gaussian
mr_in = mr*(atan((rvir-r)/dr_smooth)/!dpi+0.5)
mr_out = (atan((r-rvir)/dr_smooth)/!dpi+0.5)*Mvir
;mr_new2 =  (mr_t*(1.-exp(-p/10.)))
mr_new2 =  mr_t*(1.-1./cosh(-r/rvir))
;mr_new2 =  mr_t*(1.-1./cosh(-p/10.))
mr_new = mr_in+mr_out
plot,r,mr,yr=[1e7,1e14],xtitle='r [kpc]',ytitle='M'+sunsymbol(),/xlog,/ylog
vline,rvir,1e5,1e15,linestyle=2
oplot,r,mr_numerical,linestyle=1
oplot,r,mr_numerical_Mc,linestyle=2,color=70
hline,mr_t,1e1,1e15,linestyle=1
;hline,Mhalo,1e1,1e15,linestyle=2
;hline,Mvir,1e1,1e15,linestyle=3
;oplot,r,mr_new,color=255
;oplot,r,mr_in,color=220
;oplot,r,mr_out,color=150
;oplot,r,mr_new2,color=30


plot,r,mr,xr=[10,3000],yr=[3e10,2e12],xtitle='r [kpc]',ytitle='M'+sunsymbol(),/xlog,/ylog
hline,Mvir,1e1,1e15,linestyle=3
;oplot,r,mr_new,color=255
oplot,r,mr_numerical_Mc/7.523e12*Mvir,linestyle=2,color=70
;oplot,r,mr_in,color=220
;oplot,r,mr_out,color=150
vline,rvir-dr_smooth,1,1e15,linestyle=1
vline,rvir+dr_smooth,1,1e15,linestyle=1
;oplot,r,mr_new2,color=30


;=== find truncation correction factor and test it ===
truncate_correction = mr_numerical_Mc/mr_numerical
plot,[0],[0],/nodata,xr=[0.1,10000],yr=[1e7,1e14],xtitle='r [kpc]',ytitle='M'+sunsymbol(),/xlog,/ylog

for i=1,10 do begin
Mvir = 1.77d10
rvir = 89.d
rs = 4.6d * i
c = rvir/rs
delta_c = 200./3.*c^3./(alog(1.+c)-c/(1.+c))
H = 73. /1000. ; km/s/kpc
G = 4.3e-6 ; (m/s)^2 M^-1 kpc
rho_cri = 3.*H^2./(8.*!pi*G)
rho = delta_c / ((r/rs)*(1.+r/rs)^2.) * rho_cri
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
phi0 = G*Mhalo/rs
p = r/rs
mr2 = Mhalo*(alog(p+1.)-p/(p+1.))
oplot,r,mr2*truncate_correction,color=25*i
endfor

hline,mr_t,1e1,1e15,linestyle=1
oplot,r,mr2*truncate_correction,color=255,linestyle=2
oplot,r,mr,color=100

kpc2km = 3.08568e16
a = -G*mr/(r*kpc2km*r)
avir = -G*Mvir/(r*kpc2km*r)
erase  &  !p.multi=0
;plot,r,a,xr=[5,25],yr=[-3e-16,-1.55e-17],/xlog
plot,r,a,xr=[380,400],yr=[-1.7e-15,-1.55e-15],/xlog
oplot,r,avir,linestyle=2
vline,rvir,linestyle=1
legend,['a_en','a_vir','r_vir'],linestyle=[0,2,1],box=0
a_new = avir*(atan((r-rvir)/3.)/!dpi+0.5) + a*(atan((rvir-r)/3.)/!dpi+0.5)
oplot,r,a_new,linestyle=1,color=255


device,/close

END
