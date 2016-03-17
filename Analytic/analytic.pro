pro analytic


l = 10	; tail length in kpc
w = 0.3	; width in kpc
r = 12.2	; Galactic radius in kpc
b_imp = findgen(1000)/10.+0.1
Msub = [1e5,1e6,1e7,1e8,1e9]
nlump = [300000,30000,3000,300,30]/(4/3.*!pi*100^3.)
rs = [0.04,0.11,0.30,0.79,2.11]

set_plot,'ps'
@plot_setting
device,file='freq.eps',/color,/enc,/landscape
!p.multi=[0,2,1]

;=== frequency description ===
plot,[0],[0],/nodata,xr=[0,100],yr=[-1,3],xtitle='b [kpc]',ytitle=textoidl('LOG(N_{enc})')
Pmax = b_imp
for i=0,N_elements(rs)-1 do begin
  tmp = max([w,rs[i]])
  sub = where(b_imp gt tmp,complement=other)
; Pmax equals max([b_imp,w,rs[i]])
  Pmax[sub] = b_imp[sub]
  Pmax[other] = tmp
  Nenc = 4.*!dpi^2.*r*l*Pmax*nlump[i]
  oplot,b_imp,alog10(Nenc),linestyle=i
endfor

plot,[0],[0],/nodata,xr=[4.9,9.1],yr=[-0.1,11],xtitle=textoidl('Log(M_{subhalo}) [M'+odot+']'),ytitle=textoidl('N_{enc}')
b_imp = 0.1
for i=0,N_elements(Msub)-1 do begin
  Nenc = 4.*!dpi^2.*r*l*max([b_imp,w,rs[i]])*nlump[i]
  oplot,[alog10(Msub[i])],[Nenc],psym=5
endfor
device,/close


;=== analytic dE description ===
v_sub = 200
v_tail = 77
theta = 90*!dtor
G = 4.3e-6
v_enc = sqrt(v_tail^2.+v_sub^2.-2.*v_tail*v_sub*sin(theta))
x = findgen(2001)/100 - 10

device,file='delE.eps',/color,/enc,/landscape

b_imp = 0.1*rs
plot,[0],[0],/nodata,xr=[-10,10],yr=[-0.2,0.2],xtitle='x [kpc]',ytitle=textoidl('\DeltaE / (GM/r_s)')

for i=0,N_elements(Msub)-1 do begin
  d_enc = sqrt(b_imp[i]^2.+x^2.)
  p = d_enc/rs[i]
  mr = Msub[i]*(alog(p+1.)-p/(p+1.))
  sub = where(x lt 0)
  d_enc[sub] = -1*d_enc[sub]
  delE = G*mr/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc)
  delE[where(x eq 0)] = 0
  f_nor = G*Msub[i]/rs[i]
  oplot,x,delE/f_nor,linestyle=i
endfor


rs = rs[3]
Msub = Msub[3]
b_imp = [0,2,4,6,8,16]*rs
plot,[0],[0],/nodata,xr=[-10,10],yr=[-0.2,0.2],xtitle='x [kpc]',ytitle=textoidl('\DeltaE / (GM/r_s)')

for i=0,N_elements(b_imp)-1 do begin
  d_enc = sqrt(b_imp[i]^2.+x^2.)
  p = d_enc/rs
  mr = Msub*(alog(p+1.)-p/(p+1.))
  sub = where(x lt 0)
  d_enc[sub] = -1*d_enc[sub]
  delE = G*mr/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc)
  delE[where(x eq 0)] = 0
  f_nor = G*Msub/rs
  oplot,x,delE/f_nor,linestyle=i
endfor

device,/close

END
