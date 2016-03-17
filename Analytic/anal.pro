pro analytic


l = 10	; tail length in kpc
w = 0.3	; width in kpc
r = 12.2	; Galactic radius in kpc
b_imp = findgen(1000)/10.+0.1
Mlump = [1e5,1e6,1e7,1e8,1e9]
Nlump = [300000,30000,3000,300,30]/3.
rs = [0.04,0.11,0.30,0.79,2.11]

Nenc_e5 = b_imp
Nenc_e6 = b_imp
Nenc_e7 = b_imp
Nenc_e8 = b_imp
Nenc_e9 = b_imp
for i=0,N_elements(b_imp)-1 do begin
  Nenc_e5[i] = 4.*!dpi^2.*r*l*max([b_imp[i],w,rs[0]])*Nlump[0]
  Nenc_e6[i] = 4.*!dpi^2.*r*l*max([b_imp[i],w,rs[1]])*Nlump[1]
  Nenc_e7[i] = 4.*!dpi^2.*r*l*max([b_imp[i],w,rs[2]])*Nlump[2]
  Nenc_e8[i] = 4.*!dpi^2.*r*l*max([b_imp[i],w,rs[3]])*Nlump[3]
  Nenc_e9[i] = 4.*!dpi^2.*r*l*max([b_imp[i],w,rs[4]])*Nlump[4]
endfor

set_plot,'ps'
@plot_setting
device,file='freq.eps',/color,/enc
!p.multi=[0,2,1]
plot,[0],[0],/nodata,xr=[0,100],yr=[-1,11],xtitle='b [kpc]',ytitle=textoidl('LOG(N_{enc})')
oplot,b_imp,alog10(Nenc_e5),linestyle=5
oplot,b_imp,alog10(Nenc_e6),linestyle=4
oplot,b_imp,alog10(Nenc_e7),linestyle=3
oplot,b_imp,alog10(Nenc_e8),linestyle=2
oplot,b_imp,alog10(Nenc_e9),linestyle=1

b_imp = 0.1
Nenc = 4.*!dpi^2.*r*l*max([b_imp*w*rs])*Nlump
plot,[0],[0],/nodata,xr=[5,9],yr=[0,11],xtitle=textoidl('M_{subhalo} [M'+odot+']'),ytitle=textoidl('LOG(N_{enc})')
oplot,alog10(Mlump),alog10(Nenc),psym=5
device,/close


v_sub = 200
v_tail = 77
theta = 90*!dtor
G = 4.3e-6
v_enc = sqrt(v_tail^2.+v+sub^2.-2.*v_tail*v_sub*sin(theta))
x = findgen(2001)/100 - 10

device,file='delE.eps',/color,/enc

rs = [0.04,0.11,0.30,0.79,2.11]
Msub = [1e5,1e6,1e7,1e8,1e9]
b_imp = 0.1*rs

plot,[0],[0],/nodata,xr=[-.6,.6],yr=[-4.5,4.5],xtitle='x [kpc]',ytitle=textoidl('\DeltaE / E_0')

for i=0,N_elements(Msub)-1 do begin
  d_enc = sqrt(b_imp[i]^2.+x^2.)
  sub = where(x lt 0)
  d_enc[sub] = -1*d_enc[sub]
  delE = G*Msub[i]/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc)


d_enc = sqrt(b_imp[1]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b2 = G*Msub[1]/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc); + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[2]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b3 = G*Msub[2]/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc); + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[3]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b4 = G*Msub[3]/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc); + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[4]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b5 = G*Msub[4]/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc); + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)

f_nor1 = G*Msub[0]/rs[0]
f_nor2 = G*Msub[1]/rs[1]
f_nor3 = G*Msub[2]/rs[2]
f_nor4 = G*Msub[3]/rs[3]
f_nor5 = G*Msub[4]/rs[4]

oplot,x,delE_b1/f_nor1,
oplot,x,delE_b2/f_nor2,linestyle=1
oplot,x,delE_b3/f_nor3,linestyle=2
oplot,x,delE_b4/f_nor4,linestyle=3
oplot,x,delE_b5/f_nor5,linestyle=4


rs = 0.30
Msub = 1e7
b_imp = [0.2,2,4,6,8,16]*rs
d_enc = sqrt(b_imp[0]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b1 = G*Msub/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)

d_enc = sqrt(b_imp[1]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b2 = G*Msub/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[2]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b3 = G*Msub/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[3]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b4 = G*Msub/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
d_enc = sqrt(b_imp[4]^2.+x^2.)
sub = where(x lt 0)
d_enc[sub] = -1*d_enc[sub]
delE_b5 = G*Msub/(d_enc*v_enc)*(-v_tail*sin(theta)^2.*x/d_enc + 0.5*G*Msub/(d_enc*v_enc)*sin(theta)^4.*x^2/d_enc^2)
f_nor = G*Msub/rs

plot,x,delE_b1/f_nor,xr=[-3,3],yr=[-2.1,2.1],xtitle='x [kpc]',ytitle=textoidl('\DeltaE / E_0')
oplot,x,delE_b2/f_nor,linestyle=1
oplot,x,delE_b3/f_nor,linestyle=2
oplot,x,delE_b4/f_nor,linestyle=3
oplot,x,delE_b5/f_nor,linestyle=4

device,/close

END
