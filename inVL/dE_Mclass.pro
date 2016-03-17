pro dE_Mclass

set_plot,'ps'
@plot_setting
device,file='dE_Mclass.ps',/color,/landscape;,xsize=15,ysize=24,xoffset=0.5,yoffset=0.5

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

list = ['all','e5e6','e6e7','e7e8','e8e9','e9e10']
Nlist = N_elements(list)
Etot0 = dblarr(Nlist)
Etot_final = dblarr(Nlist)
for k=0,Nlist-1 do begin
  dir = '/media/SEADISK/LG/inVL_'+list[k]
  if k eq 0 then dir = '/media/SEADISK/LG/inVL_12790'

  chr_rdtbl,dir+'/snapshot/snap00000',0,arr,/silent
  Etot0[k] = double(arr[8,0])

  chr_rdtbl,dir+'/snapshot/snap12790',0,arr,/silent
  arr = double(arr)
  Etot_final[k] = arr[8,0]

  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = mean(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  r_tide = 0.108653
  Msat = 10000.
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
;  delE = E_total-E_total[0]
;  q = delE/epsilon
;  dJ = J-J[0]
;  scale_dJ = dJ / (r_tide/r_peri * J)
  print,k,epsilon
endfor

dE = Etot_final-Etot0
q = dE/epsilon
x = indgen(Nlist)
q = shift(q,-1)
plot,x,q,xtickname=shift(list,-1),yr=[-8,4],xtitle='Msub',ytitle='q'
plot,x,abs(q),xtickname=shift(list,-1),yr=[0,8],xtitle='Msub',ytitle='|q|'
device,/close

END
