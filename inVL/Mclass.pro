pro Mclass

set_plot,'ps'
@plot_setting
!p.charsize=1.3
!p.multi=[0,4,3]
device,file='Mclass.ps',/color,/landscape;,xsize=15,ysize=24,xoffset=0.5,yoffset=0.5

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

;list = ['all','e8e9','e9e10']
list = ['all','e5e6','e6e7','e7e8','e8e9','e9e10']
for k=1,5 do begin
  dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[k]
;  if k eq 0 then dir = '/media/SEADISK/LG/inVL_12790'
  chr_rdtbl,dir+'/snapshot/snap08440',0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  E_total = arr[8,*]
  J = arr[9,*]

  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  r_tide = 0.108653
  Msat = 10000.
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
  delE = E_total-E_total[0]
  q = delE/epsilon
  dJ = J-J[0]
  scale_dJ = dJ / (r_tide/r_peri * J)
  plot,x,y,psym=3,xr=[-20,20],yr=[-20,20],/isotropic,title=list[k]
  plot,q,scale_dJ,xr=[-2.5,2.5],yr=[-2,2],psym=3
endfor
multiplot,/reset
erase

; show a tail with close encouners
!p.multi=[0,4,4]
for k=1,5 do begin
  dir = '/media/SEADISK/LG/inVL_'+list[k]
  chr_rdtbl,dir+'/snapshot_subhalo/subhalo00000',0,arr,/silent
  Nsub = N_elements(arr[0,*])
  chr_rdtbl,dir+'/snapshot/snap12790',0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x = arr[1,*]
  y = arr[2,*]
  z = arr[3,*]
  E_total = arr[8,*]
  J = arr[9,*]

  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  r_tide = 0.108653
  Msat = 20000.
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
  delE = E_total-E_total[0]
  q = delE/epsilon
  dJ = J-J[0]
  scale_dJ = dJ / (r_tide/r_peri * J)

;  plot,x,y,psym=3,xr=[-20,20],yr=[-20,20],/isotropic,title=list[k]+' '+strtrim(Nsub,2)
  plot,x,y,psym=3,/isotropic,title=list[k]+' '+strtrim(Nsub,2)

  dir10p = dir+'_close_10p'
  dir5p = dir+'_close_5p'
  dir2p = dir+'_close_2p'
  if k eq 4 then begin
    dir10p = dir+'_close_50p'
    dir5p = dir+'_close_40p'
    dir2p = dir+'_close_30p'
  endif else if k eq 5 then begin
    dir10p = dir+'_close_20'
    dir5p = dir+'_close_15'
    dir2p = dir+'_close_10'
  endif
  chr_rdtbl,dir10p+'/snapshot/snap12790',0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x_10p = arr[1,*]
  y_10p = arr[2,*]
  z_10p = arr[3,*]
  E_total = arr[8,*]
  J = arr[9,*]
  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
  delE = E_total-E_total[0]
  q_10p = delE/epsilon
  dJ = J-J[0]
  scale_dJ_10p = dJ / (r_tide/r_peri * J)

  chr_rdtbl,dir5p+'/snapshot/snap12790',0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x_5p = arr[1,*]
  y_5p = arr[2,*]
  z_5p = arr[3,*]
  E_total = arr[8,*]
  J = arr[9,*]
  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
  delE = E_total-E_total[0]
  q_5p = delE/epsilon
  dJ = J-J[0]
  scale_dJ_5p = dJ / (r_tide/r_peri * J)

  chr_rdtbl,dir2p+'/snapshot/snap12790',0,arr,/silent
  arr = double(arr)
  t = string(arr[0,0]/1000.,f='(f6.3)')
  x_2p = arr[1,*]
  y_2p = arr[2,*]
  z_2p = arr[3,*]
  E_total = arr[8,*]
  J = arr[9,*]
  chr_rdtbl,dir+'/part_peri',0, arr,/silent
  arr = double(arr)
  t_peri = reform(arr[0,*])
  x_peri = reform(arr[1,*])
  y_peri = reform(arr[2,*])
  z_peri = reform(arr[3,*])
  r_peri = sqrt(x_peri^2.+y_peri^2.+z_peri^2.)
  p = r_peri/rs
  mr = Mhalo*(alog(p+1.)-p/(p+1.))
  epsilon = (Msat/mr)^(1./3.)*(4.3e-6*mr/r_peri)
  delE = E_total-E_total[0]
  q_2p = delE/epsilon
  dJ = J-J[0]
  scale_dJ_2p = dJ / (r_tide/r_peri * J)

  chr_rdtbl,dir10p+'/snapshot_subhalo/subhalo00000',0,arr,/silent
  N10p = ' '+strtrim(N_elements(arr[0,*]),2)
  chr_rdtbl,dir5p+'/snapshot_subhalo/subhalo00000',0,arr,/silent
  N5p = ' '+strtrim(N_elements(arr[0,*]),2)
  chr_rdtbl,dir2p+'/snapshot_subhalo/subhalo00000',0,arr,/silent
  N2p = ' '+strtrim(N_elements(arr[0,*]),2)
  plot,x_10p,y_10p,psym=3,/isotropic,title='10%'+N10p
  plot,x_5p,y_5p,psym=3,/isotropic,title='5%'+N5p
  plot,x_2p,y_2p,psym=3,/isotropic,title='2%'+N2p

  plot,q,scale_dJ,xr=[-2.5,2.5],yr=[-2,2],psym=3
  plot,q_10p,scale_dJ_10p,xr=[-2.5,2.5],yr=[-2,2],psym=3
  plot,q_5p,scale_dJ_5p,xr=[-2.5,2.5],yr=[-2,2],psym=3
  plot,q_2p,scale_dJ_2p,xr=[-2.5,2.5],yr=[-2,2],psym=3

endfor

for k=1,5 do begin
  dir = '/media/SEADISK/LG/inVL_'+list[k]
  chr_rdtbl,dir+'/snapshot_subhalo/subhalo12790',0,arr,/silent
  arr = double(arr)
  t_sub = reform(arr[0,*])
  x_sub = reform(arr[1,*])
  y_sub = reform(arr[2,*])
  z_sub = reform(arr[3,*])
  Vx_sub = reform(arr[4,*])
  Vy_sub = reform(arr[5,*])
  Vz_sub = reform(arr[6,*])
  max_delE_sub = reform(arr[8,*])
  bsize = stddev(max_delE_sub)/2. 
  plothist,max_delE_sub,xh,yh,bin=bsize,/noplot
  plothist,max_delE_sub,bin=bsize,yr=[0.8,max(yh)],/ylog
endfor
multiplot,/reset
erase







device,/close

END
