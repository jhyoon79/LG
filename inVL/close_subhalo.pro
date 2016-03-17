pro close_subhalo

set_plot,'ps'
@plot_setting
device,file='close_subhalo.ps',/color
!p.multi=[0,1,2]
list = ['e5e6','e6e7','e7e8','e8e9','e9e10','e6e7_30000']

for j=0,N_elements(list)-1 do begin
  dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[j]+'/snapshot_subhalo/'
  chr_rdtbl,dir+'subhalo08440',0,arr,/silent
  mind = double(arr[7,*])
  max_delE = double(arr[8,*])
  sub_sort = sort(max_delE)
  Nsub = N_elements(max_delE)
  
  pre_ave = 9e9
  ave = mean(max_delE)
  std = stddev(max_delE)
  while (ave ne pre_ave) do begin
    pre_ave = ave
    clip = where(abs(max_delE-ave) lt 3.*std)
    ave = mean(max_delE[clip])
    std = stddev(max_delE[clip])
  endwhile

  if j eq 0 then begin
    cri = .1;0.2
  endif else if j eq 1 then begin
    cri = .3;0.6
  endif else if j eq 2 then begin
    cri = .4;0.8
  endif else if j eq 3 then begin
    cri = .8;1.6
  endif else if j eq 4 then begin
    cri = .5;1
  endif
;  cri = ave+5*std
  Nclose = fix(Nsub/5)
  cri_E = max_delE[(reverse(sort(max_delE)))[Nclose-1]]
  cri_d = mind[(sort(mind))[Nclose-1]]
  bsize = .01
  plothist,max_delE,xh,yh,bin=bsize,/noplot

  plothist,max_delE,bin=bsize,yr=[0.7,max(yh)],/ylog
  vline,ave,ymin=0.1,ymax=9e5,color=255,/ylog
  vline,cri,ymin=0.1,ymax=9e5,linestyle=2,/ylog

  plothist,max_delE,bin=bsize,xr=[0,5],yr=[0.7,max(yh)],/ylog
  vline,ave,ymin=0.1,ymax=9e5,color=255,/ylog
  vline,cri,ymin=0.1,ymax=9e5,linestyle=2,/ylog
  
  outlier_E = where(max_delE ge cri_E,count_E)
  outlier_d = where(mind le cri_d,count_d)
  chr_rdtbl,'../frog_VLsubhalo_'+list[j]+'.dat',0,arr,/silent
  openw,1,'close_'+list[j]+'_20p_'+strtrim(Nclose,2)+'_E'
  for i=0,count_E-1 do printf,1,arr[*,outlier_E[i]],f='(6(f14.5,1x),a16,2(f13.5))'
  close,1
  openw,1,'close_'+list[j]+'_20p_'+strtrim(Nclose,2)+'_d'
  for i=0,count_d-1 do printf,1,arr[*,outlier_d[i]],f='(6(f14.5,1x),a16,2(f13.5))'
  close,1
endfor

device,/close

END
