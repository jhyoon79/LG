pro move

dir = '/media/SEADISK/LG/various_encounter/evol/'
spawn,'ls '+dir+'snapshot_xy/ > tmp'
chr_rdtbl,'tmp',0,fpart
spawn,'ls '+dir+'snapshot_subhalo_xy/ > tmp'
chr_rdtbl,'tmp',0,fsub
spawn,'rm -f tmp'
Nfile = N_elements(fpart)

!p.multi = [0,2,1]
for i=15,35,1 do begin
;for i=0,Nfile-1,2 do begin
  chr_rdtbl,dir+'snapshot_xy/'+fpart[i],0,arr,/silent
  arr = double(arr)
  tpts = arr[0,*]
  xpts = arr[1,*]
  ypts = arr[2,*]
  zpts = arr[3,*]
  chr_rdtbl,dir+'snapshot_subhalo_xy/'+fsub[i],0,arr
  arr = double(arr)
  tsub = arr[0,*]
  xsub = arr[1,*]
  ysub = arr[2,*]
  zsub = arr[3,*]
  plot,xpts,ypts,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
  oplot,xsub,ysub,psym=1
  legend,strtrim(tpts[0],2)
  plot,ypts,zpts,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
  oplot,ysub,zsub,psym=1
  wait,0.7
endfor
multiplot,/reset

END

