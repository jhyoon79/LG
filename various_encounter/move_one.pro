pro move_one

dir = '/media/SEADISK/LG/various_encounter/'
spawn,'ls '+dir+'snapshot_xy/ > tmp'
chr_rdtbl,'tmp',0,fpart
spawn,'ls '+dir+'snapshot_subhalo_xy/ > tmp'
chr_rdtbl,'tmp',0,fsub
spawn,'rm -f tmp'
Nfile = N_elements(fpart)

s = 0.0069943216
epsilon = 120.29439
r1 = [0]
r2 = [0]
r3 = [0]
!p.multi = [0,2,2]
;for i=10,40,2 do begin
for i=0,Nfile-1,2 do begin
  chr_rdtbl,dir+'snapshot_xy/'+fpart[i],0,arr,/silent
  arr = double(arr)
  tpts = arr[0,*]
  xpts = arr[1,*]
  ypts = arr[2,*]
  zpts = arr[3,*]
  dE = reform(arr[7,*])
  Etot = reform(arr[8,*])
  J = reform(arr[9,*])
  delE = Etot-Etot[0]
  q = delE/epsilon
  dJ = J-J[0]
  dJsJ = dJ / (s * J)

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

  plot,xpts,ypts,psym=3,/isotropic
  oplot,[xpts[2025]],[ypts[2025]],psym=7,color=255
  plot,q,dJsJ,psym=3,xr=[5,-5],yr=[-5,5],/isotropic

  wait,0.5
endfor
multiplot,/reset

END
