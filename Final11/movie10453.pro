pro movie10453

dir = '/media/SEADISK/LG/Final11/snapshot_10453/'
spawn,'ls '+dir+' > tmp'
chr_rdtbl,'tmp',0,fpart
spawn,'rm -f tmp'
Nfile = N_elements(fpart)

s = 0.0069943216
epsilon = 120.29439

@plot_setting
device,file='movie10453.ps',/color,yoffset=1;,/landscape
!p.multi = [0,2,1]
!p.charsize=1
for i=50,1050-1,5 do begin
  chr_rdtbl,dir+fpart[i],0,arr,/silent
  arr = double(arr)
  tpts = arr[0,*]
  xpts = arr[1,*]
  ypts = arr[2,*]
  zpts = arr[3,*]
Etotpts = reform(arr[8,*])
Jpts = reform(arr[9,*])
delE = Etotpts-Etotpts[0]
qpts = delE/epsilon
dJpts = Jpts-Jpts[0]
dJsJ = dJpts / (s * Jpts)

  xcen = xpts[0]
  ycen = ypts[0]
  zcen = zpts[0]
;  plot,xpts,ypts,psym=3,xr=[-20,20],yr=[-20,20],/isotropic,xtit='x',ytit='y'
;  plot,ypts,zpts,psym=3,xr=[-20,20],yr=[-20,20],/isotropic,xtit='y',ytit='z'
  plot,xpts-xcen,ypts-ycen,psym=3,/isotropic,xtit='x',ytit='y'
  legend,'t='+strtrim(tpts[0],2)
;  plot,ypts-ycen,zpts-zcen,psym=3,/isotropic,xtit='y',ytit='z'
  plot,qpts,dJsJ,psym=3,xr=[5,-5],yr=[-5,5]

  if i/100 eq i/100. then print,tpts[0]
wait,0.5
endfor
multiplot,/reset
device,/close

END

