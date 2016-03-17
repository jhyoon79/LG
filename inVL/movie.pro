pro movie

dir='/media/SEADISK/LG/FinalRun/inVL_inner_e5e6/'

spawn,'ls '+dir+'snapshot/ > tmp'
chr_rdtbl,'tmp',0,fname

spawn,'ls '+dir+'snapshot_subhalo/ > tmp'
chr_rdtbl,'tmp',0,fname_sub

@plot_setting
device,file='movie.ps',/color,ysize=24
!p.multi=[0,2,2]
for i=0,N_elements(fname)-1,5 do begin
  chr_rdtbl,dir+'snapshot/'+fname[i],0,arr,/silent
  t=double(arr[0,*])
  x=double(arr[1,*])
  y=double(arr[2,*])
  z=double(arr[3,*])

  chr_rdtbl,dir+'snapshot_subhalo/'+fname_sub[i],0,arr,/silent
  xsub=double(arr[1,*])
  ysub=double(arr[2,*])
  zsub=double(arr[3,*])
  mind=double(arr[7,*])
  maxE=double(arr[8,*])
  subd = where(mind lt 3)
  subE = where(maxE gt .5)

  plot,x,y,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
  oplot,xsub[subd],ysub[subd],psym=6
  plot,x-x[0],y-y[0],psym=3,xr=[-5,5],yr=[-5,5],/isotropic
  oplot,xsub[subd]-x[0],ysub[subd]-y[0],psym=6
  legend,[strtrim(min(mind[subd],sss),2)],box=0
  legend,[strtrim(maxE[subd[sss]],2)],box=0,/bottom

  plot,x,y,psym=3,xr=[-20,20],yr=[-20,20],/isotropic
  oplot,xsub[subE],ysub[subE],psym=6
  plot,x-x[0],y-y[0],psym=3,xr=[-5,5],yr=[-5,5],/isotropic
  oplot,xsub[subE]-x[0],ysub[subE]-y[0],psym=6
  legend,[strtrim(max(maxE[subE]),2)],box=0

endfor
device,/close

END
