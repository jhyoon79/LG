pro tails_bk

dir = '/media/SEADISK/LG/Final11/'
spawn,'ls '+dir+'snap* > tmp'
chr_rdtbl,'tmp',0,fname
spawn,'rm -f tmp'
Nfile = N_elements(fname)

@plot_setting
device,file='tails_bk.eps',/enc,xsize=24,ysize=18,/cmyk
loadct,13
!p.charsize=1.5
!p.multi=[0,4,3]
for i=0,11 do begin
;for i=3,3 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  tname = strmid(fname[i],30,5)
  t = string(float(tname)/1000,f='(f5.2)')+'Gyr'
    arr = double(arr)
    t_pts = (arr[0,0])[0]
    x_pts = reform(arr[1,*])
    y_pts = reform(arr[2,*])
    z_pts = reform(arr[3,*])
    vx_pts = reform(arr[4,*])
    vy_pts = reform(arr[5,*])
    vz_pts = reform(arr[6,*])
    Etot_pts = reform(arr[8,*])
    Nline = 232
    if tname eq '08770' or tname eq '09110' then Nline=290
    chr_rdtbl,dir+'part001_'+tname,0,arr,N_lines=Nline,/silent
    arr = double(arr)
    t_orbit = reform(arr[0,*])
    x_orbit = reform(arr[1,*])
    y_orbit = reform(arr[2,*])
    z_orbit = reform(arr[3,*])
  ;=== convert XYZ into RA & Dec ===
    xyz2radec,x_pts,y_pts,z_pts,ra,dec,l=l,b=b
phi = angsep(ra[0],dec[0],ra,dec)

  Nphi = 360*2
  k = (findgen(Nphi)+1)/Nphi*(360d)
  Bk = dblarr(Nphi)
  max_phi = max(abs(phi))
  for m=0,Nphi-1 do begin
    fact = k[m]*abs(phi)*!dtor
    Bk[m] = abs(1./Nphi*total(complex(cos(fact),sin(fact))))
  endfor
  k2 = 360./k
  xt = textoidl('360^{\circ} / k')
  xmin = 0.2  &  xmax = 360
  ymin = 0    &  ymax = 5
  plot,k2,Bk,xr=[xmin,xmax],yr=[ymin,ymax],/xlog

endfor
multiplot,/reset
erase
device,/close

END
