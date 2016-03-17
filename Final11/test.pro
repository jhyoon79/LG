pro test

dir = '/media/SEADISK/LG/Final11/'
spawn,'ls '+dir+'snap* > tmp'
chr_rdtbl,'tmp',0,fname
spawn,'rm -f tmp'
Nfile = N_elements(fname)

;=====================
; radec contour map
;=====================
    chr_rdtbl,dir+'00000_08770',0,arr
    arr = double(arr)
    Etot_init = reform(arr[8,*])
    Jtot_init = reform(arr[9,*])
    delE_init = Etot_init-Etot_init[0]
dJ = Jtot_init-Jtot_init[0]
epsilon = 120.29439
    q_init = delE_init/epsilon
    qmin = min(q_init);*0.9
    qmax = max(q_init);*0.9
    qcolor = round((q_init-qmin)/(qmax-qmin)*245.)+10
set_plot,'ps'
loadct,13
device,file='test_cntr.ps',xsize=15,ysize=25,/color
!p.charsize=1.2
!p.multi=[0,1,4]
multiplot
plot,q_init,dJ,xr=[max(q_init),min(q_init)],/nodata
for i=0,N_elements(dJ)-1 do oplot,[q_init[i]],[dJ[i]],color=qcolor[i],psym=3
ind = [0,2,5,10]
for ii=1,3 do begin
  i = ind[ii]
  if ii eq 0 then chr_rdtbl,'/media/SEADISK/LG/realistic/snapshot/snap08440',0,arr,/silent
  if ii gt 0 then chr_rdtbl,fname[i],0,arr,/silent
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
    xyz2radec,x_pts,y_pts,z_pts,ra_pts,dec_pts
    xyz2radec,x_orbit,y_orbit,z_orbit,ra_orbit,dec_orbit
  
    sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
    ra_orbit_sub = ra_orbit[sub_orbit]
    dec_orbit_sub = dec_orbit[sub_orbit]
  
  xmin = min(ra_pts)
  xmax = max(ra_pts)
  ymin = min(dec_pts)
  ymax = max(dec_pts)
  xt = ''  &  yt = 'Dec [deg]'
  if ii eq 3 then xt = 'RA [deg]'
  xr = xmax-xmin
  yr = ymax-ymin
  xgn = round(xr*1/0.3)
  ygn = round(yr*1/0.3)
print,'xr=',xr,xgn,'   yr=',yr,ygn
  sig = 1.5
  lev = [2*sig,3*sig,4*sig,5*sig,10*sig,15*sig,20*sig]
  pts_contour,ra_pts,dec_pts,xgn,ygn,gv,x,y
  dy = (ii-1)*0.3
  x1 = 0.15
  x2 = 0.7
  y1 = 0.7 - dy
  y2 = 0.95 - dy
  multiplot,ygap=0.02,/doxaxis
;  contour,gv,x,y,levels=lev,xtit=xt,ytit=yt,thick=1.,xr=[xmax,xmin],yr=[ymin,ymax];,pos=[x1,y1,x2,y2]
plot,ra_pts,dec_pts,/nodata,xr=[max(ra_pts),min(ra_pts)]
for k=0,N_elements(ra_pts)-1 do oplot,[ra_pts[k]],[dec_pts[k]],color=qcolor[k],psym=3
  oplot,[ra_pts[0]],dec_pts[[0]],psym=1,color=150,symsize=2,thick=7
  if ii eq 2 then arrow,235,0,232.5,0,/data
  if ii eq 3 then arrow,232,0.25,230,0.25,/data
  if ii eq 3 then arrow,244,2,244,3.5,/data

  legend,t,box=0,/bottom,charsize=1.2
endfor
multiplot,/reset
device,/close

END
