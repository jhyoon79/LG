PRO movie_e5e6

dir = '/media/SEADISK/LG/FinalRun/inVL_inner_e5e6/'
spawn,'ls '+dir+'snapshot/snap* > tmp'
chr_rdtbl,'tmp',0,fname,/silent
spawn,'ls '+dir+'snapshot_subhalo/sub* > tmp'
chr_rdtbl,'tmp',0,fsubname,/silent
spawn,'rm -f tmp'
Nf = N_elements(fname)

G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))
;chr_rdtbl,dir+'part_peri',0,arr,/silent
chr_rdtbl,'/media/SEADISK/LG/realistic/part_peri_nosub',0,arr,/silent
arr = double(arr)
t_peri = reform(arr[0,*])
x_peri = reform(arr[1,*])
y_peri = reform(arr[2,*])
z_peri = reform(arr[3,*])
Vx_peri = reform(arr[4,*])
Vy_peri = reform(arr[5,*])
Vz_peri = reform(arr[6,*])
r_peri = min(sqrt(x_peri^2.+y_peri^2.+z_peri^2.))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
Msat = 10000.
s = (Msat/mr)^(1./3.)
epsilon = s*(G*mr/r_peri)

;=== N-body ===

chr_rdtbl,dir+'snapshot/snap00000',0,arr0,/silent
Etot0 = double(reform(arr0[8,0]))
Jtot0 = double(reform(arr0[9,0]))

chr_rdtbl,'frog_VLsubhalo_inner_e5e6all.dat',0,arr,/silent
Msub = double(arr[6,*])
rs = double(arr[7,*])
chr_rdtbl,dir+'snapshot_subhalo/subhalo08440',0,arr,/silent
dmin = double(arr[7,*])
enc = where(dmin lt rs*50,Nenc)
Msub = Msub[enc]
rs = rs[enc]
help,rs,enc

@plot_setting
device,file='movie_e5e6.ps',/color,/landscape
!p.charsize=1.5
!p.multi=[0,2,2]
;for i=0,Nf-1,10 do begin
;1: 142-145 :M=432835,rs=0.03731,#=49031,q_impact=-1.8
;2: 182-185 :M=119312,rs=0.02626,#=84482,q_impact=-2
;3: 648-649 :M=929192,rs=0.04704,#=50455,q_impact=-0.6
;           :M=412268,rs=0.07541,#=42390,q_impact=-0.6
;4: 675-680 :M=534441,rs=0.04646,#=83312,q_impact=-0.6
;215-220
;645-690
;690-730
;835-850
for i=0,Nf-1,5 do begin
;for i=646,649,1 do begin
  chr_rdtbl,fname[i],0,arr,/silent
  t = string(float(arr[0,0])/1000,f='(f5.3)')
  arr = double(arr)
  x = reform(arr[1,*])
  y = reform(arr[2,*])
  z = reform(arr[3,*])
  vx = reform(arr[4,*])
  vy = reform(arr[5,*])
  vz = reform(arr[6,*])
  Etot = reform(arr[8,*])
  Jtot = reform(arr[9,*])
  dE = Etot - Etot[0]
  dJ = Jtot - Jtot[0]
  q = dE/epsilon
  dJsJ = dJ/(s*Jtot)

  q0 = (Etot - Etot0[0])/epsilon
  dJsJ0 = (Jtot-Jtot0[0])/(s*Jtot[0])

  chr_rdtbl,fsubname[i],0,arr,/silent
  arr = double(arr[*,enc])
  xsub = reform(arr[1,*])
  ysub = reform(arr[2,*])
  zsub = reform(arr[3,*])

  tmp = min(abs(q+0.6),impact)
  dsub = sqrt((x[impact]-xsub)^2.+(y[impact]-ysub)^2.+(z[impact]-zsub)^2.)
  near = where(dsub lt 5,Nnear)
help,near

  xmin = -19  &  xmax = 19
  ymin = -19  &  ymax = 19
  plot,x,y,psym=3,xtit='x [kpc]',ytit='y [kpc]',xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  oplot,[x[0]],[y[0]],psym=1
;  for j=0,Nenc-1 do begin
;    subhalo_circle = circle(xsub[j],ysub[j],rs[j])
;    oplot,subhalo_circle[0,*],subhalo_circle[1,*]
;  endfor
  legend,'t='+t+'Gyr',box=0,charsize=1

  xmin = -5  &  xmax = 5
  ymin = -5  &  ymax = 5
  zmin = -5  &  zmax = 5
  xcen = x[impact]
  ycen = y[impact]
  zcen = z[impact]
;  plot,x-xcen,y-ycen,psym=3,xtit=textoidl('\Deltax [kpc]'),ytit=textoidl('\Deltay [kpc]'),xr=[xmin,xmax],yr=[ymin,ymax],/isotropic
  plot,x-xcen,y-ycen,psym=3,xtit=textoidl('\Deltax [kpc]'),ytit=textoidl('\Deltay [kpc]'),/isotropic
  oplot,[0],[0],psym=1
  oplot,[x[impact]-xcen],[y[impact]-ycen],psym=7,color=100
  for j=0L,Nnear-1 do begin
    subhalo_circle = circle(xsub[near[j]],ysub[near[j]],rs[near[j]])
    oplot,subhalo_circle[0,*]-xcen,subhalo_circle[1,*]-ycen,color=255
  endfor
  mind = min(dsub,min_sub)
  oplot,[xsub[min_sub]]-xcen,[ysub[min_sub]]-ycen,psym=6
  legend,['Msub='+strtrim(Msub[min_sub],2),'rs='+strtrim(rs[min_sub],2),'#='+strtrim(enc[min_sub],2)],box=0,charsize=0.8
;  oplot,[xsub[4379]]-xcen,[ysub[4379]]-ycen,psym=4
;  oplot,[xsub[8438]]-xcen,[ysub[8438]]-ycen,psym=5
  oplot,[xsub[4318]]-xcen,[ysub[4318]]-ycen,psym=4
  oplot,[xsub[5146]]-xcen,[ysub[5146]]-ycen,psym=5
;  oplot,[xsub[5011]]-xcen,[ysub[5011]]-ycen,psym=5
;  oplot,[xsub[8547]]-xcen,[ysub[8547]]-ycen,psym=5

;  plot,x-xcen,z-zcen,psym=3,xtit=textoidl('\Deltax [kpc]'),ytit=textoidl('\Deltaz [kpc]'),xr=[xmin,xmax],yr=[zmin,zmax],/isotropic
  plot,x-xcen,z-zcen,psym=3,xtit=textoidl('\Deltax [kpc]'),ytit=textoidl('\Deltaz [kpc]'),/isotropic
  oplot,[0],[0],psym=1
  oplot,[x[impact]-xcen],[z[impact]-zcen],psym=7,color=100
  for j=0L,Nnear-1 do begin
    subhalo_circle = circle(xsub[near[j]],zsub[near[j]],rs[near[j]])
    oplot,subhalo_circle[0,*]-xcen,subhalo_circle[1,*]-zcen,color=255
  endfor
  oplot,[xsub[min_sub]]-xcen,[zsub[min_sub]]-zcen,psym=6
  legend,['Msub='+strtrim(Msub[min_sub],2),'rs='+strtrim(rs[min_sub],2),'#='+strtrim(enc[min_sub],2)+','+strtrim(min_sub,2)],box=0,charsize=0.8
;  oplot,[xsub[4379]]-xcen,[zsub[4379]]-zcen,psym=4
;  oplot,[xsub[8438]]-xcen,[zsub[8438]]-zcen,psym=5
  oplot,[xsub[4318]]-xcen,[zsub[4318]]-zcen,psym=4
  oplot,[xsub[5146]]-xcen,[zsub[5146]]-zcen,psym=5
;  oplot,[xsub[5011]]-xcen,[zsub[5011]]-zcen,psym=5
;  oplot,[xsub[8547]]-xcen,[zsub[8547]]-zcen,psym=5

  xt = 'q'
  yt = textoidl('\DeltaJ / sJ')
  plot,q,dJsJ,xr=[5,-5],yr=[-5,5],xtitle=xt,ytitle=yt,psym=3,/isotropic
  oplot,[q[0]],[dJsJ[0]],psym=1
  oplot,[q[impact]],[dJsJ[impact]],psym=7,color=100
  vline,0,0,linestyle=2
  hline,0,0,linestyle=2

;  xt = 'q0'
;  yt = textoidl('\DeltaJ0 / sJ')
;  plot,q0,dJsJ0,xr=[5,-5],yr=[-5,5],xtitle=xt,ytitle=yt,psym=3,/isotropic
;  oplot,[q0[0]],[dJsJ0[0]],psym=1
;  vline,0,0,linestyle=2
;  hline,0,0,linestyle=2
  if i/100 eq i/100. then print,i
endfor
multiplot,/reset
ERASE
device,/close

END
