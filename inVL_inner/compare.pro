; Compare the tail with e5-e10 subhaloes and e6-e10 subhaloes.
; FinalRun/inVL_inner & inVL_inner_e6e10
pro compare

dir = '/media/SEADISK/LG/FinalRun/'

spawn,'ls '+dir+'inVL_inner/snapshot > tmp'
chr_rdtbl,'tmp',0,fname_e5,/silent
spawn,'ls '+dir+'inVL_inner_e6e10/snapshot > tmp'
chr_rdtbl,'tmp',0,fname_e6,/silent
spawn,'ls '+dir+'inVL_inner_e7e10/snapshot > tmp'
chr_rdtbl,'tmp',0,fname_e7,/silent
spawn,'ls '+dir+'inVL_inner_e8e10/snapshot > tmp'
chr_rdtbl,'tmp',0,fname_e8,/silent
spawn,'ls '+dir+'inVL_inner_e9e10/snapshot > tmp'
chr_rdtbl,'tmp',0,fname_e9,/silent
spawn,'rm -f tmp'

Mvir = 1.77e+12
rs = 24.6d
rvir = 389.0d
Msat = 10000.
G = 4.3d-6	; kpc Msun^-1 (km/s)^2 
c = rvir/rs
Mhalo = Mvir/(alog(c+1.)-c/(c+1.))

readcol,dir+'inVL_inner/part001',t,x,y,z
r_peri = min(sqrt(x*x+y*y+z*z))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s5 = (Msat/mr)^(1./3.)
epsilon5 = s5*(G*mr/r_peri)
readcol,dir+'inVL_inner_e6e10/part001',t,x,y,z
r_peri = min(sqrt(x*x+y*y+z*z))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s6 = (Msat/mr)^(1./3.)
epsilon6 = s6*(G*mr/r_peri)
readcol,dir+'inVL_inner_e7e10/part001',t,x,y,z
r_peri = min(sqrt(x*x+y*y+z*z))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s7 = (Msat/mr)^(1./3.)
epsilon7 = s7*(G*mr/r_peri)
readcol,dir+'inVL_inner_e8e10/part001_04640',t,x,y,z
r_peri = min(sqrt(x*x+y*y+z*z))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s8 = (Msat/mr)^(1./3.)
epsilon8 = s8*(G*mr/r_peri)
readcol,dir+'inVL_inner_e9e10/part001_04640',t,x,y,z
r_peri = min(sqrt(x*x+y*y+z*z))
p = r_peri/rs
mr = Mhalo*(alog(p+1.)-p/(p+1.))
s9 = (Msat/mr)^(1./3.)
epsilon9 = s9*(G*mr/r_peri)

@plot_setting
Nfile = N_elements(fname_e5)
;for i=0,Nfile-1,50 do begin
for i=844,844 do begin
  chr_rdtbl,'snap08440',0,arr,/silent
  t56 = arr[0,0]
  x56 = double(arr[1,*])
  y56 = double(arr[2,*])
  z56 = double(arr[3,*])
  Etot56 = double(arr[8,*])
  J56 = double(arr[9,*])
  dE56 = Etot56-Etot56[0]
  q56 = dE56/epsilon5
  dJ56 = J56-J56[0]
  scale_dJ56 = dJ56 / (s5 * J56)

  chr_rdtbl,dir+'inVL_inner/snapshot/'+fname_e5[i],0,arr,/silent
  t5 = arr[0,0]
  x5 = double(arr[1,*])
  y5 = double(arr[2,*])
  z5 = double(arr[3,*])
  Etot5 = double(arr[8,*])
  J5 = double(arr[9,*])
  dE5 = Etot5-Etot5[0]
  q5 = dE5/epsilon5
  dJ5 = J5-J5[0]
  scale_dJ5 = dJ5 / (s5 * J5)

  chr_rdtbl,dir+'inVL_inner_e6e10/snapshot/'+fname_e6[i],0,arr,/silent
  t6 = arr[0,0]
  x6 = double(arr[1,*])
  y6 = double(arr[2,*])
  z6 = double(arr[3,*])
  Etot6 = double(arr[8,*])
  J6 = double(arr[9,*])
  dE6 = Etot6-Etot6[0]
  q6 = dE6/epsilon6
  dJ6 = J6-J6[0]
  scale_dJ6 = dJ6 / (s6 * J6)
  xcen = x6[0]
  ycen = y6[0]
  chr_rdtbl,dir+'inVL_inner_e7e10/snapshot/'+fname_e7[i],0,arr,/silent
  t7 = arr[0,0]
  x7 = double(arr[1,*])
  y7 = double(arr[2,*])
  z7 = double(arr[3,*])
  Etot7 = double(arr[8,*])
  J7 = double(arr[9,*])
  dE7 = Etot7-Etot7[0]
  q7 = dE7/epsilon7
  dJ7 = J7-J7[0]
  scale_dJ7 = dJ7 / (s7 * J7)
  chr_rdtbl,dir+'inVL_inner_e8e10/snapshot/'+fname_e8[i],0,arr,/silent
  t8 = arr[0,0]
  x8 = double(arr[1,*])
  y8 = double(arr[2,*])
  z8 = double(arr[3,*])
  Etot8 = double(arr[8,*])
  J8 = double(arr[9,*])
  dE8 = Etot8-Etot8[0]
  q8 = dE8/epsilon8
  dJ8 = J8-J8[0]
  scale_dJ8 = dJ8 / (s8 * J8)
  chr_rdtbl,dir+'inVL_inner_e9e10/snapshot/'+fname_e9[i],0,arr,/silent
  t9 = arr[0,0]
  x9 = double(arr[1,*])
  y9 = double(arr[2,*])
  z9 = double(arr[3,*])
  Etot9 = double(arr[8,*])
  J9 = double(arr[9,*])
  dE9 = Etot9-Etot9[0]
  q9 = dE9/epsilon9
  dJ9 = J9-J9[0]
  scale_dJ9 = dJ9 / (s9 * J9)

  xyz2radec,x56,y56,z56,ra56,dec56
  xyz2radec,x5,y5,z5,ra5,dec5
  xyz2radec,x6,y6,z6,ra6,dec6
  xyz2radec,x7,y7,z7,ra7,dec7
  xyz2radec,x8,y8,z8,ra8,dec8
  xyz2radec,x9,y9,z9,ra9,dec9

  dy = 1.5
  sym = 1.5
  device,file='compare.eps',/color,/cmyk,/enc
  plot,[0],[0],xtit='RA',ytit='Dec',/nodata,xr=[295,227],yr=[-5,16]
  oplot,ra56,dec56-dy,psym=3,color=30
  oplot,[ra56[0]],[dec56[0]-dy],psym=1,color=30,symsize=sym
  oplot,ra5,dec5,psym=3
  oplot,[ra5[0]],[dec5[0]],psym=1,symsize=sym
  oplot,ra6,dec6+1*dy,psym=3,color=70
  oplot,[ra6[0]],[dec6[0]+1*dy],psym=1,color=70,symsize=sym
  oplot,ra7,dec7+2*dy,psym=3,color=100
  oplot,[ra7[0]],[dec7[0]+2*dy],psym=1,color=100,symsize=sym
  oplot,ra8,dec8+3*dy,psym=3,color=210
  oplot,[ra8[0]],[dec8[0]+3*dy],psym=1,color=210,symsize=sym
  oplot,ra9,dec9+4*dy,psym=3,color=255
  oplot,[ra9[0]],[dec9[0]+4*dy],psym=1,color=255,symsize=sym
  legend,['e5-e10','e6-e10','e7-e10','e8-e10','e9-e10'],psym=8,color=[0,70,100,210,255],box=0,/bottom,/horizontal,charsize=1.1
  device,/close

  device,file='compare_q.eps',/color,/cmyk,/enc,ysize=15
  dy = 5
  multiplot
  plot,[0],[0],xtit='q',ytit=textoidl('\DeltaJ / sJ'),xr=[7,-7],yr=[-8,24],/nodata
  oplot,q56,scale_dJ56-dy,psym=3,color=30
  oplot,[q56[0]],[scale_dJ56[0]-dy],psym=1,color=30,symsize=sym
  oplot,q5,scale_dJ5,psym=3
  oplot,[q5[0]],[scale_dJ5[0]],psym=1,symsize=sym
  oplot,q6,scale_dJ6+1*dy,psym=3,color=70
  oplot,[q6[0]],[scale_dJ6[0]]+1*dy,psym=1,color=70,symsize=sym
  oplot,q7,scale_dJ7+2*dy,psym=3,color=100
  oplot,[q7[0]],[scale_dJ7[0]]+2*dy,psym=1,color=100,symsize=sym
  oplot,q8,scale_dJ8+3*dy,psym=3,color=210
  oplot,[q8[0]],[scale_dJ8[0]]+3*dy,psym=1,color=210,symsize=sym
  oplot,q9,scale_dJ9+4*dy,psym=3,color=255
  oplot,[q9[0]],[scale_dJ9[0]]+4*dy,psym=1,color=255,symsize=sym
  device,/close

endfor

END
