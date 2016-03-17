pro make_subhalo_VL

; LMC (Zhao 1998)
;	r = 50+-1 kpc
;	R = 127+-40 kpc
;	V = 84km/s
openr,1,'frogin'
readf,1,tmp
readf,1,tmp
readf,1,tmp
readf,1,N_halo
close,1

;if (N_halo eq 0) then goto, jump

;chr_rdtbl,'vlsubs.txt',3,arr
chr_rdtbl,'vltwosubs.txt',3,arr
arr = double(arr)
Mtidal = arr[5,*]
id = arr[0,*]
GCdistance = arr[1,*]
peakVmax = arr[2,*]
Vmax = arr[3,*]
rVmax = arr[4,*]
Mtidal = arr[5,*]
rtidal = arr[6,*]
x = arr[7,*]
y = arr[8,*]
z = arr[9,*]
Vx = arr[10,*]
Vy = arr[11,*]
Vz = arr[12,*]
rs = rVmax/2.
print,max(rs),max(rtidal)




print, 'max subhalo mass',Mtidal[0:9]
;sub_mass = where(Mtidal ge 1e7 and Mtidal lt 1e8)	;for *_10sub.ps
;sub_random = round(randomu(987,N_halo)*(N_elements(sub_mass)-1))
sub_mass = where(Mtidal ge 1e7 and Mtidal lt 1e8)	;for *_10sub.ps

if (N_elements(sub_mass) lt N_halo) then begin
;  N_halo = N_elements(sub_mass)
  print, '=== not enought number of subhalos!!!  new N_halo',N_halo
endif

;sub_random = round(randomu(987,N_halo)*(N_elements(sub_mass)-1))
;sub_selected = sub_mass[sub_random]
;sub_selected = sub_mass
;sub_mass = min(abs(3.3253e9-reform(Mtidal)),sub)
;sub_selected = sub
;print, sub,Mtidal[sub]


chr_rdtbl,'progGCdistance.txt',0,arr
progRgc = double(arr)
progRgc[where(progRgc eq 0)] = 9e9
min_Rgc = min(transpose(progRgc),dimension=2)
sub_R50 = where(min_Rgc le 50 and min_Rgc gt 0)
help,min_Rgc, sub_R50
sub_selected = sub_R50
;sub_selected = where(GCdistance lt 20+rs*40.)
Mvir_subhalo = Mtidal[sub_selected]
rs_subhalo = rs[sub_selected]
rvir_subhalo = rtidal[sub_selected]
x_subhalo = x[sub_selected]
y_subhalo = y[sub_selected]
z_subhalo = z[sub_selected]
Vx_subhalo = Vx[sub_selected]
Vy_subhalo = Vy[sub_selected]
Vz_subhalo = Vz[sub_selected]

sort_sub = reverse(sort(Mvir_subhalo))
Mvir_subhalo = Mvir_subhalo[sort_sub]
rs_subhalo = rs_subhalo[sort_sub]
rvir_subhalo = rvir_subhalo[sort_sub]
x_subhalo = x_subhalo[sort_sub]
y_subhalo = y_subhalo[sort_sub]
z_subhalo = z_subhalo[sort_sub]
Vx_subhalo = Vx_subhalo[sort_sub]
Vy_subhalo = Vy_subhalo[sort_sub]
Vz_subhalo = Vz_subhalo[sort_sub]

jump:
openw, 5, 'frog_VLsubhalo.dat'
for i=0L,N_elements(x_subhalo)-1 do $
     printf,5,x_subhalo[i],y_subhalo[i],z_subhalo[i], $
		Vx_subhalo[i],Vy_subhalo[i],Vz_subhalo[i], $
		Mvir_subhalo[i],rs_subhalo[i],rvir_subhalo[i], $
		f='(6(f14.5,1x),a16,2(f13.5))'
close, 5

openw, 5, 'frog_VLsubhalo_all.dat'
for i=0L,N_elements(x)-1 do $
     printf,5,x[i],y[i],z[i],Vx[i],Vy[i],Vz[i],Mtidal[i],rs[i],rtidal[i],$
		f='(6(f14.5,1x),a16,2(f13.5))'
close, 5


END
