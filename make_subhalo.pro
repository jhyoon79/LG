pro make_subhalo

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

if (N_halo eq 0) then goto, jump

chr_rdtbl,'vlsubs.txt',3,arr
arr = double(arr)
Mtidal = arr[5,*]
arr = arr[*,reverse(sort(Mtidal))]
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
rs = rtidal/10.

print, 'max subhalo mass',Mtidal[0:9]
;sub_mass = where(Mtidal ge 1e7 and Mtidal lt 1e8)	;for *_10sub.ps
;sub_random = round(randomu(987,N_halo)*(N_elements(sub_mass)-1))
sub_mass = where(Mtidal ge 1e7 and Mtidal lt 1e8)	;for *_10sub.ps

if (N_elements(sub_mass) lt N_halo) then begin
  N_halo = N_elements(sub_mass)
  print, '=== not enought number of subhalos!!!  new N_halo',N_halo
endif

sub_random = round(randomu(987,N_halo)*(N_elements(sub_mass)-1))
sub_selected = sub_mass[sub_random]
;sub_selected = sub_mass
;sub_mass = min(abs(3.3253e9-reform(Mtidal)),sub)
;sub_selected = sub
;print, sub,Mtidal[sub]


Mvir_subhalo = Mtidal[sub_selected]
print, Mvir_subhalo
rs_subhalo = rs[sub_selected]
rvir_subhalo = rtidal[sub_selected]
x_subhalo = x[sub_selected]
y_subhalo = y[sub_selected]
z_subhalo = z[sub_selected]
Vx_subhalo = Vx[sub_selected]
Vy_subhalo = Vy[sub_selected]
Vz_subhalo = Vz[sub_selected]

print, 'mass',Mvir_subhalo
print, 'TOTAL Mass', total(Mvir_subhalo)

@plot_setting
set_plot,'ps'
device,file='make_subhalo.ps',/color
plothist, alog10(Mtidal),bin=0.5,xtitle=textoidl('log M_{tidal}'),ytitle='#'
if (N_elements(Mvir_subhalo) ge 2) then plothist, alog10(Mvir_subhalo),bin=0.5,/fill,fcolor=255,/overplot
plot, x,y,psym=3,xtitle='x',ytitle='y'
plot, y,z,psym=3,xtitle='y',ytitle='z'
plot, x,z,psym=3,xtitle='x',ytitle='z'
device,/close

jump:
openw, 5, 'frog_subhalo.dat'
for i=0L,N_halo-1 do printf,5,x_subhalo[i],y_subhalo[i],z_subhalo[i], $
		Vx_subhalo[i],Vy_subhalo[i],Vz_subhalo[i], $
		Mvir_subhalo[i],rs_subhalo[i],rvir_subhalo[i], $
		f='(6(f14.8,1x),a15,2(f13.8))'
close, 5


END
