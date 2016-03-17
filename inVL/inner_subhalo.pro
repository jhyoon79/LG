pro inner_subhalo

dir = '/media/SEADISK/LG/FinalRun/inVL_all/snapshot_subhalo/'

spawn,'ls '+dir+' > tmp'
chr_rdtbl,'tmp',0,fname
spawn,'rm -f tmp'
Nfile = N_elements(fname)
for i=0,Nfile-150 do begin
  chr_rdtbl,dir+fname[i],0,arr
  x = double(arr[1,*])
  y = double(arr[2,*])
  z = double(arr[3,*])
  r = sqrt(x^2+y^2+z^2)
  if i eq 0 then rmin = r
  sub = where(r lt rmin,count)
  if count ge 1 then rmin[sub] = r[sub]
endfor  
forprint,rmin,textout='rmin'

END
