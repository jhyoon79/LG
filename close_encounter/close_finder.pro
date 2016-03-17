; find closest encounter
pro close_finder

dir_BI5_test = '/scratch/jhyoon/Research/LG/BI5_test/'
dir_part = '/scratch/jhyoon/Research/LG/BI5_test/snapshot_e8_1000sub/'
dir_subhalo = '/scratch/jhyoon/Research/LG/BI5_test/snapshot_subhalo_e8_1000sub/'
dir_nosub = '/scratch/jhyoon/Research/LG/realistic/snapshot_nosub/'

;spawn, 'ls '+dir_part+' > temp9'
;chr_rdtbl,'temp9',0,fname_part
;fname_part = dir_part+reform(fname_part)

spawn, 'ls '+dir_subhalo+' > temp9'
chr_rdtbl,'temp9',0,fname_subhalo,/silent
fname_subhalo = dir_subhalo+reform(fname_subhalo)

;spawn, 'ls '+dir_nosub+' > temp9'
;chr_rdtbl,'temp9',0,fname_nosub
;fname_nosub = dir_nosub+reform(fname_nosub)
spawn, 'rm -f temp9'

chr_rdtbl,'../frog_BI5subhalo_e8_1000sub.dat',0,arrBI5,/silent
Mtidal = double(arrBI5[6,*])
Rs = double(arrBI5[7,*])
Rtidal = double(arrBI5[8,*])

Nn = 1
N_fname = N_elements(fname_subhalo)
sub_all = [-999]
for i=0,N_fname-1 do begin
  chr_rdtbl,fname_subhalo[i],0,arr,/silent
  arr = double(arr)
  d_subhalo_tail = arr[8,*]
  sub = where(d_subhalo_tail lt Nn*Rs)
  sub_all = [sub_all,sub]
  if (i/100 eq i/100.) then print,fname_subhalo[i]
endfor
sub_all = sub_all[where(sub_all gt 0)]
sub_final = redundant(sub_all,/non)
openw,1,'close_'+strtrim(Nn,2)+'Rs_e8'
for i=0,N_elements(sub_final)-1 do printf,1,arrBI5[*,sub_final[i]],f='(6(f14.8,1x),a15,2(1x,f13.8))'
close,1

END
