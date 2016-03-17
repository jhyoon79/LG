pro oscil

dir = '/media/SEADISK/LG/various_encounter/'
spawn,'ls '+dir+'snapshot_xy/ > tmp'
chr_rdtbl,'tmp',0,fpart
spawn,'ls '+dir+'snapshot_xy_nosub5500/ > tmp'
chr_rdtbl,'tmp',0,fsub
spawn,'rm -f tmp'
Nfile = N_elements(fpart)

s = 0.0069943216
epsilon = 120.29439
r1 = [0]
r2 = [0]
;for i=150,150 do begin
for i=0,Nfile-1,1 do begin
  chr_rdtbl,dir+'snapshot_xy/'+fpart[i],0,arr,/silent
  arr = double(arr[*,2025])
  tpts = arr[0,*]
  xpts = arr[1,*]
  ypts = arr[2,*]
  zpts = arr[3,*]
  dE = reform(arr[7,*])
  Etot = reform(arr[8,*])
  J = reform(arr[9,*])
  delE = Etot-Etot[0]
  q = delE/epsilon
  dJ = J-J[0]
  dJsJ = dJ / (s * J)

  chr_rdtbl,dir+'snapshot_xy_nosub5500/'+fsub[i],0,arr,/silent
  arr = double(arr[*,2025])
  tnosub = arr[0,*]
  xnosub = arr[1,*]
  ynosub = arr[2,*]
  znosub = arr[3,*]
  dE1 = reform(arr[7,*])
  Etot1 = reform(arr[8,*])
  J1 = reform(arr[9,*])
  delE1 = Etot1-Etot1[0]
  q1 = delE1/epsilon
  dJ1 = J1-J1[0]
  dJsJ1 = dJ1 / (s * J1)

  r1 = [r1,sqrt(xpts^2.+ypts^2.+zpts^2.)]
  r2 = [r2,sqrt(xnosub^2.+ynosub^2.+znosub^2.)]
endfor
r1 = r1[1:*]
r2 = r2[1:*]
plot,r1
oplot,r2,color=255,linestyle=2
stop

END
