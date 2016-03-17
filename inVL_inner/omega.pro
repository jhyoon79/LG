pro omega

dir1 = '/media/SEADISK/LG/FinalRun/inVL_inner/'
dir2 = '/media/SEADISK/LG/FinalRun/inVL_inner_e6e10/'

chr_rdtbl,'frog_VLsubhalo_inner.dat',0,arr
Msub1 = double(arr[6,*])
chr_rdtbl,'frog_VLsubhalo_inner_e6e10.dat',0,arr
Msub2 = double(arr[6,*])

readcol,dir1+'part001',t1,x1,y1,z1
readcol,dir1+'part001_03400',t2,x2,y2,z2
readcol,dir1+'part001_06520',t3,x3,y3,z3
t = [t1,t2,t3]
x = [x1,x2,x3]
y = [y1,y2,y3]
z = [z1,z2,z3]
r1 = sqrt(x^2.+y^2.+z^2.)
t1 = t

dr1 = shift(r1,1) - r1
dr2 = r1 - shift(r1,-1)
dr1[0] = 0
dr2[N_elements(r1)-1] = 0
cri = dr1*dr2
sub = where(cri lt 0)
rmean = mean(r1)
apo = where(r1[sub] gt rmean,complement=peri)
t_apo1 = t1[sub[apo]]
t_peri1 = t1[sub[peri]]
r_apo1 = r1[sub[apo]]
r_peri1 = r1[sub[peri]]

readcol,dir2+'part001',t2,x2,y2,z2
r2 = sqrt(x2^2.+y2^2.+z2^2.)

dr1 = shift(r2,1) - r2
dr2 = r2 - shift(r2,-1)
dr1[0] = 0
dr2[N_elements(r2)-1] = 0
cri = dr1*dr2
sub = where(cri lt 0)
rmean = mean(r2)
apo = where(r2[sub] gt rmean,complement=peri)
t_apo2 = t2[sub[apo]]
t_peri2 = t2[sub[peri]]
r_apo2 = r2[sub[apo]]
r_peri2 = r2[sub[peri]]

dt1 = t_apo1 - shift(t_apo1,1)
dt2 = t_apo2 - shift(t_apo2,1)
tR1 = mean(dt1[where(dt1 gt 300 and dt1 lt 400)])
tR2 = mean(dt2[where(dt2 gt 300 and dt2 lt 400)])

dt1 = t_peri1 - shift(t_peri1,1)
dt2 = t_peri2 - shift(t_peri2,1)
tR3 = mean(dt1[where(dt1 gt 300 and dt1 lt 400)])
tR4 = mean(dt2[where(dt2 gt 300 and dt2 lt 400)])

print,'dOmega=',0.5*(total(Msub1)-total(Msub2))/1.77e12

print,tR1,tR2,tR3,tR4
print,(tR1-tR2)/tR2
print,(tR3-tR4)/tR4

stop
END
