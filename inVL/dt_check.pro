pro dt_check

dir_inVL100 = '/media/SEADISK/LG/inVL_allsub/dt1Myr/'
dir_inVL050 = '/media/SEADISK/LG/inVL_allsub/dt0.5Myr/'
dir_inVL025 = '/media/SEADISK/LG/inVL_allsub/dt0.25Myr/'


;=== read out data of the tail particles ===
chr_rdtbl,dir_inVL025+'part001',0,arr,/silent
arr = double(arr)
t_part025 = reform(arr[0,*])
x_part025 = reform(arr[1,*])
y_part025 = reform(arr[2,*])
z_part025 = reform(arr[3,*])
Vx_part025 = reform(arr[4,*])
Vy_part025 = reform(arr[5,*])
Vz_part025 = reform(arr[6,*])
dE_part025 = reform(arr[7,*])
Etot_part025 = reform(arr[8,*])
J_part025 = reform(arr[9,*])

chr_rdtbl,dir_inVL050+'part001',0,arr,/silent
arr = double(arr)
t_part050 = reform(arr[0,*])
x_part050 = reform(arr[1,*])
y_part050 = reform(arr[2,*])
z_part050 = reform(arr[3,*])
Vx_part050 = reform(arr[4,*])
Vy_part050 = reform(arr[5,*])
Vz_part050 = reform(arr[6,*])
dE_part050 = reform(arr[7,*])
Etot_part050 = reform(arr[8,*])
J_part050 = reform(arr[9,*])

chr_rdtbl,dir_inVL100+'part001',0,arr,/silent
arr = double(arr)
t_part = reform(arr[0,*])/1000.
x_part = reform(arr[1,*])
y_part = reform(arr[2,*])
z_part = reform(arr[3,*])
Vx_part = reform(arr[4,*])
Vy_part = reform(arr[5,*])
Vz_part = reform(arr[6,*])
dE_part = reform(arr[7,*])
Etot_part = reform(arr[8,*])
J_part = reform(arr[9,*])

chr_rdtbl,dir_inVL025+'snap12758',0,arr,/silent
arr = double(arr)
t_f025 = reform(arr[0,*])
x_f025 = reform(arr[1,*])
y_f025 = reform(arr[2,*])
z_f025 = reform(arr[3,*])
Vx_f025 = reform(arr[4,*])
Vy_f025 = reform(arr[5,*])
Vz_f025 = reform(arr[6,*])
dE_f025 = reform(arr[7,*])
Etot_f025 = reform(arr[8,*])
J_f025 = reform(arr[9,*])

chr_rdtbl,dir_inVL050+'snap12758',0,arr,/silent
arr = double(arr)
t_f050 = reform(arr[0,*])
x_f050 = reform(arr[1,*])
y_f050 = reform(arr[2,*])
z_f050 = reform(arr[3,*])
Vx_f050 = reform(arr[4,*])
Vy_f050 = reform(arr[5,*])
Vz_f050 = reform(arr[6,*])
dE_f050 = reform(arr[7,*])
Etot_f050 = reform(arr[8,*])
J_f050 = reform(arr[9,*])

chr_rdtbl,dir_inVL100+'snap12758',0,arr,/silent
arr = double(arr)
t_f = reform(arr[0,*])
x_f = reform(arr[1,*])
y_f = reform(arr[2,*])
z_f = reform(arr[3,*])
Vx_f = reform(arr[4,*])
Vy_f = reform(arr[5,*])
Vz_f = reform(arr[6,*])


set_plot,'ps'
@plot_setting
!p.multi=[0,2,3]
device,file='dt_check.ps',/color,ysize=25
plot,x_part,x_part050,psym=3,xtitle='x 1Myr',ytitle='x 0.50Myr',/isotropic
plot,Vx_part,Vx_part050,psym=3,xtitle='Vx 1Myr',ytitle='Vx 0.50Myr',/isotropic
plot,x_f,x_f050,psym=3,xtitle='x f 1Myr',ytitle='x f 0.50Myr',/isotropic
plot,x_part025,x_part050,psym=3,xtitle='x 0.25Myr',ytitle='x 0.50Myr',/isotropic
plot,Vx_part025,Vx_part050,psym=3,xtitle='Vx 0.25Myr',ytitle='Vx 0.50Myr',/isotropic
plot,x_f025,x_f050,psym=3,xtitle='x f 0.25Myr',ytitle='x f 0.50Myr',/isotropic

plot,y_part,y_part050,psym=3,xtitle='y 1Myr',ytitle='y 0.50Myr',/isotropic
plot,Vy_part,Vy_part050,psym=3,xtitle='Vy 1Myr',ytitle='Vy 0.50Myr',/isotropic
plot,y_f,y_f050,psym=3,xtitle='y f 1Myr',ytitle='y f 0.50Myr',/isotropic
plot,y_part025,y_part050,psym=3,xtitle='y 0.25Myr',ytitle='y 0.50Myr',/isotropic
plot,Vy_part025,Vy_part050,psym=3,xtitle='Vy 0.25Myr',ytitle='Vy 0.50Myr',/isotropic
plot,y_f025,y_f050,psym=3,xtitle='y f 0.25Myr',ytitle='y f 0.50Myr',/isotropic

plot,z_part,z_part050,psym=3,xtitle='z 1Myr',ytitle='z 0.50Myr',/isotropic
plot,Vz_part,Vz_part050,psym=3,xtitle='Vz 1Myr',ytitle='Vz 0.50Myr',/isotropic
plot,z_f,z_f050,psym=3,xtitle='z f 1Myr',ytitle='z f 0.50Myr',/isotropic
plot,z_part025,z_part050,psym=3,xtitle='z 0.25Myr',ytitle='z 0.50Myr',/isotropic
plot,Vz_part025,Vz_part050,psym=3,xtitle='Vz 0.25Myr',ytitle='Vz 0.50Myr',/isotropic
plot,z_f025,z_f050,psym=3,xtitle='z f 0.25Myr',ytitle='z f 0.50Myr',/isotropic

bsize = 0.001
plothist,x_part-x_part050,bin=bsize,xtitle='x 1Myr'
plothist,x_f-x_f050,bin=bsize,xtitle='x f 1Myr'
plothist,(x_f-x_f050)/x_f050,bin=bsize,xr=[-1,1],xtitle='x f 1Myr'
plothist,x_part025-x_part050,bin=bsize,xtitle='x 0.25Myr'
plothist,x_f025-x_f050,bin=bsize,xtitle='x f 0.25Myr'
plothist,(x_f025-x_f050)/x_f025,bin=bsize,xr=[-1,1],xtitle='x f 0.25Myr'


plothist,y_part-y_part050,bin=bsize,xtitle='y 1Myr'
plothist,y_f-y_f050,bin=bsize,xtitle='y f 1Myr'
plothist,(y_f-y_f050)/y_f050,bin=bsize,xtitle='y f 1Myr'
plothist,y_part025-y_part050,bin=bsize,xtitle='y 0.25Myr'
plothist,y_f025-y_f050,bin=bsize,xtitle='y f 0.25Myr'
plothist,(y_f025-y_f050)/y_f050,bin=bsize,xtitle='y f 0.25Myr'

plothist,z_part-z_part050,bin=bsize,xtitle='z 1Myr'
plothist,z_f-z_f050,bin=bsize,xtitle='z f 1Myr'
plothist,(z_f-z_f050)/z_f050,bin=bsize,xtitle='z f 1Myr'
plothist,z_part025-z_part050,bin=bsize,xtitle='z 0.25Myr'
plothist,z_f025-z_f050,bin=bsize,xtitle='z f 0.25Myr'
plothist,(z_f025-z_f050)/z_f050,bin=bsize,xtitle='z f 0.25Myr'


plot,t_part,x_part-x_part050,psym=3,xtitle='t [Gyr]',ytitle='x(1Myr) - x(0.5Myr)'
plot,t_part,y_part-y_part050,psym=3,xtitle='t [Gyr]',ytitle='y(1Myr) - y(0.5Myr)'
plot,t_part,z_part-z_part050,psym=3,xtitle='t [Gyr]',ytitle='z(1Myr) - z(0.5Myr)'
plot,t_part,(x_part-x_part050)/x_part050,psym=3,xtitle='t [Gyr]',ytitle='x(1Myr) - x(0.5Myr)',yr=[-1,1]

plot,t_part,x_part025-x_part050,psym=3,xtitle='t [Gyr]',ytitle='x(0.25Myr) - x(0.5Myr)'
plot,t_part,y_part025-y_part050,psym=3,xtitle='t [Gyr]',ytitle='y(0.25Myr) - y(0.5Myr)'
plot,t_part,z_part025-z_part050,psym=3,xtitle='t [Gyr]',ytitle='z(0.25Myr) - z(0.5Myr)'
plot,t_part,(x_part025-x_part050)/x_part025,psym=3,xtitle='t [Gyr]',ytitle='x(0.25Myr) - x(0.5Myr)',yr=[-1,1]

device,/close

END
