pro disk_pass


dir_out = '/media/SEADISK/LG/realistic/'

;=== read output of the first particle ===
chr_rdtbl, dir_out+'part001', 0, arr
arr = double(arr)
t_part1 = arr[0,*]
x_part1 = arr[1,*]
y_part1 = arr[2,*]
z_part1 = arr[3,*]
Vx_part1 = arr[4,*]
Vy_part1 = arr[5,*]
Vz_part1 = arr[6,*]
delE_part1 = arr[7,*]
Etot_part1 = arr[8,*]
delJ_part1 = arr[9,*]
r_part1 = sqrt(x_part1^2.+y_part1^2.)


plot,[0],[0],/nodata,xr=[0,25],yr=[-20,20],/isotropic
plot, r_part1, z_part1
hline,0
vline,8
vline,7.5,linestlye=1

END
