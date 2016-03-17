#!/bin/csh
rm -rf time_running

foreach run_num (1 )

echo 'run number='$run_num
if ($run_num == 1) then 
set pre = `head -45 make_subhalo_BI5.pro | grep iseed | awk '{print $3}'`
vi make_subhalo_BI5.pro << EOF > temp
:44s/$pre/$run_num
:wq
EOF
else
@ pre = $run_num - 1
vi make_subhalo_BI5.pro << EOF > temp
:44s/$pre/$run_num
:wq
EOF
endif

vi make_subhalo_BI5.pro << EOF > temp
:64s/1.e10/1.e6
:64s/1.e9/1.e6
:64s/1.e8/1.e6
:64s/1.e7/1.e6
:wq
EOF

foreach i (6 7 8 9 10) 

echo 'subhalo mass of 10^'$i' solar mass'
@ j = $i - 1
vi make_subhalo_BI5.pro << EOF > temp
:64s/1.e$j/1.e$i
:wq
EOF

set pre = `head -4 frogin | grep subhalos | awk '{print $1}'`
vi frogin << EOF > temp
:4s/$pre/1
:wq
EOF

  foreach Ns (1 10 100 1000)

if ($Ns != 1) then 
@ Ns_pre = $Ns / 10
vi frogin << EOF > temp
:4s/$Ns_pre/$Ns
:wq
EOF
endif

set Nsub = `head -4 frogin | awk '/halos/ {print $1}'`
echo 'Nsubhalos='$Nsub

head -64 make_subhalo_BI5.pro | grep Mvir_subhalo

idl make_Pal5.idl
idl make_subhalo_BI5.idl
f77 frog.f -o frog.out
./frog.out
idl Pal5_plot.idl

echo e$i\_$Nsub `cat tmp_time` >> time_running

mkdir /scratch/jhyoon/Research/LG/snapshot_e$i\_$Nsub\_run$run_num
mkdir /scratch/jhyoon/Research/LG/snapshot_e$i\_$Nsub\subhalo\_run$run_num
cp Pal5_plot.ps Pal5_plot_e$i\_$Nsub\_run$run_num.ps
cp part001 /scratch/jhyoon/Research/LG/part001_e$i\_$Nsub\_run$run_num 
cp part_peri /scratch/jhyoon/Research/LG/part_peri_e$i\_$Nsub\_run$run_num 
cp snapshot/* /scratch/jhyoon/Research/LG/snapshot_e$i\_$Nsub\_run$run_num/
cp snapshot_subhalo/* /scratch/jhyoon/Research/LG/snapshot_e$i\_$Nsub\subhalo\_run$run_num/

  end

end
end
