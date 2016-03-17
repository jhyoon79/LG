pro time_series

readcol,'time_series.dat',t,b_impact,Pheat_part,Pheat_nosub

set_plot,'ps'
@plot_setting
!p.charsize=1
device,file='time_series.ps',/color
plot,t,b_impact,xtitle='t [Myr]',ytitle='b'
plot,t,Pheat_nosub,xtitle='t [Myr]',ytitle='Pheat'
oplot,t,Pheat_part,color=255
plot,t,Pheat_part-Pheat_nosub,xtitle='t [Myr]',ytitle='Pheat'
oplot,t,b_impact/max(b_impact)*2,color=255
device,/close

END
