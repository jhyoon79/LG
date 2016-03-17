;===
PRO convert_xyz_radec,x,y,z,ra=ra,dec=dec

l = atan(y/(x+8.))*!radeg 
b = atan(z/sqrt((x+8.)^2.+y^2.))*!radeg
lcosb = l*cos(b*!dtor)
dec = asin( cos(b*!dtor)*sin((l-33.)*!dtor)*sin(62.6*!dtor)+sin(b*!dtor)*cos(62.6*!dtor) )*!radeg
ra = asin( (cos(b*!dtor)*sin((l-33.)*!dtor)*cos(62.6*!dtor)-sin(b*!dtor)*sin(62.6*!dtor))/cos(dec*!dtor) )*!radeg + 282.25

END


PRO DENnVr,x,y,z,vx,vy,vz,Etot,ra,dec,ra_proj,dec_proj,tail=tail,Npts=Npts,$
          m_phi,surf_den,surf_den_err,sigma_Vr,slope,slope_err

;=== estimate density along the orbit ===
IF (tail EQ 'leading') THEN BEGIN
  tailpart = where(Etot-Etot[0] lt 0)
ENDIF ELSE IF (tail EQ 'trailing') THEN BEGIN
  tailpart = where(Etot-Etot[0] gt 0)
ENDIF

Npts_inbin = Npts
Ntailpart = N_ELEMENTS(tailpart)
Nbin = (Ntailpart/Npts_inbin eq Ntailpart/double(Npts_inbin)) ? Ntailpart/Npts_inbin:Ntailpart/Npts_inbin+1
phi = angsep(ra[0],dec[0],ra_proj[tailpart],dec_proj[tailpart])
sort_proj = sort(phi)
al = ra_proj[tailpart[sort_proj]]
dl = dec_proj[tailpart[sort_proj]]
ra_part = ra[tailpart]
dec_part = dec[tailpart]

;=== estimate sigma_v along the orbit ===
x_sun = -8.d
v = sqrt(vx^2.+vy^2.+vz^2.)
cos_theta = ((x-x_sun)*vx+y*vy+z*vz) / (sqrt((x-x_sun)^2.+y^2.+z^2.)*abs(v))
Vradial = v*cos_theta
Vtangen = v*sqrt(1-cos_theta^2)
Vr = Vradial[tailpart]

m_phi = dblarr(Nbin)
surf_den = dblarr(Nbin)
surf_den_err = dblarr(Nbin)
sigma_Vr = dblarr(Nbin)
slope = dblarr(Nbin)
slope_err = dblarr(Nbin)
for i=0,Nbin-1 do begin
  initial = i*Npts_inbin
  final = initial + Npts_inbin-1
  if (final gt Ntailpart-1) then final = Ntailpart-1
  ai = al[initial]
  di = dl[initial]
  af = al[final]
  df = dl[final]
  N_bin = final-initial+1
  if final-initial ge 1 then begin 
    m_phi[i] = mean(phi[sort_proj[initial:final]])
    sigma_Vr[i] = stddev(Vr[sort_proj[initial:final]])
  endif else begin
    m_phi[i] = (phi[sort_proj[initial:final]])
    sigma_Vr[i] = (Vr[sort_proj[initial:final]])
  endelse
  area = angsep(ai,di,af,df)^2.
  surf_den[i] = N_bin/area
  surf_den_err[i] = sqrt(N_bin)/area

  fit = poly_fit(ra_part[initial:final],dec_part[initial:final],1,sigma=sigma)
  slope[i] = fit[1]
  slope_err[i] = sigma[1] 
endfor

END


PRO compare

  SET_PLOT,'ps'
  @plot_setting
  DEVICE,file='compare.ps',/color,ysize=20,yoffset=0.5
  !p.multi=[0,3,3]
;=== no subhalo ===
list = ['nosub','e6e7','e6e7_30000']

for k=0,N_elements(list)-1 do begin
  dir = '/media/SEADISK/LG/FinalRun/inVL_'+list[k]+'/'
  if k eq 0 then dir = '/media/SEADISK/LG/realistic/'
  chr_rdtbl,dir+'snapshot/snap08440',0,arr,/silent
;  if k eq 0 then chr_rdtbl,dir+'snapshot/snap08437',0,arr,/silent
  if k eq 0 then chr_rdtbl,'../inVL_inner/snap08440',0,arr,/silent
;  if k eq 1 then chr_rdtbl,'snap08440',0,arr,/silent
  arr = double(arr)
  t_pts = (arr[0,0])[0]
  x_pts = reform(arr[1,*])
  y_pts = reform(arr[2,*])
  z_pts = reform(arr[3,*])
  vx_pts = reform(arr[4,*])
  vy_pts = reform(arr[5,*])
  vz_pts = reform(arr[6,*])
  Etot_pts = reform(arr[8,*])

  chr_rdtbl,dir+'part001',0,arr,/silent
  arr = double(arr)
  t_orbit = reform(arr[0,*])
  x_orbit = reform(arr[1,*])
  y_orbit = reform(arr[2,*])
  z_orbit = reform(arr[3,*])

;=== convert XYZ into RA & Dec ===
  convert_xyz_radec,x_pts,y_pts,z_pts,ra=ra_pts,dec=dec_pts
  convert_xyz_radec,x_orbit,y_orbit,z_orbit,ra=ra_orbit,dec=dec_orbit

  sub_orbit = where(abs(t_orbit-t_pts) lt 120)	; orbit within 120Myr
  ra_orbit_sub = ra_orbit[sub_orbit]
  dec_orbit_sub = dec_orbit[sub_orbit]

  plot,[0],[0],/nodata,xr=[240,200],yr=[-55,8],psym=3,/isotropic,title=list[k]
 ; plot,[0],[0],/nodata,xr=[max(ra_pts),min(ra_pts)],yr=[min(dec_pts),max(dec_pts)],psym=3,/isotropic,title=list[k]
  oplot,ra_orbit_sub,dec_orbit_sub,color=100
  oplot,ra_pts,dec_pts,psym=3
  min_d_curve,ra_orbit_sub,dec_orbit_sub,ra_pts,dec_pts,x_proj=ra_proj,y_proj=dec_proj

  Npts = 150 
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='leading',Npts=Npts,$
        m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err,sigma_Vr_pts_l,slope_pts_l,slope_pts_l_err
  DENnVr,x_pts,y_pts,z_pts,vx_pts,vy_pts,vz_pts,Etot_pts,ra_pts,dec_pts,ra_proj,dec_proj,tail='trailing',Npts=Npts,$
        m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,sigma_Vr_pts_t,slope_pts_t,slope_pts_t_err

xmin = 0  &  xmax = 26;max([m_phi_pts_l,m_phi_pts_t])
if (k ge 4) then xmax = 100
case k of
  0: begin
       ymin = 1  &  ymax = 1000
     end
  1: begin
       ymin = 1  &  ymax = 1000
     end
  2: begin
       ymin = 1  &  ymax = 1000
     end
  3: begin
       ymin = 1  &  ymax = 300
     end
  4: begin
       ymin = .1  &  ymax = 800
     end
  5: begin
       ymin = .1  &  ymax = 1000
     end
  6: begin
       ymin = .1  &  ymax = 200
     end
endcase
  plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\Sigma [n/deg^2]'),/ylog
  oploterror,m_phi_pts_l,surf_den_pts_l,surf_den_pts_l_err
  oploterror,m_phi_pts_t,surf_den_pts_t,surf_den_pts_t_err,linestyle=0,color=255,errcolor=255
IF (k EQ 0) THEN  legend,['leading','trailing'],linestyle=[0,2],color=[0,255],box=0,charsize=0.8,/bottom

case k of
  0: begin
       ymin = 0  &  ymax = 20
     end
  1: begin
       ymin = 0  &  ymax = 20
     end
  2: begin
       ymin = 0  &  ymax = 20
     end
  3: begin
       ymin = 0  &  ymax = 26
     end
  4: begin
       ymin = 0  &  ymax = 41
     end
  5: begin
       ymin = 0  &  ymax = 66
     end
  6: begin
       ymin = 0  &  ymax = 20
     end
  ;plot,phi,Vradial

endcase
plot,[0],[0],/nodata,xr=[xmin,xmax],yr=[ymin,ymax],xtitle=textoidl('\Phi [deg]'),ytitle=textoidl('\sigma_{V_r} [km/s]')
print,m_phi_pts_l,sigma_Vr_pts_l
oplot,m_phi_pts_l,sigma_Vr_pts_l
oplot,m_phi_pts_t,sigma_Vr_pts_t,linestyle=0,color=255
print,'l',stddev(sigma_Vr_pts_l)
print,'t',stddev(sigma_Vr_pts_t)
endfor
multiplot,/reset
device,/close

END
