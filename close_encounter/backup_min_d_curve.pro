; Find points on the curve of which distance to initial points are minimum
; dxy : a criterion for the convergence
pro min_d_curve,x_curve,y_curve,x_pts,y_pts,x_proj=x_proj,y_proj=y_proj,dxy=dxy

Npts = N_elements(x_pts)
Ncurve = N_elements(x_curve)
x_curve = double(x_curve)
y_curve = double(y_curve)
x_pts = double(x_pts)
y_pts = double(y_pts)
x_proj = dblarr(Npts)
y_proj = dblarr(Npts)

if not keyword_set(dxy) then dxy=1d-12

;=== find minimum distance point ===
for i=0,Npts-1 do begin
  d = sqrt((x_pts[i]-x_curve)^2.+(y_pts[i]-y_curve)^2.)
  sort_sub = sort(d);dmin = min(d,min_sub)
  min_sub = sort_sub[0]
  x_curve_min1 = x_curve[sort_sub[0]]
  y_curve_min1 = y_curve[sort_sub[0]]
  dydx_min1 = (y_pts[i]-y_curve_min1)/(x_pts[i]-x_curve_min1)
  slope_compare = dydx_min1*((y_pts[i]-y_curve)/(x_pts[i]-x_curve))
  sub = where(slope_compare lt 0)

  if (sub[0] eq -1) then begin
    x_curve_min2 = x_curve[sort_sub[1]]
    y_curve_min2 = y_curve[sort_sub[1]]
  endif else begin
    dmin2 = min(d[sub],min_sub2)
    x_curve_min2 = x_curve[sub[min_sub2]]
    y_curve_min2 = y_curve[sub[min_sub2]]
  endelse

  derr = 999.
  while (derr gt dxy) do begin
    min_d1 = sqrt((x_pts[i]-x_curve_min1)^2.+(y_pts[i]-y_curve_min1)^2.)
    min_d2 = sqrt((x_pts[i]-x_curve_min2)^2.+(y_pts[i]-y_curve_min2)^2.)
    x_curve_half = (x_curve_min1+x_curve_min2)/2.d
    y_curve_half = (y_curve_min1+y_curve_min2)/2.d
    d_half = sqrt((x_pts[i]-x_curve_half)^2.+(y_pts[i]-y_curve_half)^2.)
    derr = abs(min_d2-d_half)/min_d2
;if abs(x_pts[i]-237.055) lt 0.1 then print,i,derr,x_curve_min1,x_curve_min2,x_curve_half
    if (d_half lt min_d1) then begin
      x_curve_min2 = x_curve_min1
      y_curve_min2 = y_curve_min1
      x_curve_min1 = x_curve_half
      y_curve_min1 = y_curve_half
;if abs(x_pts[i]-237.055) lt 0.1 then print,'case 1',i,derr,min_d1,min_d2,d_half
    endif else begin
      x_curve_min2 = x_curve_half
      y_curve_min2 = y_curve_half
;if abs(x_pts[i]-237.055) lt 0.1 then print,'case 2',i,derr,min_d1,min_d2,d_half
    endelse
  endwhile
  x_proj[i] = x_curve_min1
  y_proj[i] = y_curve_min1
endfor

END
