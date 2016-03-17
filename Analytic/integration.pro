PRO integration,mode

print,KEYWORD_SET(mode)
IF (KEYWORD_SET(mode) EQ 0) THEN BEGIN
  set_plot,'ps'
  @plot_setting
  device,file='integration.ps'
ENDIF


v = 200
G = 4.3e-6
Msub = [1e5,1e6,1e7,1e8,1e9]
rs = [0.04,0.11,0.30,0.79,2.11]
Msub = Msub[3]
rs = rs[3]

r = findgen(100000)/1+1/1.
dv = 2.*G*Msub/v*(-alog(r/rs+1.)/r)
plot,r,dv,/xlog







IF (KEYWORD_SET(mode) EQ 0) THEN BEGIN
  device,/close
ENDIF

stop
END
