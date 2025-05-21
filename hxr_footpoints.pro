; latest version May 21, 2025

function F0_E, ee, e0, k
  ; inputs:
  ; ee - energy in keV
  ; e0 - e0 in keV
  ; k is the spectral index for ee >> e0
  norm=1.0
  return, (1. + ee/e0) ^ (-k)
end

function inj_kappa, vel, k
  ; inputs:
  ; vel - velocity in v_t units
  ; k - spectral index
  compile_opt idl2
  norm = gamma(k + 1.) / (sqrt(!pi) * k ^ (3. / 2.) * gamma(k - 1. / 2.))
  return, norm * (1 + x * x / k) ^ (-k - 1.)
end


function density_h, h, h2 
  ; density as a function of height
  ; h heigh in cm 
  ; density is in cm^-3
 
  n0 = 1e10 ; constant density in the loop 
  ;density = n0 + 1e17 * exp(-abs(x - X0) ^ 2 / 1.38e7 ^ 2)
  ;density = n0 + 1e17 * ((x + 0.001 * X0) / (X0)) ^ (20) + 1e13 / (1. + x / (X0 / 10)) ^ 4
  h0=144*1e5 ; 144 km scale height
  density = n0 + 1e17 * exp(- h /h0) +1e15*exp(-h/h2)
  return, density
end

pro hxr_footpoints

  e = 4.8e-10
  kev = 1.6e-9 ; 1 kev in cgs units
  LogColloumn = 20.
  ;LogColloumn = 7. ; for neutral gas
  KK = 2. * !pi * e ^ 4 * LogColloumn / (kev * kev)
  h1arc=7.25e7 ; 1arcsec in cm
  h2=7.25e7*0.5
  print, KK
  ; stop
  
  restore, 'haug_brm_cross.dat'
  cs = p.cross_section ; cross-section pre-calculated to speed up
  ee = p.ee ; energy array
 
  Nz  = 249 ; number of h-points
  Nee = N_elements(ee) ; number of electron bins
  
  h = findgen(Nz)/nz*  10. * h1arc  ; 10 arcseconds
  z=reverse(h) ; reverse height
  hmax=max(h)
  dh=h[1]-h[0]
  
  NN = total(reverse(density_h(h,h2)),/cumulative)*dh
  window,0,ysize=700
  
  !P.Multi=[0,1,2]
  !P.charsize=1.5
  plot,h/h1arc,  density_h(h,h2),/ylog,ytitle=tex2idl('density, [cm$^{-3}$]'),xtitle='height [arcsec]'
  oplot,h/h1arc, density_h(h,h2/10),line=1
  
  plot,reverse(h)/h1arc, NN,/ylog,ytitle=tex2idl('Column depth, [cm$^{-2}$]'),xtitle='height [arcsec]'
  
  ; oploting the stoping depths for 10, 30, 60, keV electrons
  ; from E^2=2*kk*NN, NN_20=E^2/(2kk)
  
  oplot,[0,20],[10.,10]^2/(2.*KK),line=1
  oplot,[0,20],[30.,30]^2/(2.*KK),line=2
  oplot,[0,20],[60.,60]^2/(2.*KK),line=3

  !P.Multi=0
  !P.charsize=1
  x2jpeg,'density_h.jpeg'


  F_Eh=fltarr(Nee,Nz) 
  Ph_Eh=fltarr(Nee,Nz) ; photon spectrum

  E0=10. ;kev 
  delta=6
  Window,1
  !P.Multi=[0,1,2]
  plot_oo,ee,F0_E(ee, e0, delta),xtitle='Energy [keV]', ytitle='Flux Spectrum',$
      xrange=[3,300],/xs
  loadct,39
  For j=1,N_elements(NN)-1 do begin
    E_N=sqrt(ee^2+2.*KK*NN[j])
    F_eh[*,j]=ee*F0_E(E_N, e0, delta)/E_N
    Ph_eh[*, j] = cs ## (F_eh[*, j])*density_h(z[j],h2)*dh
    IF (j MOD 10) EQ 0 THEN oplot,ee,F_eh[*,j],color=j
  Endfor
  
;  stop
  
  ; photon spectrum
  plot_oo,ee,total(Ph_eh,2),xtitle='Energy [keV]', ytitle='Photon Flux Spectrum',xrange=[10,100],/xs
  For j=1,N_elements(NN)-1 do IF (j MOD 10) EQ 0 THEN oplot,ee,Ph_eh[*,j],color=j
  !P.Multi=0
  x2jpeg,'Flux_vs_energy.jpeg'
  
  Window,2
  !P.Multi=[0,1,2]
  plot, z/h1arc, Ph_eh[1,*]/max(Ph_eh[1,*]),ytitle=tex2idl('I(Z)')
  ; for different energies
  max_eps=ee
  max_h=fltarr(N_elements(ee))
  FOR j=1,N_elements(ee)-1 do BEGIN
    IF (j MOD 10) EQ 0 THEN oplot, z/h1arc,Ph_eh[j,*]/max(Ph_eh[j,*]),color=j
    peak=max(Ph_eh[j,*],loc)
    max_h[j]=h[N_elements(h)-loc]
  ENDFOR
  
  plot, max_eps, max_h/h1arc,ytitle=tex2idl('Footpoint heigh [arcsec]'),$
    xrange=[20,100],psym=1,xtitle='Photon Energy [keV]'
  ; for different energies
  !P.Multi=0
  
  x2jpeg,'Photon_peak.jpeg'
  
  
  print,'all plotted....  OK'
  stop

  window, 0
  !p.multi = [0, 1, 2]
  plot, x / 7.25e7, N / 1e17, /ylog, xtitle = 'distance', ytitle = 'Colomn Depth 1e17'

  plot, x / 7.25e7, dens, /ylog, xtitle = 'distance', ytitle = 'Density'
  !p.multi = 0

  stop

  ; N=499
  ; e=findgen(N)*1.00+3.5
  restore, 'D:\idl\loop_sim\haug_brm_cross.dat'
  cs = p.cross_section
  e = p.ee

  Flux = fltarr(n_elements(e), n_elements(x))
  Ph = fltarr(n_elements(e), n_elements(x))

  for j = 0, n_elements(x) - 1 do begin
    for i = 0, n_elements(e) - 1 do begin
      N_st = 3.84e17
      ; Flux(i,j)=1e10*E(i)*((1./2.)*exp(-sqrt(E(i)*E(i)+N(j)/N_st)/2.)+sqrt(E(i)*E(i)+N(j)/N_st)^(-(4+1)))
      Flux[i, j] = 1e10 * e[i] * (sqrt(e[i] * e[i] + N[j] / N_st) ^ (-(4 + 1)))
    end
    Ph[*, j] = cs ## (Flux[*, j] * dens[j])
  end

  peak = fltarr(n_elements(e))
  peak2 = fltarr(n_elements(e))
  for i = 0, n_elements(e) - 1 do begin
    pv = max(Ph[i, *], max_pos)
    peak[i] = x[max_pos]

    pv = max(Flux[i, *] * dens, max_pos)
    peak2[i] = x[max_pos]
  end

  window, 1
  !p.multi = [0, 1, 2]
  plot_oo, e, total(Ph, 2), xtitle = 'Energy, keV'
  oplot, e, cs ## ((1e5 * exp(-e / 1.3) + e ^ (-2)) * 1e20), line = 2, color = 40000
  oplot, e, Ph[*, 0], line = 1
  oplot, e, Ph[*, 50], line = 2
  oplot, e, Ph[*, 100], line = 3
  oplot, e, Ph[*, 150], line = 4
  oplot, e, Ph[*, 200], line = 5
  oplot, e, Ph[*, 250], line = 6
  plot, e, -deriv(alog(e), alog(total(Ph, 2))), thick = 3, yrange = [0, 6], /xlog
  !p.multi = 0

  window, 2
  plot, x / 7.25e7, Ph[1, *] / max(Ph[1, *]), xtitle = 'distance, arcsec', ytitle = 'Photon Flux', xrange = [0, 10]
  oplot, x / 7.25e7, Ph[2, *] / max(Ph[2, *]), line = 0
  oplot, x / 7.25e7, Ph[3, *] / max(Ph[3, *]), line = 0
  oplot, x / 7.25e7, Ph[4, *] / max(Ph[4, *]), line = 0
  oplot, x / 7.25e7, Ph[7, *] / max(Ph[7, *]), line = 0
  oplot, x / 7.25e7, Ph[15, *] / max(Ph[15, *]), line = 1
  oplot, x / 7.25e7, Ph[20, *] / max(Ph[20, *]), line = 2
  oplot, x / 7.25e7, Ph[50, *] / max(Ph[50, *]), line = 3
  oplot, x / 7.25e7, Ph[100, *] / max(Ph[100, *]), line = 4
  oplot, x / 7.25e7, Ph[300, *] / max(Ph[300, *]), line = 5

  window, 3
  ; plot,e,peak/7.25e7,/xlog
  plot, e, 20. - peak / 7.25e7, /xlog, psym = 3, yrange = [0, 4], xrange = [10, 300], xstyle = 1
  oplot, e, peak2 / 7.25e7, line = 2, color = 40000

  cs = total(Ph[*, 0 : 20], 2)
  FP = total(Ph[*, 200 : 299], 2)
  LS = total(Ph[*, 20 : 200], 2)

  window, 7
  plot_oo, e, cs, xrange = [3, 300], xstyle = 1
  oplot, e, FP, line = 1
  oplot, e, LS, line = 2, color = 20000
  oplot, e, cs + FP + LS, line = 3
  oplot, e, total(Ph[*, 50 : 100], 2), line = 0

  stop
end
