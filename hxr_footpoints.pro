; The code to produce HXR source height versus energy 
; Written with Kasia in Glasgow following the paper 
; http://dx.doi.org/10.1088/0004-637X/717/1/250
; latest version May 21, 2025

;Functions USED:
function F0_E, ee, e0, k
  ; Initial distribution of electrons N(E) injected in the corona
  ; inputs:
  ; ee - energy in keV
  ; e0 - e0 in keV
  ; k >2 is the spectral index for ee >> e0
  norm=1d36
  return, norm/E0*(k-1.)*(k-2.)/k^2*(ee/e0)/(1. + ee/e0/k)^k
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hxr_footpoints
; main code starts here...
  e = 4.8e-10
  kev = 1.6d-9 ; 1 kev in cgs units
  LogColloumn = 20.
  LogColloumn = 7. ; for neutral gas
  KK = 2. * !pi * e ^ 4 * LogColloumn / (kev * kev)
  h1arc=7.25e7 ; 1arcsec in cm
  h2=7.25e7*0.5
  print, KK
  ; stop
  
  restore, 'haug_brm_cross.dat'
  cs = p.cross_section ; cross-section pre-calculated to speed up in units 1e-27 cm^2
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


  F_Eh=dblarr(Nee,Nz) 
  Ph_Eh=dblarr(Nee,Nz) ; photon spectrum

  E0=3. ;kev 
  delta=4.
  Window,1
  !P.Multi=[0,1,2]
  plot_oo,ee,F0_E(ee, e0, delta),xtitle='Energy [keV]', ytitle='Electron Flux Spectrum',$
      xrange=[3,300],/xs
  loadct,39
  For j=0,N_elements(NN)-1 do begin
    E_N=sqrt(ee^2+2.*KK*NN[j])
    F_eh[*,j]=ee*F0_E(E_N, e0, delta)/E_N
    Ph_eh[*, j] = cs ## (F_eh[*, j])*density_h(z[j],h2)*dh/(4*!PI*1.5d13^2)*1e-27
    IF (j MOD 10) EQ 0 THEN oplot,ee,F_eh[*,j],color=j
  Endfor
  
  print, 'Total injection rate electrons/sec (above min(ee)) =', int_tabulated(ee,F0_E(ee, E0, delta))
;  stop
  
  ; photon spectrum
  plot_oo,ee,total(Ph_eh,2),xtitle='Energy [keV]', ytitle='Photon Flux Spectrum [ph/sec/cm/kev]',xrange=[3,100],/xs,$
  title='Thick target photon spectrum at 1 au'
  ;For j=1,N_elements(NN)-1 do IF (j MOD 10) EQ 0 THEN oplot,ee,Ph_eh[*,j],color=j
  !P.Multi=0
  x2jpeg,'Flux_vs_energy.jpeg'
  
  Window,2
  !P.Multi=[0,1,2]
  plot, z/h1arc, Ph_eh[1,*]/max(Ph_eh[1,*]),ytitle=tex2idl('I(h)/max(I(h))')
  ; for different energies
  max_eps=ee
  max_h=fltarr(N_elements(ee))
  FOR j=1,N_elements(ee)-1 do BEGIN
    IF (j MOD 50) EQ 0 THEN oplot, z/h1arc,Ph_eh[j,*]/max(Ph_eh[j,*]),color=j
    peak=max(Ph_eh[j,*],loc)
    max_h[j]=h[N_elements(h)-loc]
  ENDFOR
  
  plot, max_eps, max_h/h1arc,ytitle=tex2idl('HXR source height [arcsec]'),$
    xrange=[20,100],psym=1,xtitle='Photon Energy [keV]'
  ; for different energies
  !P.Multi=0
  
  x2jpeg,'Photon_peak.jpeg'
  
  
  print,'all plotted....  OK'
  
  stop

 end
 