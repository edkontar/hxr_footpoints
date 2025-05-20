function kappa, x, k
norm = gamma(k+1.)/(sqrt(!PI)*k^(3./2.)*gamma(k-1./2.))
return, norm * (1+x*x/k)^(-k-1.)
end

function density, x

e=4.8e-10
kev=1.6e-9
LogColloumn=20.
KK=2.*!PI*e^4*LogColloumn/(kev*kev)
print,kk
;stoptestet

X0=10.*7.25*1e7; 10 arcseconds
n0=1e11
density=n0+1e17*exp(-abs(x-x0)^2/1.38e7^2)
density=n0+1e17*((x+0.001*x0)/(x0))^(20)+1e13/(1.+x/(x0/10))^4

;IF (x LE 0 ) then density =n0*(x/1e9)^(1.8)
;IF (x GE x0 ) then density =n0*exp((x-X0)/0.2e8)
;IF (x GE 0 ) then density =n0*exp((x-X0)/5e8)
return, density
end


pro hxr_footpoints

x=findgen(900)*20.*7.25*1e7/900. ;60 arcseconds
N=dblarr(n_elements(x))
dens=dblarr(n_elements(x))
;for i=0, n_elements(x) do N(i)=QSIMP('density', 0., x(i),Jmax=20)
for i=0, n_elements(x)-1 do begin
N(i)=QROMB('density', 0., x(i),eps=1e-4)
dens(i)=density(x(i))
end


window,0
!P.Multi=[0,1,2]
plot,x/7.25e7,N/1e17,/ylog,xtitle='distance',ytitle='Colomn Depth 1e17'

plot,x/7.25e7,dens,/ylog,xtitle='distance',ytitle='Density'
!P.Multi=0

stop

;N=499
;e=findgen(N)*1.00+3.5
restore,'D:\idl\loop_sim\haug_brm_cross.dat'
cs=p.cross_section
e=p.ee

Flux=fltarr(n_elements(e),n_elements(x))
Ph  =fltarr(n_elements(e),n_elements(x))

for j=0,n_elements(x)-1 do begin
for i=0,n_elements(e)-1 do begin
N_st=3.84e17
;Flux(i,j)=1e10*E(i)*((1./2.)*exp(-sqrt(E(i)*E(i)+N(j)/N_st)/2.)+sqrt(E(i)*E(i)+N(j)/N_st)^(-(4+1)))
Flux(i,j)=1e10*E(i)*(sqrt(E(i)*E(i)+N(j)/N_st)^(-(4+1)))

end
Ph(*,j)=cs##(Flux(*,j)*dens(j))
end

peak=fltarr(n_elements(e))
peak2=fltarr(n_elements(e))
for i=0,n_elements(e)-1 do begin
pv=max(ph(i,*),max_pos)
peak(i)=x(max_pos)

pv=max(flux(i,*)*dens,max_pos)
peak2(i)=x(max_pos)

end

window,1
!P.Multi=[0,1,2]
plot_oo,e,total(ph,2),xtitle='Energy, keV'
oplot,e,cs##((1e5*exp(-e/1.3)+e^(-2))*1e20),line=2,color=40000
oplot,e,ph(*,0), line=1
oplot,e,ph(*,50), line=2
oplot,e,ph(*,100), line=3
oplot,e,ph(*,150), line=4
oplot,e,ph(*,200), line=5
oplot,e,ph(*,250), line=6
plot,e,-deriv(alog(e),alog(total(ph,2))),thick=3,yrange=[0,6],/xlog
!P.Multi=0

window,2
plot,x/7.25e7,ph(1,*)/max(ph(1,*)),xtitle='distance, arcsec',ytitle='Photon Flux',xrange=[0,10]
oplot,x/7.25e7,ph(2,*)/max(ph(2,*)),line=0
oplot,x/7.25e7,ph(3,*)/max(ph(3,*)),line=0
oplot,x/7.25e7,ph(4,*)/max(ph(4,*)),line=0
oplot,x/7.25e7,ph(7,*)/max(ph(7,*)),line=0
oplot,x/7.25e7,ph(15,*)/max(ph(15,*)),line=1
oplot,x/7.25e7,ph(20,*)/max(ph(20,*)),line=2
oplot,x/7.25e7,ph(50,*)/max(ph(50,*)),line=3
oplot,x/7.25e7,ph(100,*)/max(ph(100,*)),line=4
oplot,x/7.25e7,ph(300,*)/max(ph(300,*)),line=5


window,3
;plot,e,peak/7.25e7,/xlog
plot,e,20.-peak/7.25e7,/xlog,psym=3,yrange=[0,4],xrange=[10,300],xstyle=1
oplot,e,peak2/7.25e7,line=2,color=40000

CS=total(ph(*,0:20),2)
FP=total(ph(*,200:299),2)
LS=total(ph(*,20:200),2)

window,7
plot_oo,e,CS,xrange=[3,300],xstyle=1
oplot,e,FP,line=1
oplot,e,LS,line=2,color=20000
oplot,e,CS+FP+LS,line=3
oplot,e,total(ph(*,50:100),2),line=0



stop

end
