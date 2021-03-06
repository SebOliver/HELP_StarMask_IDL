
pro star_artefacts_circle, starcat, nstars,title, fit_params, limit, zero, hole_x, x_min=x_min,x_max=x_max

;----------------------------------------------------------------------
;----------------------------------------------------------------------
; star-source x-correlated catalogue structure
; title, string for labeling plot
; x_min, x_max - radial model fitting range (outputs)
; fit_params: parameters fitted to radial number density around stars
; limit - returned value estimating best limiting radius (derived from
;         fit_params and reliability thresholds)
;
;  TBD - needs to return some sort of quality metric on fit parameters
;
;----------------------------------------------------------------------

reliability=0.8    ; reliabilty target for objects at boundaries of star mask (this is a key parameter which ultimately defines how far down the radial density profile we go

; these next parameters are used to define the bin widths using in
; density estimation and are thus not critical

frac=0.1           ; first guess of fraction of maximum radius used as anulus with for background determination
snr=10.            ; SNR threshold for determining that boundary (only required for range definition)
back_accuracy=0.03 ; target accuracy for background accuracy (only required for range definition)

;----------------------------------------------------------------------
; BACKGROUND estimation
;----------------------------------------------------------------------

; maximum radius of separation star-object
rmax=max(starcat.separation)

; hardwire in first guess a annulus width as a fraction of the rmax
for i=0,1 do begin 

   dr=frac*rmax
   rmin=(1.-frac)*rmax 
   
; area of this annulus
   omega0=!pi*(rmax^2-rmin^2)
   
   
; background calculation 
   back_region=where(starcat.separation ge rmin, back_count)
   back0=back_count/nstars/omega0 ; density in sources per star per sq arc^2
   
; what we really want is an omega such that back0*omega is high enough
; to be measured to better than (say) 3% accuracy, i.e. about 200
; objects

   omega_new=back_accuracy^(-2)/(back0*nstars)
   frac_new=omega_new/omega0*frac

;   print,frac,frac_new,back0*omega0,back0*omega0*nstars,nstars
   frac=frac_new
endfor

;----------------------------------------------------------------------
; PART 0 calculating the binnned density profile
;-------------------------------------------------------------------------------

; calculating minium area of annulus
; This us used to defined the binning for density histogram
; 
; now the area that we need to measure this threshold such that excess
; density is measured above SNR given error from total counts (assumes
; Poisson error bars
 
omega_min=snr^2*reliability^2/(1-reliability)/(back0*nstars)

; and we can scale this into a derivitive in r^2 i.e. omega =
; d(pi*r^2) 

dr2=omega_min/!pi  /2

;-------------------------------------------------------------------------------
; counting sources in bins of r^2
;----------------------------------------------------------------------

h=histogram(starcat.separation^2,bin=dr2,min=0.,location=x)
; x is the radius of the bin
x=sqrt(x)
; area of the annulus
omega=!pi*dr2

; y is counts scaled to sensible units
y=h/float(nstars)/omega ; counts per star per square arc sec


;-------------------------------------------------------------------------------
;
; PART 1
;
; using an edge detection filter to determine the threshold of the
; hole
;-------------------------------------------------------------------------------


edge=convol(y,[-1,0,+1])   
peak=max(edge,peak_index)
peak_x=x[peak_index]

edge_status=0

; hardwiring in a radius in which the star itself might be found
zero=3.

; some inelegent code to make sure we have enough bins for fitting

peak_range=where(x ge zero and x le 2*peak_x-zero,npeak_range)
if npeak_range le 3 then begin 
   edge_status=1
   peak_range=where(x ge 0 and x le 2*peak_x,npeak_range)
   if npeak_range le 3 then begin 
      edge_status=2
      peak_range=indgen(5)
   endif
endif 


; then fitting this with Gaussian
Result_hole = GAUSSFIT(x[peak_range], edge[peak_range], hole_params, NTERMS=3)
hole_x=hole_params[1]
if hole_x le 0 then message,'error'



;----------------------------------------------------------------------
; 
; PART 2
;
; model fitting to the radial decline in density of artefacts
;----------------------------------------------------------------------


;----------------------------------------------------------------------
; Fitting range, hardwired for now
;----------------------------------------------------------------------

x_min=10.
x_max=rmax

range=where(x ge x_min and x le x_max, nrange)
if nrange gt 0 then begin 

;-------------------------------------------------------------------------------



;-------------------------------------------------------------------------------
; fitting this range of data with an exponential plus constant background
;-------------------------------------------------------------------------------
result = COMFIT(X[range], Y[range], [2*back0,exp(-0.036),back0], yfit=yfit, /expon,itmax=100)
fit_params=result


; using this mode fit to estimate the required threshold denisty

limit= star_artefacts_threshold(fit_params, reliability=reliability)


;----------------------------------------------------------------------
; grouping together all the plotting at the end
;----------------------------------------------------------------------

plot,x,y,xtitle='!17Separation !7[!17arc sec!7]!17',$
     ytitle='!17Master density !7[!17(arc sec)!u-2!n!7]!17',$
     title=title,psym=10,yrange=[0,back0*3]


; overplotting a flat background level

threshold_density_e_b = fit_params[2]/reliability  ; reliability is the fraction of sources that are real (the background) 
                                           ; compared to all objects

oplot,!x.crange,back0*[1,1],color=3

oplot,x[range],yfit,color=4
oplot,!x.crange,fit_params[2]*[1,1],color=4,lines=1
plotsym,0,2,/fill
oplot,[limit],[threshold_density_e_b],psym=8,colo=4

oplot,hole_params[1]*[1,1],!y.crange,color=5
oplot,zero*[1,1],!y.crange,color=6


plot,x,edge,psym=10

oplot,hole_params[1]*[1,1],!y.crange,color=5
oplot,zero*[1,1],!y.crange,color=6


endif


return
end
