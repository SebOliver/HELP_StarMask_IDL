pro star_artefacts_spike_azimuthal, starcat, nstars,title, ang_fit_params,back=back,status=status,scale

;----------------------------------------------------------------------
; working out an azimuthal profile of a difraction spike
; 
; quite a lot of this code is legacy stuff that is not required and
; should be deleted. At heart its quite simple
; 
;----------------------------------------------------------------------
; star-source x-correlated catalogue structure
; nstars: number of stars used in cross-correlation
; title, string for labeling plot
; ang_fit_params: 
; theta : angle of difraction spike
; width : width of diffraction spike
; back  : input estimate of the background level

; parameters 
;----------------------------------------------------------------------

theta=ang_fit_params[1]
width=ang_fit_params[2]*5

reliability=0.8    ; reliabilty target for objects at boundaries of star mask (this is a key parameter which 

frac=0.1           ; first guess of fraction of maximum radius used as anulus with for background determination
snr=3.            ; SNR threshold for determining that boundary (only required for range definition)
back_accuracy=0.03 ; target accuracy for background accuracy (only required for range definition)

;----------------------------------------------------------------------
; BACKGROUND estimation (same as for circle code)
;----------------------------------------------------------------------

; maximum radius of separation star-object
rmax=max(starcat.separation)


;----------------------------------------------------------------------
; THRESHOLDS
; 
; We can now look at some of the values we will need to define the
; thresholds and boundaries we will be trying to find.
; exploit some maths that probably should be written here 
; but at the moment just taken from my note book
;----------------------------------------------------------------------

; so then the target density for the boundary is (e+b)
; and excess denisty e


threshold_density_e_b = back/reliability  ; reliability is the fraction of sources that are real (the background) 
                                           ; compared to all objects


; now the area that we need to measure this threshold such that excess
; density is measured above SNR given error from total counts (assumes
; Poisson error bars
 
omega_min=snr^2*reliability^2/(1-reliability)/(back*nstars)

;----------------------------------------------------------------------
; now the code changes for diffraction spikes 
;----------------------------------------------------------------------


;----------------------------------------------------------------------
; but first we have to localise the peak and measure the width
;----------------------------------------------------------------------


;plot,starcat.theta,starcat.separation,psym=3
; angular distance from spike (wrapped to be +- 180)
dtheta=(starcat.theta-theta)
cirrange,dtheta
dtheta=dtheta-180.

; *linear) distance from spike centre r*dtheta/radians

dtheta_lin=dtheta*starcat.separation*!pi/180

; objects that are within 40 widths from spike centre
in_spike=where(abs(dtheta_lin) lt 15*width,nspike)

dtheta_bin=0.5
y=histogram(dtheta_lin[in_spike],bin=dtheta_bin,locat=x)

; area of strip and normalisation of counts
omega=dtheta_bin*sqrt(rmax^2-x^2)
y=y/omega/nstars

plot,x,y,psym=10,title='!17'+string(theta)

; fitting with Gaussian + linear slope
Result = GAUSSFIT( X, Y, A, NTERMS=5)
;print,a
oplot,x,result,colo=2



;checking we have reasonably sensible values

if abs(a[1])/ang_fit_params[2] gt 3 or (abs(a[2]-ang_fit_params[2]) gt 3) then begin 
   message,'Error',/inf
;   print,ang_fit_params
;   print,a
   stop
endif else begin
   
   dtheta_lin=dtheta_lin-a[1]
   width=5*a[2]
   a[1]=a[1]+theta
   ang_fit_params=a

endelse


;----------------------------------------------------------------------
; getting the scale of the spike, relative to the general artefact high
;----------------------------------------------------------------------


artefact=ang_fit_params[3]-back
scale= ang_fit_params[0]/artefact
;stop

return
end
