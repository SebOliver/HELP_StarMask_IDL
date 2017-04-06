;----------------------------------------------------------------------
; The HELP Star mask generation IDL code
;----------------------------------------------------------------------

;----------------------------------------------------------------------
; some generic plotting
;----------------------------------------------------------------------





starcat=mrdfits('./xmm_lss_star_match_all_120.fits',1)

;----------------------------------------------------------------------
; key parameters 
;----------------------------------------------------------------------

field_string='xmm_lss'

;----------------------------------------------------------------------
; things that we are expecting to find in starcat structure
; 
;  starcat is a cross match between a reference star catalogue
;  e.g. Gaia and a catalogue of objects which includes artefacts that
;  we are trying to mask out
;
;  star_id: an identifier for the reference star should be 1-N where N
;  is number of stars in this cross-match
;
;  separation: separation in arc sec between reference star and
;  catalogue object
;  theta: separation in angle (-180<degrees<180) between ref star and
;  catalogue object
;  dra: separation in arc sec on the sky in RA direction between ref
;  star and catalogue object
;  ddec: as dra but for Dec.
; 
;----------------------------------------------------------------------

;----------------------------------------------------------------------
; basic statistics and checks on input catalogue
;----------------------------------------------------------------------

; maximum radius of separation star-object
rmax=max(starcat.separation)

; check star ids are as they should be
star_id=starcat.star_id
uid=star_id(uniq(star_id,sort(star_id)))
if min(star_id) ne min(uid) then message,'Problem with star ids'
if max(star_id) ne max(uid) then message,'Problem with star ids'
nstars=max(star_id)



;----------------------------------------------------------------------
; Basic diagnostic plots
;----------------------------------------------------------------------
set_plot_ps,field_string+'_star_analysis_general.ps'
colours

;----------------------------------------------------------------------
; 2D histogram of the radial and transverse components plotted in
; image  (contours didn't work very well, might be better
; things you could do with kenerl density estimators, but this is not
; used in masking)
;----------------------------------------------------------------------

h2d = HIST_2D(starcat.separation^2,starcat.theta, bin1=200, bin2=4,min1=0,max1=120^2,min2=-180,max2=180)
icplot,h2d

;----------------------------------------------------------------------
; Normal  Delta RA, Delta Dec plot  
;----------------------------------------------------------------------


h2d = HIST_2D(starcat.dra,starcat.ddec, bin1=2, bin2=2,min1=-120,max1=120,min2=-120,max2=120)
icplot,h2d
contour,h2d,levels=[10,20,30,40,50]
contour,alog10(h2d)

endps

;----------------------------------------------------------------------
; this part does an angular projection which is where I (by hand)
; define the angles I want to work with 
;----------------------------------------------------------------------
set_plot_ps,field_string+'_star_spike_theta_all.ps'

h=histogram(starcat.theta,location=x)

plot,x,h,/ys,xtitle='!7h',ytitle='!17Counts'
oplot,!x.crange,median(h)*[1,1]
moms=moment(h)
oplot,!x.crange,moms[0]*[1,1]
oplot,!x.crange,(moms[0]+2*sqrt(moms[1]))*[1,1]


theta_spike=[-135.0,-90.0,-45,0,+45,+90,135,180]
width_spike=replicate(1.8,n_elements(theta_spike)) ; width is 3*sigma of a Gaussian fit to spike profile.

scale_spike=fltarr(n_elements(theta_spike))

for i=0,n_elements(theta_spike)-1 do  oplot,[1,1]*theta_spike[i],!y.crange,lines=1,colo=2


endps


;----------------------------------------------------------------------
; setting up output file to capture results
;----------------------------------------------------------------------

get_lun,unit
openw,unit,field_string+'_analysis_data.txt'


;----------------------------------------------------------------------
;----------------------------------------------------------------------
; now looping around magnitude bins
;----------------------------------------------------------------------
;----------------------------------------------------------------------


; setting up the magnitude limits and a string to be used in plot titles
mag_min=[10,6,7,8,9,10,11]
mag_max=[12,7,8,9,10,11,12]


; setting up some arrays to store values calculated in subroutine calls

limits_circle=fltarr(n_elements(mag_min))
zero_circle=fltarr(n_elements(mag_min))
hole_circle=fltarr(n_elements(mag_min))

limits_spike=fltarr(n_elements(theta_spike),n_elements(mag_min))
nstars_arr=fltarr(n_elements(mag_min))

for i=0,n_elements(mag_min)-1 do begin 

;----------------------------------------------------------------------
; reading in files (different for different magnitude ranges as we
; used a larger search radius for brighter stars)
;----------------------------------------------------------------------

   if mag_max[i] gt 9. then starcat=mrdfits('./xmm_lss_star_match_all_120.fits',1) $
   else starcat=mrdfits('./xmm_lss_star_match_r9_480.fits',1) 

;----------------------------------------------------------------------
; creating a string to use for plots
;----------------------------------------------------------------------

   title='!17'+string(mag_min[i],format='(i2)')+'!7<!17R!7<!17'+string(mag_max[i],format='(i2)')

;----------------------------------------------------------------------
; constructing the relevant sub-sample and (since this is a
; cross-matched catalogue of star-object pairs) re-constructing the
; list of unique stars
;----------------------------------------------------------------------

   samp=where(starcat.rmag le mag_max[i] and starcat.rmag ge mag_min[i] )

; tycho IDS
ty=starcat[samp].tyc2
; unique list of tycho ids
uni=uniq(ty,sort(ty))
; then count number of stars used
nstars=n_elements(uni)
nstars_arr[i]=nstars

;----------------------------------------------------------------------
; This routine defines the central hole and artefact model for each
; magnitude bin
;----------------------------------------------------------------------

set_plot_ps,field_string+'_star_analysis_circle'+string(mag_min[i],format='(i2.2)')+'-'+string(mag_max[i],format='(i2.2)')+'.ps'

star_artefacts_circle, starcat[samp],nstars,title, circle_fit_params,limit,zero,hole_x,x_min=x_min,x_max=x_max

; storing the output
limits_circle[i]=limit
zero_circle[i]=zero
hole_circle[i]=hole_x

endps

;----------------------------------------------------------------------
; doing the azimuthal modelling of the spike positions, widths and
; profiles we only do this with a relatively large bin in magnitudes
; which give enough signal-to-noise ratio to do this. All we are
; looking for is the position and width of the spike and the peak
; hight of the spike relative to the level of artefacts not in the
; spike. Hopefully this applies reasonably well at all mags.
;
; it might make more logical sense to have this outside the loop but
; as it relies on the radial profile calculated above it would
; generate lots of extra code to move it outside
;----------------------------------------------------------------------

if i eq 0 then begin 

   set_plot_ps,field_string+'_star_analaysis_spike_azimuthal.ps'+string(mag_min[i],format='(i2.2)')+'-'+string(mag_max[i],format='(i2.2)')+'.ps'

   for j=0,n_elements(theta_spike) -1 do begin 

; setting up the initial parameters
      width=width_spike[j]
      theta=theta_spike[j]
      ang_fit_params=[1,theta,width/3.,0,0]

; this is the important subroutine (note that we are using the
; circular fit to the artefacts to define the background level)

      star_artefacts_spike_azimuthal, starcat[samp], nstars,title,ang_fit_params,$
                            back=circle_fit_params[2],status=status,scale
 
; storing the returned values in arrays 
      scale_spike[j]=scale
      theta_spike[j]=ang_fit_params[1]
      width_spike[j]=ang_fit_params[2]*3.

 
   endfor
   endps


printf,unit,theta_spike,format='(6f6.1)'
printf,unit,width_spike,format='(6f6.1)'
printf,unit,scale_spike,format='(6f6.1)'

;----------------------------------------------------------------------

; closing the check on whether we are using the big bright bin
endif

;----------------------------------------------------------------------
; using the scaling of the spike amplitudes and the radial profile of
; all pairs to build a model of the radial profile of each spike and
; hence a threshold radius
;----------------------------------------------------------------------

for j=0,n_elements(theta_spike) -1 do begin 

; creating a model for the spike radial profile based on the circular model
; but scaled up by a the hight of the spike relative to artefacts
      spike_fit_params=circle_fit_params
      spike_fit_params[0]=spike_fit_params[0]*scale_spike[j]
      limits_spike[j,i]= star_artefacts_threshold(spike_fit_params, reliability=relibaility)

endfor

endfor


;----------------------------------------------------------------------
; linear fitting of magnitude trends
;----------------------------------------------------------------------

mag=(mag_min+mag_max)/2.
temp=linfit(mag,alog10(limits_circle))
faint=where (mag gt 8.)
temp2=linfit(mag[faint],alog10(hole_circle[faint]*10))

;----------------------------------------------------------------------
; closing the text file
;----------------------------------------------------------------------

close,unit
free_lun,unit

;----------------------------------------------------------------------
; summary plots
;----------------------------------------------------------------------

set_plot_ps,field_string+'_star_analysis_summary.ps'

;----------------------------------------------------------------------
; plotting the circular radius vs magnitude
;----------------------------------------------------------------------

plotsym,0,/fill
plot,mag,alog10(limits_circle),psym=8,xtitle='!17 Rmag',$
     ytitle='!17log!d10!n(r!dcircle!n/!7[!17arc sec!7]!17)',ys=16 ;,yrange=[3,10]
oplot,mag,alog10(limits_circle),psym=8,color=4
oplot,!x.crange,!x.crange*temp[1]+temp[0]

oplot,mag,alog10(hole_circle*10),psym=8,color=3
oplot,!x.crange,!x.crange*temp2[1]+temp2[0],color=3

legend,['artefacts','hole*10'],psym=[8,8],colors=[4,3],/top,/right,box=0

;----------------------------------------------------------------------
; plotting the constraints on spike radii
;----------------------------------------------------------------------


plot,mag,(limits_spike[0,*])/(limits_circle),psym=8,xtitle='!17 Rmag',$
   ytitle='!17!n(r!dspike!n/r!dcircle!n)',/ys;a,yrange=[-3,3]

for i=0,n_elements(theta_spike)-1 do  oplot,(mag_min+mag_max)/2.,(limits_spike[i,*])/(limits_circle),psym=8,color=3+i


legend,string(round(theta_spike)),colors=findgen(n_elements(theta_spike))+3,psym=8,charsize=1,/bottom,/right,box=0

endps

end
