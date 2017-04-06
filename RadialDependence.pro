;.r star_artefacts_spike
;.r star_artefacts_circle

;----------------------------------------------------------------------




set_plot_ps,'xmm_star_analysis_general.ps'
colours
starcat=mrdfits('./xmm_lss_star_match_all_120.fits',1)
;starcat=mrdfits('./xmm_lss_star_match_r9_480.fits',1)

;----------------------------------------------------------------------
; key parameters 
;----------------------------------------------------------------------

frac=0.1        ; fraction of maximum radius used as anulus with for bacground determination
reliability=0.8 ; reliabilty target for objects at boundaries of star mask
snr=3.          ; SNR threshold for determining that boundary


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

; these plots have too many points in to be useful

;plot,starcat.separation,starcat.theta,psym=3,xtitle='R',ytitle='Theta'
;plot,starcat.dra,starcat.ddec,psym=3,xtitle='DRA',ytitle='DDEC' 

;----------------------------------------------------------------------
; 2D histogram of the radial and transverse components plotted in
; image and contours (not very useful it seems)
;----------------------------------------------------------------------

h2d = HIST_2D(starcat.separation^2,starcat.theta, bin1=200, bin2=4,min1=0,max1=120^2,min2=-180,max2=180)
icplot,h2d

contour,h2d,levels=[10,20,30,40,50]
contour,alog10(h2d)

;----------------------------------------------------------------------
; this part does an angular projection
;----------------------------------------------------------------------

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

;----------------------------------------------------------------------
; Normal  Dra,Ddec plot  
;----------------------------------------------------------------------


h2d = HIST_2D(starcat.dra,starcat.ddec, bin1=2, bin2=2,min1=-120,max1=120,min2=-120,max2=120)
icplot,h2d
contour,h2d,levels=[10,20,30,40,50]
contour,alog10(h2d)


endps

;----------------------------------------------------------------------
; now looping around magnitude bins
;----------------------------------------------------------------------
get_lun,unit
openw,unit,'xmm_analysis_range.txt'
mag_min=[10,6,7,8,9,10,11]
mag_max=[12,7,8,9,10,11,12]
x_min=[100,100,50,40,30,20]
x_max=[200,200,70,60,50,40]
limits_circle=fltarr(n_elements(mag_min))
zero_circle=fltarr(n_elements(mag_min))
hole_circle=fltarr(n_elements(mag_min))

limits_spike=fltarr(n_elements(theta_spike),n_elements(mag_min))
nstars_arr=fltarr(n_elements(mag_min))

for i=0,n_elements(mag_min)-1 do begin 
;   k=k0
;   c=c0
   if mag_max[i] gt 9. then starcat=mrdfits('./xmm_lss_star_match_all_120.fits',1) $
   else starcat=mrdfits('./xmm_lss_star_match_r9_480.fits',1) 

   samp=where(starcat.rmag le mag_max[i] and starcat.rmag ge mag_min[i] )
; tycho IDS
ty=starcat[samp].tyc2
; unique list of tycho ids
uni=uniq(ty,sort(ty))
; count number of stars used
nstars=n_elements(uni)
nstars_arr[i]=nstars

;----------------------------------------------------------------------
; this is the code we are going to replace with a sub-routine
;----------------------------------------------------------------------

 
 title='!17'+string(mag_min[i],format='(i2)')+'!7<!17R!7<!17'+string(mag_max[i],format='(i2)')

h=histogram(starcat[samp].separation,bin=1,min=0.,location=x)

y=h/(2*!pi*x)/nstars


;----------------------------------------------------------------------
;----------------------------------------------------------------------

  
;----------------------------------------------------------------------
set_plot_ps,'xmm_star_analaysis_circle'+string(mag_min[i],format='(i2.2)')+'-'+string(mag_max[i],format='(i2.2)')+'.ps'
star_artefacts_circle, starcat[samp],nstars,title, circle_fit_params,limit,zero,hole_x,x_min=x_min,x_max=x_max
endps
;----------------------------------------------------------------------

printf,unit,'Range:',x_min,x_max
limits_circle[i]=limit
zero_circle[i]=zero
hole_circle[i]=hole_x

; doing the azimuthal bit of the stellar diffraction spike finding
; (only for wide range of magnitudes, should be common across all)



if i eq 0 then begin 
   set_plot_ps,'xmm_star_analaysis_spike_azimuthal.ps'+string(mag_min[i],format='(i2.2)')+'-'+string(mag_max[i],format='(i2.2)')+'.ps'
   for j=0,n_elements(theta_spike) -1 do begin 
      width=width_spike[j]
      theta=theta_spike[j]
      ang_fit_params=[1,theta,width/3.,0,0]
      star_artefacts_spike_azimuthal, starcat[samp], nstars,title,ang_fit_params,$
                            back=circle_fit_params[2],status=status,scale
      scale_spike[j]=scale
      spike_fit_params=circle_fit_params
      spike_fit_params[0]=spike_fit_params[0]*scale
      limits_spike[j,i]= star_artefacts_threshold(spike_fit_params, reliability=relibaility)

 ;     stop
      theta_spike[j]=ang_fit_params[1]
      width_spike[j]=ang_fit_params[2]*3.

   endfor
   endps
print,theta_spike,format='(6f6.1)'
print,width_spike,format='(6f6.1)'

;stop

endif

for j=0,n_elements(theta_spike) -1 do begin 
   width=2.0
   theta=theta_spike[j]
   spike_fit_params=circle_fit_params
      spike_fit_params[0]=spike_fit_params[0]*scale_spike[j]
      limits_spike[j,i]= star_artefacts_threshold(spike_fit_params, reliability=relibility)

endfor

endfor

; closing the text file
close,unit
free_lun,unit

;----------------------------------------------------------------------
; linear fitting of magnitude trends
;----------------------------------------------------------------------
mag=(mag_min+mag_max)/2.
temp=linfit(mag,alog10(limits_circle))
faint=where (mag gt 8.)
temp2=linfit(mag[faint],alog10(hole_circle[faint]*10))


;----------------------------------------------------------------------
; summary plots
;----------------------------------------------------------------------

set_plot_ps,'xmm_star_analysis_summary.ps'

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
