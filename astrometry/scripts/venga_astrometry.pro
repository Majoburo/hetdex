pro venga_astrometry, paramfile, plot=plot, quiet=quiet, pause=pause0, updateoutput=updateoutput, skyval_img=skyval_img, fluxconv=fluxconv

; ========================= MODIFY THIS PARAMETERS ==============================================
; ===============================================================================================
; INITIAL CENTER OF SEARCH BOXES AND THRESHOLD FOR FIBER FLUX TO BE CONSIDERED
  boxRA0=0.0 //centro relativo a las cooords usadas
  boxDEC0=0.0
  boxTHETA0=0.0
  fthreshold=-1e6
; INITIAL SIZE OF SEARCH BOXES, NUMBER OF GRID ELEMENTS PER AXIS AND NUMBER OF ITERATIONS FOR RECURSIVE SEARCH
  DboxRA=16.0
  DboxDEC=16.0
  DboxTHETA=0.01
  Ngrid0=5
  Nmax=3
; SIZE AND RESOLUTION OF FINAL HIGH RESOLUTION BOX
  DboxRAf=4.0
  DboxDECf=4.0
  DboxTHETAf=0.01
  lRA=0.1
  lDEC=0.1
  lTHETA=0.005
; ================================================================================================
; ================================================================================================

DEVICE, TRUE_COLOR=24, DECOMPOSED=0, RETAIN=2
t0=systime(1)
:Q
uuuuuuu

; In VENGA we run this on only one expo:qsure so Nexp is always 1
  Nexp=1

; 
  if (keyword_set(plot) eq 0) then plot=0
  if (keyword_set(quiet) eq 0) then quiet=0
  if not keyword_set(updateoutput) then updateoutput=0


; READ PARAMETER FILE
  openr, lun, paramfile, /get_lun
  imgfile=''
  errfile=''
  resfile=''
  file=''
  ptowfile=''
  coordfile=''
  fluxfile=''
  extfile=''
  VPseeing=''
  readf, lun, imgfile, errfile, resfile, file, ptowfile, coordfile, fluxfile, extfile, VPseeing
  free_lun, lun


tmp=(strsplit(file,'/', /regex, /extract))
outFile=tmp[n_elements(tmp)-1]

set_plot,'ps'
device, filename=outfile+'_astrometry.ps', xsize=20, ysize=10


; SDSS IMAGE AND IFU OBSERVATIONS SEEING FWHM

  fwhmSDSS = 1.36
  fwhmVP=float(VPseeing)

; INITIAL INPUT ERROR CHECK
file_check=(FILE_TEST(file+'pefsm.fits') and FILE_TEST(file+'pefw.fits') and FILE_TEST(imgfile) and FILE_TEST(errfile) and FILE_TEST(resfile))
if not file_check then message, 'ERROR: Input File Check Failed. Check PATH & NAMES of input files.'

; READ VP DATA

print,''
print,'--> READING VP DATA'

	outFile='tmp'
	ms_sub_tw_mkvpcube,file,ptowfile,coordfile,fluxfile,extfile,outFile	; MAKE TEMPORAL DATACUBE FOR SINGLE FRAME
	data=venga_readfits_v1(outFile,setup=setup)



wave0_vp=data.wave	
npix=(size(wave0_vp))[1]
nfib=(size(wave0_vp))[2]
ok=indgen(nfib)


; EXCLUDE FIBER OF WHICH WAVELENGTH=0

;; ; mask bad corner in Peter's data---REMOVE!!!!!!!!!
;; auxbad=[58, 159,174,188,203,217,218,232,233]
;; wave0_vp[*,auxbad]=0


badidx=where(wave0_vp eq 0.,count)

if count gt 0 then begin
	badfib=reform((array_indices(wave0_vp,badidx))[1,*])
	bad=badfib[uniq(badfib,sort(badfib))]
	print,''
	print,'FIBER '+strtrim(string(bad+1),2)+' IS EXCLUDED SINCE IT HAS WRONG VALUE FOR WAVELENGTH (=0).'
	count=0

	for i=0,nfib-1 do begin
;		if (i ne bad) then begin
            if (where(bad eq i) eq [-1]) then begin
			if count eq 0 then ok=[i] else ok=[ok,i]
			count++
		endif
	endfor
endif

;PRINT, ok


; UPDATE NUMBER OF FIBERS TO INCLUDE ONLY GOOD ONES
nfib=n_elements(ok)

wave_vp=dblarr(npix,nfib)
wave_vp=wave0_vp

flux_vp=data[ok].flux
err_vp=data[ok].err
wave0_vp=data[ok].wave	
ra_vp=data[ok].x
dec_vp=data[ok].y

;stop

; CENTER OF ROTATION OF THE IFU
cenra=mean(data.x)
cendec=mean(data.y)


; REFORMAT SO MANGA_ASTROMETRY WORKS (ALTHOUG IN THIS CASE EACH FIBER
; HAS ITS OWN WAVELENGTH ARRAY)

 flux=dblarr(nexp, nfib, npix)
 error=dblarr(nexp, nfib, npix)
 rafib0=dblarr(nexp, nfib, npix)
 decfib0=dblarr(nexp, nfib, npix)
 wave=dblarr(nexp, nfib, npix)


for i=0, nexp-1 do begin
   for j=0, nfib-1 do begin
      flux[i,j,*]=flux_vp[*,j]
      error[i,j,*]=err_vp[*,j]
      wave[i,j,*]=wave0_vp[*,j]
      rafib0[i,j,*]=ra_vp[j]
      decfib0[i,j,*]=dec_vp[j]
   endfor
endfor

  error[where(flux eq -666)]=1e8 ; BAD PIXELS

flux[where(flux ne -666)]=flux[where(flux ne -666)]*1d17
error[where(flux ne -666)]=error[where(flux ne -666)]*1d17


; ======== compute broad-band flux for each fiber by integrating over passband
  
; read filter transmision curve
  readcol, resfile, lamfil, resfil, /silent
  resfil[where(resfil le (1e-2)*resfil)]=0
  lam0fil=total(lamfil*resfil)/total(resfil)
  
  
  bbflux=dblarr(nexp, nfib)
  ebbflux=dblarr(nexp, nfib)
  
  for i=0, nexp-1 do begin
     statusline, "Integrating Broad-band Flux for Exposure "+strtrim(string(i+1),2)+" of "+strtrim(string(nexp),2) 

     for j=0, nfib-1 do begin
        
        auxflux=flux[i,j,*]
        auxerr=error[i,j,*]
        auxwave=wave[i,j,*]
        
                                ; mask sky lines and bad pixels
        nosel=where((auxwave ge 6363.78-10 and auxwave le 6363.78+10) or (auxwave ge 6300.304-10 and auxwave le 6300.304+10) or (auxwave ge 5889.953-10 and auxwave le 5889.953+10) or (auxwave ge 5577.338-10 and auxwave le 5577.338+10) or (auxflux eq -666), complement=sel) ;
        if (sel ne [-1]) then auxflux=auxflux[sel]
        if (sel ne [-1]) then auxerr=auxerr[sel]
        if (sel ne [-1]) then auxwave=auxwave[sel]
        
        interres=interpol(resfil,lamfil,auxwave)
        interres[where(interres le (1e-2)*max(interres))]=1e-6
        
        
        ftmpconv=auxwave*interres*auxflux/lam0fil^2 ; top integrand for flux
        etmpconv=auxwave*interres*auxerr/lam0fil^2  ; top integrand for error
        ftmpfilt=interres/auxwave                   ; bottom integrand
;        mtmp1=int_tabulated(auxwave,ftmpconv)
;        etmp1=sqrt(int_tabulated(auxwave,etmpconv^2))
;        mtmp2=int_tabulated(auxwave,ftmpfilt)
        mtmp1=tsum(auxwave,ftmpconv)
        etmp1=sqrt(tsum(auxwave,etmpconv^2))
        mtmp2=tsum(auxwave,ftmpfilt)
        bbflux[i,j]=(mtmp1/mtmp2) ; FILTER-CONVOLVED MANGA BROAD BAND FLUX [ergs/...]
        ebbflux[i,j]=(etmp1/mtmp2) ; ERROR IN BROAD BAND FLUX

        ; MONTE CARLO ERRORS
        Nmc=20
        auxbbflux=dblarr(Nmc)
        for l=0, Nmc-1 do begin
           ftmpconv=auxwave*interres*(auxflux+randomn(seed, n_elements(auxflux))*auxerr)/lam0fil^2 ; top i
           mtmp1=tsum(auxwave,ftmpconv)
           auxbbflux[l]=(mtmp1/mtmp2) ; FILTER-CONVOLVED MANGA BROAD BAND FLUX [ergs/...]
        endfor
        ebbflux[i,j]=stddev(auxbbflux)

                                ;   if (j eq 63) then stop
        
        if (n_elements(sel) le 2) then bbflux[i,j]=-666
        if (n_elements(sel) le 2) then ebbflux[i,j]=1d8

     endfor
  endfor
  

;; include 2% zero-point uncertainty
;  ebbflux=sqrt(ebbflux^2+(0.02*bbflux)^2) 
   
; ========== READ BROAD BAND IMAGE ===========
  
  img=mrdfits(imgfile, 0, hdr)
; eimg=mrdfits(errfile, 0, hdre)
 eimg=img*0.1+robust_sigma(img)                 ; FIX!!!! JUST ASUMING 2% ERROR WHILE I DON'T HAVE ERROR MAP
  
  
; trim image to IFU field for speed

adxy, hdr, min(ra_VP)-30.0/3600., min(dec_VP)-30.0/3600., minx, miny
adxy, hdr, max(ra_VP)+30.0/3600., max(dec_VP)+30.0/3600., maxx, maxy
;adxy, hdr, ra_VP, dec_VP, tmpx, tmpy

auxnx=(size(img))[1]
auxny=(size(img))[2]


hextract, img, hdr, imgtrim, hdrtrim, max([0,min([minx, maxx])]), min([auxnx,max([minx, maxx])]), max([0,min([miny, maxy])]), min([auxny,max([miny, maxy])]) , /SILENT

;adxy, hdrtrim, ra_VP, dec_VP, tmpx, tmpy

img=imgtrim
hdr=hdrtrim




; EXTRACT WCS INFO FROM THE REFERENCE IMAGE HEADER
  rdhd, hdr, structure=str, cmat=cmat, full=1 ; THINK IT WORKS ANYWAY FOR CDn_m ASTRONOMY CASE.
  img_ra=cmat.ra
  img_dec=cmat.dec
  nx_img=(size(img_ra))[1]
  ny_img=(size(img_ra))[2]
  px_size=abs(sxpar(hdr,'CD1_1')*3600.)
  print,'SDSS Pixel size in arcsec : ', px_size
  apr=4.235/px_size/2.
  print,'VENGA 4.235/2 arcsec aperture  = ',strtrim(string(apr),2),' pixels in ref image'
  
  
;; CONVOLVE TO MATCH SEEING 
  
  if (fwhmSDSS lt fwhmVP) then begin
     fwhmKERNEL=sqrt(fwhmVP^2-fwhmSDSS^2)/px_size
     kernel=PSF_GAUSSIAN(NPIXEL=51, FWHM=fwhmKERNEL, /NORMALIZE)
     imgconv=convolve(img, kernel)
     img=imgconv
  endif
  
  
;--------------------------------------------------
; CONVERT ADU TO PHYSICAL FLUX DENSITY [ergs s-1 cm-2 A-1]
; MODIFIED BY GUILLE TO WORK ON SDSS MOSAIC IMAGES
; DR8 MOSAICS ARE IN UNITS OF NanoMaggies
;

if not (keyword_set(skyval_img)) then skyval_img=0.      ; SKY VALUE OF THE REFERENCE IMAGE (USED FOR SKY SUBTRACTION BEFORE APERTURE PHOTOMETRY).
if not (keyword_set(fluxconv)) then fluxconv=3.631e-6    ; FLUX CONVERSION FACTOR FROM IMAGE UNITS TO Jy (3.631e-6 for SDSS DR8)

  

; CONVERT TO PHYSICAL FLUX DENSITY
  img=img-skyval_img 
  f_nu_img=img*fluxconv         ; f_nu [Jy = 10^-23 erg s-1 cm-2 Hz-1]
  ef_nu_img=eimg*fluxconv	; f_nu [Jy = 10^-23 erg s-1 cm-2 Hz-1]
  

; CONVERT TO F_LAMBDA
  c=299792458.0e10              ; c in A/s
  f_lambda_img=f_nu_img*replicate((c/float(lam0fil)^2)*1.e-23,(size(f_nu_img))[1],(size(f_nu_img))[2]) ; [ergs s-1 cm-2 A-1] 
  ef_lambda_img=ef_nu_img*replicate((c/float(lam0fil)^2)*1.e-23,(size(ef_nu_img))[1],(size(ef_nu_img))[2]) ; [ergs s-1 cm-2 A-1] 
  
  f_lambda_img=f_lambda_img*1e17
  ef_lambda_img=ef_lambda_img*1e17
  



  ; do aperture photometry at each pixel position
;  print, '====================================================='
;  print, '===== Doing Aperture Photometry on SDSS Images ======'
;  print, '====================================================='

;  f_lambda_img_sum=dblarr(nx_img, ny_img)
;  ef_lambda_img_sum=dblarr(nx_img, ny_img)

;  for i=ceil(apr), nx_img-ceil(apr) do begin
;     statusline, string(i)+" of "+string(nx_img-ceil(apr))
;     for j=ceil(apr), ny_img-ceil(apr) do begin
;        phpadu=1.
;        aper, f_lambda_img, i, j, auxflux, eauxflux, sky_img, skyerr_img, phpadu, apr, /flux, /NAN, setskyval=0., /silent
;        aper, ef_lambda_img^2, i, j, eauxflux, eeauxflux, sky_img, skyerr_img, phpadu, apr, /flux, /NAN, setskyval=0., /silent
;        f_lambda_img_sum[i,j]=auxflux
;        ef_lambda_img_sum[i,j]=sqrt(eauxflux)
;     endfor
;  endfor
 

  f_lambda_img_sum=f_lambda_img
  ef_lambda_img_sum=ef_lambda_img
 




  bbind=intarr(nfib)
  for i=0, nfib-1 do bbind[i]=where(abs(wave[0,i,*]-lam0fil) eq min(abs(wave[0,i,*]-lam0fil))) ; wavelength index corresponding to filter effective wavelength
  
  
; ===================================================================
; ======= REGISTER EACH EXPOSURE ====================================
; ===================================================================


;  window, 0, xsize=800, ysize=1200

pfinal=dblarr(Nexp,5)  
  
chi2best=fltarr(Nexp)

for k=0, nexp-1 do begin
; for k=1, 1 do begin
     
; select good fibers (only bright fibers and reject bad fibers)
     goodfib=where(bbflux[k,*] ne -666 and bbflux[k,*] ge fthreshold)
;     goodfib=where(bbflux[k,*] ne -666)
     

; fiducial rad and dec of fibers
; this is different for VENGA 
decfib=dblarr(nfib)    
rafib=dblarr(nfib)    

for i=0, nfib-1 do decfib[i]=decfib0[k,i,bbind[i]]
for i=0, nfib-1 do rafib[i]=rafib0[k,i,bbind[i]]

; RECURSIVE MINIMIZATION OVER SHRINKING 3D GRID IN RA, DEC, THETA
     pbest=[boxRA0,boxDEC0,boxTHETA0,1,0]
;     Nmax=5
     for n=0, Nmax-1 do begin     
        Ngrid=Ngrid0
        RAgrid=dblarr(Ngrid, Ngrid, Ngrid)
        DECgrid=dblarr(Ngrid, Ngrid, Ngrid)
        THETAgrid=dblarr(Ngrid, Ngrid, Ngrid)
        parr=dblarr(Ngrid, Ngrid, Ngrid, 5)
        chi2arr=dblarr(Ngrid, Ngrid, Ngrid)
        for i=0, Ngrid-1 do for j=0, Ngrid-1 do RAgrid[*,i,j]=dindgen(Ngrid)*(DboxRA/Ngrid)/2.0^n-(DboxRA*(1.-1./Ngrid)/2.)/2.0^n+pbest[0]
        for i=0, Ngrid-1 do for j=0, Ngrid-1 do DECgrid[i,*,j]=dindgen(Ngrid)*(DboxDEC/Ngrid)/2.0^n-(DboxDEC*(1.-1./Ngrid)/2.)/2.0^n+pbest[1]
        for i=0, Ngrid-1 do for j=0, Ngrid-1 do THETAgrid[i,j,*]=dindgen(Ngrid)*(DboxTHETA/Ngrid)/2.0^n-(DboxTHETA*(1.-1./Ngrid)/2.)/2.0^n+pbest[2]
        

        for i=0, Ngrid-1 do begin
           for j=0, Ngrid-1 do begin
              for l=0, Ngrid-1 do begin
                 statusline, strtrim(string((i*Ngrid+j)*Ngrid+l),2)+" of "+strtrim(string(Ngrid^3),1)
                 p0=[RAgrid[i,j,l], DECgrid[i,j,l], THETAgrid[i,j,l], 1, 0] ; Initial Guess for Parameters [dRA, dDEC, THETA, A, B]endfor

                 eval=func(p0, MAFLUX=reform(bbflux[k, goodfib]), EMAFLUX=reform(ebbflux[k, goodfib]), MARA=rafib[goodfib], MADEC=decfib[goodfib], IMG=f_lambda_img_sum, EIMG=ef_lambda_img_sum, HIMG=hdr, APR=apr, PLOT=0, ALLRA=rafib, ALLDEC=decfib, CENRA=cenra, CENDEC=cendec)  
                 

                 parr[i,j,l,*]=[p0[0], p0[1], p0[2], eval[1], eval[0]]
                 chi2arr[i,j,l]=eval[4]
              endfor
           endfor
        endfor

        indbest=array_indices(chi2arr, where(chi2arr eq min(chi2arr)))
        pbest=parr[indbest[0], indbest[1], indbest[2], *]
        print, "Best P: ", pbest
        chi2best[k]=(func(pbest, MAFLUX=reform(bbflux[k, goodfib]), EMAFLUX=reform(ebbflux[k, goodfib]), MARA=rafib[goodfib], MADEC=decfib[goodfib], IMG=f_lambda_img_sum, EIMG=ef_lambda_img_sum, HIMG=hdr, APR=apr, /PLOT, ALLRA=rafib, ALLDEC=decfib, CENRA=cenra, CENDEC=cendec))[4] 
       
;         window, 1, xsize=900, ysize=600
        !p.multi=[0,3,2,0]
        
                                ; DRA VS DDEC
        contour, chi2arr[*,*,indbest[2]], RAgrid[*,*,indbest[2]], DECgrid[*,*,indbest[2]], nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaRA'), ytitle=textoidl('\DeltaDEC')
        contour, chi2arr[*,*,indbest[2]], RAgrid[*,*,indbest[2]], DECgrid[*,*,indbest[2]], nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, RAgrid, DECgrid, psym=8, color=256
        plotsym, 0, 0.1
        oplot, RAgrid, DECgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[0]], [pbest[1]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[0]], [pbest[1]], psym=8
        
                                ; DRA VD THETA
        
        contour, reform(chi2arr[*,indbest[1], *], Ngrid, Ngrid), reform(RAgrid[*,indbest[1], *], Ngrid, Ngrid), reform(THETAgrid[*,indbest[1], *], Ngrid, Ngrid), nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaRA'), ytitle=textoidl('\theta')


        contour, reform(chi2arr[*,indbest[1], *], Ngrid, Ngrid), reform(RAgrid[*,indbest[1], *], Ngrid, Ngrid), reform(THETAgrid[*,indbest[1], *], Ngrid, Ngrid), nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, RAgrid, THETAgrid, psym=8, color=256
        plotsym, 0, 0.1
        oplot, RAgrid, THETAgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[0]], [pbest[2]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[0]], [pbest[2]], psym=8
        

                                ; DDEC VD THETA
        
        contour, reform(chi2arr[indbest[0],*, *], Ngrid, Ngrid), reform(DECgrid[indbest[0],*,*], Ngrid, Ngrid), reform(THETAgrid[indbest[0],*,*], Ngrid, Ngrid), nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaDEC'), ytitle=textoidl('\theta')


        contour, reform(chi2arr[indbest[0],*,*], Ngrid, Ngrid), reform(DECgrid[indbest[0],*,*], Ngrid, Ngrid), reform(THETAgrid[indbest[0],*,*], Ngrid, Ngrid), nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, DECgrid, THETAgrid, psym=8, color=256
        plotsym, 0, 0.1
        oplot, DECgrid, THETAgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[1]], [pbest[2]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[1]], [pbest[2]], psym=8
        
        ; 1D chi2 distributions

        auxchi2arr=dblarr(Ngrid)
        for i=0, Ngrid-1 do auxchi2arr[i]=min(chi2arr[i,*,*])

        plot, RAgrid[*,0,0], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\DeltaRA'), title=textoidl('\DeltaRA=')+strtrim(string(pbest[0], format='(f12.2)'),2)      
        oplot, [pbest[0], pbest[0]], [-1e6, 1e6], linestyle=2

        for i=0, Ngrid-1 do auxchi2arr[i]=min(chi2arr[*,i,*])
        plot, DECgrid[0,*,0], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\DeltaDEC'), title=textoidl('\DeltaDEC=')+strtrim(string(pbest[1], format='(f12.2)'),2)        
        oplot, [pbest[1], pbest[1]], [-1e6, 1e6], linestyle=2
        
        for i=0, Ngrid-1 do auxchi2arr[i]=min(chi2arr[*,*,i])
        plot, THETAgrid[0,0,*], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\theta'), title=textoidl('\theta=')+strtrim(string(pbest[2], format='(f12.2)'),2)              
        oplot, [pbest[2], pbest[2]], [-1e6, 1e6], linestyle=2
        
        
;        stop
        
     endfor
                                ; FINAL MINIMIZATION ON HIGH RESOLUTION GRID


NgridRA=ceil(DboxRAf/lRA)
NgridDEC=ceil(DboxDECf/lDEC)
NgridTHETA=ceil(DboxTHETAf/lTHETA)

; SEPARATE THESE VALUES AND MAKE MODIFYABLE PARAMETERS!!!!

        RAgrid=dblarr(NgridRA, NgridDEC, NgridTHETA)
        DECgrid=dblarr(NgridRA, NgridDEC, NgridTHETA)
        THETAgrid=dblarr(NgridRA, NgridDEC, NgridTHETA)
        parr=dblarr(NgridRA, NgridDEC, NgridTHETA, 5)
        chi2arr=dblarr(NgridRA, NgridDEC, NgridTHETA)

        for i=0, NgridDEC-1 do for j=0, NgridTHETA-1 do RAgrid[*,i,j]=dindgen(NgridRA)*(DboxRAf/NgridRA)-(DboxRAf*(1.-1./NgridRA)/2.)+pbest[0]
        for i=0, NgridRA-1 do for j=0, NgridTHETA-1 do DECgrid[i,*,j]=dindgen(NgridDEC)*(DboxDECf/NgridDEC)-(DboxDECf*(1.-1./NgridDEC)/2.)+pbest[1]
        for i=0, NgridRA-1 do for j=0, NgridDEC-1 do THETAgrid[i,j,*]=dindgen(NgridTHETA)*(DboxTHETAf/NgridTHETA)-(DboxTHETAf*(1.-1./NgridTHETA)/2.)+pbest[2]
        
        for i=0, NgridRA-1 do begin
           for j=0, NgridDEC-1 do begin
              for l=0, NgridTHETA-1 do begin
                statusline, strtrim(string((i*NgridDEC+j)*NgridTHETA+l),2)+" of "+strtrim(string(NgridRA*NgridDEC*NgridTHETA),1)
                 p0=[RAgrid[i,j,l], DECgrid[i,j,l], THETAgrid[i,j,l], 1, 0] ; Initial Guess for Parameters [dRA, dDEC, THETA, A, B]endfor

                 eval=func(p0, MAFLUX=reform(bbflux[k, goodfib]), EMAFLUX=reform(ebbflux[k, goodfib]), MARA=rafib[goodfib], MADEC=decfib[goodfib], IMG=f_lambda_img_sum, EIMG=ef_lambda_img_sum, HIMG=hdr, APR=apr, PLOT=0, ALLRA=rafib, ALLDEC=decfib, CENRA=cenra, CENDEC=cendec)  
                 

                 parr[i,j,l,*]=[p0[0], p0[1], p0[2], eval[1], eval[0]]
                 chi2arr[i,j,l]=eval[4]
              endfor
           endfor
        endfor

        indbest=array_indices(chi2arr, where(chi2arr eq min(chi2arr)))
        pbest=parr[indbest[0], indbest[1], indbest[2], *]



        print, "Best P: ", pbest
        chi2best[k]=(func(pbest, MAFLUX=reform(bbflux[k, goodfib]), EMAFLUX=reform(ebbflux[k, goodfib]), MARA=rafib[goodfib], MADEC=decfib[goodfib], IMG=f_lambda_img_sum, EIMG=ef_lambda_img_sum, HIMG=hdr, APR=apr, /PLOT, ALLRA=rafib, ALLDEC=decfib, CENRA=cenra, CENDEC=cendec))[4] 

;         window, 2, xsize=900, ysize=600
        !p.multi=[0,3,2,0]
        
                                ; DRA VS DDEC
        contour, chi2arr[*,*,indbest[2]], RAgrid[*,*,indbest[2]], DECgrid[*,*,indbest[2]], nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaRA'), ytitle=textoidl('\DeltaDEC')
        contour, chi2arr[*,*,indbest[2]], RAgrid[*,*,indbest[2]], DECgrid[*,*,indbest[2]], nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, RAgrid, DECgrid, psym=8, color=256
        plotsym, 0, 0.1
        oplot, RAgrid, DECgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[0]], [pbest[1]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[0]], [pbest[1]], psym=8
        
                                ; DRA VD THETA
        
        contour, reform(chi2arr[*,indbest[1], *], NgridRA, NgridTHETA), reform(RAgrid[*,indbest[1], *], NgridRA, NgridTHETA), reform(THETAgrid[*,indbest[1], *], NgridRA, NgridTHETA), nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaRA'), ytitle=textoidl('\theta')


        contour, reform(chi2arr[*,indbest[1], *], NgridRA, NgridTHETA), reform(RAgrid[*,indbest[1], *], NgridRA, NgridTHETA), reform(THETAgrid[*,indbest[1], *], NgridRA, NgridTHETA), nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, RAgrid, THETAgrid, psym=8, color=256
        plotsym, 0, 0.1
        oplot, RAgrid, THETAgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[0]], [pbest[2]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[0]], [pbest[2]], psym=8
        

                                ; DDEC VD THETA
        
        contour, reform(chi2arr[indbest[0],*, *], NgridDEC, NgridTHETA), reform(DECgrid[indbest[0],*,*], NgridDEC, NgridTHETA), reform(THETAgrid[indbest[0],*,*], NgridDEC, NgridTHETA), nlevels=50, /fill, charsize=1.5, xtitle=textoidl('\DeltaDEC'), ytitle=textoidl('\theta')


        contour, reform(chi2arr[indbest[0],*,*], NgridDEC, NgridTHETA), reform(DECgrid[indbest[0],*,*], NgridDEC, NgridTHETA), reform(THETAgrid[indbest[0],*,*], NgridDEC, NgridTHETA), nlevels=10, /overplot, color=cgcolor('white')
        plotsym, 0, 0.1, /fill
        oplot, DECgrid, THETAgrid, psym=8, color=256
        plotsym, 0, 1
        oplot, DECgrid, THETAgrid, psym=8
        plotsym, 0, 2, /fill
        oplot, [pbest[1]], [pbest[2]], psym=8, color=cgcolor('red')
        plotsym, 0, 2
        oplot, [pbest[1]], [pbest[2]], psym=8
        
        ; 1D chi2 distributions

        auxchi2arr=dblarr(NgridRA)
        for i=0, NgridRA-1 do auxchi2arr[i]=min(chi2arr[i,*,*])

        plot, RAgrid[*,0,0], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\DeltaRA'), title=textoidl('\DeltaRA=')+strtrim(string(pbest[0], format='(f12.2)'),2)      
        oplot, [pbest[0], pbest[0]], [-1e6, 1e6], linestyle=2

        auxchi2arr=dblarr(NgridDEC)
        for i=0, NgridDEC-1 do auxchi2arr[i]=min(chi2arr[*,i,*])
        plot, DECgrid[0,*,0], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\DeltaDEC'), title=textoidl('\DeltaDEC=')+strtrim(string(pbest[1], format='(f12.2)'),2)        
        oplot, [pbest[1], pbest[1]], [-1e6, 1e6], linestyle=2

        auxchi2arr=dblarr(NgridTHETA)
        for i=0, NgridTHETA-1 do auxchi2arr[i]=min(chi2arr[*,*,i])
        plot, THETAgrid[0,0,*], auxchi2arr, psym=-2, yrange=[0.9*min(auxchi2arr), 1.1*max(auxchi2arr)], charsize=1.5, xtitle=textoidl('\theta'), title=textoidl('\theta=')+strtrim(string(pbest[2], format='(f12.2)'),2)              
        oplot, [pbest[2], pbest[2]], [-1e6, 1e6], linestyle=2
        

 
     ; SAVE BEST FIT PARAMETERS
     pfinal[k,*]=pbest



if (keyword_set(pasue0)) then pause


if updateoutput then begin
	print,'--> GENERATING OUTPUT FILE'
		; AIRMASS CORRECTION, FLUX CALIBRATION FOR SINGLE FRAME
		tmp=(strsplit(file,'/', /regex, /extract))
		outFile=tmp[n_elements(tmp)-1]
		output=file+'pefsma.fits'

		fluxCal_snglFrame,file,ptowfile, coordfile, fluxfile, extfile, outfile, pbest, cenra, cendec, $
				FLUX=flux,FERR=ferr,COORDS=coords,HDR_S=hdr_s,HDR_W=hdr_w


		; LEAVE INFO IN THE HEADER & CREATE '*pefsma.fits', '*pefwea.fits' FILES
		h=hdr_s
		out=flux


		for i=0,1 do begin
			sxaddpar,h,'HISTORY','------ RESULTS OF fluxCal_refImg.pro -------'
			sxaddpar,h,'HISTORY',' FLUX CALIBRATION DONE:  Y (New Flux)  =  A * X (Orig. Flux) + B'
			sxaddpar,h,'A',pbest[3]
			sxaddpar,h,'B',pbest[4]*1d-17
			sxaddpar,h,'sigma_A',0
			sxaddpar,h,'sigma_B',0
			sxaddpar,h,'Delta_RA',pbest[0]
			sxaddpar,h,'Delta_DEC',pbest[1]
			sxaddpar,h,'Delta_THETA',pbest[2]
			sxaddpar,h,'HISTORY','--------------------------------------------'
			mwrfits,out,output,h,/create

			h=hdr_w
			output=file+'pefwea.fits'
			out=ferr
		endfor

		; WRITE NEW COORDINATE FILE
		postfix=outFile
		tmp=(strsplit(coordfile,'/.', /extract))
		outFile=tmp[n_elements(tmp)-2]
		coordOutFile=outFile+'_'+postfix+'_a.txt'
		makeCoordFile,coords,coordOutFile
print, '********************************************************************************
print, 'Coordinate file ', coordOutFile, ' is made.'

endif









endfor


;stop


openw, lun, outfile+'_astrometry.out', /get_lun, width=500
printf, lun, "# EXP  DRA   DDEC   THETA   A   B   chi2    FWHM"
for i=0, Nexp-1 do printf, lun, i+1, pfinal[i,0], pfinal[i,1], pfinal[i,2], pfinal[i,3], pfinal[i,4], chi2best[i], fwhmVP, format='(i,f14.2,f14.2,f14.2,f14.2,f14.2,f14.3,f14.2)'
free_lun, lun



device,/close
set_plot,'x'


print, '============================================'
print, 'Elapsed Time: '+strtrim(string(systime(1)-t0, format='(f10.2)'),2)+' seconds'
print, '============================================'
       

;stop




end


; FUNCTION THAT DOES THE FITTING AT EACH RA, DEC, THETA POSITION AND
; RETURNS A,B,ERRORS AND CHI2

function func, p, MAFLUX=maflux, EMAFLUX=emaflux, MARA=mara, MADEC=madec, IMG=img, EIMG=eimg, HIMG=himg, APR=apr, plot=plot, ALLRA=allra, ALLDEC=alldec, CENRA=cenra, CENDEC=cendec

                nfib=n_elements(MAFLUX)

                ; ROTATE COORDINATES

                dra=(mara-cenra)*3600.*cos(!pi/180.*madec)
                ddec=(madec-cendec)*3600.
                
                dra2=dra*cos(!pi/180.*p[2])+ddec*sin(!pi/180.*p[2])
                ddec2=ddec*cos(!pi/180.*p[2])-dra*(!pi/180.*p[2])
    
                mara0=mara
                madec0=madec

                mara=dra2/3600./cos(!pi/180.*madec0)+cenra
                madec=ddec2/3600.+cendec

                ; TRANSFORM TO PIXEL COORDINATES IN IMAGE
		xc=dblarr(nfib)
		yc=dblarr(nfib)
		adxy,himg,mara+(p[0]/3600.)/cos(!pi/180.*madec),madec+p[1]/3600.,xc,yc
		adxy,himg,allra+(p[0]/3600.)/cos(!pi/180.*alldec),alldec+p[1]/3600.,xcall,ycall

                ; NORMALIZE FOR PROPER PHOTOMETRY WITH APER
		tmp_scale=min(abs(img[where(img ne 0.)]))
                img=img/tmp_scale ; TEMPORAL SCALING TO > 1, SINCE 'aper.pro' MESSES UP eflux_img CALCULATION IF <1
		eimg=eimg/tmp_scale  ; TEMPORAL SCALING TO > 1, SINCE 'aper.pro' MESSES UP eflux_img CALCULATION IF <1
                
                ; APERTURE PHOTOMETRY WITH APER
                phpadu=1.
		aper,img,xc,yc,flux_img,eflux_img0,sky_img,skyerr_img,phpadu,apr,/flux,/NAN,setskyval=0.,/silent
		aper,eimg^2,xc,yc,eflux_img,eeflux_img,esky_img,eskyerr_eimg,phpadu,apr,/flux,/NAN,setskyval=0.,/silent

                ; UNDO NORMALIZATION
		img=img*tmp_scale	; SCALE BACK
		eimg=eimg*tmp_scale	; SCALE BACK

		flux_img=(flux_img*tmp_scale)
		eflux_img=sqrt(eflux_img)*tmp_scale

;;                flux_img=interpolate(img, xc, yc)
;;                eflux_img=interpolate(eimg, xc, yc)                       

                eflux_img=sqrt(eflux_img^2+(0.02*flux_img)^2) ; ADD 2% ZP ERROR
                if (where(finite(eflux_img) ne 1) ne [-1]) then eflux_img[where(finite(eflux_img) ne 1)] = 10*median(eflux_img)
                

                ; ORIGINAL MANGA FLUXES
                maflux0=maflux
                emaflux0=emaflux

                coef=linfit(maflux0, flux_img, chisq=chi2, sigma=ecoef, sdev=sqrt(emaflux0^2+eflux_img^2))


; ======REMOVE-------
;coef=[0.1, 1.0]
;stop


                ; SCALED MANGA FLUXES
                maflux=coef[1]*maflux0+coef[0]
                emaflux=coef[1]*emaflux0

                ; image size

                ; PLOT RESULTS IF KEYORD PLOT SET TO 1
                if (keyword_set(plot)) then begin
                   loadct, 0
                  ; wset, 0
                   !p.multi=[0,1,1,0]
;                   cgimage, img[floor(min(xc))-10:ceil(max(xc))+10,floor(min(yc))-10:ceil(max(yc))+10], /KEEP_ASPECT_RATIO, position=[0.1, 0.1, 0.9, 0.9], background='white', /scale, /axes, minvalue=median(img)-3*stddev(img), maxvalue=median(img)+10.*stddev(img), title='SDSS', charsize=1 , xrange=[floor(min(xc))-10,ceil(max(xc))+10], yrange=[floor(min(yc))-10,ceil(max(yc))+10]
                   
;                   px_size=abs(sxpar(himg,'CD1_1')*3600.)
;                   tvcircle, apr, xcall, ycall, color=cgcolor('red'), /data
;                   tvcircle, apr, xc, yc, color=cgcolor('yellow'), /data
                   
                   loadct, 39
                   plotsym, 0, 1, /fill
 ;                  ploterror, maflux0, flux_img, emaflux0, eflux_img, psym=8, /iso, xtitle='(FLUX VENGA)', ytitle='FLUX SDSS', title=textoidl('\chi^2=')+strtrim(string( total(((maflux-flux_img)/sqrt(emaflux^2+eflux_img^2))^2)),2), charsize=1 ;, /xlog, /ylog, type=3
                   ploterror, maflux0, flux_img, emaflux0, eflux_img, psym=8, /iso, xtitle='(FLUX VENGA)', ytitle='FLUX SDSS', title=textoidl('\chi^2=')+strtrim(string(chi2),2), charsize=1 , /xlog, /ylog, type=3, xrange=[1e-4, 1e2], yrange=[1e-4, 1e2]
                   
                oplot, [1e-6, 1e6], [1e-6, 1e6], color=150, thick=3
                oplot, 10.0^(findgen(100)/10.-4), coef[1]*10.0^(findgen(100)/10.-4)+coef[0], color=250, thick=3          
 ;                  oplot, [-1e4, 1e4], [-1e4, 1e4], color=150, thick=3
 ;                  oplot, findgen(100)-50., p[3]*(findgen(100)-50.)+p[4], color=250, thick=3 
                   
                   sauron_colormap
!p.multi=[0,1,3,0]
                  
;                   plot_velfield_ori, xc, yc, maflux0, range=[min(maflux)-stddev(maflux), max(maflux)], charsize=1, xthick=3, ythick=3, xrange=[0.9*min(xcall), 1.1*max(xcall)], yrange=[0.9*min(ycall), 1.1*max(ycall)], title='VENGA'
;                   tvcircle, apr, xcall, ycall, color=cgcolor('red'), /data
;                   tvcircle, apr, xc, yc, color=cgcolor('yellow'), /data
                   
                   plot_velfield_ori, xc, yc, maflux, range=[min(maflux)-stddev(maflux), max(maflux)], charsize=1, xthick=3, ythick=3, xrange=[0.9*min(xcall), 1.1*max(xcall)], yrange=[0.9*min(ycall), 1.1*max(ycall)], title='A*VENGA+B'
                   tvcircle, apr, xcall, ycall, color=cgcolor('red'), /data
                   tvcircle, apr, xc, yc, color=cgcolor('yellow'), /data
                   
                   plot_velfield_ori, xc, yc, flux_img, range=[min(maflux)-stddev(maflux), max(maflux)], charsize=1, xthick=3, ythick=3, xrange=[0.9*min(xcall), 1.1*max(xcall)], yrange=[0.9*min(ycall), 1.1*max(ycall)], title='SDSS'
                   tvcircle, apr, xcall, ycall, color=cgcolor('red'), /data
                   tvcircle, apr, xc, yc, color=cgcolor('yellow'), /data
                   
               
                   plot_velfield_ori, xc, yc, flux_img-maflux, range=[-0.1*max(maflux), 0.1*max(maflux)], charsize=1, xthick=3, ythick=3, xrange=[0.9*min(xcall), 1.1*max(xcall)], yrange=[0.9*min(ycall), 1.1*max(ycall)], title='SDSS-(A*VENGA+B)'
                   tvcircle, apr, xcall, ycall, color=cgcolor('red'), /data
                   tvcircle, apr, xc, yc, color=cgcolor('yellow'), /data
                   
                endif


                ; RETURN COEFFICIENTS AND CHI2

                return, [coef, ecoef, chi2]

end



; PROCEDURE TO COLLAPSE FRAME TO 1D DATACUBE
PRO ms_sub_tw_mkvpcube, file, ptowfile, coordfile, fluxfile, extfile, outFile

nframes=1

;READ PARAMETER FILE

output=outFile+'.fits'
file_delete,output,/allow_nonexistent


;READ PTOW FILE, CHECKING FOR ORDER OF WAVELENGTH SOLUTION

openr, lun, ptowfile, /get_lun
line=''
readf, lun, line
free_lun, lun

cols=n_elements(strsplit(line, /regex, /extract))
nfib=file_lines(ptowfile)
a=dblarr(cols,nfib) ; array with the wavelength solution coefficients for each fiber

openr, lun, ptowfile, /get_lun
readf, lun, a
free_lun, lun


;READ COORDINATES FILE

b=strarr(nfib)
openr, lun, coordfile, /get_lun
readf, lun, b
free_lun, lun
aux=strarr(3,nfib)
for i=0, nfib-1 do aux[*,i]=strsplit(b[i], /regex, /extract)

coords=dblarr(2,nfib)
for i=0, nfib-1 do begin get_coords, junk, instring=aux[1,i]+' '+aux[2,i] & coords[*,i]=junk & endfor
coords[0,*]=15.*coords[0,*] ;so both ra and dec are in decimal degrees

;READ FUX CALIBRATION FILE

readcol, fluxfile, junk, lamflux, factor


;READ DATA CHECKING FOR BINNING

aux=mrdfits(file+'pefsm.fits', 0, h)
nx=(size(aux))[1]
ny=(size(aux))[2]

data=dblarr(nx,ny,nframes)
weight=dblarr(nx,ny,nframes)
error=dblarr(nx,ny,nframes)
airmass=dblarr(nx,ny,nframes)
exptime=dblarr(nx,ny,nframes)

for i=0, nframes-1 do begin
    data[*,*,i]=mrdfits(file+'pefsm.fits', 0, h)
    airmass[*,*,i]=sxpar(h, 'airmass')
    exptime[*,*,i]=sxpar(h, 'EXPTIME')
    weight[*,*,i]=mrdfits(file+'pefw.fits', 0, h)
endfor
nogap=where(weight ne 0)
error[nogap]=1.0/sqrt(weight[nogap])


;CREATE WAVELENGTH MAP

;xarr=dindgen(nx+1)
;lambda=dblarr(nx+1, nfib)
;for i=0, nfib-1 do begin & for j=0, cols-1 do lambda[*,i]=lambda[*,i]+a[j,i]*xarr^j & endfor

; changed by GB to make lambda the center of the pixel
xarr=dindgen(nx+1)
lambda0=dblarr(nx+1, nfib)
lambda=dblarr(nx, nfib)
for i=0, nfib-1 do begin & for j=0, cols-1 do lambda0[*,i]=lambda0[*,i]+a[j,i]*xarr^j & endfor
for i=0, nx-1 do lambda[i,*]=(lambda0[i,*]+lambda0[i+1,*])/2. ;place lambda at center of pixel


;COMBINE DATA AND ERRORS AND CORRECT FOR AIRMASS

readcol, extfile, elam, ecoeff, comment="#"
elam=elam*10

flux=dblarr(nx, nfib)
ferr=dblarr(nx, nfib)
nval=nframes*5

for i=0, nfib-1 do begin
    statusline, "Collapsing and Calibrating Fiber: "+strtrim(string(i+1),2)
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1
    for j=0, nx-1 do begin
        val=reform(data[j,ymin:ymax,*], nval)
        wei=reform(weight[j,ymin:ymax,*], nval)
        err=reform(error[j,ymin:ymax,*], nval)
        am=reform(airmass[j,ymin:ymax,*], nval)
        exp=reform(exptime[j,ymin:ymax,*], nval)
        lam=lambda[j,i]
        linterp, elam, ecoeff, lam, ecoeffint
        good0=where(val ne -666) ; bad pixels
    if (good0 eq [-1]) then begin flux[j,i]=-666 & flux[j,i]=-666 & continue & endif
        ; CORRECT FOR AIRMASS
        val[good0]=val[good0]*10.0^(0.4*ecoeffint*am[good0])/exp[good0]
        err[good0]=err[good0]*10.0^(0.4*ecoeffint*am[good0])/exp[good0]
        ; SIGMA CLIP VALUES TO REJECT COSMIC RAYS
        medf=median(val[good0])
        mede=median(err[good0])
        ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
        good1=where(val ne -666 and val le medf+3.*mede and val ge medf-3.*mede ); bad pixels + sigma clipped pixels
        ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
        if (good1 eq [-1] or medf eq 0) then begin flux[j,i]=-666 & flux[j,i]=-666 & continue & endif
        flux[j,i]=total(val[good1]*wei[good1])/total(wei[good1])
        ferr[j,i]=sqrt(1.0/total(1.0/err^2))


    endfor
endfor



;FLUX CALIBRATE

fluxori=flux

for i=0, nfib-1 do begin
    auxl=lambda[*,i]
    linterp, lamflux, factor, auxl, factor1
    flux[*,i]=flux[*,i]*factor1
    ferr[*,i]=ferr[*,i]*factor1
endfor

flux[where(fluxori eq -666)]=-666


;WRITE MULTI EXTENSION FITS FILE

h=['COMMENT   Flux in ergs/s/cm^2/A']
mwrfits, flux, output, /create
h=['COMMENT   Error in the Flux in ergs/s/cm^2/A']
mwrfits, ferr, output
h=['COMMENT   Wavelength of the left edge of each pixel in A']
mwrfits, lambda, output
h=['COMMENT   RA, DEC in decimal degrees']
mwrfits, coords, output

;stop

END

; Read VENGA fits file and return data structure
; Optionally WRITE a fits table
;
; Usage data=venga_readfits('NGC2903_p1d1r')
;
; Version history
; 1.0 (8-4-2009): This header
;     (8-17-2009): added keyword "SETUP" for red/blue data
; To do:
; 1) Undo specres = specres*1.6/5.0 - DONE (IM 8/17/09)
; 2) Add red/blue keyword to exclude appropriate number of fibers on
;    line 60: d=d[0:nfib-2] for red, d=d[0:nfib-3] for blue - DONE (IM 8/17/09)
;
; Modified it so wavelength is the center of the pixel instead of the
; left edge. (GB, 09-02-2009)
function venga_readfits_v1,fileroot,WRITE=write, SETUP=SETUP


file=strtrim(fileroot)+'.fits'
print,'Reading '+file

;read in the data
flux    = mrdfits(file,0,hdr)
err     = mrdfits(file,1)
wave    = mrdfits(file,2) ; left edge of each pixel. *NOT* the center
coords  = mrdfits(file,3)


nw   = n_elements(flux[*,0])

;FIX THIS CHEAP MASKING !!!!!!!!!!!!!!!!!!!!!
;mask bad pixels ;
;s=where(flux lt 0d,complement=c,nbadpix)
;s=where(flux eq -666,complement=c,nbadpix)
;if (s ne [-1]) then flux[s]=median(flux)     ; NO MASKING ANY MORE
;if (s ne [-1]) then err[s]=err[s]+median(err)+median(flux)

;last and first two CCD columns have crap data.

flux= flux[2:nw  -3,*]
err = err [2:nw  -3,*]
wave= wave[2:nw  -3,*]

nw   = n_elements(flux[*,0])
nfib = (size(flux))[2]

data=create_struct('name','fileroot','x',0d,'y',0D,'flux',dblarr(nw),'err',dblarr(nw),'wave',dblarr(nw))


d=replicate(data,nfib)

;make the wavelength at the center of the pixel instead of the left edge

;wavec=dblarr(nw+1, nfib)
;for i=0, nw-1 do wavec[i,*]=(wave[i,*]+wave[i+1,*])/2.

d.flux=flux[0:nw-1,*]
d.err =err [0:nw-1,*]
;d.wave=wavec[0:nw-1,*]
d.wave=wave[0:nw-1,*]

;d.wave=wave
d.x=reform(coords[0,*])
d.y=reform(coords[1,*])

;stop

;;use SETUP keyword to get rid of junk fibers
;;for RED
;if SETUP eq 'red' then begin
;   ;Last fiber is empty: remove it.
;   d=d[0:nfib-2]
;endif
;;for BLUE
;if SETUP eq 'blue' then begin
;   ;Last two fibers are empty: remove
;   d=d[0:nfib-3]
;endif

IF KEYWORD_set(write) THEN mwrfits,d,strtrim(fileroot)+'_table.fits',/create

return,d
end


PRO fluxCal_snglFrame, file, ptowfile, coordfile, fluxfile, extfile, outfile, p, cenra, cendec, FLUX=flux,FERR=ferr,COORDS=coords, HDR_S=hdr_s,HDR_W=hdr_w

factor_brd=p[3]
zero_brd=p[4]*1d-17
Del_RA=p[0]
Del_DEC=p[1]
Del_THETA=p[2]

;READ PARAMETER FILE

output=outFile+'pefsma.fits'
file_delete,output,/allow_nonexistent

nframes=1

;READ PTOW FILE, CHECKING FOR ORDER OF WAVELENGTH SOLUTION

openr, lun, ptowfile, /get_lun
line=''
readf, lun, line
free_lun, lun

cols=n_elements(strsplit(line, /regex, /extract))
nfib=file_lines(ptowfile)
a=dblarr(cols,nfib) ; array with the wavelength solution coefficients for each fiber

openr, lun, ptowfile, /get_lun
readf, lun, a
free_lun, lun


;READ COORDINATES FILE

b=strarr(nfib)
openr, lun, coordfile, /get_lun
readf, lun, b
free_lun, lun
aux=strarr(3,nfib)
for i=0, nfib-1 do aux[*,i]=strsplit(b[i], /regex, /extract)

coords=dblarr(2,nfib)
for i=0, nfib-1 do begin get_coords, junk, instring=aux[1,i]+' '+aux[2,i] & coords[*,i]=junk & endfor
coords[0,*]=15.*coords[0,*] ;so both ra and dec are in decimal degrees

;READ FUX CALIBRATION FILE

readcol, fluxfile, junk, lamflux, factor


;READ DATA CHECKING FOR BINNING

aux=mrdfits(file+'pefsm.fits', 0, h)
nx=(size(aux))[1]
ny=(size(aux))[2]

data=dblarr(nx,ny,nframes)
weight=dblarr(nx,ny,nframes)
error=dblarr(nx,ny,nframes)
airmass=dblarr(nx,ny,nframes)
exptime=dblarr(nx,ny,nframes)

for i=0, nframes-1 do begin
    data[*,*,i]=mrdfits(file+'pefsm.fits', 0, hdr_s)
    airmass[*,*,i]=sxpar(h, 'airmass')
    exptime[*,*,i]=sxpar(h, 'EXPTIME')
    weight[*,*,i]=mrdfits(file+'pefw.fits', 0, hdr_w)
endfor
nogap=where(weight ne 0)
error[nogap]=1.0/sqrt(weight[nogap])


;CREATE WAVELENGTH MAP

xarr=dindgen(nx+1)
lambda=dblarr(nx+1, nfib)
for i=0, nfib-1 do begin & for j=0, cols-1 do lambda[*,i]=lambda[*,i]+a[j,i]*xarr^j & endfor

;COMBINE DATA AND ERRORS AND CORRECT FOR AIRMASS

readcol, extfile, elam, ecoeff, comment="#"
elam=elam*10

flux=dblarr(nx, ny)
ferr=dblarr(nx, ny)
nval=nframes*5

for i=0, nfib-1 do begin
    ;print, "Fiber: ", (i+1)
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1
    for j=0, nx-1 do begin
        val=reform(data[j,ymin:ymax,*], nval)
        wei=reform(weight[j,ymin:ymax,*], nval)
        err=reform(error[j,ymin:ymax,*], nval)
        am=reform(airmass[j,ymin:ymax,*], nval)
        exp=reform(exptime[j,ymin:ymax,*], nval)
        lam=lambda[j,i]
        linterp, elam, ecoeff, lam, ecoeffint
        good0=where(val ne -666) ; bad pixels
    if (good0 eq [-1]) then begin flux[j,ymin:ymax]=-666 & flux[j,ymin:ymax]=-666 & continue & endif
        ; CORRECT FOR AIRMASS
	;STOP
        val[good0]=val[good0]*10.0^(0.4*ecoeffint*am[good0])/exp[good0]
        err[good0]=err[good0]*10.0^(0.4*ecoeffint*am[good0])/exp[good0]
	;STOP
        ; SIGMA CLIP VALUES TO REJECT COSMIC RAYS
        medf=median(val[good0])
        mede=median(err[good0])
        ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
        good1=where(val ne -666 and val le medf+3.*mede and val ge medf-3.*mede ); bad pixels + sigma clipped pixels
        ; IF ALL DATA POINTS ARE BAD WRITE -666 AND EXIT
        if (good1 eq [-1] or medf eq 0) then begin flux[j,ymin:ymax]=-666 & flux[j,ymin:ymax]=-666 & continue & endif

        flux[j,ymin:ymax]=val
        ferr[j,ymin:ymax]=err

    endfor
endfor

bad=where(flux eq -666,complement=good)


;FIRST FLUX CALIBRATION (STANDARD STAR)
for i=0, nfib-1 do begin
    ymin=(i+1)*8-2-1
    ymax=(i+1)*8+2-1

    auxl=lambda[*,i]
    linterp, lamflux, factor, auxl, factor1
    for j=ymin,ymax do begin
	    flux[*,j]=flux[*,j]*factor1
	    ferr[*,j]=ferr[*,j]*factor1
	endfor
endfor

flux[bad]=-666  ; PRESERVE BAD PIXEL FLAG

good=where(flux ge 0,count,complement=cmpl,ncomplement=ncmpl)

	print,''
	print,'','mean(flux)','mean(err)','mean(flux/err)',format='(a10,3a20)'
	print,'BEFORE: ',mean(flux[good]),mean(ferr[good]),mean(flux[good])/mean(ferr[good]),format='(a10,2g20.4,d20.4)'	; WAVELENGTH=0 INCLUDED

;stop

; SECOND FLUX CALIBRATION (REFERENCE IMAGE)
flux[good]=flux[good]*factor_brd+zero_brd
;ONLY APPLY THE SCALING FACTOR, NOT THE ZEROPOINT
;flux[good]=flux[good]*factor_brd
ferr[good]=sqrt((ferr[good]*factor_brd)^2+(zero_brd)^2)

	print,' AFTER : ',mean(flux[good]),mean(ferr[good]),mean(flux[good])/mean(ferr[good]),format='(a10,2g20.4,d20.4)'
	print,''


                ; ROTATE COORDINATES

                dra=(coords[0,*]-cenra)*3600.*cos(!pi/180.*coords[1,*])
                ddec=(coords[1,*]-cendec)*3600.
                
                dra2=dra*cos(!pi/180.*p[2])+ddec*sin(!pi/180.*p[2])
                ddec2=ddec*cos(!pi/180.*p[2])-dra*(!pi/180.*p[2])
    
                ra0=coords[0,*]
                dec0=coords[1,*]

                ra1=dra2/3600./cos(!pi/180.*coords[1,*])+cenra
                dec1=ddec2/3600.+cendec


coords[0,*]=ra1+(p[0]/3600.)/cos(!pi/180.*dec1)
coords[1,*]=dec1+p[1]/3600.




END

;-----------------------------------------------------------------------
; WRITE NEW COORDINATE TEXT FILE (ORIGINAL NAME+'_a.txt') WITH THE SAME FORMAT AS ORIGINAL.

PRO makeCoordFile,coords,coordOutFile

	n=(size(coords))[2]
	ra=strarr(n) & dec=strarr(n)
	for i=0,n-1 do begin
		radecArr=strsplit(adstring(coords[0,i],coords[1,i],1),' ',/extract,/regex)
		ra[i]=strjoin(radecArr[0:2],':')
		dec[i]=strjoin(radecArr[3:5],':')
	endfor

	forprint,indgen(n)+1,ra,dec,textout=coordOutFile, format='(3a14)' ,/nocomment

END

