
; Run SExtractor, scamp and swarp to astrometrically calibrate,
; align, and coadd a set of LBC images.  You want to run this on
; one filter, one field at a time.  Images from multiple nights etc
; may be coadded - as long as you have some way of keeping track of
; the photometric calibration, or some standard image of the field
; to calibrate your output stack against.

; Routine originally written by Dave Sand and Jason Harris,
; with modifications and cleanup by Ben Weiner.  Some ideas/pieces 
; originated with John Moustakas.

; This routine may be redistributed freely; please retain the credit line.
; A brief description of the reduction is given in D. Sand et al
; (2009, ApJ, 704, 898) and will be described further in a paper
; by BJW; if you use these scripts, consider referencing these papers.


pro do_scampswarp, infile, preproc=preproc, sextractor=sextractor, $
  scamp=scamp, swarp=swarp,finfile=finfile,RA=RA,DEC=DEC, $
  astref_catalog=astref_catalog,astref_band=astref_band,psf_find=psf_find, $
  noweights=noweights,nofluxscale=nofluxscale,weightname=weightname, $
  coordtype=coordtype

; infile is list of input images
; finfile is the output mosaic image name
; preproc, sextractor, scamp, swarp, and psf_find control which
;   processing steps to do
; RA and DEC allow you to specify a pointing center for the output image
; astref_catalog and astref_band specify the astrometry reference
;   catalog for scamp, defaults are 'SDSS-R5' and 'i'
; noweights=1 means don't use weight maps when combining
; nofluxscale=1 means don't use the FLXSCALE header parameter 
;   computed by scamp to normalize images when combining.
;   if you use FLXSCALE, your output image should be in DN/sec
;   if you don't it should be in total DN ?
; weightname= use this image as the input weight image instead of 
;   <imagename>.weight.fits
; coordtype= if you want output in something other than EQUATORIAL



;let us assume that the current directory has your 'raw' data --
;reduced, cosmic ray rejected  and has weight maps
rawpath = './'
;rootpath='/home/bjw/'
;rootpath='/data1/bjw/egsy/reduce2/'
;data directory is where we copy the data to and modify the headers, etc
datapath = rawpath+'data/'
;sexpath = rawpath+'sex/'
;mosaicpath is where the output images are written
mosaicpath = rawpath+'swarp/'
;where to find configuration files
;confdir=rootpath+'conf/'
confdir=rawpath+'conf/'
sexparam = confdir+'lbc.output.param'
sexconv = confdir+'default.conv'
sexnnw = confdir+'default.nnw'
sexconfig = confdir+'lbcred.input.sex'
scampconfig = confdir+'lbc.scamp'
swarpconfig = confdir+'lbc.swarp'
fwhmconfig = confdir+'find_fwhm_sex.sex'
fwhmparam = confdir+'find_fwhm_sex.param'

if (not file_test(confdir,/directory)) then $
   stop, 'Couldnt find config file directory ',confdir
if (not file_test(datapath,/directory)) then begin
   print, 'Making data directory ',datapath
   spawn, 'mkdir -p '+datapath
endif
if (not file_test(datapath+'crap',/directory)) then begin
   spawn, 'mkdir -p '+datapath+'crap/'
endif
if (not file_test(mosaicpath,/directory)) then begin
   print, 'Making output mosaic directory ',mosaicpath
   spawn, 'mkdir -p '+mosaicpath
endif


; Don't absolutely need to use the weights for the astrometry pass
; and it seems to make sextractor 2.5.0 segfault.  Weights work in
; 2.4.4 though.  Set this to 0 if you run into problems.
useweightsforastrom=1

; How many iterations of scamp to use for astrometry
nscamppass=1L

splog, 'Starting scampswarp',systime(0), filename='splog.scampswarp',/append,/close


; Find out how many images we have and make a string list
; of their new-copy names with full datapath

nimages=0
tmplist = strarr(5000)
image=''
if (file_test(infile)) then begin
   openr,2,infile
   while (~EOF(2)) do begin
      readf, 2, format='(A0)',image
      newimage = datapath+file_basename(image)
      if (file_test(image)) then begin
         nimages = nimages+1
         tmplist[nimages-1] = newimage
      endif else begin
         print, 'Couldnt find image ',image
      endelse
   endwhile
   close, 2
   newimagelist = tmplist[0:nimages-1]
   print, 'Input list has ',nimages,' images'
endif else begin
   message, 'Couldnt open ',infile
endelse

;----------
; Preprocess images in the datapath directory - modify headers, make
; weight images if needed

if keyword_set(preproc) then begin
image=''
    i=0L
    if (file_test(infile)) then begin
        openr, 2, infile
        t0 = systime(1)
        while (~EOF(2)) do begin
            readf, 2, format='(A0)',image
            print, image
;num_1 = strpos(blahimage,'.fits')
;image=strmid(blahimage,num_1+5)
            newimage = datapath+file_basename(image)
            weightimage = repstr(newimage, '.fits','.weight.fits')
            if keyword_set(weightname) then begin
;               oldweightimage = 'wght.fits'
               oldweightimage = weightname
            endif else begin
               oldweightimage = repstr(image,'.fits','.weight.fits')
            endelse
            mwrfits, 0, weightimage, /create
            mwrfits, 0, newimage, /create
            bighdr = headfits(image,ext=0)
            for iext = 1L, 4L do begin
                splog, 'Processing '+image+', extension '+string(iext,format='(I0)')
                im = mrdfits(image,iext,oldhdr, /silent)
                if (file_test(oldweightimage)) then begin
                   imweight = mrdfits(oldweightimage,iext,wgholdhdr,/silent) 
                endif else begin
                   imweight = im-im + 1.0
                endelse
; fix any weight pixels that are negative
                iineg = where(imweight lt 1.0e-8)
                imweight[iineg] = 0.0
                imsize = size(im,/dim) & nx = imsize[0] & ny = imsize[1]
                orientation = -(sxpar(oldhdr, 'CROTA1') MOD 360.0)
 ; HACK: all images get same pixscale
;             if (not keyword_set(pixscale)) then pixscale = abs(sxpar(oldhdr,'CDELT1'))*3600.0D ; [arcsec/pixel]
             if (not keyword_set(pixscale)) then pixscale = sqrt(sxpar(oldhdr,'CD1_1')^2+sxpar(oldhdr,'CD1_2')^2)*3600.0D ; [arcsec/pixel]
             a = hogg_make_astr(sxpar(oldhdr,'CRVAL1'),sxpar(oldhdr,'CRVAL2'),$
               sxpar(oldhdr,'NAXIS1')*pixscale/3600.0D,sxpar(oldhdr,'NAXIS2')*pixscale/3600.0D,$
               orientation=orientation,pixscale=pixscale/3600.0D)
             a.crpix = [sxpar(oldhdr,'CRPIX1'),sxpar(oldhdr,'CRPIX2')] ; NOTE!
             mkhdr, hdr, im ; generate a basic FITS header
             gain = float(sxpar(oldhdr,'GAIN'))
             rdnoise = float(sxpar(oldhdr,'RDNOISE'))
             ;RA = float(sxpar(oldhdr,'OBJRA'))
             ;DEC = float(sxpar(oldhdr,'OBJDEC'))
;write the RA and DEC to a file for later use by the swarp portion of
;the code
             ;RADEC = RA+','+DEC
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
             sxaddpar, hdr, 'DATE_OBS', sxpar(bighdr,'DATE_OBS')
             sxaddpar, hdr, 'OBJECT', sxpar(oldhdr,'OBJECT')
             sxaddpar, hdr, 'GAIN', gain
             sxaddpar, hdr, 'RDNOISE', rdnoise
             sxaddpar, hdr, 'EXPTIME', float(sxpar(oldhdr,'EXPTIME'))
             sxaddpar, hdr, 'AIRMASS', float(sxpar(bighdr,'AIRMASS'))
             sxaddpar, hdr, 'FILTER', sxpar(bighdr,'FILTER')
             sxdelpar, hdr, 'HISTORY'
             putast, hdr, a

;sometimes, for some reason, CDELT1 and CDELT2 do not get put in
;right.  do it here
;             sxaddpar, hdr, 'CDELT1','-0.00006222'
;             sxaddpar, hdr, 'CDELT2','0.00006222'


;moustakas does a sky subtraction here.  i don't like that..
; he then rejects cosmic rays.  i've already done this.

;I need to change the weight image header here.
             splog, 'Writing '+newimage+', extension '+string(iext,format='(I0)')
             mwrfits, float(im), newimage,hdr
; Make the weightimage have 1.0 for exp time?

             splog, 'Writing '+weightimage+', extension '+string(iext,format='(I0)')
             mwrfits, float(imweight),weightimage,hdr 
         endfor

        endwhile



        splog, 'Total time to fix up the headers = ', (systime(1)-t0)/60.0, ' minutes.'
        splog, 'Total time to fix up the headers = ', (systime(1)-t0)/60.0, ' minutes.', filename='splog.scampswarp',/append,/close
endif else begin
    message, infile+' is not here!!'
endelse

close, 2

endif

;----------
; generate SExtractor catalogs to be used by scamp for astrometric matching

if keyword_set(sextractor) then begin
    rawimage=''
    t0=systime(1)
    openr, 2, infile
    while (~EOF(2)) do begin
        readf, 2, format='(A0)', rawimage
;        image = datapath+rawimage
        image = datapath+file_basename(rawimage)
        weightimage = repstr(image, '.fits','.weight.fits')
        bgimage = repstr(image, '.fits','_bg.fits')
        cat = repstr(image,'.fits','.cat')
        seg = repstr(image,'.fits','.seg.fits')
        print, 'SExtracting '+image
;moustakas's version has nthreads = 4
;        spawn, 'rm -f '+datapath+'*.cat'
        if (file_test(cat)) then begin
            spawn, 'rm -f '+cat
        endif
;        if keyword_set(noweights) then begin
        if (useweightsforastrom eq 0) then begin
           weightopts = ' -WEIGHT_TYPE NONE'
        endif else begin
           weightopts = ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightimage+' -WEIGHT_THRESH 0'
;           weightopts = ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightimage
        endelse
        spawn, 'sex '+image+' -c '+sexconfig+' -CATALOG_NAME '+cat+' -CATALOG_TYPE FITS_LDAC'+$
          ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 2.0 -ANALYSIS_THRESH 2.0 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
;          ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightimage+' -WEIGHT_THRESH 0'+$
;          ' -WEIGHT_TYPE NONE'+$
          weightopts+$
          ' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+bgimage+$
          ' -VERBOSE_TYPE NORMAL -DEBLEND_MINCONT 0.01 -BACK_TYPE AUTO'+$
          ' -MEMORY_BUFSIZE 4096 ',/sh


;        spawn, 'mv check.fits '+bgimage

    endwhile
    splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0,  'minutes'
    splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0,  'minutes', filename='splog.scampswarp',/append,/close
    close, 2

endif

;----------
; Run scamp on the .cat files made by SEx to correct astrometry
; and write new headers in the .head files

if keyword_set(scamp) then begin

; here we can look at everything in the data dir, or better try to
; restrict to the input list
;       catlist = file_search(datapath+'lbc*.cat')
; convert the list of images to a list of .cat files
       catlist = strarr(nimages)
       for i = 0,nimages-1 do begin
          catlist[i] = repstr(newimagelist[i],'.fits','.cat')
       endfor

; if your field is not in SDSS, you'll need to use some other catalog
; like USNO
       if not keyword_set(astref_catalog) then begin
          astref_catalog = 'SDSS-R5'
       endif
       astrinstru_key = 'FILTER'
       if not keyword_set(astref_band) then begin
;           astref_band = 'g'
           astref_band = 'i'
       endif
       photinstru_key = 'FILTER'
       magzero_key = 'PHOT_C'
       extinct_key = 'PHOT_K'
       photomflag_key = 'PHOTFLAG'

       t0 = systime(1)
; three iterations was too many and seems like one is fine
       for iter = 0L, nscamppass-1L do begin 
          case iter of
             0L: begin
                degree = '3'
; mosaic_type = 'FIX_FOCALPLANE'
                mosaic_type = 'LOOSE'
                position_maxerr = '2.0'
                posangle_maxerr = '3.0'
                aheader_suffix = '.ahead'
             end
             1L: begin
                degree = '3'
                mosaic_type = 'FIX_FOCALPLANE'
                position_maxerr = '0.5'
                posangle_maxerr = '2.0'
                aheader_suffix = '.head'
             end
             else: begin
                degree = '3'
                mosaic_type = 'FIX_FOCALPLANE'
                position_maxerr = '0.5'
                posangle_maxerr = '2.0'
                aheader_suffix = '.head'
             end
          endcase


          spawn, 'scamp '+strjoin(catlist,' ')+' -c '+scampconfig+' -CHECKPLOT_DEV PSC'+$
            ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
            ' -ASTREF_CATALOG '+astref_catalog+' -ASTREF_BAND '+astref_band+' -MERGEDOUTCAT_TYPE FITS_LDAC -SAVE_REFCATALOG Y'+$
            ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
            ' -CROSSID_RADIUS 4.0 -STABILITY_TYPE EXPOSURE -CDSCLIENT_EXEC aclient'+$
            ' -ASTRINSTRU_KEY '+astrinstru_key+' -NTHREADS 1 -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh

       endfor 
;copy new *.head files to *_bg.head
       if (file_test(datapath+'*_bg.head')) then begin
           file_delete, repstr(file_search(datapath+'lbc*.head'),'.head','_bg.head'),/allow_nonexistent,/quiet
       endif
       file_copy,file_search(datapath+'lbc*.head'),repstr(file_search(datapath+'lbc*.head'),'.head','_bg.head')


; Note that if scamp fails to find a photometric solution, the 
; FLXSCALE in the *.head files will be NAN.  This is a pain.
; We could try to test for this and (a) issue warnings, (b)
; set the FLXSCALE to something useful like 1/EXPTIME ?
; in the meantime, you can combine those with nofluxscale=1
; in the command line options.

;find median flux of the background images and put that in the *bg.head files
       openr, 2, infile
       rawimage=''
       numlines = file_lines(infile)
       fscale = fltarr(numlines)
       med_bg = fltarr(4)
       iimg = 0
       fscalenorm = 0.0
       while (~EOF(2)) do begin
           readf, 2, format='(A0)', rawimage
;           image = datapath+rawimage
           image = datapath+file_basename(rawimage)
           bgimage = repstr(image, '.fits','_bg.fits')
           bgheadname = repstr(bgimage, '.fits','.head')
           headnameimg = repstr(image,'.fits','.head')
           for iext = 1L, 4L do begin

;get image into array
               imback = mrdfits(bgimage,iext,yohdr,/silent)
;find median of array
               med_bg[iext-1L] = median(imback)
;find exposure time
               etime = sxpar(yohdr,'exptime')

           endfor
           truemed = median(med_bg)
           if iimg eq 0 then begin
               fscalenorm = truemed
           endif
           fscale[iimg] = truemed/fscalenorm
;read in scamp header file
           openr, 4, bgheadname

; Here it's setting the FLXSCALE in the bkgrd image .head file
; to be the ratio of the median to that of the first image.
           spawn, 'rm -f tmp.head'
           openw, 5,'tmp.head'
           headline = ''
           headlineimg = ''
           count=0ll
           while ~EOF(4) do begin
               readf, 4, format='(A0)',headline
               yo = strsplit(headline, ' ',/extract)

               if yo[0] eq 'FLXSCALE=' then begin

                   yo[1]=fscale[iimg]

               endif
               printf, 5,strjoin(yo,' ')
               count = count + 1
           endwhile
           close, 4
           spawn, 'rm -f '+bgheadname
           spawn, 'mv tmp.head '+bgheadname
           close, 5
           close, 7
           iimg = iimg+1
       endwhile
       print, fscale
       close, 2
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.', filename='splog.scampswarp',/append,/close

endif

;----------
; Run swarp to warp and combine images

if keyword_set(swarp) then begin

           if keyword_set(RA) and keyword_set(DEC) then begin
               radec = RA+','+DEC
               center_type = 'MANUAL' 
               centeropts = ' -CENTER_TYPE '+center_type+' -CENTER '+radec
           endif else begin
               print, 'Didnt input RA and DEC key words for image center, trying option ALL'
               centeropts = ' -CENTER_TYPE ALL'
           endelse

          pixelscale_type = 'MEDIAN' 
          pixel_scale = '0.224'
;          image_size = '7900,8300'
          header_only = 'N'
          suffix = ''

; Choice of combinetype is less than ideal.  WEIGHTED uses the weight
; maps but may not eliminate outliers.  AVERAGE is not robust to left
; over cosmic ray residuals etc.  MEDIAN is robust and ought to work
; well if combining a reasonable number of images taken in roughly
; the same photometric conditions, but I worry about small numbers
; of images or other problems.  It would be nice to have the option
; of trimming outliers and then doing a weighted mean or average.

;       combinetype = 'WEIGHTED' ; better to use AVERAGE
;       combinetype = 'AVERAGE' 
       combinetype = 'MEDIAN' 
       gain_keyword = 'GAIN'
       projection_err = '0.001'
       keywords = 'OBJECT,FILTER,SATURATE,RDNOISE,GAIN,EXPTIME,AIRMASS,TIME-OBS' ; keywords to copy

        mosaic_file = mosaicpath+finfile
       mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')
       spawn, 'rm -f '+mosaic_file
       spawn, 'rm -f '+mosaic_weightfile
       spawn, 'rm -f '+mosaicpath+'*resamp*fits*'
       spawn, 'rm -f '+mosaic_file+'.gz'
       spawn, 'rm -f '+mosaic_weightfile+'.gz'


          if file_test(mosaic_file,/regular) then begin
             splog, 'WARNING! Mosaic '+mosaic_file+' exists'
         endif
         spawn, 'rm -f '+datapath+'*bg*fits*'
; This moved the weight.fits files so they didn't get counted
; as part of the input list, no longer really necessary
         spawn, 'mv '+datapath+'*.weight.fits '+datapath+'crap/'

; using file_search here does bad stuff if there are any leftover
; files from other positions.  Instead use the newimagelist we
; made at the beginning.
;          imagelist = file_search(datapath+'lbc*.fits',count=numimg)
         imagelist = newimagelist
         numimg = nimages
;nimgs is the number of images and should go into the header
          exptot = 0
          for i=0, numimg-1 do begin
              currhdr = headfits(imagelist[i],exten=1)
              exptot = exptot+float(sxpar(currhdr,'EXPTIME'))
          endfor
          avgexp = exptot/numimg
          print, imagelist
; Moving the weight.fits back, could stop doing this now 
          spawn, 'mv '+datapath+'crap/*weight.fits '+datapath

; Construct list of weights from input list rather than everything
;          wghtimagelist = file_search(datapath+'lbc*.weight.fits')
          wghtimagelist = strarr(nimages)
          for i = 0,nimages-1 do begin
             wghtimagelist[i] = repstr(imagelist[i],'.fits','.weight.fits')
          endfor
          
          if keyword_set(noweights) then begin
             weightopts = ' -WEIGHT_TYPE NONE'
          endif else begin
             weightopts = ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits -WEIGHT_THRESH 0'
;             weightopts = ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits'
          endelse
; here it might have been nice to add a keyword that was 1/exptime
; that we could use
          if keyword_set(nofluxscale) then begin
             fluxscaleopts = ' -FSCALE_KEYWORD NONE -FSCALE_DEFAULT 1.0'
          endif else begin
             fluxscaleopts = ' -FSCALE_KEYWORD FLXSCALE -FSCALE_DEFAULT 1.0'
          endelse
          if (not keyword_set(coordtype)) then coordtype = 'EQUATORIAL'

          t0 = systime(1)
          splog, 'Building '+mosaic_file+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
;combine images.  ****make sure that background subtraction is
;identical to the performed by sextractor earlier.  

          spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
            ' -IMAGEOUT_NAME '+mosaic_file+' -WEIGHTOUT_NAME '+mosaic_weightfile+$
            weightopts+' -HEADER_SUFFIX .head'+fluxscaleopts+$
;            ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
            ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE '+coordtype+$
            ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+centeropts+$
            ' -PIXELSCALE_TYPE '+pixelscale_type+$
            ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE Y -GAIN_KEYWORD '+gain_keyword+' '+$
            ' -COMBINE_BUFSIZE 1024 -COPY_KEYWORDS '+keywords+' -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 4 -DELETE_TMPFILES N', /sh

; Add a constant background back to the mosaic.  This currently
; assumes that swarp normalized to e-/sec and it should do 
; something more intelligent like looking at flxscale?

; We can do a file_search on lbc*resamp.fits here and not find
; anything old because any earlier run of the script should have
; moved its resamp.fits files into mosaicpath dir.  relying on this is 
; less than ideal.
          backmean=0
          resam_im = file_search('lbc*.resamp.fits',count=numresamp)
          print, 'numresamp = ',numresamp
          for i=0, numresamp-1 do begin
              currhdr = headfits(resam_im[i])
              backmean = backmean+float(sxpar(currhdr,'BACKMEAN'))
;              print, 'backmean = ',float(sxpar(currhdr,'BACKMEAN'))
          endfor
;          backmean = backmean/(float(numresamp))
; if swarp normalizes the output into e-/sec, do the same with the 
; background we are adding back; divide by tottime*4 (because 4 chips)
;
; as a temporary hack, assume that if nofluxscale=1 the images were
; not divided by exptime, otherwise they were.
;          backmean = backmean/(float(numresamp))
          if keyword_set(nofluxscale) then begin
             backmean = backmean/(float(numresamp)*4.0)
          endif else begin
             backmean = backmean/(exptot*4.0)
          endelse
          print, 'tottime, backmean = ',exptot,backmean
; in case something went wrong with computing backmean
          if (backmean lt -1.0e6 or backmean gt +1.0e6) then backmean = 0.0

;combine background-only images.  DO NOT background subtract in swarp


          mosaicim = readfits(mosaic_file,mosaichdr)
          spawn, 'cp '+mosaic_file+' checkit.fits'
          mosaicim_wght = readfits(mosaic_weightfile ,mosaichdrwght)
          mosaicim = mosaicim+float(backmean)
; What to do where there is no good data?  If we set it to 0, then
; it puts holes in saturated stars, which is annoying.  Perhaps
; best to set it to 1e+5.
          mosaicim[where(mosaicim_wght lt 1.0e-5)] = 100000
;          mosaicim[where(mosaicim_wght lt 1.0e-5)] = 0.0
;add in the number of images to the header
          sxaddpar, mosaichdr, 'NUMIMG',strn(numimg),'Number of input images'
          sxaddpar, mosaichdr, 'TOTTIME',strn(exptot),'Total Exposure time'
          sxaddpar, mosaichdr, 'AVGTIME',strn(avgexp),'Avg Exposure time'
          sxaddpar,mosaidhdr, 'BACKMEAN',strn(backmean),'Mean background'
          writefits, mosaic_file,mosaicim,mosaichdr


          splog, 'Total time to build '+mosaic_file+' = ', (systime(1)-t0)/60.0, ' minutes.'
          splog, 'Total time to build '+mosaic_file+' = ', (systime(1)-t0)/60.0, ' minutes.', filename='splog.scampswarp',/append,/close



;mv resampled images to swarp directory (and gzip all the images?)
;
; may not want to delete older resamp*.fits files if they are
; from another pointing and we still plan to use them for something?
          spawn, 'rm -f '+mosaicpath+'*resamp*.fits*'

          spawn, 'mv lbc*.0001.resamp.fits '+mosaicpath
          spawn, 'mv lbc*.0001.resamp.weight.fits '+mosaicpath
          spawn, 'mv lbc*.0002.resamp.fits '+mosaicpath
          spawn, 'mv lbc*.0002.resamp.weight.fits '+mosaicpath
          spawn, 'mv lbc*.0003.resamp.fits '+mosaicpath
          spawn, 'mv lbc*.0003.resamp.weight.fits '+mosaicpath
          spawn, 'mv lbc*.0004.resamp.fits '+mosaicpath
          spawn, 'mv lbc*.0004.resamp.weight.fits '+mosaicpath

endif

;----------
; Get a very rough measure of the PSF of output image

if keyword_set(psf_find) then begin

;now finding psf of mosaiced image and constituent resampled images
; not doing resamp.fits yet, if we were, would have to worry about
; what is in the file_search list.
        mosaic_file = mosaicpath+finfile
        mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')
        catmosaic = repstr(mosaic_file,'.fits','_fwhm.cat')
        resamplist = file_search(mosaicpath+'*resamp*.fits*')
;        catresamp = repstr(resamplist,'.fits','_fwhm.cat')
        file_delete,catmosaic,/allow_nonexistent
;        file_delete,catresamp,/allow_nonexistent

;first find psf of the mosaicked image
;need to use weights here or it goes crazy at edges
            spawn, 'sex '+mosaic_file+' -c '+fwhmconfig+' -CATALOG_NAME '+catmosaic+' -PARAMETERS_NAME '+fwhmparam+' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+mosaic_weightfile

            spawn, 'rm -f cat.temp'
            spawn, 'grep -v "#" '+catmosaic+' > cat.temp'

            numlines = file_lines('cat.temp')
            sexflux=fltarr(numlines)
            sexfwhm=fltarr(numlines)
            sexx=fltarr(numlines)
            sexy=fltarr(numlines)
            sexx_sweet=fltarr(numlines)
            sexy_sweet=fltarr(numlines)
            sexflux_sweet=fltarr(numlines)
            sexfwhm_sweet=fltarr(numlines)
            rdfloat,'cat.temp',sexflux,col=3
            rdfloat,'cat.temp',sexfwhm,col=13
            rdfloat,'cat.temp',sexx,col=5
            rdfloat,'cat.temp',sexy,col=6
            spawn, 'rm -f cat.temp'


;find resistant mean FWHM (and sigma)
            i=0
            j=0
;pick out those points with flux 1e3 - 1e5 ?  This is lame
; since it depends on filter and exposure?

            while(i lt numlines) do begin
;                if ((sexflux[i] lt 400000) AND (sexflux[i] gt 100000)) then begin
                if ((sexflux[i] gt 1000) AND (sexflux[i] lt 100000)) then begin
                    sexflux_sweet[j]=sexflux[i]
                    sexfwhm_sweet[j]=sexfwhm[i]
                    sexx_sweet[j]=sexx[i]
                    sexy_sweet[j]=sexy[i]
                    j=j+1
                endif

                i=i+1
            endwhile

            sexflux_sweet = extrac(sexflux_sweet,0,j-1)
            sexfwhm_sweet = extrac(sexfwhm_sweet,0,j-1)
            sexx_sweet = extrac(sexx_sweet,0,j-1)
            sexy_sweet = extrac(sexy_sweet,0,j-1)

 
; my fwhm
;just take average of 20 lowest point that are greater than 0 (1.0
;pixel for good measure)
; also take median of all the points where fwhm < 8 pix (1.8")
            sexfwhm_sweet = sexfwhm_sweet(where(sexfwhm_sweet gt 1.0))
            indexfwhm = sort(sexfwhm_sweet)
            goodfwhm = sexfwhm_sweet(indexfwhm)
            if n_elements(goodfwhm) gt 20 then begin
                goodgoodfwhm = extrac(goodfwhm,0,19)
;                print, goodgoodfwhm
                finfwhm = moment(goodgoodfwhm)
                avgfwhm = finfwhm[0]
            endif else begin
                finfwhm = moment(goodfwhm)
                avgfwhm = finfwhm[0]
            endelse
            medfwhm = median( goodfwhm[where(goodfwhm lt 8)] )

            stringtitle='FWHM is: '+strn(avgfwhm)+' pix ('+strn(avgfwhm*0.224)+' arcsec)'
            medtitle='med FWHM is: '+strn(medfwhm)+' pix ('+strn(medfwhm*0.224)+' arcsec)'
            print,stringtitle
            print,medtitle
;make postscript file showing how good the fit is...
            cleanplot
            finps = repstr(mosaic_file,'.fits','_fwhm.eps')
            str_crap = 'rm -f '+finps
            spawn, str_crap
            set_plot, 'PS'
            device, filename=finps, /encapsulated
            plot, sexflux_sweet, sexfwhm_sweet,psym=4,xtitle = 'Flux (cnts)', ytitle = 'FWHM (pixels)',yrange=[0.0,7.0],title=stringtitle
            device, /close

           stringtitle='FWHM is: '+strn(avgfwhm)+' pix ('+strn(avgfwhm*0.224)+' arcsec)'
;make postscript file showing how good the fit is...
            cleanplot
            finps = repstr(mosaic_file,'.fits','_x_fwhm.eps')
            str_crap = 'rm -f '+finps
            spawn, str_crap
            set_plot, 'PS'
            device, filename=finps, /encapsulated
            plot, sexx_sweet, sexfwhm_sweet,psym=4,xtitle = 'Pixel x ', ytitle = 'FWHM (pixels)',yrange=[0.0,7.0],title=stringtitle
            device, /close

           stringtitle='FWHM is: '+strn(avgfwhm)+' pix ('+strn(avgfwhm*0.224)+' arcsec)'
;make postscript file showing how good the fit is...
            cleanplot
            finps = repstr(mosaic_file,'.fits','_y_fwhm.eps')
            str_crap = 'rm -f '+finps
            spawn, str_crap
            set_plot, 'PS'
            device, filename=finps, /encapsulated
            plot, sexy_sweet, sexfwhm_sweet,psym=4,xtitle = 'Pixel y ', ytitle = 'FWHM (pixels)',yrange=[0.0,7.0],title=stringtitle
            device, /close


;write PSF to header of image

            checkit = headfits(mosaic_file)
            sxaddpar,checkit,'fwhm_PSF',avgfwhm
            modfits,mosaic_file,0,checkit

; We could measure the psf of the individual resampled images here

; then delete the resampled images to save space
          spawn, 'rm -f '+mosaicpath+'*resamp*.fits*'



endif

splog, 'Ending scampswarp',systime(0), filename='splog.scampswarp',/append,/close


end
