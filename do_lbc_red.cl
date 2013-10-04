procedure do_lbc_red(trimnover)

# this procedure allows you to reduce a collection of LBC frames
# taken in the same filter
# by Dave Sand and Jason Harris for blue, modified by Ben Weiner for Red

# uses IRAF, mscred, and lacos_im.cl from P. van Dokkum

# steps include overscan subtraction, trimming, bias,
# flattening with twilights, scale and subtract fringe frame
# mask saturated pixels, cosmic ray removal, and construction
# of weight maps.

# This routine may be redistributed freely; please retain the credit line.
# A brief description of the reduction is given in D. Sand et al
# (2009, ApJ, 704, 898) and will be described further in a paper
# by BJW; if you use these scripts, consider referencing these papers.


bool trimnover {no,prompt="Do you want to subtract overscan and trim?"}
string trimlist {"trimlist", prompt="File containing list of images to be overscanned and trimmed"}
bool bias {no,prompt="Do you want to make a master bias to zero-subtract?"}
string biaslist {"biaslist", prompt="File containing list of bias images"}
bool biassub {no, prompt="Do you want to subtract bias from your images?"}
bool twiflat {no, prompt="Create twilight flat?"}
bool fringeframe {no, prompt="Create supersky and fringe?"}
string twilist {"twilist", prompt="File containing object images used to make twilight flat."}
string fringelist {"fringelist", prompt="File containing object images used to make supersky/fringe frame."}
bool flatdo {no,prompt="Do you want to apply twilight flat?"}
string flatdolist {"flatdolist",prompt="List of fits files to flat field"}
bool fringedo {no,prompt="Do you want to apply fringe correction?"}
string fringedolist {"fringedolist",prompt="List of fits files to de-fringe"}
bool saturcorr {no,prompt="Do you want to mask saturated and bleeding pixels?"}
string bpmlist {"bpmlist",prompt="List of bad pixel mask directories matching your flatdolist"}
bool crrays {no,prompt="Do you want to run lacosmic to remove cr rays?"}
bool addmasks {no, prompt="Do you want to combine the saturation and cr masks?"}
bool dowcs {no,prompt="Fine tune WCS?"}
string catfile {"all.cat",prompt="What is the name of the USNO catalog file to fine tune the WCS"}
bool mkweight {no,prompt="Make weight map by chopping out chip 1 of flat and normalizing"}
bool dobackup {no,prompt="Backup each step to a subdirectory?"}
string weightfile {"Flat.fits",prompt="Image to construct the weight map from"}
string databasefile {"try2.db",prompt="Name of database file to use to determine WCS initially"}

begin


####Random variables that we will use later
	string img,s,s2,nlist,tanlist,tanlistnew,tanlist_dum,tanlist_dum2,bpmfile
	string express,bpmdir,hfile,mask1,mask2,maskout,finname,thetime
	int i2,j,image_x,image_y,xdiff,ydiff,hix,hiy
	real check1,check2,wghtmean

####Check to see if the 'mscred' package is loaded.  If not, load it.

	if (!defpac("stsdas")){
		stsdas
	}
	if (!defpac("playpen")){
		playpen
	}
	if (!defpac("mscred")){
		mscred
	}

		if (!access("wcsdone")){
			mkdir("wcsdone/")
		}



######If you want, make a xtalk correction file here --- To be honest, I don't see any evidence for xtalk.
#	if (xtalk) {
#		delete ("outtalkfile")
#		unlearn("xtcoeff")
		###this is not ready for prime time, must adjust parameters
#		xtcoeff(input="@"//xtalkfile,output="outtalkfile",victim="im1,im2,im3,im4,im5,im6,im7,im8",source="im2,im1,im4,im3,im6,im5,im8,im7",smin=20000.,smax=40000.,low=2.,high=2.)
#		print("Cross-talk coefficients determined")
#	}




####If trimnover, subtract and remove overscan region (and hopefully soon xtalk correct)

	if (trimnover) {

		unlearn("ccdproc")
		#set ccdproc parameters to subtract overscan and trim, ONLY
		ccdproc.ccdtype=""
		ccdproc.bpmasks=""
		ccdproc.noproc=no
		ccdproc.merge = no
		#####We need to make a xtalk correction file and put it in here
		#ccdproc.xtalkco = yes
		ccdproc.xtalkco = no
		ccdproc.xtalkfi = "outtalkfile"
		ccdproc.oversca=yes
		ccdproc.trim=yes
		ccdproc.fixpix=no
		ccdproc.zerocor=no
		ccdproc.darkcor=no
		ccdproc.flatcor=no
		ccdproc.sflatco=no
		ccdproc.biassec="!biassec"
		ccdproc.trimsec="!trimsec"
		ccdproc.interac=no
		#ccdproc.verbose=yes
		mscred.verbose=yes
		# Don't override what the user sets
		#mscred.backup="once"
		print("Subtracting overscan and trimming images")
		#Run ccdproc
		ccdproc.images="@"//trimlist
		ccdproc("@"//trimlist)
		#save this step by copying products to directory 'trim'
		if (dobackup) {
		  if (!access("trim")){
			mkdir("trim/")
		  }
		copy("@"//trimlist, "trim/")
		}
		print("Images have been trimmed and overscan subtracted")
	}

#if you want to make a finalbias.fits frame...
	if (bias){
		print("Combining bias frames...")
		if(access("Finalbias.fits")){
			imdel("Finalbias")
		}
		zerocombine(input="@"//biaslist,output="Finalbias",combine="average",reject="avsigclip",ccdtype="",process=yes,delete=no,scale="none",statsec="*",nlow=1,nhigh=1,nkeep=2)
		print("Finalbias.fits has been made")
	}


####subtract biases from frames	
	if (biassub){
		print("Now subtracting the zero frame from the trimlist")
		#now subtract zero from images
		ccdproc.zerocor=yes
		ccdproc.merge=no
		ccdproc.zero="Finalbias.fits"
		ccdproc("@"//trimlist)
		if (dobackup) {
		  if(!access("bias")){
			mkdir("bias/")
		  }
		copy("@"//trimlist,"bias/")
		}
		print("Bias subtraction of trimlist complete")
	}	

#####Make twilight flat
	if (twiflat){
		#unlearn("flatcombine")
		flatcombine.combine="median"
		flatcombine.reject="avsigclip"
		flatcombine.scale="mode"
		flatcombine.lsigma=3
		flatcombine.hsigma=3
		flatcombine.nlow=2
		flatcombine.nhigh=2
		flatcombine.nkeep=1
		print("Making twilight flats")
		flatcombine.ccdtype = ""
		flatcombine.subsets=no
		if(access("Flat.fits")){
			imdel("Flat.fits")
		}
		flatcombine("@"//twilist,output="Flat")
		print("Flat.fits has been made")

	}


####I'm going to do the saturation correction before the  flat fielding now.  lets see what happens!

if(saturcorr){
		ccdproc.saturat=63000.
		ccdproc.bleed=63000.
		ccdproc.btrail=20.
		ccdproc.bgrow=0.
		ccdproc.fixpix=no
		ccdproc.sgrow=3
		ccdproc.fixfile="BPM"
		ccdproc.merge=no
		ccdproc.bpmasks="@"//bpmlist
		print("Masking saturated and bleeding pixels")
		list=bpmlist
		while(fscan(list,img)!=EOF){
			if (!access(img)){
				mkdir(img)
			}
			cd(img)
			imdel("*.pl")
			imdel("crmask*")
			cd("../")
		}
		ccdproc("@"//flatdolist)
		if (dobackup) {
		  if(!access("satur")){
			mkdir("satur/")
		  }
  		  copy(flatdolist,"satur/")
		  cd ("satur/")
		  imdel("@"//flatdolist)
		  delete(flatdolist)
		  cd("../")
		  copy("@"//flatdolist,"satur/")
		}
		print("Saturation masking done")
}



#####Flat field data
	if (flatdo){
		print("Flatfielding the data")
		ccdproc.flatc=yes
		#this needs some tweaking..
		ccdproc.flat="Flat"
		ccdproc.merge=no
		ccdproc("@"//flatdolist)
		if (dobackup) {
		  if(!access("flat")){
			mkdir("flat/")
		  }
		  copy(flatdolist,"flat/")
		  cd("flat/")
		  imdel("@"//flatdolist)
		  delete(flatdolist)
		  cd("../")
		  #copy flatfielded data into directory 'flat'
		  copy("@"//flatdolist, "flat")
		}
	}


# Make weight map from flat field

	if (mkweight) {
		imdel("wght_tmp.fits")
		imdel("wght.fits")
		imdel("wght*")
		copy(weightfile,"wght.fits",verb+)
		print("Making weight map")
		mscsplit("wght.fits",output="",del+,verb+)
		immean("wght_1",verb+)
		imarith("wght_1","/",real(immean.mean),"wght_1")
		immean("wght_2",verb+)
		imarith("wght_2","/",real(immean.mean),"wght_2")
		immean("wght_3",verb+)
		imarith("wght_3","/",real(immean.mean),"wght_3")
		immean("wght_4",verb+)
		imarith("wght_4","/",real(immean.mean),"wght_4")
		mscjoin("wght",output="",del+,verb+)

	}


##### Make super-sky flat by combining object frames.  Median smooth it 251x251
##### and subtract to get a high-pass spatial frequency image of the fringes.
	if (fringeframe){
		#unlearn("sflatcombine")
		sflatcombine.combine="median"
		sflatcombine.reject="ccdclip"
		sflatcombine.scale="mode"
		sflatcombine.lsigma=3
		sflatcombine.hsigma=3
		sflatcombine.rdnoise=8.0
		sflatcombine.gain=1.6
		print("Making super sky flat")
		sflatcombine.ccdtype = ""
		sflatcombine.subsets=no
		if(access("Supersky.fits")){
			imdel("Supersky.fits")
		}
		sflatcombine("@"//fringelist,output="Supersky")
		mscmedian("Supersky","Supersky_medsmooth",251,251)
		mscarith("Supersky","-","Supersky_medsmooth","Supersky_highpass")
		print("Supersky.fits, Supersky_medsmooth.fits, Supersky_highpass.fits have been made")

	}


## Do fringe removal and second flattening with a supersky here?
#  Fit scale and subtract the fringe frame, then second-flat by
#  dividing by the median-smoothed supersky.
	if (fringedo){
		print("Fringe removal and supersky-flat the data")
#		rmfringe.masks="BPM"
		rmfringe.masks=""
		rmfringe.background=""
		rmfringe("@"//fringedolist,"@"//fringedolist,"Supersky_highpass")
		ccdproc.sflatcor=yes
		ccdproc.sflat="Supersky_medsmooth"
		ccdproc.merge=no
		ccdproc("@"//fringedolist)
		if (dobackup) {
		  if(!access("fringe")){
			mkdir("fringe/")
		  }
		  copy(fringedolist,"fringe/")
		  cd("fringe/")
		  imdel("@"//fringedolist)
		  delete(fringedolist)
		  cd("../")
		  #copy superflatfielded data into directory 'fringe'
		  copy("@"//fringedolist, "fringe")
		}
	}

# Run lacos_im to remove cosmic rays.  This works on single-extension
# fits files so we mscsplit the images, and make crmask files that
# can be combined with the bad pixel masks.
# When crrays finishes, the images are still split.

if (crrays) {
   	print("Starting cosmic-ray removal")
#	nlist = flatdolist
	nlist = fringedolist
	list=nlist
	if(!access("crrays")){
		mkdir("crrays/")
	}
	lacos_im.gain=1.6
	lacos_im.readn=8.0
	lacos_im.sigclip=4.5
	lacos_im.sigfrac=0.3
	lacos_im.objlim=1.0
	lacos_im.niter=2
	while(fscan(list,img)!=EOF){
		#get root of image name
		i=strlen(img)
		i2=i-5
		s=substr(img,1,i2)
		#split the image up and run lacosmic
   	        print("Removing CRs from ",img)
		mscsplit(img,output="",del+,verb+)
		imdel("temp.fits")
		lacos_im(input=s//"_1",output="temp.fits",outmask="crmask_"//s//"_1",gain=1.6,readn=8.0)
		imdel(s//"_1")
		imrename("temp.fits",s//"_1")
		mv("crmask_"//s//"_1.fits",s//"_bpm/")
		imdel("temp.fits")
		lacos_im(input=s//"_2",output="temp.fits",outmask="crmask_"//s//"_2",gain=1.6,readn=8.0)
		imdel(s//"_2")
		imrename("temp.fits",s//"_2")
		mv("crmask_"//s//"_2.fits",s//"_bpm/")
		imdel("temp.fits")
		lacos_im(input=s//"_3",output="temp.fits",outmask="crmask_"//s//"_3",gain=1.6,readn=8.0)
		imdel(s//"_3")
		imrename("temp.fits",s//"_3")
		mv("crmask_"//s//"_3.fits",s//"_bpm/")
		imdel("temp.fits")
		lacos_im(input=s//"_4",output="temp.fits",outmask="crmask_"//s//"_4",gain=1.6,readn=8.0)
		imdel(s//"_4")
		imrename("temp.fits",s//"_4")
		mv("crmask_"//s//"_4.fits",s//"_bpm/")
		#combine the saturation mask with the lacosmic mask


	}
}

# Take the crmasks and add them to bad pixel / weight masks.
# Then rejoin the mosaic image and BPM masks using mscjoin

if (addmasks) {

#TO DO: ADD IN A STATIC BAD PIXEL MASK
#			if (access ("staticmask.pl")){
#				express="max(a,b)"
#				imexpr(express,"temppl.fits","crmask_"//img,"staticmask.pl")
#				imdel( "crmask_"//img)
#				imrename("temppl.fits","crmask_"//img)
#			}

		list=bpmlist
		print("Adding together the saturation mask and the crray mask")
		express="max(a,b)"

		while(fscan(list,bpmdir)!=EOF){
			cd(bpmdir)
			imdel("bpm*final*.pl")
			cd("../")
			j=strlen(bpmdir)
			i2=j-4
			s=substr(bpmdir,1,i2)
			cd(bpmdir)

			mask1="crmask_"//s//"_1.fits"
			mask2="bpm_LBCCHIP1.pl"
			maskout="bpm"//s//"_final_1.pl"
			imexpr(express,maskout,mask1,mask2)
			hedit(maskout,"BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			imdel("../"//maskout)
			imcopy(maskout,"../"//maskout)
			cd("../")
			hedit(s//"_1.fits","BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			fixpix(s//"_1.fits",masks="BPM",linterp="INDEF",cinterp="INDEF",verbose+,pixels-)
			cd("crrays")
			imdel(s//"_1.fits")
			cd("../")
			imcopy(s//"_1.fits","crrays/")
			cd(bpmdir)

			mask1="crmask_"//s//"_2.fits"
			mask2="bpm_LBCCHIP2.pl"
			maskout="bpm"//s//"_final_2.pl"
			imexpr(express,maskout,mask1,mask2)
			hedit(maskout,"BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			imdel("../"//maskout)
			imcopy(maskout,"../"//maskout)
			cd("../")
			hedit(s//"_2.fits","BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			fixpix(s//"_2.fits",masks="BPM",linterp="INDEF",cinterp="INDEF",verbose+,pixels-)
			cd("crrays")
			imdel(s//"_2.fits")
			cd("../")
			imcopy(s//"_2.fits","crrays/")
			cd(bpmdir)



			mask1="crmask_"//s//"_3.fits"
			mask2="bpm_LBCCHIP3.pl"
			maskout="bpm"//s//"_final_3.pl"
			imexpr(express,maskout,mask1,mask2)
			hedit(maskout,"BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			imdel("../"//maskout)
			imcopy(maskout,"../"//maskout)
			cd("../")
			hedit(s//"_3.fits","BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			fixpix(s//"_3.fits",masks="BPM",linterp="INDEF",cinterp="INDEF",verbose+,pixels-)
			cd("crrays")
			imdel(s//"_3.fits")
			cd("../")
			imcopy(s//"_3.fits","crrays/")
			cd(bpmdir)


			mask1="crmask_"//s//"_4.fits"
			mask2="bpm_LBCCHIP4.pl"
			maskout="bpm"//s//"_final_4.pl"
			imexpr(express,maskout,mask1,mask2)
			hedit(maskout,"BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			imdel("../"//maskout)
			imcopy(maskout,"../"//maskout)
			cd("../")
			hedit(s//"_4.fits","BPM",bpmdir//"/"//maskout,add-,delete-,ver-,show+,update+)
			fixpix(s//"_4.fits",masks="BPM",linterp="INDEF",cinterp="INDEF",verbose+,pixels-)
			cd("crrays")
			imdel(s//"_4.fits")
			cd("../")
			imcopy(s//"_4.fits","crrays/")
			imgets(s//"_4","EXPTIME")
			thetime = imgets.value
			mscjoin(s,del+,verb+)
			mscsplit("wght.fits",del-,verb+)
			imcopy("wght_0",s//"_wght_0")
#change the weight exposure time to match the actual image exposure time
#print("I made it this far")
print("editing weight exposure times")
			hedit(s//"_wght_0",fields="EXPTIME",value=thetime,add+,ver-,up+)
#print("but not this far")
#now make weight file for each chip and combine them!
				bpmfile = "bpm"//s//"_final_1.pl"
				imcalc("wght_1,"//bpmfile,s//"_wght_1","if im2.gt.0 then 0.0 else im1")
				hedit(s//"_wght_1",fields="EXPTIME",value=thetime,add+,ver-,up+)
				bpmfile = "bpm"//s//"_final_2.pl"
				imcalc("wght_2,"//bpmfile,s//"_wght_2","if im2.gt.0 then 0.0 else im1")
				hedit(s//"_wght_2",fields="EXPTIME",value=thetime,add+,ver-,up+)
				bpmfile = "bpm"//s//"_final_3.pl"
				imcalc("wght_3,"//bpmfile,s//"_wght_3","if im2.gt.0 then 0.0 else im1")
				hedit(s//"_wght_3",fields="EXPTIME",value=thetime,add+,ver-,up+)
				bpmfile = "bpm"//s//"_final_4.pl"
				imcalc("wght_4,"//bpmfile,s//"_wght_4","if im2.gt.0 then 0.0 else im1")
				hedit(s//"_wght_4",fields="EXPTIME",value=thetime,add+,ver-,up+)
			mscjoin(s//"_wght",output=s//".weight",del+,verb+)
			imdel("wght_*")

# print("yo, are they asking for dark images yet?")
		}

}


# BJW: I have not had much luck with msccmatch, it does find shifts
# okay but the output WCS appears not good.  Better to just skip this?
# use scamp and swarp, those work well on LBC data.

	if (dowcs){

		nlist = flatdolist
		print("Fine tuning the WCS")
			msccmatch.usebpm=yes
			msccmatch.verbose=yes
			msccmatch.nsearch=1000
			msccmatch.search=150.
			msccmatch.rsearch=0.2
			msccmatch.cbox=11
			msccmatch.maxshif=3.
			msccmatch.csig=0.07
			msccmatch.cfrac=0.5
			msccmatch.listcoo=yes
			msccmatch.nfit=50
			msccmatch.rms=0.7
			msccmatch.fitgeom="general"
  			msccmatch.reject=2.4
			msccmatch.update=yes
			msccmatch.interac=yes
			msccmatch.fit=yes
			msccmatch.accept=yes
			list=nlist
			while(fscan(list,img)!=EOF){
				print("Running msccmatch on:")
				print(img)
				mscsplit(img, output="",del-, verb+)
				i=strlen(img)
				i2=i-5
				s=substr(img,1,i2)
				i=strlen(s)
				i2 = i-2
				s2 = substr(s,1,i2)
mscsetwcs(s//"_1", databasefile)


				msccmatch(s//"_1", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=75., rsearch=0.5, cbox=11, maxshif=5., csig=0.1, cfrac=0.3, listcoo=yes, nfit=30, rms=1.5, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_1", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=0.5, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_1", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=30., rsearch=0.5, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)



				#msccmatch(s//"_2", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=100., rsearch=0.5, cbox=11, maxshif=3., csig=0.1, cfrac=0.5, listcoo=yes, nfit=30, rms=1.5, fitgeom="general", reject=4.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_2", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=75., rsearch=0.5, cbox=11, maxshif=3., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=1.0, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_2", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=75., rsearch=0.5, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=1.0, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_3", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=1.0, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=1.5, fitgeom="general", reject=4.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_3", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=30., rsearch=0.5, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_4", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=1.0, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=1.5, fitgeom="general", reject=4.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(s//"_4", coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=0.5, cbox=11, maxshif=5., csig=0.07, cfrac=0.5, listcoo=yes, nfit=30, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)
				
				#mscjoin(s,output="",del-,verb+)

				#msccmatch(img, coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=0.5, cbox=11, maxshif=5., csig=0.1, cfrac=0.5, listcoo=yes, nfit=50, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)

				#msccmatch(img, coords=catfile, outcoor="", usebpm=yes, verbose=yes, nsearch=5000, search=50., rsearch=0.5, cbox=11, maxshif=5., csig=0.1, cfrac=0.5, listcoo=yes, nfit=50, rms=0.7, fitgeom="general", reject=3.0, update=yes, interac=yes, fit=yes, accept=yes)


#now make weight file for each chip
				#bpmfile = "bpm"//s//"_final_1.pl"
				#imcalc("wght_1.fits,"//bpmfile,s//"_1_wght","if im2.gt.0 then 0.0 else im1")
				#bpmfile = "bpm"//s//"_final_2.pl"
				#imcalc("wght_2.fits,"//bpmfile,s//"_2_wght","if im2.gt.0 then 0.0 else im1")
				#bpmfile = "bpm"//s//"_final_3.pl"
				#imcalc("wght_3.fits,"//bpmfile,s//"_3_wght","if im2.gt.0 then 0.0 else im1")
				#bpmfile = "bpm"//s//"_final_4.pl"
				#imcalc("wght_4.fits,"//bpmfile,s//"_4_wght","if im2.gt.0 then 0.0 else im1")

				cd("wcsdone/")
				imdel(s//"*")
				cd("../")

				copy(s//"_1.fits","wcsdone/")
				copy(s//"_2.fits","wcsdone/")
				copy(s//"_3.fits","wcsdone/")
				copy(s//"_4.fits","wcsdone/")



			}

	}




end
