#!/bin/bash

CD=${PWD}

if [[ $(locate -b '\spm.m') == "" ]]; then
    echo "please make sure that you have SPM"    
    return 1
else
    SPMPATH=$(echo $(locate -b '\spm.m')|sed -e 's/\ /:/g'|sed -e 's|/spm.m||g')
    export MATLABPATH=${SPMPATH}:${MATLABPATH}
fi

if [[ $(which fsl) == "" ]]; then
    echo "please make sure that you have FSL"
    return 1    
fi

# make a list of run-1 fMRI experiments
ls -1 sub-*/func/swsub*1_bold.nii* > swfuncs.txt
ls -1 sub-*/func/wsub*1_bold.nii*  >  wfuncs.txt

# regress out motion effects from the run-1 by using the realignment parameters rp
regress_mot.sh
ls -1 sub-*/func/wsub*1_bold_clean.nii.gz > cleanfuncs.txt

# exclude subjects 10, 32, 76, 104, 112 with mean FD > 0.2 mm -- see mriqc/derivatives_mriqc_bold.csv
grep -v sub-010 cleanfuncs.txt | grep -v sub-032 | grep -v sub-076 | grep -v sub-104 | grep -v sub-112 > fd_cleanfuncs.txt

# make masks of fMRI too based on non-smoothed data, at 20% of the max intensity
# before cleaning / exclusion
# while read fun; do printf "." && fslmaths $fun -Tmin -thrP 20 -bin ${fun//.nii/_mask.nii}; done < wfuncs.txt && printf "\n"
# after cleaning / exclusion
while read fun; do printf "." && fslmaths ${fun} -Tmin -thrP 20 -bin ${fun//.nii/_mask.nii}; done < fd_cleanfuncs.txt && printf "\n"

# take the mean over all masks and threshold at 75% (i.e. find where 3 out of every 4 have data)
fslmerge -t fngroup $(cat wfuncs.txt|sed 's/.nii/_mask.nii/')
fslmaths fngroup -Tmean -thr 0.75 -bin fngroup

# unzip the files necessary to map the GM masks to standard space
gzip -d $(find -name c1\*nii.gz) $(find -name y_\*nii.gz)

# Map all the c1 images (GM density) to MNI space -
# if this does not run automaticall, you can just typ it))
# >> apply_deform.m
# >> exit
matlab -nodesktop -nosplash < apply_deform.m 

# Combine c1 files into one 4d group GM density file
fslmerge -t gmgroup */*/wc1*.nii*

# Threshold each volume at .2 and binarise, take group mean and threshold at .2 again
fslmaths gmgroup -thr .2 -bin -Tmean -thr .2 -bin gmgroup 

## make a matrix to apply to the structural data to be in "fMRI - MNI" space 
echo "1 0 0 -12"  > gm2fun.mat
echo "0 1 0 -14" >> gm2fun.mat
echo "0 0 1  -2" >> gm2fun.mat
echo "0 0 0   1" >> gm2fun.mat

# Resample at fMRI resolution
flirt -in gmgroup -ref fngroup -out gmgroup_3mm -applyisoxfm 3 -init gm2fun.mat

# Threshold partial volumes again, at 20% GM density
fslmaths gmgroup_3mm -thr .2 -bin gmgroup_3mm

# use the MNI 'no-cerebellum' mask (mean of MNI-prob-1mm w/o 2nd volume)
fslroi /usr/share/data/mni-structural-atlas/MNI/MNI-prob-2mm.nii.gz nocereb_1 0 1
fslroi /usr/share/data/mni-structural-atlas/MNI/MNI-prob-2mm.nii.gz nocereb_2 2 -1
fslmerge -t nocereb nocereb_1 nocereb_2
fslmaths nocereb -Tmean nocereb2
flirt -in MNI-prob-1mm_nocereb -ref fngroup -out MNI-prob-1mm_nocereb_3mm -applyisoxfm 3 -init gm2fun.mat

# the mask that we will use is the intersection of gmgroup_3mm, fngroup, and the no-cerebellum mask
fslmaths gmgroup_3mm -mas fngroup -mas MNI-prob-1mm_nocereb_3mm -bin maskgroup -odt char

################################################################################
################################################################################

# Get the BIAS repository
git clone https://github.com/amwink/bias.git
# If you don't have git you can download the link below in your browser or use wget:
# wget -c --no-check-certificate https://github.com/amwink/bias/archive/master.zip
# and unzip the file. Below we assume the bias root directory is calles 'bias'
# - again, if this does not work then typing 'fegui' should work:
# >> addpath(genpath('bias'))
# >> fegui # load the fMRI's and the groupmask, run with default options
#    -- click "use mask" and select maskgroup.nii.gz
#    -- click "add files" and select wfuncs.txt (loading and checking these takes a while!)
#    -- go to the "Settings" tab
#    -- look at the extra options, no need to change anything for now
#    -- click "estimate fECM"
# >> exit
#
# you can also run
# >> fegui_replacement.m
# in MatLab, if you are happy with the recipe in that script
# 
matlab < bias/matlab/fastECM/fegui.m

# remove NaN from ECMs
for f in $(find . -name \*fastECM\*nii\*); do
    cm="fslmaths $f -nan $f ";
    echo $cm;
    #$cm;
done > /tmp/removenan.sh
       
# Build cglm design matrix from the tsv file and the paths to fastECM images
while read wfun; do
    sub=${wfun%%/*};                                           # extract subject number
    printf "$PWD/${wfun//.nii.gz/_fastECM_fast_dyn100.nii} ";  # print the file, changing .nii.gz into _fastECM.nii
    grep $sub participants.tsv |                               # find the subject in the participants file
	awk '{print $3,$2,$6}' |                               # print gender, age, scan duration (not handedness/language)
	sed -e 's/M/-1/' |                                     # male   -> gender -1
	sed -e 's/F/1/';                                       # female -> gender +1
done < fd_cleanfuncs.txt > design.txt                                 # read from wfuncs.txt, write to design.txt  

################################################################################
################################################################################

# do spatial stats on the 4D files
DYNS=$(find $CD -name \*clean\*dyn100.nii\*|sort)
RELU=$(grep relu <<< "$DYNS")
FAST=$(grep -v relu <<< "$DYNS")

rm -f 4dstats.txt
for f in $FAST; do

    printf "%75s " $f
    # output 4D  min max rmin rmax vol_mm vol_vox mean std-dev
    fslstats $f -n -k maskgroup.nii.gz -l 0.0001 -R -r -V -M -S 

done >> 4dstats.txt

paste design.txt 4dstats.txt | sed -e 's|/home/amwink/work/data/openNeuro/WashU_120/||g' | sed -e 's|_task-rest_run-1_bold||g' | sed -e 's/.000000//g' | sed -e 's/\t/ /g' > design_data.txt 

################################################################################
################################################################################

# compute temporal stats as maps
for f in $FAST; do

    printf "%75s: \n" $f
    
    # compute range for dynamics
    cm="fslmaths $f -Tmin /tmp/tmp1"
    echo $cm;$cm;
    cm="fslmaths $f -Tmax -sub /tmp/tmp1 ${f//.nii/_range.nii}"
    echo $cm;$cm;
    
    # compute temporal mean and std for dynamics
    cm="fslmaths $f -Tmean ${f//.nii/_mean.nii}"
    echo $cm;$cm;
    cm="fslmaths $f -Tstd ${f//.nii/_std.nii}"
    echo $cm;$cm;

    # compute temporal median and interquartile range
    cm="fslmaths $f -Tmedian ${f//.nii/_median.nii}"
    echo $cm;$cm;
    cm="fslmaths $f -Tperc 25 /tmp/tmp1"
    echo $cm;$cm;
    cm="fslmaths $f -Tperc 75 -sub /tmp/tmp1 ${f//.nii/_iqr.nii}"
    echo $cm;$cm;

done

################################################################################
################################################################################

# don't do cglm while it crashes
do_cglm=0

if (( $do_cglm )); then 

    echo "cglm is fixed; running now:"
    
    # make fastECM double datatype, remove NaN (convert to 0)
    for f in $(find . -name \*fastECM\*.nii\*); do
	CM="fslmaths $f ${f//.nii.gz/.nii} -odt double"
	echo $CM
	$CM
    done
    
    # make the model: ECM differences between male/female (corrected for age and scan length)
    printf "MAIN:\n1\nINTERACTIONS:\n" > model2x1.txt
    
    # make a background for nicer output
    flirt -in $FSLDIR/data/standard/MNI152_T1_1mm -ref fngroup -out standard_3mm -applyisoxfm 3 -init gm2fun.mat
    gzip -d standard_3mm*gz maskgroup*gz $(find . -name \*fastECM\*dyn001.nii.gz)

    # call cglm
    # mkdir -p cglmdir
    cglm -m model2x1.txt -d design.txt -t standard_3mm.nii -o cglmdir -F -l u2x1 -c 1 -C auto -B maskgroup.nii -p 1
    
    # Re-zip all niftis after using SPM ( -P 8 if you have 8 CPUs )
    find -name \*.nii -type f -print0 | xargs -0 -n 1 -P 4 gzip -f9   

    echo "cglm done"
    
fi

################################################################################
################################################################################

# make group mean time average ECMs?
do_makemeans=0;

if (( $do_makemeans != 0 )); then

    # mean ECM per group
    YMALFscans=$(grep 'nii -1 -1' design_dr.txt | awk '{print $1}')
    for f in $YMALFscans; do fslmaths $f -Tmean ${f//.nii/_ymfmean.nii};printf ".";done;echo
    fslmerge -t mean_ymalef_ecm sub-???/func/*_ymfmean.nii*
    fslmaths mean_ymalef_ecm -Tmean mean_ymalef_ecm
    YMALRscans=$(echo "${YMALFscans//fast_/relu_}")
    for f in $YMALRscans; do fslmaths $f -Tmean ${f//.nii/_ymrmean.nii};printf ".";done;echo
    fslmerge -t mean_ymaler_ecm sub-???/func/*_ymrmean.nii*
    fslmaths mean_ymaler_ecm -Tmean mean_ymaler_ecm
    YFEMFscans=$(grep 'nii 1 -1' design_dr.txt | awk '{print $1}')
    for f in $YFEMFscans; do fslmaths $f -Tmean ${f//.nii/_yffmean.nii};printf ".";done;echo
    fslmerge -t mean_yfemalef_ecm sub-???/func/*_yffmean.nii*
    fslmaths mean_yfemalef_ecm -Tmean mean_yfemalef_ecm
    YFEMRscans=$(echo "${YFEMFscans//fast_/relu_}")
    for f in $YFEMRscans; do fslmaths $f -Tmean ${f//.nii/_yfrmean.nii};printf ".";done;echo
    fslmerge -t mean_yfemaler_ecm sub-???/func/*_yfrmean.nii*
    fslmaths mean_yfemaler_ecm -Tmean mean_yfemaler_ecm
    OMALFscans=$(grep 'nii -1 1' design_dr.txt | awk '{print $1}')
    for f in $OMALFscans; do fslmaths $f -Tmean ${f//.nii/_omfmean.nii};printf ".";done;echo
    fslmerge -t mean_omalef_ecm sub-???/func/*_omfmean.nii*
    fslmaths mean_omalef_ecm -Tmean mean_omalef_ecm
    OMALRscans=$(echo "${OMALFscans//fast_/relu_}")
    for f in $OMALRscans; do fslmaths $f -Tmean ${f//.nii/_omrmean.nii};printf ".";done;echo
    fslmerge -t mean_omaler_ecm sub-???/func/*_omrmean.nii* 
    fslmaths mean_omaler_ecm -Tmean mean_omaler_ecm
    OFEMFscans=$(grep 'nii 1 1' design_dr.txt | awk '{print $1}')
    for f in $OFEMFscans; do fslmaths $f -Tmean ${f//.nii/_offmean.nii};printf ".";done;echo
    fslmerge -t mean_ofemalef_ecm sub-???/func/*_offmean.nii*
    fslmaths mean_ofemalef_ecm -Tmean mean_ofemalef_ecm
    OFEMRscans=$(echo "${OFEMFscans//fast_/relu_}")
    for f in $OFEMRscans; do fslmaths $f -Tmean ${f//.nii/_ofrmean.nii};printf ".";done;echo
    fslmerge -t mean_ofemaler_ecm sub-???/func/*_ofrmean.nii*
    fslmaths mean_ofemaler_ecm -Tmean mean_ofemaler_ecm
    
    # show "fast"
    CM="fsleyes -std1mm";for f in mean*malef*gz;do CM="$CM $f -cm render3 -dr .004325 .004568 -a 75"; done; $CM
    # show "relu"
    CM="fsleyes -std1mm";for f in mean*maler*gz;do CM="$CM $f -cm render3 -dr .004325 .004568 -a 75"; done; $CM
    
    FMEANscans=$(cat design_dr.txt | awk '{print $1}')
    RMEANscans=$(echo "${FMEANscans//fast_/relu_}")
    fslmerge -t fmeanscan $FMEANscans
    fslmaths fmeanscan -Tstd  fstdscan
    fslmaths fmeanscan -Tmean fmeanscan
    fslmerge -t rmeanscan $RMEANscans
    fslmaths rmeanscan -Tstd  rstdscan
    fslmaths rmeanscan -Tmean rmeanscan
    
    # to make a 2x2 stack of maps, do something like this:
    # 
    # make a vertical stack of 2:
    # convert -append mean_ymale.png mean_omale.png col1.png
    # make another vertical stack:
    # convert -append mean_yfemale.png mean_ofemale.png col2.png
    # append them horizontally:
    # convert +append col1.png col2.png dual_regression.png 

fi

################################################################################
################################################################################

# do melodic on dynamic ECM and input fMRI?
do_melodic=0

if (( $do_melodic )); then

    # prepare randomise and melodic
    mkdir -p randomisedir
    while read f; do ls $PWD/${f};done < wfuncs.txt > randomisedir/gicafmri.txt
    rm -rf randomisedir/gicafmri && mkdir -p randomisedir/gicafmri
    while read f; do ls $PWD/${f//bold.nii.gz/bold_fastECM_relu_dyn100.nii.gz};done < wfuncs.txt > randomisedir/gicafast.txt
    rm -rf randomisedir/gicafast && mkdir -p randomisedir/gicafast
    while read f; do ls $PWD/${f//bold.nii.gz/bold_fastECM_relu_dyn100.nii.gz};done < wfuncs.txt > randomisedir/gicarelu.txt
    rm -rf randomisedir/gicarelu && mkdir -p randomisedir/gicarelu
    
    # run melodic
    for (( seedval=0; seedval<10; seedval++ )); do
	seedstr=$(printf "%04d" $seedval)
	
	#mkdir -p randomisedir/gicarelu/${seedstr}
	#cm="melodic -i randomisedir/gicarelu.txt -o randomisedir/gicarelu/${seedstr} -m maskgroup -d 25 -n 25 --nobet --report --Oall --bgimage=standard_3mm --seed=${seedval} -v"
	#echo $cm
	#time $cm &> randomisedir/gicarelu.log &
	
	mkdir -p randomisedir/gicafast/${seedstr}
	cm="melodic -i randomisedir/gicafast.txt -o randomisedir/gicafast/${seedstr} -m maskgroup -d 25 -n 25 --nobet --report --Oall --bgimage=standard_3mm --seed=${seedval} -v"
	echo $cm
	#time $cm &> randomisedir/gicarelu.log &
	
	mkdir -p randomisedir/gicafmri/${seedstr}
	cm="melodic -i randomisedir/gicafmri.txt -o randomisedir/gicafmri/${seedstr} -m maskgroup -d 25 -n 25 --nobet --report --Oall --bgimage=standard_3mm --seed=${seedval} -v"
	echo $cm
	#time $cm &> randomisedir/gicarelu.log &    
    done > melodic_commands.sh
    time parallel -j 4 < melodic_commands.sh &

fi

################################################################################
################################################################################

# do randomise on single ECMs?
do_randomise=0

if (( $do_randomise != 0 )); then
    
    # call randomise
    mkdir -p randomisedir/wfunc
    fslmerge -t randomisedir/ecm4d_wfunc  $(cat  wfuncs.txt | sed -e 's/.nii.gz/_fastECM.nii/')
    randomise -i randomisedir/ecm4d_wfunc -o randomisedir/wfunc/wfunc -D -m maskgroup.nii -d glm_ecm.mat -t glm_ecm.con -f glm_ecm.fts -e glm_ecm.grp -T -x -R --seed=4 -n 10000
    mkdir -p randomisedir/swfunc
    fslmerge -t randomisedir/ecm4d_swfunc $(cat swfuncs.txt | sed -e 's/.nii.gz/_fastECM.nii/')
    randomise -i randomisedir/ecm4d_swfunc -o randomisedir/swfunc/swfunc -D -m maskgroup.nii -d glm_ecm.mat -t glm_ecm.con -f glm_ecm.fts -e glm_ecm.grp -T -x -R --seed=4 -n 10000

done

################################################################################
################################################################################

# do PALM on single ECMs?
do_palm=0

if (( $do_palm != 0 )); then

    # add PALM to the path
    export PATH=/usr/local/palm:${PATH}
    
    # 
    fslmerge -t randomisedir/ecm4d_wfunc_fast $(cat  wfuncs.txt | sed -e 's/.nii.gz/_fastECM_fast_dyn001/')
    #gzip -d randomisedir/ecm4d_wfunc_fast.nii.gz 
    #gzip -d maskgroup.nii.gz
    palm -i randomisedir/ecm4d_wfunc_fast.nii.gz -o randomisedir/wfunc_fast/wfunc_fast -m maskgroup.nii.gz \
	 -d ECM_sex_binage.mat -t ECM_sex_binage.con -f ECM_sex_binage.fts -eb ECM_sex_binage.grp \
	 -T -demean -noniiclass -corrcon -logp -fdr -accel tail -n 500
    #-d ECM_sexonly.mat -t ECM_sexonly.con -f ECM_sexonly.fts -eb ECM_sexonly.grp \
    #gzip -9 randomisedir/ecm4d_wfunc_fast.nii
    #gzip -9 maskgroup.nii
    
    fslmerge -t randomisedir/ecm4d_wfunc_relu $(cat wfuncs.txt | sed -e 's/.nii.gz/_fastECM_relu_dyn001/')
    #gzip -d randomisedir/ecm4d_wfunc_relu.nii.gz 
    #gzip -d maskgroup.nii.gz
    palm -i randomisedir/ecm4d_wfunc_relu.nii.gz -o randomisedir/wfunc_relu/wfunc_relu -m maskgroup.nii.gz \
	 -d ECM_sex_binage.mat -t ECM_sex_binage.con -f ECM_sex_binage.fts -eb ECM_sex_binage.grp \
	 -T -demean -noniiclass -corrcon -logp -fdr -accel tail -n 500
    #-d ECM_sexonly.mat -t ECM_sexonly.con -f ECM_sexonly.fts -eb ECM_sexonly.grp \
    #gzip -9 randomisedir/ecm4d_wfunc_relu.nii
    #gzip -9 maskgroup.nii
 
    # View FAST results
    CM="fsleyes -std1mm"
    CM="${CM} randomisedir/wfunc_fast/wfunc_fast_tfce_tstat_fwep_c3.nii.gz -cm yellow -dr 1.301 2 -a 75" # O > Y
    CM="${CM} randomisedir/wfunc_fast/wfunc_fast_tfce_tstat_fwep_c4.nii.gz -cm green  -dr 1.301 2 -a 75" # Y > O
    CM="${CM} randomisedir/wfunc_fast/wfunc_fast_tfce_tstat_fwep_c1.nii.gz -cm red    -dr 1.301 2 -a 75" # F > M
    CM="${CM} randomisedir/wfunc_fast/wfunc_fast_tfce_tstat_fwep_c2.nii.gz -cm blue   -dr 1.301 2 -a 75" # M > F
    echo $CM
    $CM 2> /dev/null
    
    # View RELU results
    CM="fsleyes -std1mm"
    CM="${CM} randomisedir/wfunc_relu/wfunc_relu_tfce_tstat_fwep_c3.nii.gz -cm yellow -dr 1.301 2 -a 75" # O > Y
    CM="${CM} randomisedir/wfunc_relu/wfunc_relu_tfce_tstat_fwep_c4.nii.gz -cm green  -dr 1.301 2 -a 75" # Y > O
    CM="${CM} randomisedir/wfunc_relu/wfunc_relu_tfce_tstat_fwep_c1.nii.gz -cm red    -dr 1.301 2 -a 75" # F > M
    CM="${CM} randomisedir/wfunc_relu/wfunc_relu_tfce_tstat_fwep_c2.nii.gz -cm blue   -dr 1.301 2 -a 75" # M > F
    echo $CM
    $CM 2> /dev/null

done

################################################################################
################################################################################

# do dual regression on dynamic ECMs?
do_dr=1

if (( $do_dr != 0 )); then
 
    # make a design with the binary contrast used for dual regression
    a=$(cat design.txt|awk '{print $1,$2}')
    echo "$a" > /tmp/a.txt
    b=$(grep evg ECM_sex_binage115.fsf | grep '.2)' | awk '{print $3}')
    echo "$b" > /tmp/b.txt
    paste /tmp/a.txt /tmp/b.txt | tr "\t" " " > design_dr.txt
    
    # use FSL's Glm GUI to make designs
    # Glm &
    
    # dual regression of dynamic ECM using Smith's masks -- age dichotomised
    ls sub-???/func/*clean*fast_dyn100.nii* > alldyn.txt 
    dual_regression_maja_mask_twosided PNAS_Smith09_bm10_3mm 1 ECM_sex_binage115.mat ECM_sex_binage115.con 10000 dual_regression_sex_binage115 $(cat alldyn.txt)
    for (( rsn=0; rsn<10; rsn++ )); do
	rsnstr=$(printf "%04d" $rsn);
	echo "component $rsn"
	for sign in n p; do
	    printf "sign: $sign -- "
	    fslmerge -t /tmp/tmpstats dual_regression_sex_binage115_twosided/dr_stage3_ic${rsnstr}${sign}_tfce_corrp_tstat?.nii.gz
	    fslmaths    /tmp/tmpstats -Tmax  /tmp/tmpstats
	    statvol=$(fslstats /tmp/tmpstats -l .95 -V)
	    statvol=${statvol%%\ *}
	    if [[ ${statvol} -gt 1 ]]; then
		echo "$statvol voxels found" 
		CM="fsleyes -std1mm dual_regression_sex_binage115_twosided/${sign}ICMASK${rsn}.nii.gz -cm pink -dr .99 99 -a 50"
		CM="${CM} dual_regression_sex_binage115_twosided/dr_stage3_ic${rsnstr}${sign}_tfce_corrp_tstat3.nii.gz -cm yellow -dr .95 1 -a 75" # O > Y
		CM="${CM} dual_regression_sex_binage115_twosided/dr_stage3_ic${rsnstr}${sign}_tfce_corrp_tstat4.nii.gz -cm green  -dr .95 1 -a 75" # Y > O
		CM="${CM} dual_regression_sex_binage115_twosided/dr_stage3_ic${rsnstr}${sign}_tfce_corrp_tstat1.nii.gz -cm red    -dr .95 1 -a 75" # F > M
		CM="${CM} dual_regression_sex_binage115_twosided/dr_stage3_ic${rsnstr}${sign}_tfce_corrp_tstat2.nii.gz -cm blue   -dr .95 1 -a 75" # M > F
		echo $CM
		$CM 2> /dev/null
	    else
		echo
	    fi
	done
    done

    # to make a 2x3 stack of maps, do something like this:
    # 
    # make a vertical stack of 3:
    # convert -append dyn100_dr115_bin_ic0_pos_medvis_vox231330.png dyn100_dr115_bin_ic2_pos_latvis_vox052228.png dyn100_dr115_bin_ic3_pos_dmn_vox242237.png col1.png
    # make another vertical stack:
    # convert -append dyn100_dr115_bin_ic5_pos_mot_vox183244.png dyn100_dr115_bin_ic7_pos_execcontr_vox184241.png dyn100_dr115_bin_ic9_pos_execcontr_vox391741.png col2.png
    # append them horizontally:
    # convert +append col1.png col2.png dual_regression.png 
    
fi
