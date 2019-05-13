#!/bin/bash

CD=$PWD
CNTR=1

for f in sub-*/func; do

    cd $f # go into subject dir

    # true filepath
    tf=`ls wsub-???_task-rest_run-1_bold.nii.gz`
    f2=${tf:1} # minus the w

    # true path of filtered_func_data
    td=${CD}/${f}/${tf}
    printf "%03d: %s\n" $CNTR $td

    # don't run more than 5 fsl_glm at one time
    NUMGLM=$(top -n1|grep fsl_glm|wc -l) 
    while [[ $NUMGLM -gt 4 ]];do
	sleep .5;
	NUMGLM=$(top -n1|grep fsl_glm|wc -l) 
	printf "$NUMGLM GLMs running\r";
    done

    # process the next one
    if [[ ! -f ${tf//.nii/_clean.nii} ]]; then

	echo "${tf//.nii/_clean.nii} not found, regressing out motion now"
	WD=${PWD};
cat <<EOF > /tmp/tmp${CNTR}.sh
cd $WD
sed rp_${f2%.nii*}.txt -e 's/^/1 /' > motion_design_run-1.txt    
fsl_glm -i ${tf} -d motion_design_run-1.txt --out=betas --out_res=resid
fslroi betas mean 0 1    
fslmaths mean -add resid.nii.gz ${tf%.nii*}_clean -odt float
rm -f betas.nii.gz resid.nii.gz mean.nii.gz
EOF

    chmod +x /tmp/tmp${CNTR}.sh
    /tmp/tmp${CNTR}.sh &
    
    fi

    cd ../.. # leave subject dir
    let CNTR=$CNTR+1

done
