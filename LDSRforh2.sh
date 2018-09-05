#!/bin/bash

#derived from saige_chunking_script.sh and /net/dumbo/home/zhowei/projects/UKBIOBANK/SAIGE/testInCloud11012017/step2/step2_SPATests.sh

filePrefix=$1 # *.rda
chrom=$2 #1-22
outputPrefix=$3

if [ ! -d "${outputPrefix}.output" ]; then
    mkdir ${outputPrefix}.output
fi

pids=() #array of process ids

#chunk 00-31
for chunk in $(seq -f %02g 00 31)
do
    sleep 60s #pause because of I/O errors when all jobs start at once
    Rscript /net/snowwhite/home/bwolford/scripts/SAIGE_step2.R \
	--bgenFileIndex=/net/fantasia/home/goncalo/UKBioBank/download/geno2018/EGAD00010001474/ukb_imp_chr${chrom}_v3.bgen.bgi \
	--bgenFile=/net/fantasia/home/goncalo/UKBioBank/download/geno2018/EGAD00010001474/ukb_imp_chr${chrom}_v3.bgen \
        --sampleFile=/net/hunt/zhowei/project/imbalancedCaseCtrlMixedModel/Rpackage_SPAGMMAT/gz_quantitaive_08022017/UKbgen/ukbgen.sample \
	--idstoIncludeFile=/net/hunt/disk2/bwolford/UKBB_CAD_2/UK10K_imputation/splitUK10KVariants_nofilter/chr_${chrom}_${chunk}  \
        --GMMATmodelFile=${filePrefix}.rda \
        --varianceRatioFile=${filePrefix}.varianceRatio.txt \
        --SAIGEOutputFile=./${outputPrefix}.output/${outputPrefix}.${chrom}.${chunk}.SAIGE.txt \
	--filetype BGEN &
  pids+=($!) #process id
done

rc=0
for i in $(seq 0 31); do
  p=${pids[i]}
  wait $p 
  wait_rc=$? #return code of subprocess
  rc=$(( rc + wait_rc ))

  [[ $wait_rc != 0 ]] && echo "[$(date)] Trait ${trait} chrom ${chrom} chunk ${i} FAILED!  (PID ${p}, exit status: ${wait_rc})"
done

[[ $rc != 0 ]] && exit 1 || exit 0
