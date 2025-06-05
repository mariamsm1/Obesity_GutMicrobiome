#!/bin/bash

plink=/ludc/Tools/Software/Plink/plink.2.0_a_5.8


genotypefile="../update_ids/output_updated_ids_round2/all_chr_MOS_SCAPIS_BMI_545"
phenofile="gmm.txt"
covfile="cov.txt"

output_dir="output_gmm"

#phenotype column names (skipping the first two columns: FID and IID)
colnames=$(head -n 1 ${phenofile} | awk '{for(i=3;i<=NF;i++) print $i}')

#looping through each phenotype column to run the association analysis

for colname in $colnames
do
  csubmit.sh -b $plink/plink2 \
    -a "--pfile ${genotypefile} \
    --pheno ${phenofile} \
    --covar ${covfile} \
    --pheno-name ${colname} \
    --glm \
    --covar-variance-standardize \
    --threads 4 \
    --out ${output_dir}/lin_reg_${colname}" \
    -m 10000 \
    -c 4 \
    -p ../condor_out
  
  # sleep 1
done
