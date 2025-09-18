#!/bin/bash

plink=/ludc/Tools/Software/Plink/plink.2.0_a_5.8
genotypefile="../update_ids/output_updated_ids_round2/all_chr_MOS_SCAPIS_BMI_545"
phenofile="pheno_sh_rich_invsimp.txt"
covfile="cov.txt"

#shannon_nds and richness are calculated from non dowsized data. shannon_ds and richness_ds are calculated from downsized data
colname="inverse_simpson" #change this to "shannon_nds" or "shannon_ds" if the association is with shannon and to "richness" or "richness_ds" if the association is with richness
#change the name of the output file according to the above

csubmit.sh -b $plink/plink2 -a "--pfile ${genotypefile} --pheno ${phenofile} --covar ${covfile} --pheno-name ${colname} --glm --covar-variance-standardize --threads 4 --out output_sh_rich_invsimp/lin_reg_invsimpson" -m 10000 -c 4 -p ../condor_out

