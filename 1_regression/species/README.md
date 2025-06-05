# Mariam Miari
# 250604

## Description of the scripts:
1. species_bmi.Rmd : association between CLR-transformed species relative abundance as exposure and BMI as outcome. This is run for the analysis of all individuals, sex-stratified, and sex-interaction.
2. src_species_bmi.R : source code for (1) above.
3. bmi_species.Rmd : association between BMI as exposure and CLR-transformed species relative abundance as outcome. This is run for the analysis of all individuals, sex-stratified, and for sex-interaction.
4. bmi_species_PA.Rmd : association between BMI as exposure and species presence/absence as outcome. This is run for the analysis of all individuals, sex-stratified, and for sex-interaction. 
5. bmi_species_pres.Rmd : association between BMI as exposure and species presence (discarding all zero abundances) as outcome. This is run for the analysis of all individuals, sex-stratified, and for sex-interaction.
6. sens1_ppi_adj.Rmd : association between BMI as exposure and CLR-transformed species relative abundance as outcome, additionally adjusting for PPI use. This is run for the analysis of all individuals, sex-stratified, and for sex-interaction.
7. sens2_antibiotics_adj.Rmd : association between BMI as exposure and CLR-transformed species relative abundance as outcome, additionally adjusting for antibiotics use. This is run for the analysis of all individuals,sex-stratified, and for sex-interaction.    
8. sens3_antibiotics_excl.Rmd : association between BMI as exposure and CLR-transformed species relative abundance as outcome. Individuals taking antibiotics during the previous 6 months were excluded. This is run for the analysis of all individuals, sex-stratified, and for sex-interaction.
9. get_most_sig_families.R : using fisher's test to get the families that are over-represented among all associations, positive associations and negative associations.
