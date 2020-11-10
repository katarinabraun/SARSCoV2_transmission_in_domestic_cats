# Data derived
To accompany the manuscript entitled "Transmission of SARS-CoV-2 in domestic cats imposes a narrow bottleneck"

Contains VCF files for each cat, timepoint, and replicate. 

The `.vcf` files are output directly from the SARSquencer pipeline and the `.vcf.recode.vcf` files are these same files but with variants in primer-binding-sites dropped. 

- To do this, we filtered the VCFs created by the SARSquencer pipeline using the `filterSNV.sh` script, which identifies primer-binding regions using the `2019-nCoV-tokyo.bed` file. 

The `data_cleaning` jupyter notebook takes these files in to clean them and outputs csvs to the `cleaned` directory as well as intersection csvs for the replicates. 

### Contact info 

Please feel free to reach out to Katarina Braun (kmbraun2@wisc.edu or kmbraun2@gmail.com) or Gage Moreno (gkmoreno@wisc.edu) with questions. 
