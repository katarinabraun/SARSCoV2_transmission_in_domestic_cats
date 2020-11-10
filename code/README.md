# Code
To accompany the manuscript entitled "Transmission of SARS-CoV-2 in domestic cats imposes a narrow bottleneck"

This directory contains all the scripts that were used to process the raw data. The pipeline to process the raw fastqs (pair, merge, map, call variants) can be found in the `SARSquencer` directory. All Jupyter Notebook files (file extension = `.ipynb`) were written to process the derived data (primarily the VCFs as well as the nucleotide diversity estimates) and to generate the figures that can be found in the manuscript. All other files are scripts that the Jupyter Notebooks pulls in (i.e. snpgenie.pl is used to generate nucleotide diversity estimates). 

Here's a brief overview of each of the Jupyter Notebooks, in the order that they are easiest to run:  
    - `data_cleaning.ipynb` cleans the raw-recoded VCF files and generates simplified CSV files that downstream scripts use. This notebook also compares replicate files and generates *intersection* CSV files, which only contain the variants that are found in independent replicates.   
    - `variants_over_time.ipynb` pulls in the intersection CSVs created above and plots all of these per annotation type, cat, and transmission pair. This notebook generates **figure 3** as well as **supplemental figures 6 and 7**.  
    - `plot_all_variants_by_genome_position_and_frequency.ipynb` also pulls in the intersection CSVs to create **figure 2A**. This notebook also looks at the amount of shared variation across cats. 
    - `SFS.ipynb` pulls in the intersection CSVs to generate **figure 2B** and includes the stats to compare cat SNP frequency spectrums (SFS) to the neutral expectation as well as to compare index cats' SFS to contact cats' SFS. 
    - `bottleneck_estimates.ipynb` generates input files for the [Sobel-Leonard/Koelle](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5487570/) exact beta-binomial method. It writes these files to `data_derived/beta_binomial/*-input.txt` and then uses the script (`Bottleneck_size_estimation_exact.r`) to generate Nb estimates as well as to plot them as TV/JT(McCrone) plots -- **figure 4**. 
    - `diversity_estimates.ipynb` uses Chase Nelson's `snpgenie.pl` script to generate nucleotide estimates for each cat, timepoint, and replicate. This script also statistically compares nucleotide diversity estimates in index and contact cats and plots π and πN/πS estimates -- **supplemental figure 8**. 

### Contact info 

Please feel free to reach out to Katarina Braun (kmbraun2@wisc.edu or kmbraun2@gmail.com) or Gage Moreno (gkmoreno@wisc.edu) with questions. 
