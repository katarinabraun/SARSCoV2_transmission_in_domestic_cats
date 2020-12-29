# SARSCoV2_transmission_in_domestic_cats
To accompany the manuscript entitled *Transmission of SARS-CoV-2 in domestic cats imposes a narrow bottleneck*. 

## Overview of the repository 

There are five primary directories in this repository: 
1. `data_raw` contains a link to the raw fastqs on SRA. The raw fastqs are not included within this repository because the files are too large. 
2. `data_derived` contains all output VCF files from the SARSquencer pipeline (more below) as well as the cleaned version of these VCFs, nucleotide diversity estimates in the `diversity` directory, and input files for the beta-binomial method in the `beta_binomial` directory. 
    - The VCF files included here have also been modified to discard any variants that fall in ARTIC v3 primer-binding regions and therefore have the file extension `*.vcf.recode.vcf`
    - To do this, we filtered the VCFs created by the SARSquencer pipeline using the `filterSNV.sh` script, which identifies primer-binding regions using the `2019-nCoV-tokyo.bed` file. 
3. `code` contains all the scripts that were used to process the raw data. The pipeline to process the raw fastqs (pair, merge, map, call variants) can be found in the `SARSquencer` directory. All Jupyter Notebook files (file extension = `.ipynb`) were written to process the derived data (primarily the VCFs as well as the nucleotide diversity estimates) and to generate the figures that can be found in the manuscript. All other files are scripts that the Jupyter Notebooks pull in (i.e. snpgenie.pl is used to generate nucleotide diversity estimates). 

    Here's a brief overview of each of the Jupyter Notebooks, in the order that they are easiest to run:  
    - `data_cleaning.ipynb` cleans the raw VCF files and generates simplified CSV files that downstream scripts use. This notebook also compares replicate files and generates *intersection* CSV files, which only contain the variants that are found in independent replicates.   
    - `variants_over_time.ipynb` pulls in the intersection CSVs created above and plots all of these per annotation type, cat, and transmission pair. This notebook generates **figure 3** as well as **supplemental figures 6 and 7**.  
    - `plot_all_variants_by_genome_position_and_frequency.ipynb` pulls in the intersection CSVs to create **figure 2A**. This notebook also looks at the amount of shared variation across cats. 
    - `SFS.ipynb` pulls in the intersection CSVs to generate **figure 2B** and includes the stats to compare cat SNP frequency spectrums (SFS) to the neutral expectation as well as to compare index cats' SFS to contact cats' SFS. 
    - `bottleneck_estimates.ipynb` generates input files for the [Sobel-Leonard/Koelle](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5487570/) exact beta-binomial method. It writes these files to `data_derived/beta_binomial/*-input.txt` and then uses the script (`Bottleneck_size_estimation_exact.r`) to generate Nb estimates as well as to plot them as TV/[JT McCrone](https://elifesciences.org/articles/35962) plots -- **figure 4**. 
    - `diversity_estimates.ipynb` uses Chase Nelson's `snpgenie.pl` script to generate nucleotide diversity estimates for each cat, timepoint, and replicate. This script also statistically compares nucleotide diversity estimates in index and contact cats and plots π and πN/πS estimates -- **supplemental figure 8**. 

4. The `figures` directory contains all of the vector files that the above jupyter notebooks generate. In addition, the final Adobe Illustrator (formatted versions of the SVG files) can be found in the `illustrator_files` sub-directory. 
- All figures apart from **figure 1** and **supplemental figures 2 and 3** can be re-created using the above Jupyter Notebooks. 
- **Figure 1** was created by hand in Adobe Illustrator (the .ai file is included in the figures directory)
- **Supplemental figures 2 and 3** were created using simple samtools command line prompts (below) and were then edited (e.g. colors, font sizes) for clarity in Adobe Illustrator. 

    ```bash 
    samtools sort SAMPLE.aln.sam > SAMPLE.sorted.bam
    samtools index SAMPLE.sorted.bam
    samtools depth -d 0 SAMPLE.sorted.bam > SAMPLE.sorted.coverage.tsv
    echo -e "Reference\tPosition\tDepth" | cat - SAMPLE.sorted.coverage.tsv
    SAMPLE.headers.sorted.coverage.tsv
    ```

5. `MW219695.1` contains the reference GTF and fasta. 

### Abstract 

The evolutionary mechanisms by which SARS-CoV-2 viruses adapt to mammalian hosts and, potentially, undergo antigenic evolution depend on the ways genetic variation is generated and selected within and between individual hosts. Using domestic cats as a model, we show that SARS-CoV-2 consensus sequences remain largely unchanged over time within hosts, while dynamic sub-consensus diversity reveals processes of genetic drift and weak purifying selection. We further identify a notable variant at amino acid position 655 in Spike (H655Y), which was previously shown to confer escape from human monoclonal antibodies. This variant arises rapidly and persists at intermediate frequencies in index cats. It also becomes fixed following transmission in two of three pairs. These dynamics suggest this site may be under positive selection in this system and illustrate how a variant can quickly arise and become fixed in parallel across multiple transmission pairs. Transmission of SARS-CoV-2 in cats involved a narrow bottleneck, with new infections founded by fewer than ten viruses. In RNA virus evolution, stochastic processes like narrow transmission bottlenecks and genetic drift typically act to constrain the overall pace of adaptive evolution. Our data suggest that here, positive selection in index cats followed by a narrow transmission bottleneck may have instead accelerated the fixation of S H655Y, a potentially beneficial SARS-CoV-2 variant. Overall, our study suggests species- and context-specific adaptations are likely to continue to emerge. This underscores the importance of continued genomic surveillance for new SARS-CoV-2 variants as well as heightened scrutiny for signatures of SARS-CoV-2 positive selection in humans and mammalian model systems.

### Contact info 

Please feel free to reach out to Katarina Braun (kmbraun2@wisc.edu or kmbraun2@gmail.com) or Gage Moreno (gkmoreno@wisc.edu) with questions. 
