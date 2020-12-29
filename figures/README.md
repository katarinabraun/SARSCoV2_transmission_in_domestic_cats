# Figures
To accompany the manuscript entitled "Transmission of SARS-CoV-2 in domestic cats imposes a narrow bottleneck"

The `figures` directory contains all of the vector files that the jupyter notebooks scripts output.   
In addition, the final Adobe Illustrator (formatted versions of the SVG files) can be found in the `illustrator_files` sub-directory. 


- All figures apart from **figure 1**, **supplemental figure 5**, and **supplemental figures 6 and 7** can be re-created using the above Jupyter Notebooks. 
- **Figure 1** was created by hand in Adobe Illustrator (the .ai file is included in the figures directory)
- **Supplemental figure 5** was prepared as described in Hodcroft et al. (2020), with different mutations targeted for the S:655 mutation {33269368}. Briefly: sequences with a mutation at nucleotide position 23525 (corresponding to a change at the 655 position in the spike glycoprotein) were selected from all available sequences on GISAID as of 29th December 2020. These sequences were included as the 'focal' set for a Nextstrain phylogenetic analysis, to which 'context' sequences were added, with the most genetically similar sequences given priority.
- **Supplemental figures 6 and 7** were created using a simple samtools command line prompts (below) and were then edited (e.g. colors, font sizes) for clarity in Adobe Illustrator. 

    ```bash 
    samtools sort SAMPLE.aln.sam > SAMPLE.sorted.bam
    samtools index SAMPLE.sorted.bam
    samtools depth -d 0 SAMPLE.sorted.bam > SAMPLE.sorted.coverage.tsv
    echo -e "Reference\tPosition\tDepth" | cat - SAMPLE.sorted.coverage.tsv
    SAMPLE.headers.sorted.coverage.tsv
    ```
### Contact info 

Please feel free to reach out to Katarina Braun (kmbraun2@wisc.edu or kmbraun2@gmail.com) or Gage Moreno (gkmoreno@wisc.edu) with questions. 
