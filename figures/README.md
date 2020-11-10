# Figures
To accompany the manuscript entitled "Transmission of SARS-CoV-2 in domestic cats imposes a narrow bottleneck"

The `figures` directory contains all of the vector files that the jupyter notebooks scripts output.   
In addition, the final Adobe Illustrator (formatted versions of the SVG files) can be found in the `illustrator_files` sub-directory. 


- All figures apart from **figure 1** and **supplemental figures 2 and 3** can be re-created using the above Jupyter Notebooks. 
- **Figure 1** was created by hand in Adobe Illustrator (the .ai file is included in the figures directory)
- **Supplemental figures 2 and 3** were created using a simple samtools command line prompts (below) and were then edited (e.g. colors, font sizes) for clarity in Adobe Illustrator. 

    ```bash 
    samtools sort SAMPLE.aln.sam > SAMPLE.sorted.bam
    samtools index SAMPLE.sorted.bam
    samtools depth -d 0 SAMPLE.sorted.bam > SAMPLE.sorted.coverage.tsv
    echo -e "Reference\tPosition\tDepth" | cat - SAMPLE.sorted.coverage.tsv
    SAMPLE.headers.sorted.coverage.tsv
    ```
### Contact info 

Please feel free to reach out to Katarina Braun (kmbraun2@wisc.edu or kmbraun2@gmail.com) or Gage Moreno (gkmoreno@wisc.edu) with questions. 
