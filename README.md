# SGRP5-genomics
Genomics pipelines

<p align="center">
  <img src="https://github.com/nicolo-tellini/SGRP5-genomics/blob/main/SGRP5org.png" alt="SGRP5org"/>
</p>

- fastQC evalution
  - fastQC [v.0.11.9](https://github.com/s-andrews/FastQC/releases/tag/v0.11.9)
  - MultiQC [v.1.12](https://github.com/ewels/MultiQC/releases/tag/v1.12)
- *spp*. compositions [sppComp-SGRP5](https://github.com/nicolo-tellini/sppComb) 
- phylogeny [phylo-SGRP5]
  - bwa mem [v.0.7.17-r1198-dirty]
  - samtools [v.1.14] (fixmate,sort,markdup,view,index) [workflow FASTQtoBAM](http://www.htslib.org/workflow/fastq.html)
  - bcftools [v.1.15.1] (mpileup,call,view,annotate,merge)
  - bgzip/tabix [v.1.11] 
  - PLINK [v1.90b6.21] 
  - R packages: 
    - SNPRelate [v.1.22.0]
    - gdsfmt [v.1.24.1]
    - ape [v.5.6-2]
    - treeio [v.1.12.0]
    - ggtree [v.2.2.4]
