
# Kostyrko Multiomics Analysis for Methylation EPIC and RNAseq
# Description
This repository contains a set of R scripts specifically designed for the data analysis performed in Kostyrko et al. The repository provides a set of R scripts for analyzing Methylation EPIC and RNAseq (correlation and survival analysis) used in the original paper. This respository is a resource for researchers looking to replicate or build upon the findings of Kostyrko et al.
# Features
## Methylation EPIC Analysis & RNAseq Correlation: 
  + This repo is tailored to handle Methylation EPIC data and RNAseq correlation analysis as performed in the Kostyrko et al study.
  + Survival Analysis & Miscellaneous RNAseq Analysis: 
# Data
The dataset used for the Kostyrko et al study can be accessed from the GEO SuperSeries under the accession number GSE198289.
### This SuperSeries is composed of the following SubSeries:
  + GSE198289	UHRF1 is a mediator of KRAS driven oncogenesis in lung adenocarcinoma [RNA-seq]
  + GSE198446	UHRF1 is a mediator of KRAS driven oncogenesis in lung adenocarcinoma [epic_methyl]
  + GSE209923	UHRF1 is a mediator of KRAS driven oncogenesis in lung adenocarcinoma [shRNA]
# Minimum Requirements
  * Download and install R from here. R (3.6.2) 
  * Download the dataset from GEO SuperSeries GSE198289.
  * R packages
 
    + STAR (2.5.3a)
    + edgeR (3.28.1) 
    + Minfi (1.32)
    + EpiDISH (2.2.2)
    + limma package(3.42.2)
    + missMethyl (1.20.4)
    + DMRcate (2.0.7)
    + survminer (0.4.6.999) and survival (3.1.11) 
    + DGCA (1.0.2 )
    + methylGSA (1.4.9) , Cytoscape ( 3.8.2) clueGO (2.5.9)

## Additional public datasets

| LUAD, survival, correlations | RNAseq Counts      | Broad GDAC ( Firehose ) | [https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz](https://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/gdac.broadinstitute.org_LUAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz) |
| ---------------------------- | ------------------ | ----------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Pathway                      | database           | MSigDB                  | [https://www.gsea-msigdb.org/gsea/msigdb](https://www.gsea-msigdb.org/gsea/msigdb)                                                                                                                                                                                                                                                                                                                                                     |
| EPIC annotation              | Probe annotations  | Bioconductor            | [https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylationEPICanno.ilm10b4.hg19.html](https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylationEPICanno.ilm10b4.hg19.html)                                                                                                                                                                                                 |
| TSG annotation               | Annotations        | UTHealth                | [https://bioinfo.uth.edu/TSGene/download.cgi?csrt=9609930864565957852](https://bioinfo.uth.edu/TSGene/download.cgi?csrt=9609930864565957852)                                                                                                                                                                                                                                                                                           |
| Probe Filter                 | Annotations        | github                  | [https://github.com/sirselim/illumina450k_filtering](https://github.com/sirselim/illumina450k_filtering)                                                                                                                                                                                                                                                                                                                               |
| ERV                          | Genomic annotation | UCSC genome             | [https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doMainPage=1&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk](https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doMainPage=1&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk)                                                                                                                                                                                                                 |
| RNAseq index                 | Genomic index      | GENCODE                 | [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz)                                                                                                                                                                                                                     |
| RNAseq GTF                   | Gene Annotation    | GENCODE                 | [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz)                                                                                                                                                                                                                                   |
