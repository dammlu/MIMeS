# Functional metagenomics and supervised machine learning reveal microbial signatures of deep-sea mining disturbance in abyssal sediments

## Description
This repository contains the bioinformatic pipeline and R scripts 
used in Damm et al. (2026)

The pipeline covers:
- Quality trimming and read processing (BBDuk, clumpify)
- Metagenomic assembly (MEGAHIT)
- Gene prediction and rRNA removal (Prodigal, SortMeRNA)
- Gene catalog construction at 95% and 99% nucleotide identity (MMseqs2)
- Taxonomic profiling (phyloFlash (EMIRGE, SPAdes))
- Abundance estimation and CLR transformation (kallisto, compositions)
- Feature selection (Boruta, r2VIM, RFE, PIMP), and supervised machine learning (Random Forest, Support Vector Machine)
- Statistical analyses

## Repository structure

/Metagenomic_pipelines/ --> Bash scripts for the bioinformatic pipeline
/Feature_Selection_Supervised_Machine_Learning_pipeline/ --> R scripts for feature selection, machine learning, statistics, and visualization

## Requirements
- Software: BBMap, MEGAHIT, Prodigal, SortMeRNA, MMseqs2, 
  kallisto, phyloFlash
- R version: 4.3.3
- Key R packages: caret, kernlab, Boruta, pomona, vegan (see 
  manuscript Methods for full list and versions)

## Data availability
Raw sequences are available via NCBI SRA under accession: will be submitted later

## Citation
Damm et al. (2026). Functional metagenomics and supervised machine 
learning reveal microbial signatures of deep-sea mining disturbance 
in abyssal sediments. Nature Communications. [DOI will be submitted later]

## Contact
Lukas Damm — ldamm@rptu.de
