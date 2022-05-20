# A Reproducible Workflow for Project "Studying _Salmonella_ Gene Expression Dynamics in Response to Novobiocin" that was performed as a part of studies in [Bioinformatics institute](https://bioinf.me/en) on program [Bioinformatics for Biologists](https://bioinf.me/en/education#!/tab/40660730-1) during spring semester 2022.

## Intoduction

#### In this project we analysed publicly available [transcriptomic data](https://www.sciencedirect.com/science/article/pii/S2352340920301918#sec1) 
#### obtained in time-series experiments on **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S** treated 
#### with different concentrations of antibiotic **novobiocin**.

## Background

## Aim
#### _To investigate effects of novobiocin on gene expression in **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S**._

## Project Objectives
#### 1. align reads to reference genome;
#### 2. obtain transcript count matrix;
#### 3. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) analysis of calculated expression matrix;
#### 4. weighted gene co-ecspression network analysis using the [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/); 
#### 5. [Wopper](https://wopper.ba.itb.cnr.it) based visualisation of moduls of interest based on results of WGCNA analysis;
#### 6. develop a reproducible workflow for transcriptomic data analysis based on [Snakemake pipeline](https://snakemake.readthedocs.io/en/stable/).

## Data Description
#### Data is publicly available and can be found in paper of [Natalia Gogoleva _et. al_, 2020](). 

## Snakemake Pipeline Desciption
#### Snakemake pipeline consist of ... number of rules:
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) analysis of samples (SnakeMake Wrapper *v 1.5.0/bio/fastqc*)
* [Hisat2](http://daehwankimlab.github.io/hisat2/) Index (*v 2.2.0*)
* [Hisat2](http://daehwankimlab.github.io/hisat2/) Align (Snakemake wrapper *v 1.5.0/bio/hisat2/align*)
* [Samtools](http://daehwankimlab.github.io/hisat2/) sort (Snakemake wrapper *v1.5.0/bio/samtools/sort*)
* [StringTie](https://ccb.jhu.edu/software/stringtie/) assembly (*v 2.2.1)
* Matrix creation by [prepDE.py3](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3) (Python3 script)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (R Markdown)
* Genes' names converter (Jupyter Notebook)
* 


## Results

## Usefull References

## Authors and Acknowledgements
#### Students
* [Valeriia Ladyhina](https://github.com/ValeriiaLadyhina)
* [Semyon Kupriyanov](https://github.com/immbiochem)
#### Supervisor
* [Alexandr Tkachenko](https://github.com/castrofiber)