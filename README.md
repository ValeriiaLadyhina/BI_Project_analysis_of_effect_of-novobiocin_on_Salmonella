# A Reproducible Workflow for Project "Studying Salmonella Gene Expression Dynamics in Response to Novobiocin" that was performed as a part of studies in [Bioinformatics institute](https://bioinf.me/en) during spring semester 2022.

## Intoduction

#### In this project we analysed publicly available [transcriptomic data](https://www.sciencedirect.com/science/article/pii/S2352340920301918#sec1) 
#### obtained in time-series experiments on **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S** treated 
#### with different concentrations of antibiotic **novobiocin**.

## Background

## Aim
#### _To investigate effects of novobiocin on gene expression in **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S**._

## Project objectives
https://www.sciencedirect.com/science/article/pii/S2352340920301918#sec1#### 1. align reads to reference genome;
#### 2. obtain transcript count matrix;
#### 3. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) analysis of calculated expression matrix;
#### 4. weighted gene co-ecspression network analysis using the [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/); 
#### 5. [Wopper](https://wopper.ba.itb.cnr.it) based visualisation of moduls of interest based on results of WGCNA analysis;
#### 6. develop a reproducible workflow for transcroptomic data analysis based on [snakemake pipeline](https://snakemake.readthedocs.io/en/stable/).

## Data description
#### Data is publicly available and can be found in paper of [Natalia Gogoleva _et. al_, 2020](). 

## Snakemake pipeline


## Results