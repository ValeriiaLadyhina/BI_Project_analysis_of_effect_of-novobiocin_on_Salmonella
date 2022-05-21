# A Reproducible Snakemake-based Workflow for Project "Studying _Salmonella_ Gene Expression Dynamics in Response to Novobiocin"
## Project was performed as a part of studies in [Bioinformatics institute](https://bioinf.me/en) on program [Bioinformatics for Biologists](https://bioinf.me/en/education#!/tab/40660730-1) during spring semester 2022.
<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Salmonella.png" width="300" height="300">

## README content
* [Introduction](#Introduction)
* [Background](#Background)
* [Aim](#Aim)
* [Project Objectives](#Project-Objectives)
* [Data Description](#Data-Description)
* [Snakemake Pipeline Description](#Snakemake-Pipeline-Description)
* [Installation of Snakemake](#Installation-of-Snakemake)
* [Start running Snakemake pipeline](#Start-Snakemake-pipeline)
  * [Few errors solutions](#Some-mistakes-solutions,-suggestions-and-explanations)
* [Results](#Results)
* [Useful References](#Useful-References)
* [Author and Acknowledgements](#Author-and-Acknowledgements)
* [Feedback contacts](#Feedback)


## Introduction

#### In this project we analysed publicly available [transcriptomic data](https://www.sciencedirect.com/science/article/pii/S2352340920301918#sec1)  obtained in time-series experiments on **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S** treated with different concentrations of antibiotic **novobiocin**.

## Background
#### DNA superspiralization is an important mechanism of gene regulation in bacteria and superspiralization alteration, for instance, in response to antibiotics can change gene expression levels. Antibiotics of aminocoumarin class, like novobiocin, are able to change both superspiralization of genomic DNA and expression of genes responding to superspiralization levels. The object of the study is _Salmonella enterica_, bacterium largely resistant to novobiocin and DNA superspiralization alterations in general. Here we aim to analyze gene expression changes in _Salmonella_ cultures grown on media with varying concentrations of antibiotics and to identify clusters of coexpressing genes.

## Aim
#### _To investigate effects of novobiocin on gene expression in **_Salmonella enterica_ subsp. _enterica_ serovar Typhimurium str. 14028S**._

## Project Objectives
#### 1. align reads to reference genome;
#### 2. obtain transcript count matrix;
#### 3. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) analysis of the count matrix;
#### 4. weighted gene co-expression network analysis using the [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/); 
#### 5. [Wopper](https://wopper.ba.itb.cnr.it) based visualization of modules of interest based on results of WGCNA analysis;
#### 6. develop a reproducible workflow for transcriptomic data analysis based on [Snakemake pipeline](https://snakemake.readthedocs.io/en/stable/).

## Data Description
#### Data is publicly available and can be found in the paper of [Natalia Gogoleva _et. al_, 2020](https://www.sciencedirect.com/science/article/pii/S2352340920301918?via%3Dihub). 
#### The provided dataset consists of RNA-seq reads obtained from samples of _S. enterica_ cultures treated with novobiocin at concentrations of 100 and 500 μg per mL, and untreated cultures. These reads were already filtered using [BBDuk (v. 37.23)](http://jgi.doe.gov/data-and-tools/bb-tools/). In our project we performed analysis only on reads from untreated culture and with addition of 500 μg per mL of novobiocin at three different timepoints: 10, 20 and 60 min.

## Snakemake Pipeline Desciption
#### Snakemake pipeline consists of 15 rules:
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) analysis of samples (SnakeMake Wrapper *v 1.5.0/bio/fastqc*)
* [Hisat2](http://daehwankimlab.github.io/hisat2/) Index (*v 2.2.0*)
* [Hisat2](http://daehwankimlab.github.io/hisat2/) Align (Snakemake wrapper *v 1.5.0/bio/hisat2/align*)
* [Samtools](http://daehwankimlab.github.io/hisat2/) sort (Snakemake wrapper *v1.5.0/bio/samtools/sort*)
* [StringTie](https://ccb.jhu.edu/software/stringtie/) assembly (*v 2.2.1)
* Matrix creation by [prepDE.py3](https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3) (Python3 script)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) (R Markdown)
* Genes' names converter neg _pos (Python script)
* [KEGG](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html)_[WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) (R Markdown)
* Genes' names converter modules (Python script)
* KEGG modules (R Markdown)
* Genes' names converter hubs (Python script)
* KEGG hub (R Markdown)
* Reverse path (Python script)
* Final Report (R Markdown)

##### R markdown files will automatically install the following packages during the run of Snakemake pipeline:
*__DESeq2, dplyr, ggplot2, gplots,clusterProfiler, WGCNA, knitr, KEGGREST,stringr__*

##### The rest of the used software will be installed during the creation of snakemake environment. All installed tools can be found in file config.yaml.

## Installation of Snakemake
#### Full installation description you can find on the [webpage of Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
#### 1) Install a Conda-based Python3 distribution
#### 2) Install [Mamba](https://github.com/mamba-org/mamba)
##### Download miniforge
```markdown
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-MacOSX-x86_64.sh"
bash Mambaforge-pypy3-MacOSX-x86_64.sh
```
or
```markdown
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-pypy3-MacOSX-x86_64.sh"
bash Mambaforge-pypy3-MacOSX-x86_64.sh
```
##### Install mamba
```markdown
conda install -n base -c conda-forge mamba
```
#### 3) Create snakemake environment
```markdown
conda activate base
mamba env create --name snakemake --file environment.yaml
```
- To delete environment
```markdown
conda env remove -n snakemake
```
#### 4) Activate snakemake environment
```markdown
conda activate snakemake
```
- To deactivate snakemake environment
```markdown
conda deactivate
```
##### Snakemake help
```markdown
snakemake --help
```
#### 5) Dowmload repositoty

#### 6) Go to folder snakemake
```markdown
cd Path/to/folder/snakemake
```
## Start Snakemake pipeline
#### Data
##### RNA seq data
1) Download needed RNA seq data sets to folder "snakemake/data/reads".
2) Make sure that you have all needed files to run pipeline in folder "snakemake/data:
   * [Reference genome](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/022/165/GCF_000022165.1_ASM2216v1/GCF_000022165.1_ASM2216v1_genomic.fna.gz) 
   * Reference annotation of 
   * pheno.csv # file with description of classes of experimental data
   * File sorted_p_adj_genes.txt
   * genefile.txt

```markdown
snakemake --use-conda --cores all
```
####  Some mistakes solutions, suggestions and explanations
##### Possible mistakes and their solutions
##### LockException
```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to 
create the same files in the following directory:
```
##### LockException Solution
```markdown
snakemake --unlock argument
```
##### Suggestions
```markdown
--cores all # use all possible cores
--cores N   # use N number of cores
```
```markdown
snakemake --use-conda --cores all -F # force to rerun fully Snakemake pipeline, otherwise 
                                     # it will run only that parts that did not run before
```
##### Explanations
*This is how the working pipeline starting its processes. In this part Snakemake describes which process will be performed and how many jobs per process (count).
```
Building DAG of jobs...
Using shell: /bin/bash
Job stats:
job                              count    min threads    max threads
-----------------------------  -------  -------------  -------------
KEGG_WGCNA                           1              1              1
KEGG_hub                             1              1              1
KEGG_modules                         1              1              1
all                                  1              1              1
deseq                                1              1              1
fastqc                              18              1              1
gene_name_converter_hubs             1              1              1
gene_name_converter_modules          1              1              1
gene_names_convereter_neg_pos        1              1              1
hisat2_align                        18              2              2
hisat2_index                         1              1              1
matrix_creation                      1              1              1
reverse_path                         1              1              1
samtools_sort                       18              8              8
stringtie_assembly                  18              8              8
total                               83              1              8

```

## Results
#### Description of project results can be found in presentation in folder results

### Additional analysis of the AT composition
##### Analysis of the AT-composition of the studied nodal genes; presented in a separate directory "AT_composition". The course of the analysis is written in the AT-composition notebook. Also in the directory is an additional file sequence.txt required for analysis (preprocessed file obtained from the annotation).

### The result of Snakemake is Final_report.html file and opened in your browser maps of pathes of interest obtained via KEGG. In the results folder there is example of such file.

## Useful References
* [Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research 2021, 10:33](https://doi.org/10.12688/f1000research.29032.2)
* [Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019)](https://www.nature.com/articles/s41587-019-0201-4)
* [Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.](https://doi.org/10.1186/s13059-014-0550-8)
* [Puccio, S., Grillo, G., Licciulli, F., Severgnini, M., Liuni, S., Bicciato, S., Peano, C. (2017). WoPPER: Web server for Position Related data analysis of gene Expression in Prokaryotes. Nucleic Acids Research, 45(W1), W109–W115.](https://academic.oup.com/nar/article/45/W1/W109/3782601)
* [Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834](https://www.degruyter.com/document/doi/10.2202/1544-6115.1128/html)
* [Twelve years of SAMtools and BCFtools Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=false)
* [Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M Transcriptome assembly from long-read RNA-seq alignments with StringTie2, Genome Biology 20, 278 (2019), doi:10.1186/s13059-019-1910-1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1)
* [Muskhelishvili G, Forquet R, Reverchon S, Meyer S, Nasser W. Coherent Domains of Transcription Coordinate Gene Expression During Bacterial Growth and Adaptation. Microorganisms. 2019 Dec 13;7(12):694. doi: 10.3390/microorganisms7120694. PMID: 31847191; PMCID: PMC6956064.](https://pubmed.ncbi.nlm.nih.gov/31847191/)


## Authors and Acknowledgements
#### Students
* [Valeriia Ladyhina](https://github.com/ValeriiaLadyhina)
* [Semyon Kupriyanov](https://github.com/immbiochem)
#### Supervisor
* [Alexandr Tkachenko](https://github.com/castrofiber)

## Feedback
If you have any questions, suggestions or complains please approach authors via email: [valeriia.ladyhina@slu.se](https://internt.slu.se/cv-originalen/valeriia-ladyhina/)
