# A Reproducible Snakemake-based Workflow for Project "Studying _Salmonella_ Gene Expression Dynamics in Response to Novobiocin"
## Project was performed as a part of studies in [Bioinformatics institute](https://bioinf.me/en) on program [Bioinformatics for Biologists](https://bioinf.me/en/education#!/tab/40660730-1) during spring semester 2022.

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

## Instalation of Snakemake
#### Full installation description you can find on the [webpage of Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
#### 1) Install a Conda-based Python3 distribution
#### 2) Install [Mamba](https://github.com/mamba-org/mamba)
```markdown
conda install -n base -c conda-forge mamba
```
#### 3) Create snakemake environment
```markdown
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
##### To delete environment
```markdown
conda env remove -n snakemake
```
#### 4) Activate snakemake environment
```markdown
conda activate snakemake
```
##### To deactivate snakemake environment
```markdown
conda deactivate
```
##### Snakemake help
```markdown
snakemake --help
```
#### 5) Clone [repositoty]() to folder snakemake


#### 6) Install pandoc *v 2.12
```markdown
conda install pandoc
```
#### 7) Go to folder snakemake
```markdown
cd Path/to/folder/snakemake
```
#### 8) Start snakemake 
```markdown
snakemake --use-conda --cores all
```
#### 9) Some mistakes solutions, suggestions and explanations
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
```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 10
Rules claiming more threads will be scaled down.
Job stats:
job                   count    min threads    max threads
------------------  -------  -------------  -------------
all                       1              1              1
deseq                     1              1              1
fastqc                   18              1              1
hisat2_align             18              2              2
hisat2_index              1              1              1
matrix_creation           1              1              1
samtools_sort            18              8              8
stringtie_assembly       18              8              8
total                    76              1              8

Select jobs to execute...

```
## Repository map
* 
## Results

## Usefull References
*  [Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 2; peer review: 2 approved]. F1000Research 2021, 10:33](https://doi.org/10.12688/f1000research.29032.2)

## Authors and Acknowledgements
#### Students
* [Valeriia Ladyhina](https://github.com/ValeriiaLadyhina)
* [Semyon Kupriyanov](https://github.com/immbiochem)
#### Supervisor
* [Alexandr Tkachenko](https://github.com/castrofiber)

## Feedback
If you have any questions, suggestions or complains please approach authors via email: [valeriia.ladyhina@slu.se](https://internt.slu.se/cv-originalen/valeriia-ladyhina/)