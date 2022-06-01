# A Reproducible Snakemake-based Workflow for Project "Studying _Salmonella_ Gene Expression Dynamics in Response to Novobiocin"
## Project was performed as a part of studies in [Bioinformatics institute](https://bioinf.me/en) on program [Bioinformatics for Biologists](https://bioinf.me/en/education#!/tab/40660730-1) during spring semester 2022.
<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/extra/Salmonella.png" width="300" height="300">

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
#### 3) Download [repository](https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella) and go to folder of snakemake
```markdown
cd Path/to/folder/snakemake
```
#### 4) Create snakemake environment
```markdown
conda activate base
mamba env create --name snakemake --file environment.yaml
```
- To delete environment
```markdown
conda env remove -n snakemake
```
#### 5) Activate snakemake environment
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
##### CreateCondaEnvironmentException
```
CreateCondaEnvironmentException:
Your conda installation is not configured to use strict channel priorities. This is however crucial for having robust and correct environments (for details, see https://conda-forge.org/docs/user/tipsandtricks.html). Please configure strict priorities by executing 'conda config --set channel_priority strict'.
```
##### CreateCondaEnvironmentException solution
```markdown

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

### DESeq2
Genes for which log2FoldChange<1.5 and FDR<0.05 were considered insignificant; a total of 789 significant genes were identified.
370 positive genes and 419 negative genes were found:
<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Volcano_plot.png" width="400" height="400">

A large number of differentially expressed genes does not give a clear idea of the nature of the changes under the action of an antibiotic:
<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/DEGs_Matrix_plot.png" width="400" height="400">

For positive and negative differentially expressed genes, a KEGG enrichment analysis was performed.
Positive DEGs:
- C5-Branched dibasic acid metabolism 
- Valine, leucine and isoleucine biosynthesis
- Two-component system 
- 2-Oxocarboxylic acid metabolism 
- Ascorbate and aldarate metabolism 

Negative DEGs:
- Flagellar assembly 
- Oxidative phosphorylation 
- Glycolysis / Gluconeogenesis 
- Bacterial chemotaxis 
- Starch and sucrose metabolism 

### WGCNA

#### Samples clustering:

Biological replicas are best grouped with each other, with the exception of sample 25. 
The processed samples at the time points of 10 and 20 minutes and the control samples at the same time points were grouped into two separate groups. 
The control and treated samples were separated into a separate group at 60 minutes. 
The data demonstrate an increase in the similarity of gene expression profiles of treated and untreated samples at 60 minutes.

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Samples_clustering.png" width="400" height="400">

*The image shows hierarchical clustering of samples (two variables are highlighted: Treated: 0 - white, 1 - black; Time: 10 - pink, 20 - blue, 60 - green).*

#### Modules building:
After building the network, 10 connected modules and 1 zero module were defined. Number of genes in modules: 154, 056, 937, 927, 462, 425, 294, 179, 176, 72, 47.
**Module-clustering dendrogram:**

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Module-clustering_dendrogram.png" width="400" height="400">

**Modules-Traits relationships:**

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Modules-Traits_relationships.png" width="400" height="400">

*Modules associated with the “Treated” variable: brown, turquoise, black, yellow. Modules associated with the “Time” variable: green, blue, red.*

**Modules clustering:**

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Modules_clustering.png" width="400" height="400">

*Modules black, brown, yellow with a strong positive correlation with the “Treated” variable were combined into one group.*

For significant modules, KEGG enrichment analysis was carried out.
**"Treated"-associated modules**

Negative module:
- Glycolysis / Gluconeogenesis 
- Amino sugar and nucleotide sugar metabolism 
- ABC transporters 
- Porphyrin metabolism 
- Histidine metabolism 
- beta-Lactam resistance 

Positive modules:
- Cationic antimicrobial peptide (CAMP) resistance 
- Two-component system   
- Lipopolysaccharide biosynthesis  
- C5-Branched dibasic acid metabolism 
- 2-Oxocarboxylic acid metabolism 
- Two-component system 
- Ascorbate and aldarate metabolism 
- Bacterial secretion system 

**"Time"-associated modules**

Positive module:
- Sulfur metabolism 
- Lysine degradation 
- Arginine and proline metabolism 
- Microbial metabolism in diverse environments 

Negative modules:
- Sulfur relay system  
- Ribosome   
- Aminoacyl-tRNA biosynthesis


#### Hub genes analysis

The functionality of the R language was used to identify the hub genes in the module (hub genes were considered to have a correlation coefficient with the eigengene of more than 0.9).
We obtained hub genes for modules associated with novobiocin and time. Number of hub genes for different modules: black - 42; brown - 167; yellow - 108; turquoise - 284; red - 94; blue - 147; green - 160.
For hub genes from modules, KEGG enrichment analysis was carried out.

**Pathways were found that were consistently found both among the hub genes of the modules and among the differentially expressed genes:**

- C5-Branched dibasic acid metabolism;
- 2-Oxocarboxylic acid metabolism;
- Two-component system;
- Glycolysis / Gluconeogenesis.

These pathways are the most sensitive pathways to novobiocin treatment; they were used to visualize the dynamics of their expression.

### Analysis of the dynamics of significant pathways

The visualization of the dynamics of the hub genes pathways was carried out. As a generalization of gene expression, the first principal component was used for a subset of genes included in the module and included in a significant biological pathway.

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Two_component_system.png" width="200" height="200">

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/C5_Branched_dibasic_acid_metabolism.png" width="200" height="200">

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/2-Oxocarboxylic_acid_metabolism.png" width="200" height="200">

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/Glycolysis_Gluconeogenesis.png" width="200" height="200">

### Mapping modules per chromosome

The modules associated with "Treated" and "Time" were mapped to the chromosome using the WoPPER online tool: https://wopper.ba.itb.cnr.it/WoPPER#!/View/v1nnzudodwda9keleutqqrw8znqo 

Time modules: | Treated modules:
--------------|-----------------
<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/WoPPER_Time.png" width="200" height="200"> | <img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/WoPPER_Treated.png" width="200" height="200">


Treated modules:

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/WoPPER_Treated.png" width="200" height="200">

### AT-composition

The topological state of DNA influences its affinity for some DNA binding proteins, especially in DNA sequences that have a high A + T base content.
We compared the AT composition of the hub gene sequences associated with novobiocin treatment with the AT composition of the hub genes associated with time.

The AT composition was calculated as the sum of A and T divided by the length of the sequence. This value was calculated for each gene from the list of hub genes. 
The results for the “Treated” and “Time” modules were combined. The distributions for the values were significantly different from normal (p-value = 2.2e-16 and p-value = 1.273e-12 for the "Treated" and "Time" modules, respectively), so the nonparametric Mann-Whitney test was used. 
The groups were significantly different from each other (p-value = 6.823e-11), although the difference between the means was small (0.496 and 0.470 for "Time" and "Treated", respectively).

The slight difference in the AT composition can be explained by the fact that it is important for supercoiling for the promoter regions of genes, so we studied it in the initial segments of the sequence using the sliding window method.
A "window" of 30 nucleotides long (the approximate length of a promoter) was iteratively shifted one nucleotide from the beginning of the gene sequence; 
The AT composition was calculated within the limits of the window according to the above method.
For each module, a sequence of values of the AT composition averaged for the hub genes of the module within the "window" is calculated.

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/AT_composition.png" width="400" height="400">

*The graph shows the frequencies of A and T, calculated within a 30-nucleotide "window". The X-axis shows a section from 1 to 300 nucleotides of the gene under study. The Y-axis plots the average frequency of A and T (AT-composition) for the hub genes of the module in a given position of the "window". Blue shows the frequencies for modules that are positively correlated with the "Treated" variable; the module shown in blue is negatively correlated with the "Treated" variable; modules in purple are not correlated with the "Treated" variable (associated with the "Time" variable).*

<img src="https://github.com/ValeriiaLadyhina/BI_Project_analysis_of_effect_of-novobiocin_on_Salmonella/blob/main/Results/Images/AT_mean_composition.png" width="400" height="400">

*The negatively correlated module showed frequencies close to the "Time" modules, so we focused on the positive "Treated" modules. For positive modules "Treated" and for modules "Time" the values were averaged. Green indicates A and T frequencies for positive "Treated" modules, blue indicates frequencies for "Time" modules. The X and Y values correspond to the previous plot.* 

## Discussion

Biological pathways may have different dynamics over time; Two-component system genes increase their expression by 20 minutes with a subsequent decrease; 
C5-Branched dibasic acid metabolism and 2-Oxocarboxylic acid metabolism genes increase expression over time; 
Glycolysis/Gluconeogenesis genes increase their expression over time in control samples, however, when treated with novobiocin, they do not show any noticeable dynamics.
Two-component system is a system for perceiving changes in the environment; this system can also be associated with a change in supercoiling, which also perceives a large number of stress stimuli (pH, osmotic composition, etc.). 
The lack of dynamics in the expression of glycolysis pathway genes may be associated with the bacteriostatic effect of novobiocin (these processes are associated with anabolic pathways); 
C5-Branched dibasic acid metabolism and 2-Oxocarboxylic acid metabolism may be involved in some of the signaling pathways associated with changes in supercoiling.

The co-expression modules show a diffuse distribution along the length of the chromosome, which may correspond to the influence of a systemic process, such as a change in supercoiling.
Genes in co-expressed modules are located at significant distances from each other (more than 1Mp). This can be explained by the topology of the chromosome, but does not exclude the influence of DNA supercoiling on the activation of modules.
The initial regions of genes sensitive to novobiocin are characterized by a richer AT composition, which suggests sensitivity to changes in DNA supercoiling.


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
