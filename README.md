# sc_tools
this a repository used for collection single cell tools.

## scRNAseq tools
### piplines
##### A flexible cross-platform single-cell data processing pipeline(NC2022)
(https://www.nature.com/articles/s41467-022-34681-z) 
We introduce UniverSC (https://github.com/minoda-lab/universc), a universal single-cell RNA-seq data processing tool that supports any unique molecular identifier-based platform. Our command-line tool, docker image, and containerised graphical application enables consistent and comprehensive integration, comparison, and evaluation across data generated from a wide range of platforms. We also provide a cross-platform application to run UniverSC via a graphical user interface, available for macOS, Windows, and Linux Ubuntu, negating one of the bottlenecks with single-cell RNA-seq analysis that is data processing for researchers who are not bioinformatically proficient.

##### Multi-level cellular and functional annotation of single-cell transcriptomes using scPipeline
https://www.nature.com/articles/s42003-022-04093-2#code-availability (Communications Biology 2022)
Here, we introduce scPipeline, a single-cell data analysis toolbox that builds on existing methods and offers modular workflows for multi-level cellular annotation and user-friendly analysis reports. Advances to scRNA-seq annotation include: (i) co-dependency index (CDI)-based differential expression, (ii) cluster resolution optimization using a marker-specificity criterion, (iii) marker-based cell-type annotation with Miko scoring, and (iv) gene program discovery using scale-free shared nearest neighbor network (SSN) analysis. Both unsupervised and supervised procedures were validated using a diverse collection of scRNA-seq datasets and illustrative examples of cellular transcriptomic annotation of developmental and immunological scRNA-seq atlases are provided herein. Overall, scPipeline offers a flexible computational framework for in-depth scRNA-seq analysis.
In tandem to scPipeline, we developed the scMiko R package that comprises a collection of functions for application-specific scRNA-seq analysis and generation of scPipeline analytic reports. 
###### access
The scMiko R package (https://github.com/NMikolajewicz/scMiko) and scPipeline scripts (https://github.com/NMikolajewicz/scPipeline) are available on GitHub. Step-by-step tutorials and documentation are also provided at https://nmikolajewicz.github.io/scMiko/.

### methods
##### Online single-cell data integration through projecting heterogeneous datasets into a common cell-embedding space(NC 2022)
https://www.nature.com/articles/s41467-022-33758-z
 SCALEX, a deep-learning method that integrates single-cell data by projecting cells into a batch-invariant, common cell-embedding space in a truly online manner. SCALEX substantially outperforms online iNMF and other state-of-the-art non-online integration methods on benchmark single-cell datasets of diverse modalities.The online data integration capacity and superior performance makes SCALEX particularly appropriate for large-scale single-cell applications to build upon previous scientific insights.
###### access
 https://github.com/jsxlei/SCALEX
 
##### scJoint integrates atlas-scale single-cell RNA-seq and ATAC-seq data with transfer learning(2022 NBT)
(https://www.nature.com/articles/s41587-021-01161-6)
Here, we present scJoint, a transfer learning method to integrate atlas-scale, heterogeneous collections of scRNA-seq and scATAC-seq data. scJoint leverages information from annotated scRNA-seq data in a semisupervised framework and uses a neural network to simultaneously train labeled and unlabeled data, allowing label transfer and joint visualization in an integrative framework. 
###### access
https://github.com/SydneyBioX/scJoint 

##### Deep learning of cross-species single-cell landscapes identifies conserved regulatory programs underlying cell types(2022 NG)
(https://www.nature.com/articles/s41588-022-01197-7)
Nvwa (the name of a mother god in ancient Chinese legend),a deep-learning-based strategy, to predict gene expression from DNA sequence in individual cells. Finally, we interpreted the cell-type-specific sequence rules and identified regulatory programs conserved across species. 
###### access
https://github.com/JiaqiLiZju/Nvwa/

##### Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data (2022 NC)
(https://www.nature.com/articles/s41467-022-28803-w)
ScType is a computational platform, which enables data-driven, fully-automated and ultra-fast cell-type identification based solely on given scRNA-seq data, combined with our comprehensive cell marker database as background information.
###### access
R source-code of the ScType algorithm is available at https://github.com/IanevskiAleksandr/sc-type/. ScType is also freely available as an interactive web-tool at http://sctype.app. 
 
##### CellWalker integrates single-cell and bulk data to resolve regulatory elements across cell types in complex tissues(GB 2021)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02279-1
CellWalker, a method that integrates single-cell open chromatin (scATAC-seq) data with gene expression (RNA-seq) and other data types using a network model that simultaneously improves cell labeling in noisy scATAC-seq and annotates cell type-specific regulatory elements in bulk data. 
###### access
https://github.com/PollardLab/CellWalker 
 
##### BRIE2: computational identification of splicing phenotypes from single-cell transcriptomic experiments(GB 2021)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02461-5#availability-of-data-and-materials
However, the intrinsic coverage limitations of scRNA-seq technologies make it challenging to associate specific splicing events to cell-level phenotypes. BRIE2 is a scalable computational method that resolves these issues by regressing single-cell transcriptomic data against cell-level features. 
###### access
https://pypi.org/project/brie/ 

##### scSTEM: clustering pseudotime ordered single-cell data(GB 2021)
We develop scSTEM, single-cell STEM, a method for clustering dynamic profiles of genes in trajectories inferred from pseudotime ordering of single-cell RNA-seq (scRNA-seq) data. s
###### access
https://github.com/alexQiSong/scSTEM
 
##### Specific splice junction detection in single cells with SICILIAN(GB 2021)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02434-8
Precise splice junction calls are currently unavailable in scRNA-seq pipelines such as the 10x Chromium platform but are critical for understanding single-cell biology. Here, we introduce SICILIAN, a new method that assigns statistical confidence to splice junctions from a spliced aligner to improve precision. SICILIAN is a general method that can be applied to bulk or single-cell data, but has particular utility for single-cell analysis due to that data’s unique challenges and opportunities for discovery. SICILIAN’s precise splice detection achieves high accuracy on simulated data, improves concordance between matched single-cell and bulk datasets, and increases agreement between biological replicates.
###### access
 https://github.com/salzmanlab/SICILIAN
 
##### pipeComp, a general framework for the evaluation of computational pipelines, reveals performant single cell RNA-seq preprocessing tools(GB 2020) 
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02136-7
We present pipeComp (), a flexible R framework for pipeline comparison handling interactions between analysis steps and relying on multi-level evaluation metrics.
###### access
https://github.com/plger/pipeComp 
 
##### Phitest for analyzing the homogeneity of single-cell populations(Bioinformatics 2020) 
https://academic.oup.com/bioinformatics/article/38/9/2639/6541626?searchresult=1
We propose a bioinformatics tool named Phitest to analyze the homogeneity of single-cell populations. Phitest is able to distinguish between homogeneous and heterogeneous cell populations, providing an objective and automatic method to optimize the performance of single-cell clustering analysis.
###### access
https://github.com/Vivianstats/PhitestR

### tools
##### scIMC: a platform for benchmarking comparison and visualization analysis of scRNA-seq data imputation methods(NAR 2022)
(https://academic.oup.com/nar/article/50/9/4877/6582166?searchresult=1)
scIMC (single-cell Imputation Methods Comparison platform), the first online platform that integrates all available state-of-the-art imputation methods for benchmarking comparison and visualization analysis.
###### access
https://server.wei-group.net/scIMC/

##### TraSig: inferring cell-cell interactions from pseudotime ordering of scRNA-Seq data(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02629-7
TraSig, a computational method for improving the inference of cell-cell interactions in scRNA-Seq studies that utilizes the dynamic information to identify significant ligand-receptor pairs with similar trajectories, which in turn are used to score interacting cell clusters. 
###### access
https://github.com/doraadong/TraSig.
 
##### scDART: integrating unmatched scRNA-seq and scATAC-seq data and learning cross-modality relationship simultaneously(GB 2022) 
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02706-x
 scDART, a deep learning framework that integrates scRNA-seq and scATAC-seq data and learns cross-modalities relationships simultaneously. Specifically, the design of scDART allows it to preserve cell trajectories in continuous cell populations and can be applied to trajectory inference on integrated data.
###### access
https://github.com/PeterZZQ/scDART

##### scAI: an unsupervised approach for the integrative analysis of parallel single-cell transcriptomic and epigenomic profiles
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1932-8
Here, we present a single-cell aggregation and integration (scAI) approach to integrate transcriptomic and epigenomic profiles (i.e., chromatin accessibility or DNA methylation) that are derived from the same cells. Unlike existing integration methods , scAI takes into consideration the extremely sparse and near-binary nature of single-cell epigenomic data. Through iterative learning in an unsupervised manner, scAI aggregates epigenomic data in subgroups of cells that exhibit similar gene expression and epigenomic profiles.Those similar cells are computed through learning a cell-cell similarity matrix simultaneously from both transcriptomic and aggregated epigenomic data using a unified matrix factorization model. 
###### access
MATLAB package: https://github.com/amsszlh/scAI  and R package: https://github.com/sqjin/scAI

##### Interfacing Seurat with the R tidy universe (Bioinformatics 2021)
https://academic.oup.com/bioinformatics/article/37/22/4100/6283576?searchresult=1
To provide Seurat with a tidyverse-oriented interface without compromising efficiency, we developed tidyseurat, a lightweight adapter to the tidyverse. Tidyseurat displays cell information as a tibble abstraction, allowing intuitively interfacing Seurat with dplyr, tidyr, ggplot2 and plotly packages powering efficient data manipulation, integration and visualization. Iterative analyses on data subsets are enabled by interfacing with the popular nest-map framework.
###### access
https://github.com/stemangiola/tidyseurat

##### ShinyCell: simple and sharable visualization of single-cell gene expression data(Bioinformatics 2021)
https://academic.oup.com/bioinformatics/article/37/19/3374/6198103?searchresult=1
In order to help address this we have developed ShinyCell, an R package that converts single-cell RNA sequencing datasets into explorable and shareable interactive interfaces. These interfaces can be easily customized in order to maximize their usability and can be easily uploaded to online platforms to facilitate wider access to published data.
###### access
https://github.com/SGDDNB/ShinyCell

##### V-SVA: an R Shiny application for detecting and annotating hidden sources of variation in single-cell RNA-seq data(Bioinformatics 2020)
https://academic.oup.com/bioinformatics/article/36/11/3582/5771333?searchresult=1
we developed an R Shiny application [Visual Surrogate Variable Analysis (V-SVA)] that provides a web-browser interface for the identification and annotation of hidden sources of variation in scRNA-seq data. This interactive framework includes tools for discovery of genes associated with detected sources of variation, gene annotation using publicly available databases and gene sets, and data visualization using dimension reduction methods.
###### access
https://github.com/nlawlor/V-SVA

##### ScisorWiz: visualizing differential isoform expression in single-cell long-read data(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/13/3474/6590641?searchresult=1
Here, we present ScisorWiz, a streamlined tool to visualize isoform expression differences across single-cell clusters in an informative and easily communicable manner. ScisorWiz achieves this with an easy, fast and reliable method of visualizing differential isoform expression data across multiple clusters and is executable from the command line with the R language
###### access
http://github.com/ans4013/ScisorWiz

##### Sciviewer enables interactive visual interrogation of single-cell RNA-Seq data from the Python programming environment(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/37/21/3961/6380550?searchresult=1
Here, we present the Single-cell Interactive Viewer (Sciviewer), a tool that overcomes this limitation by allowing interactive visual interrogation of embeddings from within Python. Beyond differential expression analysis of user-selected cells, Sciviewer implements a novel method to identify genes varying locally along any user-specified direction on the embedding. Sciviewer enables rapid and flexible iteration between interactive and programmatic modes of scRNA-Seq exploration, illustrating a useful approach for analyzing high-dimensional data.
###### access
https://github.com/colabobio/sciviewer

##### Cellsnp-lite: an efficient tool for genotyping single cells(Bioinformatics 2021) 
https://academic.oup.com/bioinformatics/article/37/23/4569/6272512?searchresult=1 
Here, we introduce a software, cellsnp-lite, implemented in C/C++ and based on well-supported package htslib, for genotyping in single-cell sequencing data for both droplet and well-based platforms.
###### access
https://github.com/single-cell-genetics/cellsnp-lite

##### scTPA: a web tool for single-cell transcriptome analysis of pathway activation signatures(Bioinformatics 2020) 
https://academic.oup.com/bioinformatics/article/36/14/4217/5841657?searchresult=1
scTPA, a web-based platform for pathway-based analysis of single-cell RNA-seq data in human and mouse. scTPA incorporates four widely-used gene set enrichment methods to estimate the pathway activation scores of single cells based on a collection of available biological pathways with different functional and taxonomic classifications. 
###### access
http://sctpa.bio-data.cn/sctpa


### database
##### EmAtlas: a comprehensive atlas for exploring spatiotemporal activation in mammalian embryogenesis(NAR 2022)
EmAtlas, which collects the most comprehensive multi-omics data and provides multi-scale tools to explore spatiotemporal activation during mammalian embryogenesis. EmAtlas contains data on multiple types of gene expression, chromatin accessibility, DNA methylation, nucleosome occupancy, histone modifications, and transcription factors, which displays the complete spatiotemporal landscape in mouse and human across several time points, involving gametogenesis, preimplantation, even fetus and neonate, and each tissue involves various cell types. To characterize signatures involved in the tissue, cell, genome, gene and protein levels during mammalian embryogenesis, analysis tools on these five scales were developed.
###### access
http://bioinfor.imu.edu.cn/ematlas

##### HTCA: a database with an in-depth characterization of the single-cell human transcriptome(NAR 2022)
(https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac791/6709243?searchresult=1)
HTCA, an interactive database constructed based on ∼2.3 million high-quality cells from ∼3000 scRNA-seq samples and comprised in-depth phenotype profiles of 19 healthy adult and matching fetal tissues. HTCA provides a one-stop interactive query to gene signatures, transcription factor (TF) activities, TF motifs, receptor–ligand interactions, enriched gene ontology (GO) terms, etc. across cell types in adult and fetal tissues. At the same time, HTCA encompasses single-cell splicing variant profiles of 16 adult and fetal tissues, spatial transcriptomics profiles of 11 adult and fetal tissues, and single-cell ATAC-sequencing (scATAC-seq) profiles of 27 adult and fetal tissues. Besides, HTCA provides online analysis tools to perform major steps in a typical scRNA-seq analysis.
###### access
http://www.htcatlas.org/





## single-cell epigenomic tools
### tools
##### EpiScanpy: integrated single-cell epigenomic analysis(NC 2021)
(https://www.nature.com/articles/s41467-021-25131-3)
EpiScanpy is a toolkit for the analysis of single-cell epigenomic data, namely single-cell DNA methylation and single-cell ATAC-seq data. To address the modality specific challenges from epigenomics data, epiScanpy quantifies the epigenome using multiple feature space constructions and builds a nearest neighbour graph using epigenomic distance between cells. EpiScanpy makes the many existing scRNA-seq workflows from scanpy available to large-scale single-cell data from other -omics modalities, including methods for common clustering, dimension reduction, cell type identification and trajectory learning techniques, as well as an atlas integration tool for scATAC-seq datasets. The toolkit also features numerous useful downstream functions, such as differential methylation and differential openness calling, mapping epigenomic features of interest to their nearest gene, or constructing gene activity matrices using chromatin openness. We successfully benchmark epiScanpy against other scATAC-seq analysis tools and show its outperformance at discriminating cell types.
###### access
https://github.com/colomemaria/epiScanpy, https://colomemaria.github.io/episcanpy_doc

##### scMethBank: a database for single-cell whole genome DNA methylation maps(NAR 2022)
https://academic.oup.com/nar/article/50/D1/D380/6376025?searchresult=1#325783260
scMethBank, the first open access and comprehensive database dedicated to the collection, integration, analysis and visualization of single-cell DNA methylation data and metadata.Current release of scMethBank includes processed single-cell bisulfite sequencing data and curated metadata of 8328 samples derived from 15 public single-cell datasets, involving two species (human and mouse), 29 cell types and two diseases. 
###### access
 https://ngdc.cncb.ac.cn/methbank/scm/

##### The Neuroscience Multi-Omic Archive: a BRAIN Initiative resource for single-cell transcriptomic and epigenomic data from the mammalian brain(NAR 2022) 
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac962/6786191?searchresult=1
Neuroscience Multi-Omic Archive (NeMO Archive; nemoarchive.org), which serves as the primary repository for genomics data from the BRAIN Initiative. Working closely with other BRAIN Initiative researchers, we have organized these data into a continually expanding, curated repository, which contains transcriptomic and epigenomic data from over 50 million brain cells, including single-cell genomic data from all of the major regions of the adult and prenatal human and mouse brains, as well as substantial single-cell genomic data from non-human primates.
###### access
Data resources described in this manuscript are available in the NeMO Archive (https://nemoarchive.org/) and NeMO Analytics (https://nemoanalytics.org/).

##### SCAFE: a software suite for analysis of transcribed cis-regulatory elements in single cells(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/22/5126/6730725?searchresult=1
We developed SCAFE—Single-Cell Analysis of Five-prime Ends—a software suite that processes sc-end5-seq data to de novo identify TSS clusters based on multiple logistic regression. It annotates tCREs based on the identified TSS clusters and generates a tCRE-by-cell count matrix for downstream analyses. The software suite consists of a set of flexible tools that could either be run independently or as pre-configured workflows.
###### access
https://github.com/chung-lab/SCAFE

##### powerEQTL: an R package and shiny application for sample size and power calculation of bulk tissue and single-cell eQTL analysis(Bioinformatics 2021)
https://academic.oup.com/bioinformatics/article/37/22/4269/6278296?searchresult=1
Here, we presented an R package called powerEQTL with flexible functions to estimate power, minimal sample size or detectable minor allele frequency for both bulk tissue and single-cell eQTL analysis. A user-friendly, program-free web application is also provided, allowing users to calculate and visualize the parameters interactively.
###### access
powerEQTL R package: https://cran.r-project.org/web/packages/powerEQTL/
powerEQTL shiny application:https://bwhbioinfo.shinyapps.io/powerEQTL/

##### scGAD: single-cell gene associating domain scores for exploratory analysis of scHi-C data(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/14/3642/6598798?searchresult=1
We present single-cell gene associating domain (scGAD) scores as a dimension reduction and exploratory analysis tool for scHi-C data. scGAD enables summarization at the gene unit while accounting for inherent gene-level genomic biases. Low-dimensional projections with scGAD capture clustering of cells based on their 3D structures. Significant chromatin interactions within and between cell types can be identified with scGAD. We further show that scGAD facilitates the integration of scHi-C data with other single-cell data modalities by enabling its projection onto reference low-dimensional embeddings. This multi-modal data integration provides an automated and refined cell-type annotation for scHi-C data.
###### access
https://sshen82.github.io/BandNorm/articles/scGAD-tutorial.html

### methods
##### MIRA: joint regulatory modeling of multimodal expression and chromatin accessibility in single cells(NM 2022)
(https://www.nature.com/articles/s41592-022-01595-z)
Rigorously comparing gene expression and chromatin accessibility in the same single cells could illuminate the logic of how coupling or decoupling of these mechanisms regulates fate commitment. Here we present MIRA, probabilistic multimodal models for integrated regulatory analysis, a comprehensive methodology that systematically contrasts transcription and accessibility to infer the regulatory circuitry driving cells along cell state trajectories. MIRA leverages topic modeling of cell states and regulatory potential modeling of individual gene loci. MIRA thereby represents cell states in an efficient and interpretable latent space, infers high-fidelity cell state trees, determines key regulators of fate decisions at branch points and exposes the variable influence of local accessibility on transcription at distinct loci. Applied to epidermal differentiation and embryonic brain development from two different multimodal platforms, MIRA revealed that early developmental genes were tightly regulated by local chromatin landscape whereas terminal fate genes were titrated without requiring extensive chromatin remodeling.
###### access
MIRA is available as a Python package at https://github.com/cistrome/MIRA. Frankencell, a Python program we developed to generate synthetic differentiation trajectories for benchmarking, is available at https://github.com/AllenWLynch/frankencell-dynverse.

##### scBasset: sequence-based modeling of single-cell ATAC-seq using convolutional neural networks(NM 2022)
(https://www.nature.com/articles/s41592-022-01562-8)
Here we introduce scBasset, a sequence-based convolutional neural network method to model scATAC data. We show that by leveraging the DNA sequence information underlying accessibility peaks and the expressiveness of a neural network model, scBasset achieves state-of-the-art performance across a variety of tasks on scATAC and single-cell multiome datasets, including cell clustering, scATAC profile denoising, data integration across assays and transcription factor activity inference.
###### access
https://github.com/calico/scBasset

##### Predictive modeling of single-cell DNA methylome data enhances integration with transcriptome data(2020 genome research)
(https://genome.cshlp.org/content/31/1/101.full?sid=7bc26163-4a2c-40cd-a46c-48a3b62f6d49)
MAPLE, a computational framework that learns the association between DNA methylation and expression using both gene- and cell-dependent statistical features. Using multiple data sets generated with different experimental protocols, we show that using predicted gene activity values significantly improves several analysis tasks, including clustering, cell type identification, and integration with transcriptome data. 
###### access
https://github.com/tanlabcode/MAPLE.1.0

##### scGWAS: landscape of trait-cell type associations by integrating single-cell transcriptomics-wide and genome-wide association studies(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02785-w
we propose scGWAS (scRNA-seq assisted GWAS analysis) to investigate the transcriptional changes of genetic variants in specific cell type contexts by leveraging a wide variety of gene–gene relationships in the human genome. scGWAS can not only identify the genetically mediated associations between cell types and traits but also construct the biological networks that are overrepresented with disease risk genes and transcriptionally active genes in a cell type. As shown below, scGWAS utilizes only the average gene expression for each cell type, which makes it scalable to large scRNA-seq datasets. 
###### access
The scGWAS package is available at GitHub

##### Regulatory analysis of single cell multiome gene expression and chromatin accessibility data with scREG(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02682-2
We develop scREG, a dimension reduction methodology, based on the concept of cis-regulatory potential, for single cell multiome data. This concept is further used for the construction of subpopulation-specific cis-regulatory networks. The capability of inferring useful regulatory network is demonstrated by the two-fold increment on network inference accuracy compared to the Pearson correlation-based method and the 27-fold enrichment of GWAS variants for inflammatory bowel disease in the cis-regulatory elements. The R package scREG provides comprehensive functions for single cell multiome data analysis.
###### access
https://github.com/Durenlab/RegNMF 

##### CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02583-w
We introduce CellPhy, a maximum likelihood framework for inferring phylogenetic trees from somatic single-cell single-nucleotide variants. CellPhy leverages a finite-site Markov genotype model with 16 diploid states and considers amplification error and allelic dropout.
###### access
https://github.com/amkozlov/cellphy

##### SCONCE: a method for profiling copy number alterations in cancer evolution using single-cell whole genome sequencing(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/7/1801/6515610?searchresult=1
We present a theoretical framework to use tumor evolutionary history to accurately call CNAs in a principled manner. To model the tumor evolutionary process and account for technical noise from low coverage single-cell whole genome sequencing data, we developed SCONCE, a method based on a Hidden Markov Model to analyze read depth data from tumor cells using matched normal cells as negative controls. Using a combination of public data sets and simulations, we show SCONCE accurately decodes copy number profiles, and provides a useful tool for understanding tumor evolution.
###### access
https://github.com/NielsenBerkeleyLab/sconce

##### SCONCE: a method for profiling copy number alterations in cancer evolution using single-cell whole genome sequencing(Bioinformatics 2022)
We present a theoretical framework to use tumor evolutionary history to accurately call CNAs in a principled manner. To model the tumor evolutionary process and account for technical noise from low coverage single-cell whole genome sequencing data, we developed SCONCE, a method based on a Hidden Markov Model to analyze read depth data from tumor cells using matched normal cells as negative controls.
###### access
https://github.com/NielsenBerkeleyLab/sconce

##### SimSCSnTree: a simulator of single-cell DNA sequencing data(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/10/2912/6551250?searchresult=1
We report on a new single-cell DNA sequence simulator, SimSCSnTree, which generates an evolutionary tree of cells and evolves single nucleotide variants (SNVs) and copy number aberrations (CNAs) along its branches. Data generated by the simulator can be used to benchmark tools for single-cell genomic analyses, particularly in cancer where SNVs and CNAs are ubiquitous.
###### access
https://github.com/compbiofan/SimSCSnTree.git 



## ST tools
### methods
##### Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data with SpaTalk(2022 NC)
(https://www.nature.com/articles/s41467-022-32111-8)
Spatially resolved transcriptomics provides genetic information in space toward elucidation of the spatial architecture in intact organs and the spatially resolved cell-cell communications mediating tissue homeostasis, development, and disease. To facilitate inference of spatially resolved cell-cell communications, we here present SpaTalk, which relies on a graph network and knowledge graph to model and score the ligand-receptor-target signaling network between spatially proximal cells by dissecting cell-type composition through a non-negative linear model and spatial mapping between single-cell transcriptomic and spatially resolved transcriptomic data. The benchmarked performance of SpaTalk on public single-cell spatial transcriptomic datasets is superior to that of existing inference methods. Then we apply SpaTalk to STARmap, Slide-seq, and 10X Visium data, revealing the in-depth communicative mechanisms underlying normal and disease tissues with spatial structure. SpaTalk can uncover spatially resolved cell-cell communications for single-cell and spot-based spatially resolved transcriptomic data universally, providing valuable insights into spatial inter-cellular tissue dynamics.
###### access
https://github.com/ZJUFanLab/SpaTalk

##### Region-specific denoising identifies spatial co-expression patterns and intra-tissue heterogeneity in spatially resolved transcriptomics data(2022 NC)
(https://www.nature.com/articles/s41467-022-34567-0)
a computational tool called Missing-value Imputation and in silico region detection for Spatially resolved Transcriptomics (MIST). MIST detects tissue regions based on their molecular content by maintaining neighboring spots that are both molecularly similar and physically adjacent. Assuming each detected molecular region has a limited number of cell types, MIST then denoises the missing values by approximating a low-rank gene expression matrix through nuclear-norm minimization algorithm
###### access
https://github.com/linhuawang/MIST.git

##### Elucidating tumor heterogeneity from spatially resolved transcriptomics data by multi-view graph collaborative learning(2022 NC)
Here, we propose stMVC, a multi-view graph collaborative-learning model that integrates histology, gene expression, spatial location, and biological contexts in analyzing SRT data by attention. Specifically, stMVC adopting semi-supervised graph attention autoencoder separately learns view-specific representations of histological-similarity-graph or spatial-location-graph, and then simultaneously integrates two-view graphs for robust representations through attention under semi-supervision of biological contexts. stMVC outperforms other tools in detecting tissue structure, inferring trajectory relationships, and denoising on benchmark slices of human cortex. Particularly, stMVC identifies disease-related cell-states and their transition cell-states in breast cancer study, which are further validated by the functional and survival analysis of independent clinical data.
###### access
https://github.com/cmzuo11/stMVC

##### DestVI identifies continuums of cell types in spatial transcriptomics data(2022 NBT)
https://www.nature.com/articles/s41587-022-01272-8
Here we introduce DestVI, a Bayesian model for multi-resolution deconvolution of cell types in ST data. Unlike other methods, DestVI learns both discrete cell-type-specific profiles and continuous sub-cell-type latent variations using a conditional deep generative model.

###### access
https://scvi-tools.org


##### Cell2location maps fine-grained cell types in spatial transcriptomics(2022 NBT)
(https://www.nature.com/articles/s41587-021-01139-4)
Here we present сell2location, a Bayesian model that can resolve fine-grained cell types in spatial transcriptomic data and create comprehensive cellular maps of diverse tissues. Cell2location accounts for technical sources of variation and borrows statistical strength across locations, thereby enabling the integration of single-cell and spatial transcriptomics with higher sensitivity and resolution than existing tools.
###### access
https://github.com/BayraktarLab/cell2location/

##### STRIDE: accurately decomposing and integrating spatial transcriptomics using single-cell RNA sequencing (NAR 2022) 
(https://academic.oup.com/nar/article/50/7/e42/6543547?searchresult=1)
STRIDE a computational method to decompose cell types from spatial mixtures by leveraging topic profiles trained from single-cell transcriptomics. STRIDE accurately estimated the cell-type proportions and showed balanced specificity and sensitivity compared to existing methods.
###### access
https://github.com/wanglabtongji/STRIDE

##### Explainable multiview framework for dissecting spatial relationships from highly multiplexed data(GB 2022) 
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02663-5
MISTy, a flexible, scalable, and explainable machine learning framework for extracting relationships from any spatial omics data, from dozens to thousands of measured markers. MISTy builds multiple views focusing on different spatial or functional contexts to dissect different effects. 
###### access
the source code of mistyR is publicly available from Bioconductor

##### Phiclust: a clusterability measure for single-cell transcriptomics reveals phenotypic subpopulations(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02590-x
 phiclust (ϕclust), a clusterability measure derived from random matrix theory that can be used to identify cell clusters with non-random substructure, testably leading to the discovery of previously overlooked phenotypes.

##### BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02734-7
Spatial transcriptomic studies are reaching single-cell spatial resolution, with data often collected from multiple tissue sections. Here, we present a computational method, BASS, that enables multi-scale and multi-sample analysis for single-cell resolution spatial transcriptomics. BASS performs cell type clustering at the single-cell scale and spatial domain detection at the tissue regional scale, with the two tasks carried out simultaneously within a Bayesian hierarchical modeling framework. 
###### access
https://www.starmapresources.org/data 


### tools
##### Giotto: a toolbox for integrative analysis and visualization of spatial expression data(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2 
Giotto, a comprehensive and open-source toolbox for spatial data analysis and visualization. The analysis module provides end-to-end analysis by implementing a wide range of algorithms for characterizing tissue composition, spatial expression patterns, and cellular interactions.The visualization module allows users to interactively visualize analysis outputs and imaging features.  
###### access
https://rubd.github.io/Giotto_site/

##### Squidpy: a scalable framework for spatial omics analysis(NM 2022)
https://pubmed.ncbi.nlm.nih.gov/35102346/
Here, we present Squidpy, a Python framework that brings together tools from omics and image analysis to enable scalable description of spatial molecular data, such as transcriptome or multivariate proteins. Squidpy provides efficient infrastructure and numerous analysis methods that allow to efficiently store, manipulate and interactively visualize spatial omics data. Squidpy is extensible and can be interfaced with a variety of already existing libraries for the scalable analysis of spatial omics data.
###### access
https://github.com/theislab/squidpy


##### RNALocate v2.0: an updated resource for RNA subcellular localization with increased coverage and annotation(NAR 2022)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8728251/
RNALocate v2.0 is a repository of integrated experimentally validated information on subcellular localization of RNA through manual curation of the literature and five other resources, along with analyses of 35 datasets from the Gene Expression Omnibus (GEO) under a common framework It also supports three RNA subcellular localization prediction tools: DM3Loc, iLoc-lncRNA and iLoc-mRNA. In total, RNALocate v2.0 integrates more than 213 000 RNA subcellular localization entries at 171 locations across 104 species. 
###### access
http://www.rnalocate.org/ 

##### Stitching and registering highly multiplexed whole-slide images of tissues and tumors using ASHLAR(Bioinformatics 2022) 
https://academic.oup.com/bioinformatics/article/38/19/4613/6668278?searchresult=1
We describe the development and testing of ASHLAR, a Python tool for coordinated stitching and registration of 103 or more individual multiplexed images to generate accurate whole-slide mosaics. ASHLAR reads image formats from most commercial microscopes and slide scanners, and we show that it performs better than existing open-source and commercial software. ASHLAR outputs standard OME-TIFF images that are ready for analysis by other open-source tools and recently developed image analysis pipelines.
###### access
https://github.com/labsyspharm/ashlar

##### 3D-Cell-Annotator: an open-source active surface tool for single-cell segmentation in 3D microscopy images(Bioinformatics 2020) 
https://academic.oup.com/bioinformatics/article/36/9/2948/5709038?searchresult=1
Here, we propose 3D-Cell-Annotator to solve this task using 3D active surfaces together with shape descriptors as prior information in a semi-automated fashion. The software uses the convenient 3D interface of the widely used Medical Imaging Interaction Toolkit (MITK). Results on 3D biological structures (e.g. spheroids, organoids and embryos) show that the precision of the segmentation reaches the level of a human expert.
###### access
http://www.3d-cell-annotator.org/

##### cytomapper: an R/Bioconductor package for visualization of highly multiplexed imaging data(Bioinformatics 2020) 
https://academic.oup.com/bioinformatics/article/36/24/5706/6050702?searchresult=1
Here, we describe cytomapper, a computational tool written in R, that enables visualization of pixel- and cell-level information obtained by multiplexed imaging.
###### access
https://www.bioconductor.org/packages/release/bioc/html/cytomapper.html

##### spatialGE: quantification and visualization of the tumor microenvironment heterogeneity using spatial transcriptomics(Bioinformatics 2022) 
https://academic.oup.com/bioinformatics/article/38/9/2645/6544582?searchresult=1
Hence, we have developed spatialGE, a software that provides visualizations and quantification of the tumor microenvironment heterogeneity through gene expression surfaces, spatial heterogeneity statistics that can be compared against clinical information, spot-level cell deconvolution and spatially informed clustering, all using a new data object to store data and resulting analyses simultaneously.
###### access
https://github.com/FridleyLab/spatialGE

##### scFeatures: multi-view representations of single-cell and spatial data for disease outcome prediction(Bioinformatics 2022) 
https://academic.oup.com/bioinformatics/article/38/20/4745/6678979?searchresult=1
Here, we present scFeatures, an approach that creates interpretable cellular and molecular representations of single-cell and spatial data at the sample level
###### access
https://github.com/SydneyBioX/scFeatures

##### PySpacell: A Python Package for Spatial Analysis of Cell Images(2020)
https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23955
We developed a toolbox that allows one to test for the presence of a spatial effect in microscopy images of adherent cells and estimate the spatial scale of this effect. The proposed Python module can be used for any light microscopy images of cells as well as other types of single-cell data such as in situ transcriptomics or metabolomics. The input format of our package matches standard output formats from image analysis tools such as CellProfiler, Fiji, or Icy and thus makes our toolbox easy and straightforward to use, yet offering a powerful statistical approach for a wide range of applications
###### access
https://pypi.org/project/pySpacell/

##### Spacemake: processing and analysis of large-scale spatial transcriptomics data(Gigascience 2022)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9295369/
Here, we present spacemake, a modular, robust, and scalable spatial transcriptomics pipeline built in Snakemake and Python. Spacemake is designed to handle all major spatial transcriptomics data sets and can be readily configured for other technologies. It can process and analyze several samples in parallel, even if they stem from different experimental methods.Spacemake's unified framework enables reproducible data processing from raw sequencing data to automatically generated downstream analysis reports. 
###### access
https://github.com/rajewsky-lab/spacemake

##### Region-specific denoising identifies spatial co-expression patterns and intra-tissue heterogeneity in spatially resolved transcriptomics data(NC 2022)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9663444/
Spatially resolved transcriptomics is a relatively new technique that maps transcriptional information within a tissue. Analysis of these datasets is challenging because gene expression values are highly sparse due to dropout events, and there is a lack of tools to facilitate in silico detection and annotation of regions based on their molecular content. Therefore, we develop a computational tool for detecting molecular regions and region-based Missing value Imputation for Spatially Transcriptomics (MIST).
###### access
https://github.com/linhuawang/MIST.git


## single cell multi-omics tools
### methods
##### Integrated analysis of multimodal single-cell data with structural similarity(NAR 2022) 
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac781/6709246?searchresult=1
SAILERX is a deep learning framework, for efficient, robust, and flexible analysis of multi-modal single-cell data.Distinct from existing methods, SAILERX can handle both parallel scRNA-seq and scATAC-seq multiome data, single modal scATAC-seq data, and a hybrid of these two types of data.
###### access
https://github.com/uci-cbcl/SAILERX

##### PyLiger: scalable single-cell multi-omic data integration in Python(Bioinformatics 2022) 
https://academic.oup.com/bioinformatics/article/38/10/2946/6561542?searchresult=1
PyLiger, a Python package for integrating single-cell multi-omic datasets. PyLiger offers faster performance than the previous R implementation (2–5× speedup), interoperability with AnnData format, flexible on-disk or in-memory analysis capability and new functionality for gene ontology enrichment analysis. The on-disk capability enables analysis of arbitrarily large single-cell datasets using fixed memory.
###### access
https://github.com/welch-lab/pyliger

### tools
##### webSCST: an interactive web application for single-cell RNA-sequencing data and spatial transcriptomic data integration(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/article/38/13/3488/6590646?searchresult=1
webSCST is a freely available spatial transcriptomic database and web service that aims to conveniently allow users to obtain predicted spatial information using scRNA-seq data. Users with raw scRNA-seq data can utilize webSCST by quality-controlling the scRNA-seq data, matching with spatial datasets in the webSCST database and finally obtaining the integrated results through four popular integrating tools. In the future, we plan to add more publicly available published spatial sequencing data to our database and add newly developed integration tools to maintain webSCST as an up-to-date resource.
###### access
webSCST browsers:http://www.webscst.com
webSCST R package:https://github.com/swsoyee/webSCST

##### High-performance single-cell gene regulatory network inference at scale: the Inferelator 3.0(Bioinformatics 2022)
In this work, we present the Inferelator 3.0, which has been significantly updated to integrate data from distinct cell types to learn context-specific regulatory networks and aggregate them into a shared regulatory network, while retaining the functionality of the previous versions. The Inferelator is able to integrate the largest single-cell datasets and learn cell-type-specific gene regulatory networks.
###### access
https://github.com/flatironinstitute/inferelator



### database
##### DISCO: a database of Deeply Integrated human Single-Cell Omics data(NAR 2022) 
(https://academic.oup.com/nar/article/50/D1/D596/6430491?searchresult=1)
DISCO , a database of Deeply Integrated Single-Cell Omics data. The current release of DISCO integrates more than 18 million cells from 4593 samples, covering 107 tissues/cell lines/organoids, 158 diseases, and 20 platforms.DISCO provides three online tools, namely Online FastIntegration, Online CELLiD, and CellMapper, for users to integrate, annotate, and project uploaded single-cell RNA-seq data onto a selected atlas. Collectively, DISCO is a versatile platform for users to explore published single-cell data and efficiently perform integrated analysis with their own data.
[FastIntegration, a fast and high-capacity version of Seurat Integration.CELLiD, an atlas guided automatic cell type identification tool.]
###### access
https://www.immunesinglecell.org/

##### EDomics: a comprehensive and comparative multi-omics database for animal evo-devo(NAR 2022) 
(https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac944/6786194#378500188）
EDomics a comparative multi-omics database for animal evo-devo containing comprehensive genomes, bulk transcriptomes, and single-cell data across 40 representative species, many of which are generally used as model organisms for animal evo-devo study.
[EDomics provides a systematic view of genomic/transcriptomic information from various aspects, including genome assembly statistics, gene features and families, transcription factors, transposable elements, and gene expressional profiles/networks. It also exhibits spatiotemporal gene expression profiles at a single-cell level, such as cell atlas, cell markers, and spatial-map information. Moreover, EDomics provides highly valuable, customized datasets/resources for evo-devo research, including gene family expansion/contraction, inferred core gene repertoires, macrosynteny analysis for karyotype evolution, and cell type evolution analysis.]
###### access
http://edomics.qnlm.ac

##### BIOMEX: an interactive workflow for (single cell) omics data interpretation and visualization(NAR 2020)
(https://academic.oup.com/nar/article/48/W1/W385/5835814?searchresult=1)
BIOMEX, a browser-based software. BIOMEX integrates state-of-the-art statistical tools and field-tested algorithms into a flexible but well-defined workflow that accommodates metabolomics, transcriptomics, proteomics, mass cytometry and single cell data from different platforms and organisms.BIOMEX guides the user through omics-tailored analyses, such as data pretreatment and normalization, dimensionality reduction, differential and enrichment analysis, pathway mapping, clustering, marker analysis, trajectory inference, meta-analysis and others.
###### access
https://carmelietlab.sites.vib.be/en/biomex


## other tools
##### Gene set enrichment analysis for genome-wide DNA methylation data(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02388-x
GOmeth and GOregion are new methods for performing unbiased gene set testing following differential methylation analysis. 
###### access
Both methods are publicly available in the missMethyl Bioconductor R package.
The GitHub repository associated with the analysis website is at: https://github.com/Oshlack/methyl-geneset-testing.

##### Integrative analysis of 3604 GWAS reveals multiple novel cell type-specific regulatory associations(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02560-3
To gain insight into potential functional mechanisms underlying GWAS associations, we developed FORGE2, which is an updated version of the FORGE web tool. FORGE2 uses an expanded atlas of cell type-specific regulatory element annotations, including DNase I hotspots, five histone mark categories and 15 hidden Markov model (HMM) chromatin states, to identify tissue- and cell type-specific signals.
###### access
https://forge2.altiusinstitute.org/

##### scMAGeCK links genotypes with multiple phenotypes in single-cell CRISPR screens(GB 2020)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1928-4
scMAGeCK, a computational framework to identify genomic elements associated with multiple expression-based phenotypes in CRISPR/Cas9 functional screening that uses single-cell RNA-seq as readout. scMAGeCK outperforms existing methods, identifies genes and enhancers with known and novel functions in cell proliferation, and enables an unbiased construction of genotype-phenotype network.Collectively, scMAGeCK is a novel tool to study genotype-phenotype relationships at a single-cell level.
###### access

##### Single-cell diploid Hi-C reveals the role of spatial aggregations in complex rearrangements and KMT2A fusions in leukemia(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02740-9
By analyzing published single-cell diploid Hi-C datasets, we find partner genes fused in leukemia exhibit smaller spatial distances than those fused in solid tumor and control gene pairs. Intriguingly, multiple partners tend to colocalize with KMT2A in the same cell. 

##### scTAM-seq enables targeted high-confidence analysis of DNA methylation in single cells(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02796-7
Single-cell DNA methylation profiling currently suffers from excessive noise and/or limited cellular throughput. We developed scTAM-seq, a targeted bisulfite-free method for profiling up to 650 CpGs in up to 10,000 cells per experiment, with a dropout rate as low as 7%. 

##### GSEApy: a comprehensive package for performing gene set enrichment analysis in Python(Bioinformatics 2022)
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac757/6847088?searchresult=1
We present a package (GSEApy) that performs GSEA in either the command line or Python environment. GSEApy uses a Rust implementation to enable it to calculate the same enrichment statistic as GSEA for a collection of pathways. 
###### access
https://github.com/zqfang/GSEApy

##### HiCRep.py: fast comparison of Hi-C contact matrices in Python(Bioinformatics 2021)
https://academic.oup.com/bioinformatics/article/37/18/2996/6133255?searchresult=1
We introduce a Python implementation of HiCRep and demonstrate that it is much faster and consumes much less memory than the existing R implementation. Furthermore, we give examples of HiCRep’s ability to accurately distinguish replicates from non-replicates and to reveal cell type structure among collections of Hi-C data.
###### access
https://github.com/Noble-Lab/hicrep

##### PaintSHOP enables the interactive design of transcriptome- and genome-scale oligonucleotide FISH experiments(NM 2022)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8349872/
Here, we introduce Paint Server and Homology Optimization Pipeline (PaintSHOP), an interactive platform for the design of oligo FISH experiments. PaintSHOP enables researchers to identify probes for their experimental targets efficiently, to incorporate additional necessary sequences such as primer pairs, and to easily generate files documenting library design. PaintSHOP democratizes and standardizes the process of designing complex probe sets for the oligo FISH community.
###### access
The source code for the PaintSHOP web application is available as Supplementary Software 1 and at https://github.com/beliveau-lab/PaintSHOP. The source code for the Homology Optimization Pipeline is available as Supplementary Software 2 and at https://github.com/beliveau-lab/PaintSHOP_pipeline.
























