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

### tools
##### Giotto: a toolbox for integrative analysis and visualization of spatial expression data(GB 2022)
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2 
Giotto, a comprehensive and open-source toolbox for spatial data analysis and visualization. The analysis module provides end-to-end analysis by implementing a wide range of algorithms for characterizing tissue composition, spatial expression patterns, and cellular interactions.The visualization module allows users to interactively visualize analysis outputs and imaging features.  
###### access
 https://rubd.github.io/Giotto_site/
 
 

## single cell multi-omics tools
### methods
##### Integrated analysis of multimodal single-cell data with structural similarity(NAR 2022) 
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkac781/6709246?searchresult=1
SAILERX is a deep learning framework, for efficient, robust, and flexible analysis of multi-modal single-cell data.Distinct from existing methods, SAILERX can handle both parallel scRNA-seq and scATAC-seq multiome data, single modal scATAC-seq data, and a hybrid of these two types of data.
###### access
https://github.com/uci-cbcl/SAILERX


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














