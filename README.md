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

### data integration methods
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
 
 

## single-cell epigenomic tools
##### EpiScanpy: integrated single-cell epigenomic analysis(NC 2021)
(https://www.nature.com/articles/s41467-021-25131-3)
EpiScanpy is a toolkit for the analysis of single-cell epigenomic data, namely single-cell DNA methylation and single-cell ATAC-seq data. To address the modality specific challenges from epigenomics data, epiScanpy quantifies the epigenome using multiple feature space constructions and builds a nearest neighbour graph using epigenomic distance between cells. EpiScanpy makes the many existing scRNA-seq workflows from scanpy available to large-scale single-cell data from other -omics modalities, including methods for common clustering, dimension reduction, cell type identification and trajectory learning techniques, as well as an atlas integration tool for scATAC-seq datasets. The toolkit also features numerous useful downstream functions, such as differential methylation and differential openness calling, mapping epigenomic features of interest to their nearest gene, or constructing gene activity matrices using chromatin openness. We successfully benchmark epiScanpy against other scATAC-seq analysis tools and show its outperformance at discriminating cell types.

###### access
https://github.com/colomemaria/epiScanpy, https://colomemaria.github.io/episcanpy_doc


##### MIRA: joint regulatory modeling of multimodal expression and chromatin accessibility in single cells(NM 2022)
(https://www.nature.com/articles/s41592-022-01595-z)
Rigorously comparing gene expression and chromatin accessibility in the same single cells could illuminate the logic of how coupling or decoupling of these mechanisms regulates fate commitment. Here we present MIRA, probabilistic multimodal models for integrated regulatory analysis, a comprehensive methodology that systematically contrasts transcription and accessibility to infer the regulatory circuitry driving cells along cell state trajectories. MIRA leverages topic modeling of cell states and regulatory potential modeling of individual gene loci. MIRA thereby represents cell states in an efficient and interpretable latent space, infers high-fidelity cell state trees, determines key regulators of fate decisions at branch points and exposes the variable influence of local accessibility on transcription at distinct loci. Applied to epidermal differentiation and embryonic brain development from two different multimodal platforms, MIRA revealed that early developmental genes were tightly regulated by local chromatin landscape whereas terminal fate genes were titrated without requiring extensive chromatin remodeling.
###### access
MIRA is available as a Python package at https://github.com/cistrome/MIRA. Frankencell, a Python program we developed to generate synthetic differentiation trajectories for benchmarking, is available at https://github.com/AllenWLynch/frankencell-dynverse.


## ST tools
##### Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data with SpaTalk(2022 NC)
(https://www.nature.com/articles/s41467-022-32111-8)
Spatially resolved transcriptomics provides genetic information in space toward elucidation of the spatial architecture in intact organs and the spatially resolved cell-cell communications mediating tissue homeostasis, development, and disease. To facilitate inference of spatially resolved cell-cell communications, we here present SpaTalk, which relies on a graph network and knowledge graph to model and score the ligand-receptor-target signaling network between spatially proximal cells by dissecting cell-type composition through a non-negative linear model and spatial mapping between single-cell transcriptomic and spatially resolved transcriptomic data. The benchmarked performance of SpaTalk on public single-cell spatial transcriptomic datasets is superior to that of existing inference methods. Then we apply SpaTalk to STARmap, Slide-seq, and 10X Visium data, revealing the in-depth communicative mechanisms underlying normal and disease tissues with spatial structure. SpaTalk can uncover spatially resolved cell-cell communications for single-cell and spot-based spatially resolved transcriptomic data universally, providing valuable insights into spatial inter-cellular tissue dynamics.
###### access
https://github.com/ZJUFanLab/SpaTalk



















