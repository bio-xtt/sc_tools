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
Online single-cell data integration through projecting heterogeneous datasets into a common cell-embedding space(NC 2022)
https://www.nature.com/articles/s41467-022-33758-z
 SCALEX, a deep-learning method that integrates single-cell data by projecting cells into a batch-invariant, common cell-embedding space in a truly online manner. SCALEX substantially outperforms online iNMF and other state-of-the-art non-online integration methods on benchmark single-cell datasets of diverse modalities.The online data integration capacity and superior performance makes SCALEX particularly appropriate for large-scale single-cell applications to build upon previous scientific insights.
###### access
 https://github.com/jsxlei/SCALEX




























