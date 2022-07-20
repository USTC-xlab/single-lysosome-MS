## Single-lysosome mass spectrometry (SLMS)
&emsp; Lysosomes are critical for cellular metabolism and heterogeneously involved in various cellular
processes such as endocytosis, autophagy and senescence. The ability to measure lysosomal 
metabolic heterogeneity is essential for understanding its physiological roles. We therefore 
built a single-lysosome mass spectrometry (SLMS) platform integrating lysosomal patch-clamp 
recording and induced nanoESI/MS that enabled concurrent metabolic and electrophysiological 
profiling of individual enlarged lysosomes.<br/>

<div align = center> 
<img src="https://user-images.githubusercontent.com/62707821/114822069-197aee80-9df4-11eb-8b4a-81fc353a10e6.png" width = "400" height = "300" />
</div>

## Main R packages

* **xcms** (version 3.9.1)
  * Framework for processing and visualization of single-spectra mass spectral data.  
* **Seurat** (version 3.2.0)
  * Identify the main lysosome types and visualize the clustering information in 2D space by t-SNE.
* **monocle** (version 2.16.0)
  * Analyze lysosome trajectories by sorting lysosomes in pseudo-time trajectories.
* **FactoMineR** (version 2.3), **factoextra** (version 1.0.7)
  * Perform principal components analysis (PCA) to check if there is a batch effect in data. 
* **pheatmap** (version 1.0.12)
  * Plot heatmaps and perform hierarchical cluster analysis (HCA) to check if there is a batch effect in data.


## Description of files in this repository
###  Scripts
1. **find peaks.r**
   * Preprocess the original data and get matrices with row names (m/z) and column names (sample names) after preprocessing.
2. **statistical analysis.R**
   * Normalize the intensity of each lysosome, check the batch effect and differential analysis.
3. **seurat.R**
   * Identify the main lysosome types and visualize the clustering information in 2D space by t-SNE.
4. **monocle.R**
   * Analyze lysosome trajectories by sorting lysosomes in pseudo-time trajectories.
5. **plot_peaks.R**
   * Map the average mass spectra for each subpopulation of lysosomes.


## Data repository
&emsp; The mass spectrometric raw data of single lysosomes have been deposited to the **MassIVE 
database** under accession code **[MSV000087208](https://massive.ucsd.edu/ProteoSAFe/QueryMSV?id=MSV000087208)**
and can be visualized with Thermo Xcalibur Qual Browser.


## Cite
**[Zhu, H., Li, Q., Liao, T. et al. Metabolomic profiling of single enlarged lysosomes. Nat Methods 18, 788–798 (2021).](https://www.nature.com/articles/s41592-021-01182-8)**
