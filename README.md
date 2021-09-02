This repository contains the scripts used for the following study:

## Reconstruction of the Mouse Olfactory Glomerular Map using Single-Cell Technologies

### I-Hao Wang<sup>1</sup>, Evan Murray<sup>2,10</sup>, Greg Andrews <sup>3,10</sup>, Hao-Ching Jiang<sup>1,10</sup>, Sung Jin Park<sup>1,10</sup>, Elisa Donnard<sup>3,10</sup>, Violeta Durán-Laforet<sup>4</sup>, Daniel M. Bear<sup>5,6</sup>, Travis E. Faust<sup>4</sup>, Manuel Garber<sup>1,3</sup>, Christina E. Baer<sup>7</sup>, Dorothy P. Schafer<sup>4</sup>,Zhiping Weng<sup>3</sup>, Fei Chen<sup>2,8</sup>, Evan Z. Macosko<sup>2,9</sup>, Paul L. Greer<sup>1,*</sup>

<sup>1</sup>Program in Molecular Medicine, University of Massachusetts Medical School, Worcester, MA, USA, 
<sup>2</sup>Broad Institute of Harvard and MIT, Cambridge, MA, USA, 
<sup>3</sup>Program in Bioinformatics and Integrative Biology, University of Massachusetts Medical School, Worcester, MA, USA, 
<sup>4</sup>Department of Neurobiology and Brudnick Neuropsychiatric Research Institute, University of Massachusetts Medical School, Worcester, MA, USA
<sup>5</sup>Department of Psychology, Stanford University, Palo Alto, CA, USA 
<sup>6</sup>Wu Tsai Neurosciences Institute, Stanford University, Palo Alto, CA, USA 
<sup>7</sup>Sanderson Center for Optical Imaging and Department of Microbiology and Physiological Systems, University of Massachusetts Medical School, Worcester, MA, USA
<sup>8</sup>Department of Stem Cell and Regenerative Biology, Harvard University, Cambridge, MA, USA, 
<sup>9</sup>Department of Psychiatry, Massachusetts General Hospital, Boston, MA, USA, 
<sup>10</sup>These authors contributed equally, 
<sup>*</sup>Corresponding author.


The raw data can be downloaded at [GSE169021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169021) (data will be available after publication)

Note:
1. the scripts start with "bk" refers to the bulk RNAseq analysis
2. the scripts start with "sc" refers to the single cell RNAseq analysis
3. the scripts start with "ss" refers to the Slide-seq V2 analysis
4. the scripts start with "map" refers to the reconstructed map analysis
5. for the Slide-seq V2 analysis, there were two sets of data used in the this study and one copy of the script provided since the same data processing procedue was performed
6. `Figures.R` contains all the codes for generating figures in this study
7. `data/map_pre.rds` contains all the information of reconstructed map
8. Due to the size, some of the intermediate files are not included in this repository, but it can be generated from the raw data provided in GEO
9. Contact information: I-Hao Wang (I-Hao.Wang@umassmed.edu), Dr. Paul Greer (Paul.Greer@umassmed.edu)





```R
sessionInfo()
R version 4.0.4 (2021-02-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 10.16

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] EnhancedVolcano_1.8.0       ggtree_3.1.2                ggnewscale_0.4.5            ape_5.5                    
 [5] ComplexHeatmap_2.6.2        monocle3_0.2.3.0            SingleCellExperiment_1.12.0 SeuratWrappers_0.3.0       
 [9] GO.db_3.12.1                ggtext_0.1.1                rstatix_0.7.0               gridExtra_2.3              
[13] scales_1.1.1                ggrepel_0.9.1.9999          reshape2_1.4.4              pheatmap_1.0.12            
[17] RColorBrewer_1.1-2          BioVenn_1.1.3               DESeq2_1.30.1               SummarizedExperiment_1.20.0
[21] MatrixGenerics_1.2.1        matrixStats_0.60.0          GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
[25] openxlsx_4.2.4              reticulate_1.20             harmony_1.0                 Rcpp_1.0.7                 
[29] org.Mm.eg.db_3.12.0         AnnotationDbi_1.52.0        IRanges_2.24.1              S4Vectors_0.28.1           
[33] Biobase_2.50.0              BiocGenerics_0.36.1         clusterProfiler_3.18.1      biomaRt_2.46.3             
[37] SeuratObject_4.0.2          Seurat_4.0.3                viridis_0.6.1               viridisLite_0.4.0          
[41] ggpubr_0.4.0                ggplot2_3.3.5               dplyr_1.0.7                 patchwork_1.1.1            
[45] jpeg_0.1-8.1               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3         scattermore_0.7        R.methodsS3_1.8.1      tidyr_1.1.3            bit64_4.0.5           
  [6] knitr_1.33             R.utils_2.10.1         irlba_2.3.3            DelayedArray_0.16.3    data.table_1.14.0     
 [11] rpart_4.1-15           RCurl_1.98-1.3         generics_0.1.0         cowplot_1.1.1          RSQLite_2.2.7         
 [16] shadowtext_0.0.8       RANN_2.6.1             future_1.21.0          bit_4.0.4              enrichplot_1.10.2     
 [21] spatstat.data_2.1-0    xml2_1.3.2             httpuv_1.6.1           assertthat_0.2.1       xfun_0.24             
 [26] hms_1.1.0              promises_1.2.0.1       fansi_0.5.0            progress_1.2.2         dbplyr_2.1.1          
 [31] readxl_1.3.1           igraph_1.2.6           DBI_1.1.1              geneplotter_1.68.0     htmlwidgets_1.5.3     
 [36] spatstat.geom_2.2-2    purrr_0.3.4            ellipsis_0.3.2         backports_1.2.1        annotate_1.68.0       
 [41] deldir_0.2-10          vctrs_0.3.8            Cairo_1.5-12.2         remotes_2.4.0          ROCR_1.0-11           
 [46] abind_1.4-5            cachem_1.0.5           withr_2.4.2            ggforce_0.3.3          checkmate_2.0.0       
 [51] sctransform_0.3.2      treeio_1.17.1.992      prettyunits_1.1.1      goftest_1.2-2          svglite_2.0.0         
 [56] cluster_2.1.2          DOSE_3.16.0            lazyeval_0.2.2         crayon_1.4.1           genefilter_1.72.1     
 [61] pkgconfig_2.0.3        labeling_0.4.2         tweenr_1.0.2           vipor_0.4.5            nlme_3.1-152          
 [66] nnet_7.3-16            rlang_0.4.11           globals_0.14.0         lifecycle_1.0.0        miniUI_0.1.1.1        
 [71] downloader_0.4         extrafontdb_1.0        BiocFileCache_1.14.0   rsvd_1.0.5             ggrastr_0.2.3         
 [76] cellranger_1.1.0       polyclip_1.10-0        lmtest_0.9-38          aplot_0.0.6            Matrix_1.3-4          
 [81] carData_3.0-4          zoo_1.8-9              beeswarm_0.4.0         base64enc_0.1-3        GlobalOptions_0.1.2   
 [86] ggridges_0.5.3         rjson_0.2.20           png_0.1-7              bitops_1.0-7           R.oo_1.24.0           
 [91] KernSmooth_2.23-20     blob_1.2.1             shape_1.4.6            stringr_1.4.0          qvalue_2.22.0         
 [96] parallelly_1.27.0      ggsignif_0.6.2         memoise_2.0.0          magrittr_2.0.1         plyr_1.8.6            
[101] ica_1.0-2              zlibbioc_1.36.0        compiler_4.0.4         scatterpie_0.1.6       ash_1.0-15            
[106] clue_0.3-59            plotrix_3.8-1          fitdistrplus_1.1-5     XVector_0.30.0         listenv_0.8.0         
[111] pbapply_1.4-3          htmlTable_2.2.1        Formula_1.2-4          MASS_7.3-54            mgcv_1.8-36           
[116] tidyselect_1.1.1       stringi_1.7.3          forcats_0.5.1          proj4_1.0-10.1         GOSemSim_2.16.1       
[121] askpass_1.1            locfit_1.5-9.4         latticeExtra_0.6-29    fastmatch_1.1-0        tools_4.0.4           
[126] future.apply_1.7.0     rio_0.5.27             circlize_0.4.13        rstudioapi_0.13        foreign_0.8-81        
[131] farver_2.1.0           Rtsne_0.15             ggraph_2.0.5           digest_0.6.27          rvcheck_0.1.8         
[136] BiocManager_1.30.16    shiny_1.6.0            gridtext_0.1.4         car_3.0-10             broom_0.7.7           
[141] ggalt_0.4.0            later_1.2.0            RcppAnnoy_0.0.19       httr_1.4.2             colorspace_2.0-2      
[146] XML_3.99-0.6           tensor_1.5             splines_4.0.4          uwot_0.1.10            tidytree_0.3.4        
[151] spatstat.utils_2.2-0   graphlayouts_0.7.1     plotly_4.9.4.1         systemfonts_1.0.2      xtable_1.8-4          
[156] jsonlite_1.7.2         tidygraph_1.2.0        R6_2.5.1               Hmisc_4.5-0            pillar_1.6.2          
[161] htmltools_0.5.1.1      mime_0.11              glue_1.4.2             fastmap_1.1.0          BiocParallel_1.24.1   
[166] codetools_0.2-18       maps_3.3.0             fgsea_1.16.0           utf8_1.2.2             lattice_0.20-44       
[171] spatstat.sparse_2.0-0  tibble_3.1.3           ggbeeswarm_0.6.0       curl_4.3.2             leiden_0.3.9          
[176] Rttf2pt1_1.3.8         zip_2.2.0              openssl_1.4.4          survival_3.2-11        munsell_0.5.0         
[181] GetoptLong_1.0.5       DO.db_2.9              GenomeInfoDbData_1.2.4 haven_2.4.1            gtable_0.3.0          
[186] extrafont_0.17         spatstat.core_2.3-0
```
