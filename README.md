This repository contains the scripts used for the following study:

## Reconstruction of the Mouse Olfactory Glomerular Map using Single-Cell Technologies

### I-Hao Wang<sup>1</sup>, Evan Murray<sup>2,8</sup>, Greg Andrews <sup>3,8</sup>, Elisa Donnard<sup>3,8</sup>, Hao-Ching Jiang<sup>1</sup>, Sung Jin Park<sup>1</sup>, Daniel M. Bear<sup>4,5</sup>, Manuel Garber<sup>1,3</sup>, Zhiping Weng<sup>3</sup>, Fei Chen<sup>2,6</sup>, Evan Z. Macosko<sup>2,7</sup>, Paul L. Greer<sup>1,*</sup>

<sup>1</sup>Program in Molecular Medicine, University of Massachusetts Medical School, Worcester, MA, USA, 
<sup>2</sup>Broad Institute of Harvard and MIT, Cambridge, MA, USA, 
<sup>3</sup>Program in Bioinformatics and Integrative Biology, University of Massachusetts Medical School, Worcester, MA, USA, 
<sup>4</sup>Department of Psychology, Stanford University, 
<sup>5</sup>Wu Tsai Neurosciences Institute, Stanford University, 
<sup>6</sup>Department of Stem Cell and Regenerative Biology, Harvard University, Cambridge, MA, USA, 
<sup>7</sup>Department of Psychiatry, Massachusetts General Hospital, Boston, MA, USA, 
<sup>8</sup>These authors contributed equally, 
<sup>*</sup>Corresponding author.
  

Note:
1. the scripts start with "bk" refers to the bulk RNAseq analysis
2. the scripts start with "sc" refers to the single cell RNAseq analysis
3. the scripts start with "ss" refers to the Slide-seq V2 analysis
4. the scripts start with "map" refers to the reconstructed map analysis
5. for the Slide-seq V2 analysis, there were two sets of data used in the this study and one copy of the script provided since the same data processing procedue was performed.
6. `Figures.R` contains all the codes for generating figures in this study





```R
sessionInfo()
R version 4.0.4 (2021-02-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reticulate_1.18             spdep_1.1-5                 spData_0.3.8                sp_1.4-5                   
 [5] MLmetrics_1.1.1             jpeg_0.1-8.1                oce_1.3-0                   sf_0.9-7                   
 [9] gsw_1.0-5                   testthat_3.0.2              circlize_0.4.12             ComplexHeatmap_2.6.2       
[13] viridis_0.5.1               viridisLite_0.3.0           forcats_0.5.1               Hmisc_4.5-0                
[17] Formula_1.2-4               survival_3.2-7              lattice_0.20-41             enrichplot_1.10.2          
[21] EnhancedVolcano_1.8.0       ggrepel_0.9.1               rstatix_0.7.0               ggpubr_0.4.0               
[25] gridExtra_2.3               ggtext_0.1.1                scales_1.1.1                dplyr_1.0.5                
[29] harmony_1.0                 Rcpp_1.0.6                  SeuratObject_4.0.0          Seurat_4.0.0               
[33] protr_1.6-2                 seqinr_4.2-5                msa_1.22.0                  Biostrings_2.58.0          
[37] XVector_0.30.0              reshape2_1.4.4              patchwork_1.1.1             ggplot2_3.3.3              
[41] clusterProfiler_3.18.1      GO.db_3.12.1                org.Mm.eg.db_3.12.0         AnnotationDbi_1.52.0       
[45] biomaRt_2.46.3              BioVenn_1.1.1               pheatmap_1.0.12             RColorBrewer_1.1-2         
[49] DESeq2_1.30.1               SummarizedExperiment_1.20.0 Biobase_2.50.0              MatrixGenerics_1.2.1       
[53] matrixStats_0.58.0          GenomicRanges_1.42.0        GenomeInfoDb_1.26.2         IRanges_2.24.1             
[57] S4Vectors_0.28.1            BiocGenerics_0.36.0        

loaded via a namespace (and not attached):
  [1] ica_1.0-2              svglite_2.0.0          class_7.3-18           foreach_1.5.1          lmtest_0.9-38         
  [6] crayon_1.4.1           MASS_7.3-53.1          nlme_3.1-152           backports_1.2.1        GOSemSim_2.16.1       
 [11] rlang_0.4.10           ROCR_1.0-11            readxl_1.3.1           irlba_2.3.3            extrafontdb_1.0       
 [16] extrafont_0.17         BiocParallel_1.24.1    rjson_0.2.20           bit64_4.0.5            glue_1.4.2            
 [21] sctransform_0.3.2      vipor_0.4.5            classInt_0.4-3         DOSE_3.16.0            haven_2.3.1           
 [26] tidyselect_1.1.0       rio_0.5.26             fitdistrplus_1.1-3     XML_3.99-0.5           tidyr_1.1.3           
 [31] zoo_1.8-8              proj4_1.0-10.1         xtable_1.8-4           magrittr_2.0.1         cli_2.3.1             
 [36] zlibbioc_1.36.0        rstudioapi_0.13        miniUI_0.1.1.1         rpart_4.1-15           fastmatch_1.1-0       
 [41] maps_3.3.0             shiny_1.6.0            xfun_0.21              askpass_1.1            clue_0.3-58           
 [46] cluster_2.1.1          tidygraph_1.2.0        expm_0.999-6           tibble_3.1.0           listenv_0.8.0         
 [51] png_0.1-7              future_1.21.0          withr_2.4.1            bitops_1.0-6           ggforce_0.3.3         
 [56] plyr_1.8.6             cellranger_1.1.0       e1071_1.7-4            coda_0.19-4            pillar_1.5.1          
 [61] GlobalOptions_0.1.2    cachem_1.0.4           raster_3.4-5           GetoptLong_1.0.5       gmodels_2.18.1        
 [66] vctrs_0.3.6            ellipsis_0.3.1         generics_0.1.0         tools_4.0.4            foreign_0.8-81        
 [71] beeswarm_0.3.1         munsell_0.5.0          tweenr_1.0.1           fgsea_1.16.0           DelayedArray_0.16.2   
 [76] fastmap_1.1.0          compiler_4.0.4         abind_1.4-5            httpuv_1.5.5           plotly_4.9.3          
 [81] GenomeInfoDbData_1.2.4 deldir_0.2-10          utf8_1.1.4             later_1.1.0.1          BiocFileCache_1.14.0  
 [86] jsonlite_1.7.2         pbapply_1.4-3          carData_3.0-4          genefilter_1.72.1      lazyeval_0.2.2        
 [91] LearnBayes_2.15.1      promises_1.2.0.1       spatstat_1.64-1        car_3.0-10             doParallel_1.0.16     
 [96] latticeExtra_0.6-29    goftest_1.2-2          spatstat.utils_2.0-0   checkmate_2.0.0        openxlsx_4.2.3        
[101] ash_1.0-15             cowplot_1.1.1          Rtsne_0.15             downloader_0.4         uwot_0.1.10           
[106] igraph_1.2.6           plotrix_3.8-1          systemfonts_1.0.1      htmltools_0.5.1.1      memoise_2.0.0         
[111] locfit_1.5-9.4         graphlayouts_0.7.1     digest_0.6.27          assertthat_0.2.1       mime_0.10             
[116] rappdirs_0.3.3         Rttf2pt1_1.3.8         units_0.7-0            RSQLite_2.2.3          future.apply_1.7.0    
[121] data.table_1.14.0      blob_1.2.1             splines_4.0.4          labeling_0.4.2         Cairo_1.5-12.2        
[126] gridtext_0.1.4         RCurl_1.98-1.2         broom_0.7.5            hms_1.0.0              colorspace_2.0-0      
[131] base64enc_0.1-3        BiocManager_1.30.10    ggbeeswarm_0.6.0       shape_1.4.5            ggrastr_0.2.3         
[136] nnet_7.3-15            RANN_2.6.1             fansi_0.4.2            parallelly_1.23.0      R6_2.5.0              
[141] ggridges_0.5.3         lifecycle_1.0.0        zip_2.1.1              curl_4.3               ggsignif_0.6.1        
[146] gdata_2.18.0           leiden_0.3.7           DO.db_2.9              Matrix_1.3-2           qvalue_2.22.0         
[151] RcppAnnoy_0.0.18       iterators_1.0.13       stringr_1.4.0          htmlwidgets_1.5.3      polyclip_1.10-0       
[156] purrr_0.3.4            shadowtext_0.0.7       mgcv_1.8-34            globals_0.14.0         openssl_1.4.3         
[161] htmlTable_2.1.0        codetools_0.2-18       gtools_3.8.2           prettyunits_1.1.1      dbplyr_2.1.0          
[166] gtable_0.3.0           DBI_1.1.1              tensor_1.5             httr_1.4.2             KernSmooth_2.23-18    
[171] stringi_1.5.3          progress_1.2.2         farver_2.1.0           annotate_1.68.0        xml2_1.3.2            
[176] rvcheck_0.1.8          boot_1.3-27            ggalt_0.4.0            ade4_1.7-16            geneplotter_1.68.0    
[181] scattermore_0.7        bit_4.0.4              scatterpie_0.1.5       spatstat.data_2.0-0    ggraph_2.0.5          
[186] pkgconfig_2.0.3        knitr_1.31
```
