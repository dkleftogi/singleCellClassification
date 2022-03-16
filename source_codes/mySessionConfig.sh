R version 4.1.3 (2022-03-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.6.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base 

other attached packages:
 [1] minpack.lm_1.2-1            ggExtra_0.9                 matlab_1.0.2                pmpp_0.1.1                 
 [5] pracma_2.3.3                ggbiplot_0.55               scales_1.1.1                dendextend_1.15.2          
 [9] corrplot_0.89               caret_6.0-88                lattice_0.20-45             ggridges_0.5.3             
[13] flowDensity_1.26.0          pheatmap_1.0.12             Rphenograph_0.99.1          ConsensusClusterPlus_1.56.0
[17] ggrepel_0.9.1               limma_3.48.3                reshape2_1.4.4              matrixStats_0.59.0         
[21] ggraph_2.0.5                tidygraph_1.2.0             RColorBrewer_1.1-2          readr_1.4.0                
[25] survminer_0.4.9             ggpubr_0.4.0                glmnet_4.1-3                Matrix_1.4-0               
[29] survival_3.2-13             gridExtra_2.3               ggplot2_3.3.5               Rtsne_0.15                 
[33] data.table_1.14.0           FlowSOM_2.0.0               igraph_1.2.6                flowCore_2.4.0             
[37] tidyr_1.1.3                 magrittr_2.0.1              dplyr_1.0.7                 plyr_1.8.6                 
[41] readxl_1.3.1  

loaded via a namespace (and not attached):
  [1] utf8_1.2.1           tidyselect_1.1.1     miscTools_0.6-26     pROC_1.17.0.1        aws.signature_0.6.0 
  [6] munsell_0.5.0        codetools_0.2-18     miniUI_0.1.1.1       Rwave_2.4-8          withr_2.4.2         
 [11] colorspace_2.0-1     flowViz_1.56.0       Biobase_2.52.0       knitr_1.33           stats4_4.1.3        
 [16] ggsignif_0.6.3       Rdpack_2.1.3         splancs_2.01-42      RSEIS_3.9-3          KMsurv_0.1-5        
 [21] polyclip_1.10-0      farver_2.1.0         flowWorkspace_4.4.0  vctrs_0.3.8          generics_0.1.0      
 [26] ipred_0.9-11         xfun_0.24            R6_2.5.0             graphlayouts_0.7.2   fields_12.3         
 [31] bitops_1.0-7         promises_1.2.0.1     nnet_7.3-17          RFOC_3.4-6           rgeos_0.5-5         
 [36] gtable_0.3.0         spam_2.6-0           sandwich_3.0-1       RProtoBufLib_2.4.0   timeDate_3043.102   
 [41] rlang_0.4.11         splines_4.1.3        rstatix_0.7.0        ModelMetrics_1.2.2.2 hexbin_1.28.2       
 [46] broom_0.7.7          yaml_2.2.1           abind_1.4-5          backports_1.2.1      httpuv_1.6.5        
 [51] IDPmisc_1.1.20       RBGL_1.68.0          tools_4.1.3          lava_1.6.9           ellipsis_0.3.2      
 [56] gplots_3.1.1         BiocGenerics_0.38.0  Rcpp_1.0.7           base64enc_0.1-3      zlibbioc_1.38.0     
 [61] purrr_0.3.4          rpart_4.1.16         MBA_0.0-9            viridis_0.6.1        S4Vectors_0.30.0    
 [66] zoo_1.8-9            haven_2.4.1          cluster_2.1.2        colorRamps_2.3       ncdfFlow_2.38.0     
 [71] scattermore_0.7      openxlsx_4.2.4       RPMG_2.2-3           lmtest_0.9-39        RANN_2.6.1          
 [76] ggnewscale_0.4.5     mime_0.10            hms_1.1.0            xtable_1.8-4         XML_3.99-0.6        
 [81] rio_0.5.27           jpeg_0.1-8.1         shape_1.4.6          ggcyto_1.20.0        bdsmatrix_1.3-4     
 [86] compiler_4.1.3       GEOmap_2.4-4         tibble_3.1.2         maps_3.3.0           KernSmooth_2.23-20  
 [91] crayon_1.4.1         htmltools_0.5.2      minqa_1.2.4          ggpointdensity_0.1.0 later_1.3.0         
 [96] Formula_1.2-4        RcppParallel_5.1.4   lubridate_1.7.10     aws.s3_0.3.21        DBI_1.1.1           
[101] tweenr_1.0.2         MASS_7.3-55          car_3.0-10           rbibutils_2.2.7      parallel_4.1.3      
[106] dotCall64_1.0-1      gower_0.2.2          forcats_0.5.1        pkgconfig_2.0.3      km.ci_0.5-2         
[111] foreign_0.8-82       sp_1.4-5             recipes_0.1.16       xml2_1.3.2           foreach_1.5.1       
[116] plm_2.4-3            prodlim_2019.11.13   stringr_1.4.0        digest_0.6.27        graph_1.70.0        
[121] cellranger_1.1.0     survMisc_0.5.5       maxLik_1.5-2         curl_4.3.1           shiny_1.7.1         
[126] gtools_3.9.2         lifecycle_1.0.0      nlme_3.1-155         jsonlite_1.7.2       carData_3.0-4       
[131] viridisLite_0.4.0    fansi_0.5.0          pillar_1.6.1         fastmap_1.1.0        httr_1.4.2          
[136] glue_1.4.2           zip_2.2.0            png_0.1-7            iterators_1.0.13     Rgraphviz_2.36.0    
[141] ggforce_0.3.3        class_7.3-20         stringi_1.6.2        moments_0.14         CytoML_2.4.0        
[146] latticeExtra_0.6-29  caTools_1.18.2       cytolib_2.4.0       