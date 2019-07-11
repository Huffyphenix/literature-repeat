########## get CTL(cytotoxic T lymphocytes) infiltration level was estimated as the average expression level of  CD8A, CD8B, GZMA,GZMB and  PRF1 .
########## Reference:Signatures of T cell dysfunction and exclusion predict cancer immunotherapy response

basic_path <- file.path("/home/huff/project")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")

# load expression data -----
expr <- readr::read_rds(file.path(TCGA_path, "pancan33_expr.rds.gz"))
ctl_standard_genes <- c("CD8A", "CD8B", "GZMA","GZMB", "PRF1")

# estimation of CTL levels by the average expression level of  CD8A, CD8B, GZMA,GZMB and  PRF1------
expr %>%
  dplyr::mutate(CTL_from_CD8A_CD8B_GZMA_GZMB_PRF1 = purrr::map(expr,.f=function(.x){
    .x %>%
      dplyr::filter(symbol %in% ctl_standard_genes) %>%
      dplyr::select(-entrez_id) %>%
      tidyr::gather(-symbol,key="sample",value="exp") %>%
      tidyr::spread(key="symbol",value="exp") %>%
      dplyr::mutate(`CTL(average_exp)` = (CD8A+CD8B+GZMA+GZMB+PRF1)/5)
  })) %>%
  dplyr::select(-expr) -> CTL_estimated_from_CD8A_CD8B_GZMA_GZMB_PRF1

CTL_estimated_from_CD8A_CD8B_GZMA_GZMB_PRF1 %>%
  readr::write_rds("/home/huff/project/data/TCGA/CTL_level_estimated/CTL_estimated_from_CD8A_CD8B_GZMA_GZMB_PRF1.rds.gz",compress = "gz")
