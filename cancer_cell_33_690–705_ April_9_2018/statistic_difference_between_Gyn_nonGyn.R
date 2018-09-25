library(magrittr)

# load data path ----------------------------------------------------------
tcga_path <- "/data/shiny-data/GSCALite/TCGA/snv"
cosmic_path <- "/project/huff/huff/data/COSMIC/cancer-related-genes"

# load data ---------------------------------------------------------------

mutation_data <- readr::read_rds(file.path(tcga_path,".rds_snv_all_gene_snv_count.rds.gz"))
oncogenes <- readr::read_tsv(file.path(cosmic_path,"Tier1+Tier2_Census_allWed Apr 11 03_19_58 2018.tsv"))
  

# functions ---------------------------------------------------------------

filter_gene_list <- function(.x, gene_list) {
  .x %>%
    dplyr::filter(symbol %in% gene_list)
}


# filter data -------------------------------------------------------------
mutation_data %>%
  dplyr::mutate(filter_snv = purrr::map(mut_count,filter_gene_list, gene_list = oncogenes$`Gene Symbol`)) %>%
  dplyr::select(-mut_count) -> muta_of_genelist


# dealing with data -------------------------------------------------------

muta_of_genelist %>%
  dplyr::mutate(num_of_cancer = n) %>%
  dplyr::select(-n) -> muta_of_genelist

## For Pan-Gyn ----

pan.Gyn <- c("CESC","UCS","UCEC","BRCA","OV")
muta_of_genelist %>%
  dplyr::filter(cancer_types %in% pan.Gyn) %>%
  dplyr::mutate(num_of_panGyn = sum(num_of_cancer)) %>%
  tidyr::unnest() %>%
  dplyr::group_by(symbol) %>% 
  dplyr::mutate(mute_sum_in_panGyn = sum(sm_count)) %>%
  dplyr::select(symbol,mute_sum_in_panGyn,num_of_panGyn) %>%
  dplyr::ungroup() %>% unique()  -> panGyn_mutate

## For non Pan-Gyn ----

muta_of_genelist %>%
  dplyr::filter(! cancer_types %in% pan.Gyn) %>%
  dplyr::mutate(num_of_non_Gyn = sum(num_of_cancer)) %>%
  tidyr::unnest() %>%
  dplyr::group_by(symbol) %>% 
  dplyr::mutate(mute_sum_in_non_Gyn = sum(sm_count)) %>%
  dplyr::select(symbol,mute_sum_in_non_Gyn,num_of_non_Gyn) %>%
  dplyr::ungroup() %>% unique() -> non_Gyn_mutate

## Combine data ----
fn_fisher <- function(.x){
  .x %>%
    as.data.frame() %>%
    .[1,] %>%
    unlist() %>%
    matrix(nr=2) -> .matrix
  fisher.test(.matrix) %>%
    matrix() %>%
    .[1,1] %>%
    unlist() %>%
    p.adjust(method = "bonferroni")
}

panGyn_mutate %>%
  dplyr::inner_join(non_Gyn_mutate,by="symbol") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(pvalue=purrr::map(data,fn_fisher)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::filter(pvalue < 0.01)
