#################### get mutation burden ######################


# data path ---------------------------------------------------------------

tcga_path <- "/data/shiny-data/GSCALite/TCGA/snv"


# load data ---------------------------------------------------------------

snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))


# analysis ----------------------------------------------------------------
cancers_26 <- c("SKCM","LUSC","LUAD","BLCA","COAD","STAD","UCEC","HNSC","READ","ACC","LIHC","CHOL","UCS","PAAD",
                "BRCA","KIRC","KICH","KIRP","UVM","TGCT","PRAD","LGG","GBM","OV","PCPG","THCA")

fn_sample_mutation <- function(.x){
  .x %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(key = barcode, value = count, -symbol) %>% 
    dplyr::mutate(samples = ifelse(count > 0, 1, 0)) %>% 
    dplyr::group_by(barcode) %>% 
    dplyr::summarise(sm_count = sum(count)) -> .d_count
}

fn_high_mutate_per <- function(.x){
  .x %>%
    dplyr::filter(sm_count >= 192) -> .y      # mutation burden cutoff: 192
  nrow(.y)/nrow(.x) -> .out
}


snv %>%
  dplyr::filter(cancer_types %in% cancers_26) %>%
  dplyr::mutate(snv_count = purrr::map(snv,fn_sample_mutation)) %>%
  dplyr::select(-snv) -> cancer_26_count


# plot --------------------------------------------------------------------
# repeat the results of literatures:
# Burden of Nonsynonymous Mutations among TCGA Cancers and Candidate Immune Checkpoint Inhibitor Responses
# Published OnlineFirst May 18, 2016; DOI: 10.1158/0008-5472.CAN-16-0170

library(ggplot2)

cancer_26_count %>%
  dplyr::mutate(mutation_load = purrr::map(snv_count,fn_high_mutate_per)) %>%
  dplyr::select(-snv_count) %>%
  tidyr::unnest() %>%
  dplyr::arrange(dplyr::desc((mutation_load))) -> cancer_rank

cancer_26_count %>%
  tidyr::unnest() %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  ggplot(aes(x=cancer_types,y=sm_count)) +
  geom_violin(trim=FALSE,fill= "gray") +
  geom_boxplot(width=0.1) +
  theme_classic() +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  geom_hline(yintercept = log2(192),color="red")


# classify the TCGA samples with mutation burden cutoff: 192 --------------

cancer_26_count %>%
  tidyr::unnest() %>%
  dplyr::mutate(mutation_status = ifelse(sm_count >= 192, "high_muation_burden","low_mutation_burden")) -> cancer_26_mutation_burden_classification

out_path <- "/project/huff/huff/data/TCGA"
cancer_26_mutation_burden_classification %>% 
  tidyr::nest(-cancer_types) %>%
  readr::write_rds(file.path(out_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz"),compress = "gz")
