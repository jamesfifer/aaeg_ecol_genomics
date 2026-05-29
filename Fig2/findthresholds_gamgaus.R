


required_packages <- c("ggplot2", "dplyr", "stringr", "furrr", "plotly", "betareg","tidyr","purrr","mgcv", "future","future.apply","car")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

lapply(required_packages, library, character.only = TRUE)



library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(mgcv)
library(future)
library(future.apply)
library(car)
plan(multicore, workers = 20)     # adjust cores


setwd("/tscc/projects/ps-roselab/jfifer/SahelII/Intermediate/")
args <- commandArgs(trailingOnly = TRUE)
RES <- args[1]
data <- args[2]
window <- as.numeric(args[3])
#window=100000
#RES="gamgausregresults_prefonly_100kb_AafW_AafE_combined.f4rationocovariate.RDS"
#data=read.table(file="prefonly_100kb_AafW_AafE_combined.f4ratio", header=T)


thresholds <- seq(0.60, 0.999, by = 0.001)

data=read.table(file=data,header=T)
data$rollmeanSTART=as.numeric(data$rollmeanSTART)
data$rollmeanEND=as.numeric(data$rollmeanEND)
behave=read.table(file="all_behavior.txt", header=T)
colnames(behave)[1]="Pop"

bams=data.frame(read.table(file=paste0("./","sahelii_meta.txt"), header = T)) # list of bam files
admix=read.table(file=paste0("./","K3_KPI6_allloci_md15_pop_admix.txt"), fill=T, header=T)

run_one_threshold <- function(percen) {
  #percen=thresholds[i]
  print(percen)
  results=readRDS(file=RES)
  print(RES)
  head(results)
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  chr123=results
  chr123 <- chr123 %>%
    mutate(CHR = case_when(
      CHR == "NC_035107.1" ~ "1",
      CHR == "NC_035108.1" ~ "2",
      CHR == "NC_035109.1" ~ "3",
      TRUE ~ CHR
    ))
  chr123 <- chr123[!is.na(chr123$p_adj) & is.finite(chr123$p_adj), ]
  chr123$CHR=as.numeric(chr123$CHR)
  chr123$p_adj=as.numeric(chr123$p_adj)
  chr123$snp=as.character(row.names(chr123))
  chr123$logp=-log10(chr123$p_adj)


  pvalue_90th_percentile <- quantile(chr123$R2, percen, na.rm = TRUE)
  chr123$above_90th <- ifelse(chr123$R2 > pvalue_90th_percentile, TRUE, FALSE)

  pval95th=subset(chr123,above_90th==TRUE)
  pval95th$rollmeanEND=pval95th$rollmeanSTART+window

  pval95th <- pval95th %>%
    arrange(CHR, rollmeanSTART) %>%
    group_by(CHR) %>%
    mutate(
      diff = rollmeanSTART - lag(rollmeanSTART, default = first(rollmeanSTART)),
      block_id = cumsum(if_else(diff > 10000, 1, 0)) #each 10kb start gets it own block now
    ) %>%
    ungroup()
  block_summary <- pval95th %>%
    group_by(CHR, block_id) %>%
    summarise(
      start = min(rollmeanSTART),
      end = max(rollmeanEND),
      length_kb = (end - start + 10000) / 1000,
      n_loci = n(),
      avg_p = mean(p_value, na.rm = TRUE),
      #avg_beta = mean(beta_coef, na.rm = TRUE),
      .groups = "drop"
    )
  block_summary <- block_summary %>%
    mutate(CHR = as.character(CHR)) %>%
    mutate(CHR = case_when(
      CHR == "1" ~ "NC_035107.1",
      CHR == "2" ~ "NC_035108.1",
      CHR == "3" ~ "NC_035109.1",
      TRUE ~ CHR
    ))
  out <- paste(strsplit(RES, "_")[[1]][1:5], collapse = "_")
#out
  outputname=paste0(percen,"percentile",out,".f4ratio.bed")
print(outputname)
    write.table(block_summary %>% dplyr::select(CHR,start,end), file=outputname, col.names =F, row.names = F,quote=F, sep='\t')
  lassobed=read.table(file=outputname, sep='\t', header=F)
    colnames(lassobed)=c("CHR","rollmeanSTART","rollmeanEND")
  lassobed$rollmeanSTART=as.numeric(lassobed$rollmeanSTART)
  lassobed$rollmeanEND=as.numeric(lassobed$rollmeanEND)

  head(lassobed)

  lassobed <- lassobed %>%
    mutate(region_label = paste0("region_", row_number()))

  annotated_data <- map_dfr(1:nrow(lassobed), function(i) {
    region_row <- lassobed[i, , drop = FALSE]  # avoid naming conflict
    region_row
    matched <- data %>%
      dplyr::filter(
        CHR == region_row$CHR,
        rollmeanSTART >= region_row$rollmeanSTART,
        rollmeanEND <= region_row$rollmeanEND
      ) %>%
      mutate(
        region_label = region_row$region_label,
        beta = as.numeric(beta)
      )

    if (nrow(matched) == 0) return(NULL)

    matched %>%
      group_by(Pop) %>%
      summarise(
        mean_beta = mean(beta, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        region = region_row$region_label,
        CHR = region_row$CHR,
        rollmeanSTART = region_row$rollmeanSTART,
        rollmeanEND = region_row$rollmeanEND
      )
  })

  combined_data <- purrr::map_dfr(1:nrow(lassobed), function(i) {
    region_row <- lassobed[i, , drop = FALSE]
    data %>%
      dplyr::filter(
        CHR == region_row$CHR,
        rollmeanSTART >= region_row$rollmeanSTART,
        rollmeanEND <= region_row$rollmeanEND
      )
  })

  allregions_summary <- combined_data %>%
    group_by(Pop) %>%
    summarise(
      allregions = mean(as.numeric(beta), na.rm = TRUE),
      .groups = "drop"
    )

  wide_df <- annotated_data %>%
    dplyr::select(Pop, region, mean_beta) %>%
    pivot_wider(
      id_cols = Pop,
      names_from = region,
      values_from = mean_beta,
      values_fill = list(mean_beta = 0)
    ) %>%
    left_join(allregions_summary, by = "Pop")

meandf=merge(wide_df,behave, by = "Pop")
merged_df <- merge(admix, meandf, by = "Pop") #
merged_df$scaledpref <- (merged_df$pref + 1) / 2
full_model=gam(logit(scaledpref) ~ s(allregions, bs="ts"), family=gaussian(),select=F, data=merged_df, method="REML")
full_modelncv=gam(logit(scaledpref) ~ s(allregions, bs="ts"), family=gaussian(),select=F, data=merged_df, method="NCV")
#full_model=betareg(scaledpref ~ allregions, data = merged_df)
#full_model=lm(pref ~ allregions +aaa, data = merged_df)
sum=summary(full_model)
print(sum)
#rsquareds[[i]]=sum[["r.sq"]]
merged_df$fit <- plogis(predict(full_model, merged_df))
merged_df$ncvfit <- plogis(predict(full_modelncv, merged_df))

r2     <- cor(merged_df$scaledpref, merged_df$fit)^2
  r2_ncv <- cor(merged_df$scaledpref, merged_df$ncvfit)^2

  return(data.frame(
    percentile = percen,
    r2 = r2,
    r2_ncv = r2_ncv
  ))
}

results_list <- future_lapply(thresholds, run_one_threshold)
df <- dplyr::bind_rows(results_list)

windowinkb=window/1000
write.table(
  df,
  file = paste0("gamgausspercentiles_prefonly_", windowinkb,
                "kbAafW_AafE_nocovariate_combined.f4ratio.txt"),
  row.names = FALSE
)





