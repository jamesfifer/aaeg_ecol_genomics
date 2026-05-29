options(repos = c(CRAN = "https://cloud.r-project.org"))


required_packages <- c("ggplot2", "dplyr", "stringr", "furrr", "plotly", "betareg")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) {
  install.packages(new_packages)
}

lapply(required_packages, library, character.only = TRUE)



library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)
library(furrr)
library(betareg)
library(mgcv)
library(car)

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
admixfile <- args[2]
minpop <- as.integer(args[3])


#filename="prefonly_AafW_AafE_combined.f4ratio"
#admixfile="K3_KPI6_allloci_md15_pop_admix.txt"
setwd("/tscc/projects/ps-roselab/jfifer/SahelII/Intermediate/")

df <- read.table(filename, header = TRUE, fill = TRUE)
admix <- read.table(admixfile, header = TRUE, fill = TRUE)

# Load and merge data
behave <- read.table("all_behavior.txt", header = TRUE)
df$Pop <- gsub("\\.\\/", "", df$Pop)
predata <- merge(behave, df)
print(nrow(predata))
data <- merge(predata,admix)
print(nrow(data))
perform_betareg <- function(group_data) {
  if (nrow(group_data) < minpop || sd(group_data$beta, na.rm = TRUE) == 0) {
    return(data.frame(
      CHR = unique(group_data$CHR),
      rollmeanSTART = unique(group_data$rollmeanSTART),
      p_value = NA,
        edf = NA,
      R2 = NA
    ))
  }

  group_data$pref_t <- (group_data$pref + 1) / 2

  fit <- tryCatch(
gam(logit(pref_t) ~ s(beta, bs="ts"), family = gaussian(), select=F, data = group_data, method="REML"),
#betareg(pref_t ~ beta + aaa, data = group_data),
    error = function(e) NA
  )

  if (inherits(fit, "gam")) {
sum=summary(fit)
    R2 <- tryCatch(sum[["r.sq"]], error = function(e) NA)
    edf=tryCatch(sum[["edf"]], error = function(e) NA)
    pval=tryCatch(sum[["p.table"]][4], error = function(e) NA)
  } else {
    pval <- NA
    R2 <- NA
    edf <- NA
npops <- nrow(group_data)
  }

  data.frame(
    CHR = unique(group_data$CHR),
    rollmeanSTART = unique(group_data$rollmeanSTART),
    p_value = pval,
        edf=edf,
    R2 = R2,
npops = nrow(group_data)
  )
}

plan(multicore, workers = future::availableCores())

grouped <- split(data, list(data$CHR, data$rollmeanSTART), drop = TRUE)

#Run in parallel
results <- future_map_dfr(grouped, perform_betareg, .progress = TRUE)

saveRDS(results, file = paste0("gamgausregresults_", filename, "minpop_",minpop ,"nocovariate.RDS"))

