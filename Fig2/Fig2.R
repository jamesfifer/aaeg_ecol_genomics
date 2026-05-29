library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)

  library(plotly)
  library(scales)
  library(data.table)
library(zoo)
library(boot)
  library(purrr)
library(mgcv)
library(extrafont)

extrafont::font_import(prompt = FALSE)
loadfonts(device = "win")


theme_set(theme_classic(base_family = "Arial"))
  

  ##############################################
  ###########Frequency based f4ratio############
  ##############################################
  
  library(data.table)

  #Sahel
  #Sum
  noprefAafW=c("DBE","FTK","KAF","NGY","SBR","SCS","TCD","TWA","VLG","RBT","KLK","DMS","TBA","DHR","LIG","KLD","BKJ","RNR","GNS","GLB") #Hmm DBE &TWA has a tiny amount of AafE
  
  #OHInani SHM outgroup for (only really relevant for west guys) 
  files<- list.files( pattern = ".*OHInani.*.BKKa.*.SHMna_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
  files<- list.files( pattern = ".*OHInani.*.BKKa.*(STL|OGD|NGO|PKT).*.SHMna_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
  files<- list.files( pattern = ".*OHInani.*.BKKa.*(DBE|FTK|KAF|NGY|SBR|SCS|TCD|TWA|VLG|RBT|KLK|DMS|TBA|DHR|LIG|KLD|BKJ|RNR|GNS|GLB).*.SHMna_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
  
  
  noprefAafE=c("BUN","KAR","KIC","LUA","RABDOM", "SKU","YAO")
  
  #SHM vs BKK with masc outgroup (can be used for east and west, but harder to interpret for west)
  files<- list.files( pattern = ".*SHMna.*.BKKa.*.masc_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
  AafE=c("RAB","GND","KWA","ABK","KBO","KAK","ENT","VMB","FCV","LPV","LBV", "CPV")
  
    files<- list.files( pattern = ".*SHMna.*.BKKa.*(BUN|KAR|KIC|LUA|RABDOM|SKU|YAO).*.masc_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
  
    files<- list.files( pattern = ".*SHMna.*.BKKa.*(RAB|FCV).*.masc_verilyfilt_subset.10kb.sums.f4ratio", full.names = TRUE)
    
  
  
  
  files  
windows=c("1000","500","100","50","10")
for (window in 1:length(windows)){
  Window=as.numeric(windows[window])
  print(Window)
  combine=list()
for (f4output in 1:length(files)){
  filename=files[f4output]
  matches <- regmatches(filename, gregexpr("(?<=\\.)[A-Z]+(?=_chr)", filename, perl = TRUE))[[1]]
  matches <- regmatches(filename, gregexpr("([^.]+)_verilyfilt_subset", filename))[[1]]
  Xpop=sub("_verilyfilt_subset", "", matches)[3]
  print(Xpop)
  if (Xpop %in% c("BKK", "brom", "masc", "BTT", "SHM")) {
    message("Skipping population: ", Xpop)
    next
  }
  
  if (!(Xpop %in% AafE )) {
    message("Skipping population: ", Xpop)
    next
  }
  df=fread(file=filename)
  colnames(df)=c("CHR","START",
                 "END","Num","Denom")
  df=subset(df, CHR=="NC_035107.1" | CHR=="NC_035108.1" | CHR== "NC_035109.1") 
  df[df$Num == ".", c("Num", "Denom")] <- NA
  df[df$Denom == ".", c("Num", "Denom")] <- NA
  #if using f4ratiov2
  df_clean=df
df_clean$Num=as.numeric(df_clean$Num)
df_clean$Denom=as.numeric(df_clean$Denom)

  df_clean[, beta := rollapply(Num, width = Window, FUN = mean, na.rm = TRUE, fill = NA, align = "center") /
       rollapply(Denom, width = Window, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
     by = CHR]
  
  setorder(df_clean, CHR, START)
  df_clean[, rollmeanSTART := rollapply(START, Window, FUN = function(x) x[1], fill = NA, align = "center"), by = .(CHR)]
  df_clean[, rollmeanEND   := rollapply(END,   Window, FUN = function(x) x[Window], fill = NA, align = "center"), by = .(CHR)]
  #remove centromeres
  filtered_data <- df_clean %>%
    dplyr::filter(!(CHR == "NC_035107.1" & START >= 145000000 & START <= 155000000)) %>%
    dplyr::filter(!(CHR == "NC_035108.1" & START >= 227000000 & START <= 232000000)) %>%
    dplyr::filter(!(CHR == "NC_035109.1" & START >= 196000000 & START <= 201000000))
  filtered_data$Pop=Xpop
  key=paste0(Xpop)
  
  combine[[key]]=filtered_data
  average=mean(filtered_data$beta,na.rm=TRUE)
  print("Printing average:")
  print(average)
}    
  ##Combine everything
  combined_df <- do.call(rbind, combine)

  # what about then only extract windows that are between -1 and 1 for beta? (ending up doing this because everything >1 and <-1 are tree violations so hard to interpret)
  combined_df <- combined_df[combined_df$beta >= -1 & combined_df$beta <= 1, ]
 
  #Take the ones that have prefdata
  allbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
  
  #setting ancestry <0 to 0 because otherwise generalist's with practically 0 hs will have more variation for model to pull from
  combined_df$beta=ifelse(combined_df$beta<0,0, combined_df$beta)
  #get rid of NAs!
  combined_df <- combined_df[!is.na(combined_df$beta), ]
  
  combined_df%>%
    group_by(Pop) %>%
    summarise(combined_df = mean(beta, na.rm = TRUE))
  
  
  combined_df_pref=subset(combined_df, (combined_df$Pop) %in% allbehave$Population) #This should be 35 pops
  unique(combined_df_pref$Pop)
  actualwindow=as.numeric(Window)*10
  #write.table(combined_df_pref, file=paste0("prefonly_",actualwindow,"kb_OHInani_BKKab_SHMna.sums.f4ratio"), quote=F, row.names=F)
  write.table(combined_df_pref, file=paste0("prefonly_AafE_",actualwindow,"kb_OHInani_BKKab_SHMna.sums.f4ratio"), quote=F, row.names=F)
  
  #write.table(combined_df, file="allpops_10mb_OHInani_BKKab_SHMna.sums.f4ratio", quote=F, row.names=F)
}
#
combined_df_prefAaE=read.table(file="prefonly_AafE_500kb_OHInani_BKKab_SHMna.sums.f4ratio", header=T)
combined_df_prefAaboth=read.table(file="prefonly_500kb_OHInani_BKKab_SHMna.sums.f4ratio",header=T)
combined_df_prefAaW=subset(combined_df_prefAaboth, !(Pop %in% AafE))
unique(combined_df_prefAaW$Pop)
combined_df_prefall=rbind(combined_df_prefAaE,combined_df_prefAaW)
write.table(combined_df_prefall, file="prefonly_500kb_AafW_AafE_combined.f4ratio", quote=F, row.names=F)
#

#write.table(combined_df_pref, file="prefonly_10mb_OHInani_SKU_masc.sums.f4ratio", quote=F, row.names=F)
#write.table(combined_df_pref, file="prefonly_10mb_SHMna_BKKab_masc.sums.f4ratio", quote=F, row.names=F)
write.table(combined_df, file="allpops_10mb_SHMna_BKKab_masc.sums.f4ratio", quote=F, row.names=F)

noprefAafW=c("DBE","FTK","KAF","NGY","SBR","SCS","TCD","TWA","VLG","RBT","KLK","DMS","TBA","DHR","LIG","KLD","BKJ","RNR","GNS","GLB") #Hmm DBE &TWA has a tiny amount of AafE
noprefAafE=c("BUN","KAR","KIC","LUA","RABDOM", "SKU","YAO")

combined_df_noprefAaW=subset(combined_df, (combined_df$Pop) %in% noprefAafW)
unique(combined_df_noprefAaW$Pop)
write.table(combined_df_noprefAaW, file="noprefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio", quote=F, row.names=F)

combined_df_noprefAaE=subset(combined_df, combined_df$Pop %in% noprefAafE)
unique(combined_df_noprefAaE$Pop)
write.table(combined_df_noprefAaE, file="noprefonly_AafE_10mb_SHMna_BKKab_masc.sums.f4ratio", quote=F, row.names=F)



AafE=c("RAB","GND","KWA","ABK","KBO","KAK","ENT","VMB","FCV","LPV","LBV", "CPV")
combined_df_prefAaW=subset(combined_df_pref, !((combined_df_pref$Pop) %in% AafE))
unique(combined_df_prefAaW$Pop)

write.table(combined_df_prefAaW, file="prefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio", quote=F, row.names=F)
#write.table(combined_df_prefAaW, file="prefonly_AafW_10mb_OHInani_SKU_masc.sums.f4ratio", quote=F, row.names=F)


##Plotting 
#the pref only 

prefonly_AafW_AafE_combined.f4ratio=read.table(file="prefonly_AafW_AafE_combined.f4ratio", header=T)
#the no pref
combined_df_noprefAaE=read.table(file="noprefonly_AafE_10mb_SHMna_BKKab_masc.sums.f4ratio", header=T)
combined_df_noprefAaW=read.table(file="noprefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio",header=T)
noprefAafW=c("DBE","FTK","KAF","NGY","SBR","SCS","TCD","TWA","VLG","RBT","KLK","DMS","TBA","DHR","LIG","KLD","BKJ","RNR","GNS","GLB") #Hmm DBE &TWA has a tiny amount of AafE
noprefAafE=c("BUN","KAR","KIC","LUA","RABDOM", "SKU","YAO")
unique(combined_df_noprefAaW$Pop)
combined_df_all=rbind(combined_df_noprefAaE,combined_df_noprefAaW,prefonly_AafW_AafE_combined.f4ratio)
write.table(combined_df_all,file="sahel2_all_10mb.sums.f4ratio.txt", quote=F, row.names=F)
combined_df_noprefall=rbind(combined_df_noprefAaE,combined_df_noprefAaW)


combined_df=prefonly_AafW_AafE_combined.f4ratio
combined_df=combined_df_noprefall
#right now I made it so the mean line is always the mean of the pref data, but could change that to all?
library(ggrastr)
combined_df=as.data.table(combined_df)
combined_df[, rollmid := (rollmeanSTART + rollmeanEND) / 2]
prefonly_AafW_AafE_combined.f4ratio=as.data.table(prefonly_AafW_AafE_combined.f4ratio)
prefonly_AafW_AafE_combined.f4ratio[, rollmid := (rollmeanSTART + rollmeanEND) / 2]

mean_df <- aggregate(beta ~ rollmid + CHR, data = prefonly_AafW_AafE_combined.f4ratio, FUN = mean)
lassobed=read.table(file="0.992percentilegamgausregresults_prefonly_AafW_AafE_combined.f4rationocovariate.RDS.f4ratio.bed") #this is the gam
colnames(lassobed)[1]="CHR"
plot_list <- list()
alphabetical_pops <- sort(unique(combined_df$Pop))
alphabetical_pops
for (i in 1:length(alphabetical_pops)) {
  POP=alphabetical_pops[i]
  print(POP)
  subcombined_df=subset(combined_df, Pop==POP)
  p1=ggplot(subcombined_df, aes(x = rollmid / 1e6, y = beta, color = Pop)) +
    geom_rect(data = lassobed,
              aes(xmin = V2/1e6, xmax = V3/1e6, ymin = 0, ymax = 1), 
              color = "lightgrey", fill = "lightgrey", linewidth=1,inherit.aes = FALSE,group = "CHR")+
    rasterize(geom_point(data = mean_df, aes(x = rollmid / 1e6, y = beta), color = "black", size = 0.3),dpi=400) +
    rasterize(geom_point(size=0.25),dpi=400)+
    labs(title = paste(POP), x = "Position (Mb)", y = "Human-specialist ancestry") +
    scale_y_continuous(breaks = c(0,0.5, 1)) +  
    facet_wrap(~CHR, scales = "free_x",nrow=3) +  
    theme_classic() + theme(legend.title = element_blank(),strip.background = element_blank()) + ylim(0,1)+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(
      legend.position = "none",         
      strip.text = element_blank(),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  p1
  
  plot_list[[i]]=p1
}
#
library(patchwork)

blank_plot <- ggplot() + theme_void()

for (start in seq(1, length(plot_list), by = 8)) {
  end <- min(start + 7, length(plot_list))
  page_plots <- plot_list[start:end]
  
  npad <- 8 - length(page_plots)
  if (npad > 0) {
    page_plots <- c(page_plots, rep(list(blank_plot), npad))
  }
  
  pdf_name <- paste0("NoPrefF4Traceplots_", start, "_to_", end, ".pdf")
  
  pdf(pdf_name, width = 12, height = 6)  # adjust width/height as needed
  print(wrap_plots(page_plots, ncol = 4))
  dev.off()
}

#######################################################
#####Creating predictions from the nopref dataset######
combined_df_noprefAaE=read.table(file="noprefonly_AafE_10mb_SHMna_BKKab_masc.sums.f4ratio", header=T)
combined_df_noprefAaW=read.table(file="noprefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio",header=T)
noprefAafW=c("DBE","FTK","KAF","NGY","SBR","SCS","TCD","TWA","VLG","RBT","KLK","DMS","TBA","DHR","LIG","KLD","BKJ","RNR","GNS","GLB") #Hmm DBE &TWA has a tiny amount of AafE
noprefAafE=c("BUN","KAR","KIC","LUA","RABDOM", "SKU","YAO")



unique(combined_df_noprefAaW$Pop)
combined_df_noprefall=rbind(combined_df_noprefAaE,combined_df_noprefAaW)
head(combined_df_noprefall)
lassobed=read.table(file="0.992percentilegamgausregresults_prefonly_AafW_AafE_combined.f4rationocovariate.RDS.f4ratio.bed") #this is the gam

colnames(lassobed)=c("CHR","rollmeanSTART","rollmeanEND")


#create a region label to identify unique regions
lassobed <- lassobed %>%
  mutate(region_label = paste0("region_", row_number()))

head(combined_df_noprefall)
#For each region in lassobed, annotate matching windows in `data`
annotated_data <- map_dfr(1:nrow(lassobed), function(i) {
  region_row <- lassobed[i, , drop = FALSE]  
  region_row
  matched <- combined_df_noprefall %>%
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

#creating the all regions column
combined_data <- purrr::map_dfr(1:nrow(lassobed), function(i) {
  region_row <- lassobed[i, , drop = FALSE]
  combined_df_noprefall %>%
    dplyr::filter(
      CHR == region_row$CHR,
      rollmeanSTART >= region_row$rollmeanSTART,
      rollmeanEND <= region_row$rollmeanEND
    )
})

#Compute overall mean beta per Pop (weighted by window count)
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

nopref_meandf <- subset(wide_df, Pop %in% c(noprefAafE, noprefAafW))

#get the full model from further below when u run with the pref pops
full_model=readRDS(file="GAMGaus_GenoPrefPred_Model.RDS")

summary(full_model)
nopref_meandf$pred_pref <- plogis(predict(full_model, nopref_meandf))

nopref_meandf$unscaled_pred_pref=(nopref_meandf$pred_pref * 2)-1

write.table(nopref_meandf,file="0.992percentileAafW_AafE_noprefpops_prefpredictionsgamgaus.txt")

p1 <- ggplot(nopref_meandf, aes(x = allregions, y = unscaled_pred_pref , label = as.factor(Pop))) +
   geom_label_repel(
        size = 2.5,
        min.segment.length = 0,
        box.padding = 0.0,
        label.padding = 0.15,
        point.padding = 0.0,
        fill = "white",       
        label.size = 0.01,      
        segment.size = 0.35,
        segment.color = "black",
    max.overlaps = Inf)+
labs(
  x="Human-specialist ancestry (predictive windows)",
  y="Preference index (predicted)") +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )+xlim(0,1)+ylim(-1,1)
p1


cairo_pdf(
  file = paste0("predictedpref_0.992percentilenoprefonly_AafW_AafE_labelsonly_gamgaussplot.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p1)
dev.off()





##Run behavef4.R

#trying different percentiles
######
#Note: this takes a while so better to just run findthresholds_gamgaus.R on cluster
######
library(purrr)
 library(tidyr)
thresholds <- seq(0.60, 0.999, by = 0.001)
#RES="gamgausregresults_prefonly_1000kb_AafW_AafE_combined.f4rationocovariate.RDS"
#RES="gamgausregresults_prefonly_100kb_AafW_AafE_combined.f4rationocovariate.RDS"
RES="gambetaregresults_prefonly_10000kb_AafW_AafE_combined.f4ratiominpop_21nocovariate.RDS"



#data=read.table(file="prefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio", header=T)
#data=read.table(file="prefonly_AafW_AafE_combined.f4ratio", header=T)
data=read.table(file="prefonly_100kb_AafW_AafE_combined.f4ratio", header=T)
data=read.table(file="prefonly_10000kb_AafW_AafE_combined.f4ratio", header=T)

data$rollmeanSTART=as.numeric(data$rollmeanSTART)
data$rollmeanEND=as.numeric(data$rollmeanEND)
behave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
colnames(behave)[1]="Pop"

norelpath="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/Sahelii_norel/"
bams=data.frame(read.table(file=paste0(norelpath,"../../sahelii_meta.txt"), header = T)) # list of bam files
admix=read.table(file=paste0(norelpath,"K3_KPI6_allloci_md15_pop_admix.txt"), fill=T, header=T)


#A better alt is just to run findf4thresh.sh on the cluster
window=100000 #note true window is "window"*1,000, because input data is already in 10kb
rsquareds=list()
rsquaredsncv=list()
for (i in 1:length(thresholds)){
  percen=thresholds[i]
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
  pval95th$rollmeanEND=pval95th$rollmeanSTART+(window*1000)
  
  pval95th <- pval95th %>%
    arrange(CHR, rollmeanSTART) %>%
    group_by(CHR) %>%
    mutate(
      diff = rollmeanSTART - lag(rollmeanSTART, default = first(rollmeanSTART)),
      block_id = cumsum(if_else(diff > 10000, 1, 0))  
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
out
  outputname=paste0(percen,"percentile",out,window,"kb",".f4ratio.bed")
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
    region_row <- lassobed[i, , drop = FALSE]  
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
    matched
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
full_model=gam(scaledpref ~ s(allregions, bs="ts"), family=betar(),select=F, data=merged_df, method="REML")
full_modelncv=gam(scaledpref ~ s(allregions, bs="ts"), family=betar(),select=F, data=merged_df, method="NCV")

sum=summary(full_model)
print(sum)
merged_df$fit <- plogis(predict(full_model, merged_df))
merged_df$ncvfit <- plogis(predict(full_modelncv, merged_df))
rsquareds[[i]]=cor(merged_df$scaledpref, merged_df$fit)^2
sumncv=summary(full_modelncv)
print(sumncv)
rsquaredsncv[[i]]=cor(merged_df$scaledpref, merged_df$ncvfit)^2
names(rsquareds)[i]=paste0("Percentile_",percen)
names(rsquaredsncv)[i]=paste0("Percentile_",percen)
}


r2_vals <- unlist(rsquareds)
r2_ncvvals <-unlist(rsquaredsncv)
percentiles <- as.numeric(gsub("Percentile_", "", names(r2_vals)))
df=data.frame(cbind(r2_vals,r2_ncvvals,percentiles))

# Plot
p1=ggplot(data=df, aes(x=percentiles,y=r2_vals))+
  geom_point(fill = "white", size = 1, stroke = 0.25, color = "black", show.legend = FALSE, shape=21)+
#  ylim(0.5,0.85)+
 # ylim(0.7,0.8)+
  
  xlab("Percentile")+ylab("R²")+
  theme_classic(base_size = 18)
p1
write.table(df,file="gamgausspercentiles_prefonly_AafW_AafE_covariate_combined.f4ratio.txt") #0.999 best, also .999 with ncv
write.table(df,file="gamgausspercentiles_prefonly_AafW_AafE_nocovariate_combined.f4ratio.txt") #0.992 best, also .992 with ncv
write.table(df,file="gamgausspercentiles_prefonly_1000kbAafW_AafE_nocovariate_combined.f4ratio.txt") 
write.table(df,file="gamgausspercentiles_prefonly_100kbAafW_AafE_nocovariate_combined.f4ratio.txt") 


df=read.table(file="gamgausspercentiles_prefonly_1000kbAafW_AafE_nocovariate_combined.f4ratio.txt") #0.991
df=read.table(file="gamgausspercentiles_prefonly_100kbAafW_AafE_nocovariate_combined.f4ratio.txt", header=T) #0.999
df=read.table(file="gamgausspercentiles_prefonly_500kbAafW_AafE_nocovariate_combined.f4ratio.txt", header=T) #0.997
df=read.table(file="gamgausspercentiles_prefonly_5000kbAafW_AafE_nocovariate_combined.f4ratio.txt", header=T) #0.998


#pdf("betaregpercentiles_prefonly_AafW_AafE_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed
#pdf("betaregpercentiles_prefonly_AafW_AafE_nocovariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed
#pdf("gampercentiles_prefonly_AafW_AafE_nocovariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed
#pdf("gampercentiles_ncv_prefonly_AafW_AafE_nocovariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed

pdf("gamgausspercentiles_ncv_prefonly_AafW_AafE_nocovariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed
pdf("gamgausspercentiles_prefonly_AafW_AafE_nocovariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed

pdf("gamgausspercentiles_ncv_prefonly_AafW_AafE_covariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed
pdf("gamgausspercentiles_prefonly_AafW_AafE_covariate_combined.f4ratio_plot.pdf", width = 8, height = 4)  # adjust width and height as needed


print(p1)
dev.off()



#load in results here


# results=readRDS(file="gamgausregresults_prefonly_AafW_AafE_combined.f4ratiocovariate.RDS"
# )
results=readRDS(file="gamgausregresults_prefonly_10000kb_AafW_AafE_combined.f4rationocovariate.RDS"
) 

results=readRDS(file="gamgausregresults_prefonly_5000kb_AafW_AafE_combined.f4rationocovariate.RDS"
)
results=readRDS(file="gamgausregresults_prefonly_10000kb_AafW_AafE_combined.f4ratiominpop_30nocovariate.RDS"
) 


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
chr1=subset(chr123, CHR=="1")



# 
# percen=0.969
# percen=0.994
 percen=0.992


 
 
pvalue_90th_percentile <- quantile(chr123$R2, percen, na.rm = TRUE)
chr123$above_90th <- ifelse(chr123$R2 > pvalue_90th_percentile, TRUE, FALSE)



p1=ggplot(chr123, aes(x = rollmeanSTART / 1e6, y = R2)) +
  geom_point(color="black", size=1) +
  geom_point(data = subset(chr123, above_90th == TRUE), 
             aes(x = rollmeanSTART/1e6, y = R2), color = "red", size = 1) +
scale_x_continuous(name = "Position (Mb)") +
  facet_wrap(~CHR, scales = "free_x",nrow=3) +  
  theme_classic()+theme(strip.text = element_blank(), strip.background = element_blank())+
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )+
ylab("R²")
p1




library(ggplot2)
library(ggrastr)

fig2path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig2/"

#Rasterize points so I can actually manipulate in inkscape 
p1_raster <- ggplot(chr123, aes(x = rollmeanSTART / 1e6, y = R2)) +
  rasterise(geom_point(color="black", size=1), dpi=300) +
  rasterise(geom_point(data = subset(chr123, above_90th == TRUE),
                       aes(x = rollmeanSTART/1e6, y = R2),
                       color="red", size=1), dpi=300) +
  facet_wrap(~CHR, scales="free_x", nrow=3) +
  scale_x_continuous(name="Position (Mb)") +
  ylab("R²") +
  theme_classic() +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        text = element_text(family="Arial"),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14))

p1_raster
# ggsave(
#   paste0(fig2path, "0.992percentileprefonly_AafW_AafE_nocovariate_gamgaus.pdf"),
#   plot = p1_raster,
#   device = cairo_pdf,
#   width = 8,
#   height = 4
# )

# ggsave(
#      paste0(fig2path, "0.999percentileprefonly_AafW_AafE_covariate_gamgaus.pdf"),
#      plot = p1_raster,
#      device = cairo_pdf,
#      width = 8,
#      height = 4
#    )


ggsave(
  paste0(fig2path, "0.998percentileprefonly_500kbAafW_AafE_covariate_gamgaus.pdf"),
  plot = p1_raster,
  device = cairo_pdf,
  width = 8,
  height = 4
)


################################ MODEL TESTING #########################################
########################################################################################

###Extracting f4ratio values at the "significant" regions and looking at relationship with preference
#masc0
library(car)
library(dplyr)
library(tidyr)
library(purrr)

#data=read.table(file="prefonly_AafW_10mb_OHInani_BKKab_SHMna.sums.f4ratio", header=T)
data=read.table(file="prefonly_AafW_AafE_combined.f4ratio", header=T)
data$rollmeanSTART=as.numeric(data$rollmeanSTART)
data$rollmeanEND=as.numeric(data$rollmeanEND)

head(data)

behave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
colnames(behave)[1]="Pop"

lassobed=read.table(file="0.992percentilegamgausregresults_prefonly_AafW_AafE_combined.f4rationocovariate.RDS.f4ratio.bed") #this is the gam used

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
meandf$scaledpref <- (meandf$pref + 1) / 2

library(mgcv)

full_model=lm(pref ~ allregions, data = meandf)
summary(full_model)
AIC(full_model)
plot(full_model)             # for lm
plot(residuals(full_model))  # for betareg

plot(meandf$pref,full_model$fitted.values)

full_model=betareg(scaledpref ~ allregions, data = meandf)
summary(full_model)
AIC(full_model)
plot(full_model)
plot(meandf$scaledpref,full_model$fitted.values)

AIC(full_model)
region_vars <- grep("^region", names(meandf), value = TRUE)

top_vars <- head(sort(apply(meandf[region_vars], 2, var), decreasing = TRUE), 20)
top_vars
formula <- as.formula(paste("pref ~", paste(names(top_vars), collapse = " + ")))


#GAM
library(mgcv)
full_model <- gam(scaledpref ~ s(allregions, bs="ts"), family = betar(), select=F,
                  data= meandf, method="REML")
sum=summary(full_model)
sum
AIC(full_model)

full_model <- gam(logit(scaledpref) ~ s(allregions, bs="ts"), family = gaussian(), select=F,
                  data= meandf, method="REML")
saveRDS(full_model,file="GAMGaus_GenoPrefPred0.992_Model.RDS")
full_model=readRDS(file="GAMGaus_GenoPrefPred0.992_Model.RDS")
sum=summary(full_model)
sum
AIC(full_model)


# #loocv (actually doesnt matter because I already did this with the thresholding)
# full_modelcv <- gam(scaledpref ~ s(allregions, bs="ts"), family = betar(), select=F,
#                   data= meandf, method="NCV")
# sum=summary(full_modelcv)
# sum
# AIC(full_model)

pred <- predict(full_model, type = "link", se.fit = TRUE)
meandf$fit <- plogis(pred$fit)
link_upper <- pred$fit + 1.96 * pred$se.fit
link_lower <- pred$fit - 1.96 * pred$se.fit

meandf$upper <- plogis(link_upper)
meandf$lower <- plogis(link_lower)

plot(meandf$scaledpref,meandf$fit)
#unscale
meandf$unscaledfit=(meandf$fit*2)-1
meandf$unscaled_lower=(meandf$lower*2)-1
meandf$unscaled_upper=(meandf$upper*2)-1


r2_cv_report <- cor(meandf$scaledpref, meandf$fit)^2
r2_cv_report
r_squared <- round(r2_cv_report, 2)

annotation_text1 <- expression(italic("p") ~ "< 0.001")
annotation_text2 <- bquote( R^2 * "=" * .(r_squared))
# Now plot
#Add AafE/W
AafE=c("RAB","GND","KWA","ABK","KBO","KAK","ENT","VMB","FCV","LPV","LBV", "CPV")
meandf$AafE=ifelse(meandf$Pop %in% AafE, "E","W")
#Add batch 
#add batch
meandf
dat<-read.csv('C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/behavior/2021_big_olfac_experiment.csv')
dat$Pop<-as.character(dat$Pop)
dat$Pop[dat$Pop=='U52']<-'ZIK'
dat<-dat[!is.na(dat$N),]
old<-dat
dat<-read.csv('C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/behavior/2021_big_olfac_experiment.csv')
dat$Pop<-as.character(dat$Pop)
dat$Pop[dat$Pop=='U52']<-'ZIK'
dat<-dat[!is.na(dat$N),]
dat <- dat %>% filter(Pop != "NGO")
#
meandf$batch=as.numeric(ifelse(meandf$Pop %in% dat$Pop, "1","0"))
meandf$batch=as.numeric(ifelse(meandf$Pop =="ZIK", "0",meandf$batch))

#install.packages("extrafont")
library(extrafont)

extrafont::font_import(prompt = FALSE)
loadfonts(device = "win")
#fonts()
theme_set(theme_classic(base_family = "Arial"))

write.table(meandf,file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/introgression/0.992percentileAafW_AafE_prefpopsonly_prefpredictionsgamgaus.txt")



fig2path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig2/"

p1 <- ggplot(meandf, aes(allregions, pref, label=Pop)) +
  geom_line(aes(y = unscaledfit), color = "black", linetype = "dashed", size = 0.8) +
  geom_ribbon(aes(ymin = unscaled_lower, ymax = unscaled_upper), fill = "grey80", alpha = 0.5)+  # CI band
  scale_fill_manual(values=c("E" = "#07539dff", "W" = "#b9d7e6ff" ))+
    geom_label_repel( 
      size = 2.5, 
      min.segment.length = 0, 
      box.padding = 0.0, 
      label.padding = 0.15, 
      point.padding = 0.0, 
      fill = "white",  
      label.size = 0.01, 
      segment.size = 0.35, 
      segment.color = "black", 
      max.overlaps = Inf)+
  scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  labs(
    x = "Human-specialist ancestry \n (predictive regions)",
    y = "Preference index"
  ) +
  annotate("text", x = .04, y = .9, label = "italic('p') < 0.001", parse = TRUE, size = 6) +
  annotate("text", x = .04, y = .8, 
           label = paste0("R^2 == ", signif(r_squared, 3)), parse = TRUE, size = 6) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  ) +
  ylim(-1, 1)
p1

fig2path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/212023_Fifer_Sahel_Genomics/Manuscript/Fig2/"


cairo_pdf(
  file = paste0(fig2path, "pref_0.992percentileprefonly_AafW_AafE_labelsonly_gamgaussplot.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p1)
dev.off()

#Model selection gam

region_vars <- grep("^region", names(meandf), value = TRUE)

top_vars <- head(sort(apply(meandf[region_vars], 2, var), decreasing = TRUE), 20)
top_vars
smooth_terms <- paste0("s(", names(top_vars), ", bs='ts')", collapse = " + ")

# Combine into full model formula
formula <- as.formula(paste("scaledpref ~", smooth_terms))
 
formula
# Fit model
full_model2 <- gam(
  formula = formula,
  family = betar(),
  select = FALSE,
  data = meandf,
  method = "REML"
)
sum=summary(full_model2)
sum


full_model3=gam(scaledpref ~ s(region_1, bs="ts")+s(region_7, bs="ts")+
                  s(region_6, bs="ts"), family = betar(), select=F,data= meandf, method="REML")
sum=summary(full_model3)
sum
AIC(full_model3)

pred <- predict(full_model3, type = "link", se.fit = TRUE)
meandf$fit <- plogis(pred$fit)

r2_cv_report <- cor(meandf$scaledpref, meandf$fit)^2
r2_cv_report
r_squared <- round(r2_cv_report, 2)




full_model4=gam(scaledpref ~ s(region_4, bs="ts")+s(region_7, bs="ts")
                 , family = betar(), select=F,data= meandf, method="NCV")
sum=summary(full_model4)
sum
AIC(full_model4)

full_model5=gam(scaledpref ~ s(region_3, bs="ts")+s(region_7, bs="ts")
                , family = betar(), select=F,data= meandf, method="NCV")
sum=summary(full_model5)
sum
AIC(full_model5)

full_model6=gam(scaledpref ~ s(region_7, bs="ts")
                , family = betar(), select=F,data= meandf, method="NCV")
sum=summary(full_model6)
sum
AIC(full_model5)

#full_model 3 does the best with the NCV

