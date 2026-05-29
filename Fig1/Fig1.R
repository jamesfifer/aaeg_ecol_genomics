library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(readxl)
library(stringr)
library(terra)
library(geodata)
library(ggh4x)
library(tidyverse)
library(sf)
library(betareg)
library(car)

library(extrafont)

extrafont::font_import(prompt = FALSE)
loadfonts(device = "win")


fonts()
theme_set(theme_classic(base_family = "Arial"))


#Obtain GPS info
geo_data=read.csv(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental/SahelII_coordinates.csv")



### Plotting preference whole genome ancestry relationship
#
norelpath="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/Sahelii_norel/"
fig1path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig1/"

bams=data.frame(read.table(file=paste0(norelpath,"../../sahelii_meta.txt"), header = T)) # list of bam files

admix=read.table(file=paste0(norelpath,"K3_KPI6_allloci_md15.all.qopt"), fill=T, header=T)

newcombined1=merge(bams,admix, by.x="Adapter",by.y="adapter")
#behave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
behave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/all_behavior_cis.txt", header=T)


colnames(behave)[1]="Pop"
merged_df <- merge(behave, newcombined1, by = "Pop",by.y="Location") #
merged_df$pref=as.numeric(merged_df$pref)



meandf <- merged_df %>%
  group_by(Pop) %>%
  summarize(across(matches("Invaded|pref|upf|dwnf"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"))

# colnames(meandf)
# model.1=lm(mean_pref ~mean_Invaded, data = meandf)
# mod_summary=summary(model.1)
# AIC(model.1)
# mod_summary
# r_squared <- round(mod_summary$r.squared, 2)
# r_squared
# 
# model.1=betareg(scaled_pref ~mean_Invaded, data = meandf)
# mod_summary=summary(model.1)
# AIC(model.1)
# mod_summary
# r_squared <- round(mod_summary$pseudo.r.squared, 2)
# r_squared



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
meandf$batch=as.numeric(ifelse(meandf$Pop =="CPV", "0",meandf$batch))


#GAM Preference 
meandf$scaled_pref=(meandf$mean_pref +1)/2
#Also made a GAM Gauss vers for predictions for other projects
meandf$logitscaled_pref=logit(meandf$scaled_pref)
gaus_model <- gam(logitscaled_pref ~ s(mean_Invaded, bs="ts"), family = gaussian(), select=F,
                  data= meandf, method="REML")
pred <- predict(gaus_model, type = "response", se.fit = TRUE)
cor(meandf$logitscaled_pref,pred$fit)^2
saveRDS(gaus_model,file="GAMGaus_WholeGenoPrefPred_Model.RDS")


#scale pref for GAM betareg (used this for Sahel II manuscript)
library(mgcv)
full_model <- gam(scaled_pref ~ s(mean_Invaded, bs="ts"), family = betar(), select=F,
                  data= meandf, method="REML")
sum=summary(full_model)
sum
AIC(full_model)
pred <- predict(full_model, type = "response", se.fit = TRUE)

meandf$fit <- pred$fit
meandf$se.fit <- pred$se.fit

#upper and lower 95% confidence limits
meandf$upper <- meandf$fit + 1.96 * meandf$se.fit
meandf$lower <- meandf$fit - 1.96 * meandf$se.fit

#back-transform predictions to original (-1 to 1) scale
meandf$fit_unscaled   <- meandf$fit   * 2 - 1
meandf$upper_unscaled <- meandf$upper * 2 - 1
meandf$lower_unscaled <- meandf$lower * 2 - 1


r2_cv_report <- cor(meandf$scaled_pref, meandf$fit)^2
r2_cv_report
r_squared <- round(r2_cv_report, 2)

annotation_text1 <- "italic('p') < 0.001"   # now a character
annotation_text2 <- paste0("R^2 == ", signif(r_squared, 3))


# Now plot
p1 <- ggplot(meandf, aes(mean_Invaded, mean_pref, label=Pop)) +
  geom_ribbon(aes(ymin = lower_unscaled, ymax = upper_unscaled), fill = "grey80", alpha = 0.5) +  # CI band
  geom_line(aes(y = fit_unscaled), color = "black", linetype = "dashed", size = 0.8) +
   geom_point(aes(shape = factor(batch)), 
              fill = "white", size = 5, stroke = 0.25, color = "black", show.legend = FALSE) +
 # geom_label_repel( 
 #    size = 2.5, 
 #    min.segment.length = 0, 
 #    box.padding = 0.0, 
 #    label.padding = 0.15, 
 #    point.padding = 0.0, 
 #    fill = "white", # box outline color 
 #    label.size = 0.01, # thickness of outline 
 #    segment.size = 0.35, 
 #    segment.color = "black", 
 #    max.overlaps = Inf)+
   scale_shape_manual(values = c("0" = 21, "1" = 24))  +
   labs(
     x = "Human-specialist ancestry",
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

#pdf(paste0(norelpath,"pref_genomewide_aaa_dotsonly_plot.pdf"), width = 6, height = 6)  # adjust width and height as needed
#pdf(paste0(norelpath,"pref_genomewide_aaa_labelsonly_plot.pdf"), width = 6, height = 6)  # adjust width and height as needed
#print(p1)
cairo_pdf(
      file = paste0(fig1path, "pref_genomewide_aaa_dotsonly_gamplot.pdf"),
       width = 7,
       height = 7,
       family = "Arial" )
print(p1)
dev.off()

#loocv
full_model <- gam(scaled_pref ~ s(mean_Invaded, bs="ts"), family = betar(), select=F,
                  data= meandf, method="NCV")
sum=summary(full_model)
sum
pred <- predict(full_model, type = "response", se.fit = TRUE)

meandf$fit <- pred$fit
meandf$se.fit <- pred$se.fit

#upper and lower 95% confidence limits
meandf$upper <- meandf$fit + 1.96 * meandf$se.fit
meandf$lower <- meandf$fit - 1.96 * meandf$se.fit

#back-transform predictions to original (-1 to 1) scale
meandf$fit_unscaled   <- meandf$fit   * 2 - 1
meandf$upper_unscaled <- meandf$upper * 2 - 1
meandf$lower_unscaled <- meandf$lower * 2 - 1


r2_cv_report <- cor(meandf$scaled_pref, meandf$fit)^2
r2_cv_report
r_squared <- round(r2_cv_report, 2)

annotation_text1 <- "italic('p') < 0.001"   # now a character
annotation_text2 <- paste0("R^2 == ", signif(r_squared, 3))


# Now plot
p1 <- ggplot(meandf, aes(mean_Invaded, mean_pref, label=Pop)) +
  geom_ribbon(aes(ymin = lower_unscaled, ymax = upper_unscaled), fill = "grey80", alpha = 0.5) +  # CI band
  geom_line(aes(y = fit_unscaled), color = "black", linetype = "dashed", size = 0.8) +
  #geom_point(aes(shape = factor(batch)), 
  #           fill = "white", size = 5, stroke = 0.25, color = "black", show.legend = FALSE) +
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
    x = "Human-specialist ancestry",
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

cairo_pdf(
  file = paste0(fig1path, "pref_genomewide_aaa_labelsonlyloocv_gamplot.pdf"),
  width = 7,
  height = 7,
  family = "Arial" )
print(p1)
dev.off()




###### MAP for preference ####
geo_data=read.csv(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental/SahelII_coordinates.csv")
behave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/all_behavior_cis.txt", header=T)
merged=merge(geo_data, behave, by.x ="Location", by.y="Population")

#adding batch 
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
merged$batch=as.numeric(ifelse(merged$Location %in% dat$Pop, "1","0"))
merged$batch=as.numeric(ifelse(merged$Location =="ZIK", "0",merged$batch))
merged$batch=as.numeric(ifelse(merged$Location =="CPV", "0",merged$batch))

###creating map

#geographic extent of our data
max_lat <- 40
min_lat <- -45
max_lon <- 50
min_lon <- -30

#boundaries in a single extent object
geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))

#Download data with geodata's world function to use for our base map
world_map <- world(resolution = 3,
                   path = "data/")


my_map <- crop(x = world_map, y = geographic_extent)
pdf("C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig1/Preferences_Map_Plotv2.pdf", width = 7, height = 7)
plot(my_map, axes = T, col = "grey95",xlim = c(-30, 50),ylim=c(-33,30))
#Adding points 

blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
colors <- blue_to_red_palette(100)[as.numeric(cut(merged$pref, breaks = seq(-1, 1, length.out = 101)))]
colors

merged$jitter_x=merged$Longitude
merged$jitter_y=merged$Latitude  

library(dplyr)
merged <- merged %>%
  mutate(jitter_x = ifelse(Location =="STL" | Location =="NGO"| Location =="MIN", jitter_x - 4, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "THI", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "KAN", jitter_x - 1, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="OHI", jitter_x + 1, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="SYL", jitter_y - 1, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="STL"|Location =="KAN"|Location =="ZDR"|Location =="MAR"|Location =="DRI"|Location =="LPF"|Location =="SEV" |Location =="OHI", jitter_y + 2, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="NGO"|Location =="PKT"|Location =="KED"|Location =="BTT"|Location =="KUM"|Location =="KIN"|Location =="BOA"|Location =="KAD"|Location =="LPV"|Location =="RAB"|Location =="ZIK", jitter_y - 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="MIN"|Location =="ABK", jitter_y - 5, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="KAN"|Location =="MAR"|Location =="LBV"|Location =="ENT", jitter_x - 2, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="KED"|Location =="BOA"|Location =="SHM"|Location =="ABK", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="BTT"|Location =="KIN"|Location =="KBO"|Location =="GND", jitter_x + 3, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="BMK"|Location =="SNM"|Location =="KAK"|Location =="KWA"|Location =="GND"|Location =="KAY", jitter_y + 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="SEV", jitter_y + 3, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="DRI"|Location =="OGD"|Location =="BNF"|Location =="FCV"|Location =="VMB", jitter_x + 2, jitter_x))




jitter_x=merged$jitter_x
jitter_y=merged$jitter_y

#pies
for (i in 1:nrow(merged)) {
  jitter_x=merged$jitter_x[i]
  jitter_y=merged$jitter_y[i]
  
  point_shape <- ifelse(merged$batch[i] == 0, 21, 24) 
  
  segments(
    x0 = merged$Longitude[i], y0 = merged$Latitude[i], 
    x1 = jitter_x, y1 = jitter_y, 
    col = "black",lwd=1
  )
  points(
    x = jitter_x, y = jitter_y,  
    col = "black", 
    pch = point_shape,
    bg = colors[i],
    cex = 3.2
  )
  
  text(
    x = jitter_x, 
    y = jitter_y,  
    labels = merged$Location[i], 
    cex = 0.6, 
    col = "black"
  )
}


x_left <- -29
x_right <- -27
y_bottom <- -30
y_top <- 5
n_colors <- 100

color_steps <- blue_to_red_palette(n_colors)
y_steps <- seq(y_bottom, y_top, length.out = n_colors + 1)

for (j in 1:n_colors) {
  rect(
    xleft = x_left, xright = x_right,
    ybottom = y_steps[j], ytop = y_steps[j + 1],
    col = color_steps[j], border = NA
  )
}

rect(x_left, y_bottom, x_right, y_top, border = "black", lwd = 1)

labels <- c(-1, -0.5, 0, 0.5, 1)
label_y <- y_bottom + (labels + 1) / 2 * (y_top - y_bottom)

text(
  x = x_right + 1.5,
  y = label_y,
  labels = labels,
  cex = 0.7,
  adj = 0
)

segments(x_right, label_y, x_right + 1, label_y, lwd = 1)

text(
  x = (x_left + x_right) / 2,
  y = y_top + 1.5,
  labels = "Preference",
  cex = 0.8,
  font = 2
)


dev.off()

#########Preferences Plot ###############


##Plot pref with error bars
#Add a country grouping 
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

obs_sf <- st_as_sf(merged, coords = c("Longitude", "Latitude"), crs = 4326)

#spatial join to match points to countries
obs_with_country <- st_join(obs_sf, world[, c("admin", "iso_a3")])
#colors
blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
colors <- blue_to_red_palette(100)[as.numeric(cut(merged$pref, breaks = seq(-1, 1, length.out = 101)))]
colors


location_order <- obs_with_country %>%
  group_by(Location) %>%
  summarise(mean_pref_loc = mean(pref, na.rm = TRUE)) %>%
  arrange(mean_pref_loc) %>%  
  pull(Location)

obs_with_country$Location <- factor(obs_with_country$Location, levels = location_order)

coords <- st_coordinates(obs_with_country$geometry)

obs_with_country$Longitude <- coords[, "X"]
obs_with_country$Latitude  <- coords[, "Y"]

if("Longitude" %in% colnames(obs_with_country)){
  admin_order <- obs_with_country %>%
    group_by(admin) %>%
    summarise(admin_long = mean(Longitude, na.rm = TRUE)) %>%
    arrange(admin_long) %>%
    pull(admin)
  
  obs_with_country$admin <- factor(obs_with_country$admin, levels = admin_order)
}



p1=ggplot(na.omit(obs_with_country), aes(x=Location, y=pref, fill=pref)) +
  geom_point(aes(fill=as.numeric(pref), shape= factor(batch)), size=6.5, stroke=.25, show.legend = F, color="black") +
  scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  
  geom_errorbar(aes(ymin=dwnf, ymax=upf), width=.1) + 
  scale_fill_gradientn(colors = blue_to_red_palette(100)) +  
  xlab("") +
  ylab("Preference index") +
  ylim(-1, 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5, size = 10),
    text = element_text(size = 15),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  ) +
  facet_nested(~admin, nest_line = element_line(linetype = 1), scales="free_x", space="free") +
  theme(
    strip.background = element_blank(),
    ggh4x.facet.nestline = element_line(colour = "black")
  )

pdf(paste0(fig1path,"PreferencesOnlyPlot.pdf"), width = 10, height = 3)  # adjust width and height as needed

print(p1)
dev.off()
