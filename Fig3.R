setwd("C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sf,raster, rgdal,terra, geodata,predicts,grDevices, RColorBrewer,dplyr,viridis,plotrix,ggplot2,MASS, StepBeta,betareg,car,rgeos)
library(sf,raster, rgdal,terra, geodata,predicts,grDevices, RColorBrewer,dplyr,viridis,plotrix,ggplot2,MASS)
library(corrplot)
library(mgcv)

allbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
allcisbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/all_behavior_cis.txt", header=T)
allcisbehave=allcisbehave %>% subset(Population!="NGO")
behaveoldNGO=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/behavior/ngoold_behavior_cis.txt", header=T)
behave=rbind(allcisbehave, behaveoldNGO)

colnames(behave)[1]="Pop"
geo_data=read.csv(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental/SahelII_coordinates.csv")

merged_df <- merge(behave, geo_data, by = "Pop",by.y="Location") #Excludes T51 and ORL
merged_df$pref=as.numeric(merged_df$pref)
#Remove ZIK (lab colony) and CPV (dont understand its demographic history) from ecological models
merged_df=merged_df %>% subset(!(Pop %in% c("ZIK", "CPV")))
#Should be 38 pops
merged_df <- merged_df %>%
  mutate(
    Long = Longitude,
    Lat = Latitude
  ) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

merged_df <- merged_df %>%
  st_transform(crs = 3395) %>%
  mutate(
    buffer_5km = st_buffer(geometry, dist = 5000),
    buffer_10km = st_buffer(geometry, dist = 10000),
    buffer_20km = st_buffer(geometry, dist = 20000),
    buffer_50km = st_buffer(geometry, dist = 50000),
    buffer_100km = st_buffer(geometry, dist = 100000)
  ) %>%
  st_transform(crs = 4326)


#Loading farming practice data
raster_data <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
#^Just using this to convert st
targetCRS <- st_crs(raster_data)
#FS <- st_read("FS_Lev_1/FS_Lev_1.shp") #data from 2000
FS <- st_read("FS_Lev_2/FS_Lev_2.shp") #data from 2015 #this matches pop density
head(FS)
names(FS) # check attribute columns


#FS$newsys <- sub("^[0-9]+\\.\\s*", "", FS$SYSTEM)
FS$newsys <- sub("^[0-9]+\\.\\s*", "", FS$LEV1_DESC)

FSWGS84<- st_transform(FS, targetCRS)



# Filter for agro-pastoral polygons
# Ensure geometry column remains an sf-compatible object
#2000 
# FS_ag <- FSWGS84 %>% 
#   filter(newsys == "Agro-pastoral millet/sorghum") %>%
#   mutate(geometry = st_make_valid(geometry))
# #2015
FS_ag <- FSWGS84 %>% 
  filter(LEV1_DESC == "2. Agropastoral farming system") %>%
  mutate(geometry = st_make_valid(geometry))

FS_ps<- FSWGS84 %>% 
  filter(LEV1_DESC == "8. Pastoral farming system") %>%
  mutate(geometry = st_make_valid(geometry))



library(dplyr)
library(sf)

# Transform, validate, and process all buffer columns
buffer_columns <- c("buffer_5km", "buffer_10km", "buffer_20km", "buffer_50km", "buffer_100km")

# Ensure all buffers are transformed to match FS_ag CRS and made valid
reference_objects <- list(FS_ag, FS_ps)

for (ref_obj in reference_objects) {
  merged_df[buffer_columns] <- lapply(
    merged_df[buffer_columns],
    function(x) {
      x %>%
        st_transform(crs = st_crs(ref_obj)) %>%
        st_make_valid()
    }
  )
}


# Ensure FS_ag geometries are valid
FS_ag <- FS_ag %>%
  mutate(geometry = st_make_valid(geometry))
FS_ps <- FS_ps %>%
  mutate(geometry = st_make_valid(geometry))

# Iterate through buffer columns
for (buffer_col in buffer_columns) {
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(
      # Intersection area with FS_ag
      !!paste0("intersection_ag_", buffer_col) := {
        intersection_ag <- st_intersection(get(buffer_col), FS_ag$geometry)
        if (length(intersection_ag) > 0) as.numeric(st_area(intersection_ag)) else 0
      },
      # Intersection area with FS_ps
      !!paste0("intersection_ps_", buffer_col) := {
        intersection_ps <- st_intersection(get(buffer_col), FS_ps$geometry)
        if (length(intersection_ps) > 0) as.numeric(st_area(intersection_ps)) else 0
      },
      #buffer area
      !!paste0("buffer_area_", buffer_col) := as.numeric(st_area(get(buffer_col))),
      #Proportion of FS_ag intersection
      !!paste0("proportion_ag_", buffer_col) := {
        intersection_area_ag <- get(paste0("intersection_ag_", buffer_col))
        buffer_area <- get(paste0("buffer_area_", buffer_col))
        if (buffer_area > 0) intersection_area_ag / buffer_area else 0
      },
      #Proportion of FS_ps intersection
      !!paste0("proportion_ps_", buffer_col) := {
        intersection_area_ps <- get(paste0("intersection_ps_", buffer_col))
        buffer_area <- get(paste0("buffer_area_", buffer_col))
        if (buffer_area > 0) intersection_area_ps / buffer_area else 0
      }
    ) %>%
    ungroup()
}

#adding the binary yes/no agropastoral

st_crs(merged_df) <- st_crs(FSWGS84)
FSWGS84 <- st_make_valid(FSWGS84)

# Perform spatial join to extract SYSTEM values
extracted_values <- st_join(merged_df, FSWGS84)
extracted_values
# Add SYSTEM values to obs_data
#For 2000 data
# obs_dataPref$SYSTEM <- extracted_values$SYSTEM #for 2000 data
# obs_dataPref$Agropastoral=ifelse(obs_dataPref$SYSTEM=="11. Agro-pastoral millet/sorghum", "1","0") #2000
#for 2015 data
merged_df$SYSTEM <- extracted_values$LEV1_DESC #for 2015 data
merged_df$Agropastoral=ifelse(merged_df$SYSTEM=="2. Agropastoral farming system", "1","0") #2015
merged_df$Agropastoral=as.numeric(merged_df$Agropastoral)

merged_df$Pastoral=ifelse(merged_df$SYSTEM=="8. Pastoral farming system", "1","0") #2015
merged_df$Pastoral=as.numeric(merged_df$Pastoral)
library(terra)
chelsa_path <- "chelsadata/chelsav2/GLOBAL/climatologies/1981-2010/bio/"

chelsa_files <- list.files(chelsa_path, pattern = "\\.tif$", full.names = TRUE)
chelsa_files

coords <- merged_df[, c("Long", "Lat")]
points <- terra::vect(coords)  #
terra::crs(points) <- "EPSG:4326"  #

for (i in 1:length(chelsa_files)){
  r1 <- terra::rast(chelsa_files[i])
  ext1 <- terra::ext(r1)
  print(i)
  print(chelsa_files[i])
  print(ext1)
}
#21:24 has different dimensions

chelsa_data <- rast(chelsa_files[-c(21:24)])

# Extract the climate data from CHELSA rasters
chelsa_extract <- terra::extract(chelsa_data, points)
chelsa_extract <- chelsa_extract[, colSums(is.na(chelsa_extract)) == 0]
chelsa_extract$Pop=merged_df$Pop
alldata=merge(chelsa_extract,merged_df, by="Pop")
#adding urban
raster_data <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
alldata$dens20<-terra::extract(raster_data,alldata[, c("Long", "Lat")],buffer=20000,fun=mean)
alldata_no_geometry <- alldata[, !names(alldata) %in% c("geometry","buffer_5km","buffer_10km","buffer_20km","buffer_50km","buffer_100km")]
alldata$pref_scaled <- (alldata$pref + 1) / 2





###################################################################################################
############################# MODEL SELECTION ######################################################
####################################################################################################
#GAMs, starting with all parameters, not doing double penality (see EnvModelingPref for details) and using ts
#because better for spatial data
colnames(alldata)=gsub("1981-2010","1981.2010",colnames(alldata))
#get rid of the ngd stuff because they are all 365 for all environments 
#also get rid of scd because doesnt vary 
alldata <- alldata %>% dplyr::select(-contains("ngd"))
alldata<- alldata %>% dplyr::select(-contains("scd"))

alldata <- alldata %>%
  mutate(
    across(starts_with("proportion_ag_buffer_"), 
           ~ round(., 6))
  )

alldata <- alldata %>%
  mutate(
    across(starts_with("proportion_ps_buffer_"), 
           ~ round(., 6))
  )

#saveRDS(alldata,file="alldata.RDS")
alldata=readRDS(file="alldata.RDS")


mod_5km = gam(pref_scaled ~ s(proportion_ag_buffer_5km, k=5) + s(Lat, Long),family = betar(), data = alldata)
mod_10km = gam(pref_scaled ~ s(proportion_ag_buffer_10km, k=5) + s(Lat, Long),family = betar(), data = alldata)
mod_20km = gam(pref_scaled ~ s(proportion_ag_buffer_20km) + s(Lat, Long),family = betar(), data = alldata)
mod_50km = gam(pref_scaled ~ s(proportion_ag_buffer_50km)+ s(Lat, Long) ,family = betar(), data = alldata)
mod_100km = gam(pref_scaled ~ s(proportion_ag_buffer_100km)+ s(Lat, Long) ,family = betar(), data = alldata)
mod_agbinary = gam(pref_scaled ~ Agropastoral + s(Lat, Long), family = betar(), data = alldata)

mod_5km = gam(pref_scaled ~  s(proportion_ps_buffer_5km,k=5) + s(Lat, Long),family = betar(), data = alldata)
mod_10km = gam(pref_scaled ~  s(proportion_ps_buffer_10km,k=5) +  s(Lat, Long),family = betar(), data = alldata)
mod_20km = gam(pref_scaled ~ s(proportion_ps_buffer_20km,k=5) + s(Lat, Long),family = betar(), data = alldata)
mod_50km = gam(pref_scaled ~ s(proportion_ps_buffer_50km,k=5) + s(Lat, Long) ,family = betar(), data = alldata)
mod_100km = gam(pref_scaled ~  s(proportion_ps_buffer_100km,k=5) + s(Lat, Long) ,family = betar(), data = alldata)
mod_agbinary = gam(pref_scaled ~  Pastoral+s(Lat, Long), family = betar(), data = alldata)


AIC(mod_5km, mod_agbinary, mod_20km, mod_50km, mod_100km) #can do this when select=F for gam (default)
summary(mod_5km)
summary(mod_10km)
summary(mod_20km)
summary(mod_50km)
summary(mod_100km)
summary(mod_agbinary)

#5km explains the most variance and is tied with 100km for lowest aic


cor=cor(alldata[c(3:60,75,90)])
corrplot(cor(alldata[c(3:60,75,90)]), method = "circle",tl.cex=.3)

full_model <- gam(pref_scaled ~ s(CHELSA_ai_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio1_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio10_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio11_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio12_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio13_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio14_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio16_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio17_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio18_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio19_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio2_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio3_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio4_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio5_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio6_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio7_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio8_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio9_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_cmi_max_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_cmi_mean_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_cmi_min_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_cmi_range_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gdd0_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gdd10_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gdd5_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gsl_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gsp_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_gst_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_hurs_max_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_hurs_mean_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_hurs_min_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_hurs_range_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_npp_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_pet_penman_max_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_pet_penman_mean_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_pet_penman_min_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_pet_penman_range_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_rsds_1981.2010_max_V.2.1,bs="ts") + 
                    s(CHELSA_rsds_1981.2010_mean_V.2.1,bs="ts") + 
                    s(CHELSA_rsds_1981.2010_min_V.2.1,bs="ts") + 
                    s(CHELSA_rsds_1981.2010_range_V.2.1,bs="ts") + 
                    s(CHELSA_vpd_max_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_vpd_mean_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_vpd_min_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_vpd_range_1981.2010_V.2.1,bs="ts")+ s(dens20,bs="ts") + s(proportion_ag_buffer_100km, bs="ts"),family = betar(),select=F,
                  data = alldata, method ="REML")

summary(full_model)

#Best model per AIC and R2
full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                    s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                    s(dens20, bs="ts"), family = betar(), select=F,
                  data= alldata, method="REML")
summary(full_model)
AIC(full_model)
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq  p-value    
# s(CHELSA_bio15_1981.2010_V.2.1)          1.1295      9 76.177 3.82e-06 ***
#   s(CHELSA_bio17_1981.2010_V.2.1)          0.6494      9  2.484 0.046624 *  
#   s(CHELSA_bio19_1981.2010_V.2.1)          0.8915      9  7.432 0.003767 ** 
#   s(CHELSA_bio9_1981.2010_V.2.1)           0.6155      9  1.588 0.092115 .  
# s(CHELSA_npp_1981.2010_V.2.1)            1.0414      9 14.007 0.000913 ***
#   s(CHELSA_pet_penman_max_1981.2010_V.2.1) 0.8132      9  4.224 0.031922 *  
#   s(dens20)                                0.8753      9  6.567 0.005582 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.801   Deviance explained =   84%
# -REML = -28.797  Scale est. = 1         n = 38
# > AIC(full_model)
# [1] -64.35939

alldata$fit=full_model$fitted.values


r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.8379018
r_squared <- round(r2_cv_report, 2)
saveRDS(full_model,file="GAMbetar_EcoPref_Model.RDS")
#
plot(full_model, select = 2, shade = TRUE)

p <- ggplot(alldata, aes(x = pref_scaled)) +
  geom_point(aes(y = fit), color = "grey20", alpha = 0.5, size = 1.5) +
  geom_label_repel(
    aes(y = fit, label = Pop),
    size = 2.5, min.segment.length = 0,
    box.padding = 0.0, label.padding = 0.15, point.padding = 0.0,
    fill = "white", label.size = 0.01,
    segment.size = 0.35, segment.color = "black",
    max.overlaps = Inf
  ) +
  xlab("Observed Pref") + ylab("EcoPredicted Pref")
  theme_classic()
p
cbind(alldata$Pop,alldata$pref_scaled,alldata$fit)




########## plotting partial effects + partial residuals
library(ggplot2)
library(ggrepel)
library(car)

# get smooth term names
term_names <- attr(terms(full_model), "term.labels")

#predict
term_pred <- predict(full_model, type = "terms", se.fit = TRUE)

plots=list()
for (term in term_names) {
  
  x_var <- gsub("^s\\(|,.*|\\)$", "", term)
  x_vals <- alldata[[x_var]]
  
  #partial fit of this term on logit scale
  partial_fit <- term_pred$fit[, paste0("s(", term, ")")]
  
  #logit partial residuals 
  partial_resid <- residuals(full_model, type = "working") + partial_fit
  
  df <- data.frame(
    x = x_vals,
    y = logit(alldata$pref_scaled),  #not gonna plot this y gonna plot residuals instead
    z = alldata$Pop,
    partial_fit = partial_fit,
    partial_resid = partial_resid
  )
  
  edf_val <- round(summary(full_model)$s.table[paste0("s(", term, ")"), "edf"], 2)
  
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = partial_fit), color = "red", size = 1.2) +
    geom_point(aes(y = partial_resid), color = "grey20", alpha = 0.5, size = 1.5) +
    geom_label_repel(
      aes(y = partial_resid, label = z),
      size = 2.5, min.segment.length = 0,
      box.padding = 0.0, label.padding = 0.15, point.padding = 0.0,
      fill = "white", label.size = 0.01,
      segment.size = 0.35, segment.color = "black",
      max.overlaps = Inf
    ) +
    labs(
      title = paste0(term, " \n(EDF = ", edf_val, ")"),
      y = "Partial residuals (logit scale)"
    ) +
    xlab("") +
    theme_classic()
  plots[[term]]=p
  
  
  cairo_pdf(
    file = paste0("GAM_partial_", term, "_plot.pdf"),
    width = 6, height = 6, family = "Arial"
  )
  print(p)
  dev.off()
}

library(patchwork)


for (start in seq(1, length(plots), by = 8)) {
  end <- min(start + 7, length(plots))
  page_plots <- plots[start:end]
  
  pdf_name <- paste0("AllPartialResidualEffectPlots", start, "_to_", end, ".pdf")
  pdf(pdf_name, width = 12, height = 12)  # adjust width/height as needed
  print(wrap_plots(page_plots, ncol = 3))
  dev.off()
}






###########
full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                    s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                    s(dens20, bs="ts"), family = betar(), select=F,
                  data= alldata, method="NCV")
summary(full_model)
AIC(full_model)


#Some other models I tried:
full_model <- gam(pref_scaled ~  
                   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                   s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
                   s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                   s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                   s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                   s(dens20, bs="ts"), family = betar(), select=F,
                 data= alldata, method="REML")
summary(full_model)
AIC(full_model)


full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                    s(dens20, bs="ts"), family = betar(), select=F,
                  data= alldata, method="REML")
summary(full_model)
AIC(full_model)

full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts"), family = betar(), select=F,
                  data= alldata, method="REML")
summary(full_model)
AIC(full_model)

  
full_model <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                         s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                         s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                         s(CHELSA_pet_penman_mean_1981.2010_V.2.1, bs="ts")+
                         s(dens20, bs="ts"), family = betar(), select=F,
                       data= alldata, method="REML")
summary(full_model)
AIC(full_model)

#Testing the best model from the predicted pref+actual pref dataset (see bioclimpredictedPref.R)

sel_model1 <- gam(pref_scaled ~ 
                   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                   s(CHELSA_vpd_max_1981.2010_V.2.1,bs="ts") + 
                   s(dens20,bs="ts") ,family = betar(),select=F,
                 data = alldata, method ="REML")
summary(sel_model1)
AIC(sel_model1)

alldata$fit=sel_model1$fitted.values


r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.74
r_squared <- round(r2_cv_report, 2)

#testing rose 2020 model

sel_model <- gam(pref_scaled ~ 
                   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                   s(CHELSA_bio18_1981.2010_V.2.1,bs="ts") + 
                   s(dens20,bs="ts") ,family = betar(),select=F,
                 data = alldata, method ="REML")
summary(sel_model)
AIC(sel_model)

sel_model <- gam(pref_scaled ~ 
                   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                   s(dens20,bs="ts") ,family = betar(),select=F,
                 data = alldata, method ="REML")
summary(sel_model)
AIC(sel_model)
alldata$fit=sel_model$fitted.values


r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.74
r_squared <- round(r2_cv_report, 2)


plot(alldata$CHELSA_bio19_1981.2010_V.2.1,alldata$pref_scaled,xlab="npp", ylab="Scaled Preference") 

plot(full_model)
AIC(full_model)




alldata$predicted=as.numeric(mod$fitted.values)
sum=summary(mod)
sum$r.sq
r_squared <- round(sum$r.sq, 2)
r_squared
annotation_text1 <- "italic('p') < 0.001"   # now a character
annotation_text2 <- paste0("R^2 == ", signif(r_squared, 3))
p1 <- ggplot(alldata, aes(x = pref_scaled, y = predicted, label = as.factor(Pop))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size=1) +
  geom_point(shape=21, size=5, stroke=.25, show.legend = F,color="black", fill="white") +
  #scale_fill_manual(values=c("E" = "#07539dff", "W" = "#b9d7e6ff" ))+
  # geom_label_repel(
  #      size = 2.5,
  #      min.segment.length = 0,
  #      box.padding = 0.0,
  #      label.padding = 0.15,
  #      point.padding = 0.0,
  #      fill = "white",       # box outline color
  #      label.size = 0.01,      # thickness of outline
  #      segment.size = 0.35,
  #      segment.color = "black",
  #  max.overlaps = Inf)+
labs(
  x="Preference index (real)",
  y="Preference index (predicted)") +
  annotate("text", x = 0, y = .9, label = annotation_text1, hjust = 0, vjust = 0, size = 6,parse=T)+
  annotate("text", x = 0, y = .75, label = annotation_text2, hjust = 0, vjust = 0, size = 6,parse=T)+
  
  theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 18)
  )+ylim(0,1)+xlim(0,1)
p1


#https://noamross.github.io/gams-in-r-course/reference-and-resources
concurvity(mod, full = TRUE)
concurvity(mod, full = FALSE)

#Here small p-values woudl indicate that residuals are not randomly distributed. This often means there are not enough basis functions.
gam.check(mod)

plot(mod$fitted.values, mod$y)
mod$pref_scaled
#It's often useful to plot the standard errors of a partial effect term 
#combined with the standard errors of the model intercept. This is because confidence intervals at the mean value of a variable can be very tiny, and don't reflect overall uncertainty in our model. Using the seWithMean argument adds in this uncertainty.
plot(mod, seWithMean = TRUE)




#method=NCV with no nei specified just does loocv 
##NOTE BUT IT ALSO FORCES GAUSSIAN, FAMILY=BETAR NOT AVAILABLE FOR LOOCV! THIS IS BAD SINCE PREFERENCE HAS BEEN ALTERED
#TO BE USED FOR BETAR

# full_model <- gam(pref_scaled ~  
#                     s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
#                     s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
#                     s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
#                     s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
#                     s(dens20, bs="ts")+
#                     s(proportion_ag_buffer_5km, bs="ts"), family = betar(), select=F,
#                   data= alldata, method="NCV")
# summary(full_model)

#what happens when you control for lat/lon? Answer: lose bio15,17,9
full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                    s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                                        s(dens20, bs="ts")+
                    Long*Lat, family = betar(), select=F,
                  data= alldata, method="REML")
summary(full_model)
AIC(full_model)


##CV using NCV
#  construct neighborhood structure

construct_nei <- function(df, radius_km) {
  n <- nrow(df)
  coords <- cbind(df$Long, df$Lat)
  distances <- distm(coords) / 1000  # km
  
  nei <- list(k = integer(0), m = integer(0), i = integer(0), mi = integer(0))
  
  block_start <- 1
  for (i in 1:n) {
    neighbors <- which(distances[i, ] <= radius_km & distances[i, ] > 0)
    
    if (length(neighbors) > 0) {
      nei$k <- c(nei$k, neighbors)
      nei$i <- c(nei$i, neighbors)  
      block_end <- block_start + length(neighbors) - 1
      nei$m <- c(nei$m, block_end)
      nei$mi <- c(nei$mi, block_end)
      block_start <- block_end + 1
    }
  }
  return(nei)
}

#
#deciding what radius to use 
library(gstat)
selected_model2 <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                         s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                         s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                         s(CHELSA_pet_penman_mean_1981.2010_V.2.1, bs="ts")+
                         s(dens20, bs="ts"), family = betar(), select=F,
                       data= alldata, method="REML")
summary(selected_model2)

selected_model2 <-gam(pref_scaled ~   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                                          s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                                           s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                                            s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                                                                 s(dens20, bs="ts"), family = betar(), select=F,
                                        data= alldata, method="REML")
summary(selected_model2)

library(ncf)

z <- resid(selected_model2)

corres <- correlog(
  x = alldata$Long, 
  y = alldata$Lat, 
  z = z,
  increment = 100,   # km bin width
  resamp = 0,
  latlon = TRUE
)




cairo_pdf(
  file = paste0("./", "GAMbetar_ecopref_correlogram.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
plot(corres$mean.of.class, corres$cor, type='b',
     xlab='Distance (km)', ylab='Correlation', main='Correlogram')
dev.off()

#find first distance where correlation <= 0 (or <= 0.1 for conservative)
idx <- which(corres$cor <= 0)[1]
radius_km <- if(!is.na(idx)) corres$mean.of.class[idx] else max(corres$mean.of.class)


library(spdep)
library(geosphere)

res <- resid(selected_model2)

#great-circle distance matrix (km)
coords <- cbind(alldata$Long, alldata$Lat)
D <- distm(coords) / 1000

radius_km <- radius_km
radius_km
nb <- dnearneigh(coords, d1 = 0, d2 = radius_km, longlat = TRUE)

#convrt to weights list (row-standardized)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Moran's I test for spatial autocorrelation
mtest <- moran.test(res, lw, zero.policy = TRUE)
mtest


#34km radius neighborhoods for CV, 

radius_km <- 34
nei <- construct_nei(alldata, radius_km)

selected_model1 <-gam(pref_scaled ~   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                        s(CHELSA_bio17_1981.2010_V.2.1, bs="ts") +
                        s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                        s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                        s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                        s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                        s(dens20, bs="ts"), family = betar(), select=F,
                      data= alldata, method="REML")
#Checking fit:
x=fitted(selected_model1)
alldata$fit=x
df=alldata %>% dplyr::select(Pop, fit)
summary(selected_model1)
selected_model1$smooth
sp_values <- selected_model1$sp  
#I set select=F, but in theory that should be automatic when providing sp values, not sure if this package does that automatically tho
selected_model1 <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1) + 
                         s(CHELSA_bio19_1981.2010_V.2.1) +
                         s(CHELSA_npp_1981.2010_V.2.1)+
                         s(CHELSA_pet_penman_mean_1981.2010_V.2.1)+
                         s(CHELSA_bio17_1981.2010_V.2.1)+
                         s(dens20), select=F,
                       data= alldata, method="NCV", nei=nei, sp=sp_values)


selected_model1 <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1) +
                         s(CHELSA_bio17_1981.2010_V.2.1)+
                         s(CHELSA_bio19_1981.2010_V.2.1) +
                         s(CHELSA_bio9_1981.2010_V.2.1)+
                         s(CHELSA_npp_1981.2010_V.2.1)+
                         s(CHELSA_pet_penman_max_1981.2010_V.2.1)+
                         s(dens20), select=F,
                       data= alldata, method="NCV", nei=nei, sp=sp_values)

summary(selected_model1)
AIC(selected_model1)

# approximate significance of smooth terms:
#   edf Ref.df      F  p-value    
# s(CHELSA_bio15_1981.2010_V.2.1)           1.035  1.069 13.704 0.000885 ***
#   s(CHELSA_bio19_1981.2010_V.2.1)           1.013  1.026  7.135 0.012979 *  
#   s(CHELSA_npp_1981.2010_V.2.1)             1.002  1.005 10.087 0.003337 ** 
#   s(CHELSA_pet_penman_mean_1981.2010_V.2.1) 1.005  1.009  3.607 0.066273 .  
# s(CHELSA_bio17_1981.2010_V.2.1)           1.006  1.013  2.896 0.103951    
# s(dens20)                                 1.006  1.012  5.338 0.027647 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

alldata$fit=selected_model1$fitted.values


r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.83
r_squared <- round(r2_cv_report, 2)





#npp is Net primaryproductivity, highly correlated with moisture index
#pet_penman_maxMaximummonthly potentialevapotranspiration



#Create a raster of preference predictions 
selected_model2 <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                         s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                         s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                         s(CHELSA_pet_penman_mean_1981.2010_V.2.1, bs="ts")+
                         s(dens20, bs="ts"), family = betar(), select=F,
                       data= alldata, method="REML")
summary(selected_model2)

selected_model2 <-gam(pref_scaled ~   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                        s(CHELSA_bio17_1981.2010_V.2.1, bs="ts")+
                        s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                        s(CHELSA_bio9_1981.2010_V.2.1, bs="ts")+
                        s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                        s(CHELSA_pet_penman_max_1981.2010_V.2.1, bs="ts")+
                        s(dens20, bs="ts"), family = betar(), select=F,
                      data= alldata, method="REML")


bio15 <- rast(paste0(chelsa_path,"CHELSA_bio15_1981-2010_V.2.1.tif"))
bio19 <- rast(paste0(chelsa_path,"CHELSA_bio19_1981-2010_V.2.1.tif"))
npp <- rast(paste0(chelsa_path,"CHELSA_npp_1981-2010_V.2.1.tif"))
penman=rast(paste0(chelsa_path,"CHELSA_pet_penman_max_1981-2010_V.2.1.tif"))
bio9=rast(paste0(chelsa_path,"CHELSA_bio9_1981-2010_V.2.1.tif"))
bio17=rast(paste0(chelsa_path,"CHELSA_bio17_1981-2010_V.2.1.tif"))
dens20 <- rast("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
buffered_raster <- focal(dens20, w=9, fun=mean, na.policy="omit") #w of 9 means 6 cell radius which sould be about 20km with 2.5res

#I dont want to make the chelsa variables so coarse so going to extrapolate dens20 instead
#default is bilinear which is fine for continuous data
dens20_resampled <- resample(dens20, bio15)


predictors <-c(bio15, bio19, dens20_resampled,npp,penman)
names(predictors) <- c("CHELSA_bio15_1981.2010_V.2.1", 
                       "CHELSA_bio19_1981.2010_V.2.1", 
                       "dens20","CHELSA_npp_1981.2010_V.2.1", "CHELSA_pet_penman_mean_1981.2010_V.2.1")
# Predict values over the raster
predicted_raster <- predict(predictors, selected_model2, type = "response")

predictors <-c(bio15, bio19, npp, bio9, bio17, dens20_resampled,penman)
names(predictors) <- c("CHELSA_bio15_1981.2010_V.2.1", 
                       "CHELSA_bio19_1981.2010_V.2.1", "CHELSA_npp_1981.2010_V.2.1",
                       "CHELSA_bio9_1981.2010_V.2.1","CHELSA_bio17_1981.2010_V.2.1",
                       "dens20", "CHELSA_pet_penman_max_1981.2010_V.2.1")
# Predict values over the raster
predicted_raster <- predict(predictors, selected_model2, type = "response")




# # Save the output raster
# writeRaster(predicted_raster, "bio15_bio19_dens20_npp_penman_ag_model_values.tif", overwrite=TRUE)
# 
# predicted_raster_pref <- (predicted_raster * 2) - 1
# writeRaster(predicted_raster_pref, "bio15_bio19_dens20_npp_penman_ag_model_prefvalues.tif", overwrite=TRUE)
 writeRaster(predicted_raster, "bio15_bio19_bio9_bio17_dens20_npp_penman_model_values.tif", overwrite=TRUE)
# 
 predicted_raster_pref <- (predicted_raster * 2) - 1
 writeRaster(predicted_raster_pref, "bio15_bio19_bio9_bio17_dens20_npp_penman_model_prefvalues.tif", overwrite=TRUE)

 ##############################################
 ##########Plotting populations with overlays#############
 
 ###Making a clim effect/dry season intensity index raster:
 

 # Stack climate rasters with names matching the model's variable names
 clim_stack <- c(bio15, bio17, bio19, bio9, npp, penman)
 names(clim_stack) <- c(
   "CHELSA_bio15_1981.2010_V.2.1",
   "CHELSA_bio17_1981.2010_V.2.1",
   "CHELSA_bio19_1981.2010_V.2.1",
   "CHELSA_bio9_1981.2010_V.2.1",
   "CHELSA_npp_1981.2010_V.2.1",
   "CHELSA_pet_penman_max_1981.2010_V.2.1"
 )
 
 # Add population density (already resampled to same grid)
 full_stack <- c(clim_stack, dens20_resampled)
 names(full_stack)[7] <- "dens20"
 
 predict_clim_effect <- function(model, newdata, ...) {
   terms_out <- predict(model, newdata = newdata, type = "terms", exclude = "batch")
   term_names <- paste0("s(", c(
     "CHELSA_bio15_1981.2010_V.2.1",
     "CHELSA_bio17_1981.2010_V.2.1",
     "CHELSA_bio19_1981.2010_V.2.1",
     "CHELSA_bio9_1981.2010_V.2.1",
     "CHELSA_npp_1981.2010_V.2.1",
     "CHELSA_pet_penman_max_1981.2010_V.2.1"
   ), ")")
   rowSums(terms_out[, term_names, drop=FALSE])
 }
 
 clim_effect_rast <- terra::predict(full_stack,
                                    model = selected_model2,
                                    fun = predict_clim_effect,
                                    na.rm = TRUE)
 
 writeRaster(clim_effect_rast, "clim_effect_index.tif", overwrite=TRUE)
 
 ##Map for clim_effect_rast
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
 
 ###creating map
 

 max_lat <- 40
 min_lat <- -45
 max_lon <- 50
 min_lon <- -30
 
 geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))
 
 world_map <- world(resolution = 3,
                    path = "data/")
 
 # Crop the map to our area of interest
 my_map <- crop(x = world_map, y = geographic_extent)
 #plot(clim_effect_rast)
 
 crs(my_map) <- "EPSG:4326"  
 pdf("ActualPref_DrySeasonIntensity_Map_Plotv1.pdf", width = 7, height = 7)
 
 plot(my_map, axes = T, col = "grey95",xlim = c(-30, 50),ylim=c(-33,30))
 plot(clim_effect_rast,axis=T, alpha=0.7,add=T)
 
 blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
 colors <- blue_to_red_palette(100)[as.numeric(cut(merged$pref, breaks = seq(-1, 1, length.out = 101)))]
 colors
 
 merged$jitter_x=merged$Longitude
 merged$jitter_y=merged$Latitude  
 
 merged=merged %>% subset(!(Location %in% c("ZIK", "CPV")))
 
 library(dplyr)
 merged <- merged %>%
   mutate(jitter_x = ifelse(Location =="STL" | Location =="NGO"| Location =="MIN", jitter_x - 4, jitter_x))%>%
   mutate(jitter_x = ifelse(Location == "THI", jitter_x - 3, jitter_x))%>%
   mutate(jitter_x = ifelse(Location == "KAN", jitter_x - 1, jitter_x))%>%
   mutate(jitter_x = ifelse(Location =="OHI", jitter_x + 1, jitter_x))%>%
   mutate(jitter_y = ifelse(Location =="SYL", jitter_y - 1, jitter_y))%>%
   mutate(jitter_y = ifelse(Location =="STL"|Location =="KAN"|Location =="ZDR"|Location =="MAR"|Location =="DRI"|Location =="LPF"|Location =="SEV" |Location =="OHI", jitter_y + 2, jitter_y))%>%
   mutate(jitter_y = ifelse(Location =="NGO"|Location =="PKT"|Location =="KED"|Location =="BTT"|Location =="KUM"|Location =="KIN"|Location =="BOA"|Location =="KAD"|Location =="LPV"|Location =="RAB", jitter_y - 3, jitter_y))%>%
   mutate(jitter_y = ifelse(Location =="MIN"|Location =="ABK", jitter_y - 5, jitter_y))%>%
   mutate(jitter_x = ifelse(Location =="KAN"|Location =="MAR"|Location =="LBV"|Location =="ENT", jitter_x - 2, jitter_x))%>%
   mutate(jitter_x = ifelse(Location =="KED"|Location =="BOA"|Location =="SHM"|Location =="ABK", jitter_x - 3, jitter_x))%>%
   mutate(jitter_x = ifelse(Location =="BTT"|Location =="KIN"|Location =="KBO"|Location =="GND", jitter_x + 3, jitter_x))%>%
   mutate(jitter_y = ifelse(Location =="BMK"|Location =="SNM"|Location =="KAK"|Location =="KWA"|Location =="GND"|Location =="KAY", jitter_y + 3, jitter_y))%>%
   mutate(jitter_y = ifelse(Location =="SEV", jitter_y + 3, jitter_y))%>%
   mutate(jitter_x = ifelse(Location =="DRI"|Location =="OGD"|Location =="BNF"|Location =="FCV"|Location =="VMB", jitter_x + 2, jitter_x))
 
 
 
 
 jitter_x=merged$jitter_x
 jitter_y=merged$jitter_y
 
 
 
 
 for (i in 1:nrow(merged)) {
   jitter_x=merged$jitter_x[i]
   jitter_y=merged$jitter_y[i]
   
   point_shape <- ifelse(merged$batch[i] == 0, 21, 24)  # 22 = square, 21 = circle
   
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
   
   # Add the label
   text(
     x = jitter_x, 
     y = jitter_y,  
     labels = merged$Location[i], 
     cex = 0.6, 
     col = "black"
   )
 }
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

###creating map

max_lat <- 40
min_lat <- -45
max_lon <- 50
min_lon <- -30

geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))

world_map <- world(resolution = 3,
                   path = "data/")

my_map <- crop(x = world_map, y = geographic_extent)

chelsa_path <- "chelsadata/chelsav2/GLOBAL/climatologies/1981-2010/bio/"






# Load raster layers (replace with actual file paths)
bio15 <- rast(paste0(chelsa_path,"CHELSA_bio15_1981-2010_V.2.1.tif"))
bio19 <- rast(paste0(chelsa_path,"CHELSA_bio19_1981-2010_V.2.1.tif"))
npp <- rast(paste0(chelsa_path,"CHELSA_npp_1981-2010_V.2.1.tif"))
petpenmax_max <- rast(paste0(chelsa_path,"CHELSA_pet_penman_max_1981-2010_V.2.1.tif"))

dens20 <- rast("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
#buffered_raster <- focal(dens20, w=matrix(1, nrow=8, ncol=8), fun=mean, na.policy="omit")
buffered_raster <- focal(dens20, w=9, fun=mean, na.policy="omit") #w of 9 means 6 cell radius which sould be about 20km with 2.5res
raster_data <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
merged$dens20<-terra::extract(raster_data,alldata[, c("Long", "Lat")],buffer=20000,fun=mean)


ext(bio15)
ext(dens20)

#oldres
res(bio15)
#newres
res(bio15_resampled)
res(dens20)
#I dont want to make the chelsa variables so coarse so going to extrapolate dens20 instead
#default is bilinear which is fine for continuous data
dens20_resampled <- resample(dens20, bio19)


# Plot the base map

pdf("ActualPref_bio15_Map_Plotv1.pdf", width = 7, height = 7)

#### FARMING PRACTICE MAP ####
#Alt do grey everything except for millet
# Plot the base map
color_palette<-rep("grey",length(unique(FSWGS84$newsys)))
system_colors <- setNames(color_palette, unique(FSWGS84$newsys))
system_colors["Agropastoral farming system"]<-"#9C8F55"
plot(my_map, axes = T, col = "grey95",xlim = c(-30, 50),ylim=c(-33,30))


#plot(buffered_raster,axis=T, alpha=0.7,add=T,
#col = colorRampPalette(c("black", "red", "orange"))(100))
#plot(agro_rast, axis=T, alpha=0.7, add=T)
plot(bio15,axis=T, alpha=0.7,add=T)
plot(bio19,axis=T, alpha=0.7,add=T)
plot(npp,axis=T, alpha=0.7,add=T)
plot(petpenmax_max,axis=T, alpha=0.7,add=T)

 # Overlay FS geometry only for Agropastoral
 if (any(FS$newsys == "Agropastoral farming system")) {
   plot(
     st_geometry(FS[FS$newsys == "Agropastoral farming system", ]), 
     col = system_colors[FS$newsys[FS$newsys == "Agropastoral farming system"]], 
     add = TRUE
   )
 }




 #pdf("C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig1/Preferences_Map_Plotv1.pdf", width = 7, height = 7)



blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
colors <- blue_to_red_palette(100)[as.numeric(cut(merged$pref, breaks = seq(-1, 1, length.out = 101)))]
colors

merged$jitter_x=merged$Longitude
merged$jitter_y=merged$Latitude  

merged=merged %>% subset(!(Location %in% c("ZIK", "CPV")))

library(dplyr)
merged <- merged %>%
  mutate(jitter_x = ifelse(Location =="STL" | Location =="NGO"| Location =="MIN", jitter_x - 4, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "THI", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "KAN", jitter_x - 1, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="OHI", jitter_x + 1, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="SYL", jitter_y - 1, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="STL"|Location =="KAN"|Location =="ZDR"|Location =="MAR"|Location =="DRI"|Location =="LPF"|Location =="SEV" |Location =="OHI", jitter_y + 2, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="NGO"|Location =="PKT"|Location =="KED"|Location =="BTT"|Location =="KUM"|Location =="KIN"|Location =="BOA"|Location =="KAD"|Location =="LPV"|Location =="RAB", jitter_y - 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="MIN"|Location =="ABK", jitter_y - 5, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="KAN"|Location =="MAR"|Location =="LBV"|Location =="ENT", jitter_x - 2, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="KED"|Location =="BOA"|Location =="SHM"|Location =="ABK", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="BTT"|Location =="KIN"|Location =="KBO"|Location =="GND", jitter_x + 3, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="BMK"|Location =="SNM"|Location =="KAK"|Location =="KWA"|Location =="GND"|Location =="KAY", jitter_y + 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="SEV", jitter_y + 3, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="DRI"|Location =="OGD"|Location =="BNF"|Location =="FCV"|Location =="VMB", jitter_x + 2, jitter_x))




jitter_x=merged$jitter_x
jitter_y=merged$jitter_y




# # Add pie charts with jitter based on proximity and avoid overlaps
for (i in 1:nrow(merged)) {
  jitter_x=merged$jitter_x[i]
  jitter_y=merged$jitter_y[i]
  
  point_shape <- ifelse(merged$batch[i] == 0, 21, 24)  # 22 = square, 21 = circle
  
  # Draw a line from the original position to the jittered position
  segments(
    x0 = merged$Longitude[i], y0 = merged$Latitude[i], 
    x1 = jitter_x, y1 = jitter_y, 
    col = "black",lwd=1
  )
  # points(
  #   x = jitter_x, y = jitter_y,  # Jittered points
  #   col = "black", 
  #   pch = point_shape,
  #   bg = colors[i],
  #   cex = 3.2
  # )
  points(
    x = jitter_x, y = jitter_y,  # Jittered points
    col = "black", 
    pch = point_shape,
    bg = colors[i],
    cex = scale_cex[i]
  )
  
  # Add the label
  text(
    x = jitter_x, 
    y = jitter_y,  # Slightly above the pie chart
    labels = merged$Location[i], 
    cex = 0.6, 
    col = "black"
  )
}
dev.off()

plot(my_map, axes = T, col = "grey95",xlim = c(-30, 50),ylim=c(-33,30))

scale_cex <- scales::rescale(sqrt(merged$dens20), to = c(2.5, 6))
# Define representative values to show in the legend
#legend_vals <- quantile(merged$dens20, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
legend_vals <- c(10, 100, 1000)

legend_sizes <- scales::rescale(sqrt(legend_vals), to = c(2.5, 6))

blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
colors <- blue_to_red_palette(100)[as.numeric(cut(merged$pref, breaks = seq(-1, 1, length.out = 101)))]
colors

merged$jitter_x=merged$Longitude
merged$jitter_y=merged$Latitude  

merged=merged %>% subset(!(Location %in% c("ZIK", "CPV")))

library(dplyr)
merged <- merged %>%
  mutate(jitter_x = ifelse(Location =="STL" | Location =="NGO"| Location =="MIN", jitter_x - 4, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "THI", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location == "KAN", jitter_x - 1, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="OHI", jitter_x + 1, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="SYL", jitter_y - 1, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="STL"|Location =="KAN"|Location =="ZDR"|Location =="MAR"|Location =="DRI"|Location =="LPF"|Location =="SEV" |Location =="OHI", jitter_y + 2, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="NGO"|Location =="PKT"|Location =="KED"|Location =="BTT"|Location =="KUM"|Location =="KIN"|Location =="BOA"|Location =="KAD"|Location =="LPV"|Location =="RAB", jitter_y - 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="MIN"|Location =="ABK", jitter_y - 5, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="KAN"|Location =="MAR"|Location =="LBV"|Location =="ENT", jitter_x - 2, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="KED"|Location =="BOA"|Location =="SHM"|Location =="ABK", jitter_x - 3, jitter_x))%>%
  mutate(jitter_x = ifelse(Location =="BTT"|Location =="KIN"|Location =="KBO"|Location =="GND", jitter_x + 3, jitter_x))%>%
  mutate(jitter_y = ifelse(Location =="BMK"|Location =="SNM"|Location =="KAK"|Location =="KWA"|Location =="GND"|Location =="KAY", jitter_y + 3, jitter_y))%>%
  mutate(jitter_y = ifelse(Location =="SEV", jitter_y + 3, jitter_y))%>%
  mutate(jitter_x = ifelse(Location =="DRI"|Location =="OGD"|Location =="BNF"|Location =="FCV"|Location =="VMB", jitter_x + 2, jitter_x))




jitter_x=merged$jitter_x
jitter_y=merged$jitter_y




# # Add pie charts with jitter based on proximity and avoid overlaps
for (i in 1:nrow(merged)) {
  jitter_x=merged$jitter_x[i]
  jitter_y=merged$jitter_y[i]
  
  point_shape <- ifelse(merged$batch[i] == 0, 21, 24)  # 22 = square, 21 = circle
  
  # Draw a line from the original position to the jittered position
  segments(
    x0 = merged$Longitude[i], y0 = merged$Latitude[i], 
    x1 = jitter_x, y1 = jitter_y, 
    col = "black",lwd=1
  )
  points(
    x = jitter_x, y = jitter_y,  # Jittered points
    col = "black", 
    pch = point_shape,
    bg = colors[i],
    cex = scale_cex[i]
  )
  
  # Add the label
  text(
    x = jitter_x, 
    y = jitter_y,  # Slightly above the pie chart
    labels = merged$Location[i], 
    cex = 0.6, 
    col = "black"
  )
}


usr <- par("usr")  # (xmin, xmax, ymin, ymax)
legend(
  x = usr[1] + 0.2 * (usr[2] - usr[1]),  # 20% from the left edge
  y = usr[3] + 0.5 * (usr[4] - usr[3]), # a little above bottom
  legend = legend_vals,
  title = "Human Density per km",
  pt.cex = legend_sizes,
  pch = 21,
  pt.bg = "black",
  col = "black",
  bty = "n",
  y.intersp = 2.5,
  cex = 0.8
)




