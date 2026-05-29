
setwd("C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sf,raster, rgdal,terra, geodata,predicts,grDevices, RColorBrewer,dplyr,viridis,plotrix,ggplot2,MASS, StepBeta,betareg,car,rgeos)
library(sf,raster, rgdal,terra, geodata,predicts,grDevices, RColorBrewer,dplyr,viridis,plotrix,ggplot2,MASS, StepBeta,betareg,car,rgeos)


allbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/all_behavior.txt", header=T)
allcisbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Admixture/all_behavior_cis.txt", header=T)
allcisbehave=allcisbehave %>% subset(Population!="NGO")
behaveoldNGO=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/behavior/ngoold_behavior_cis.txt", header=T)
prebehave=rbind(allcisbehave, behaveoldNGO)
colnames(prebehave)[1]="Pop"
predbehave=read.table(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/introgression/0.992percentileAafW_AafE_noprefpops_prefpredictionsgamgaus.txt")
predbehave1=data.frame(cbind(predbehave$Pop,predbehave$unscaled_pred_pref))
colnames(predbehave1)=c("Pop","pref")
predbehave1$pref=as.numeric(predbehave1$pref)
behave=bind_rows(prebehave,predbehave1)

geo_data=read.csv(file="C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental/SahelII_coordinates.csv")

merged_df <- merge(behave, geo_data, by = "Pop",by.y="Location") #Excludes T51 and ORL and RABDOM

merged_df$pref=as.numeric(merged_df$pref)
#Remove ZIK (lab colony) and CPV (dont understand its demographic history) from ecological models
merged_df=merged_df %>% subset(!(Pop %in% c("ZIK", "CPV","RABDOM")))
#Should be 64 pops
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
names(FS) 

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

buffer_columns <- c("buffer_5km", "buffer_10km", "buffer_20km", "buffer_50km", "buffer_100km")

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


FS_ag <- FS_ag %>%
  mutate(geometry = st_make_valid(geometry))
FS_ps <- FS_ps %>%
  mutate(geometry = st_make_valid(geometry))

for (buffer_col in buffer_columns) {
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(
      !!paste0("intersection_ag_", buffer_col) := {
        intersection_ag <- st_intersection(get(buffer_col), FS_ag$geometry)
        if (length(intersection_ag) > 0) as.numeric(st_area(intersection_ag)) else 0
      },
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
points <- terra::vect(coords)  # Create the vector without specifying CRS
terra::crs(points) <- "EPSG:4326"  # Set the coordinate reference system

for (i in 1:length(chelsa_files)){
  r1 <- terra::rast(chelsa_files[i])
  ext1 <- terra::ext(r1)
  print(i)
  print(chelsa_files[i])
  print(ext1)
}

chelsa_data <- rast(chelsa_files[-c(21:24)])

# Extract the climate data from CHELSA rasters
chelsa_extract <- extract(chelsa_data, points)
chelsa_extract <- chelsa_extract[, colSums(is.na(chelsa_extract)) == 0]
chelsa_extract$Pop=merged_df$Pop
alldata=merge(chelsa_extract,merged_df, by="Pop")
#adding urban
raster_data <- raster("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
alldata$dens20<-terra::extract(raster_data,alldata[, c("Long", "Lat")],buffer=20000,fun=mean)
alldata_no_geometry <- alldata[, !names(alldata) %in% c("geometry","buffer_5km","buffer_10km","buffer_20km","buffer_50km","buffer_100km")]



###################################################################################################
############################# MODEL SELECTION ######################################################
####################################################################################################

library(corrplot)
library(mgcv)

alldata$pref_scaled <- (alldata$pref + 1) / 2

#Remove NGY
alldata=subset(alldata,Pop !="NGY")
alldata=subset(alldata,Pop !="RABDOM")

# Round all proportion_ag_buffer_* columns to 6 decimal places otherwise will get variance that isnt real (e.g. 0.99999999999)
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


mod_5km = gam(pref_scaled ~ s(proportion_ag_buffer_5km, k=5) + s(Lat, Long),family = betar(), data = alldata)
mod_10km = gam(pref_scaled ~ s(proportion_ag_buffer_10km,k=5) + s(Lat, Long),family = betar(), data = alldata)
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



#GAMs, starting with all parameters, not doing double penality (see EnvModelingPref for details) and using ts
#because better for spatial data
colnames(alldata)=gsub("1981-2010","1981.2010",colnames(alldata))
#get rid of the ngd stuff because they are all 365 for all environments 
#also get rid of scd because doesnt vary 
alldata <- alldata %>% dplyr::select(-contains("ngd"))
alldata<- alldata %>% dplyr::select(-contains("scd"))

#saveRDS(alldata,file="alldata_withpred.RDS")
alldata=readRDS(file="alldata_withpred.RDS")

cor=cor(alldata[c(3:60,91:92)])
corrplot(cor(alldata[c(3:60,91:92)]), method = "circle",tl.cex=.3)

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
                    s(CHELSA_vpd_range_1981.2010_V.2.1,bs="ts")+ s(dens20,bs="ts") + s(proportion_ag_buffer_50km, bs="ts"),family = betar(),select=F,
                  data = alldata, method ="REML")

summary(full_model)



full_model <- gam(pref_scaled ~ 
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                     s(CHELSA_pet_penman_range_1981.2010_V.2.1)+
                                    s(CHELSA_bio4_1981.2010_V.2.1,bs="ts") + 
                                        s(CHELSA_npp_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_rsds_1981.2010_range_V.2.1)+
                    s(dens20,bs="ts"),family = betar(),select=F,
                  data = alldata, method ="REML")
summary(full_model)
AIC(full_model)


sel_model <- gam(pref_scaled ~ 
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_vpd_max_1981.2010_V.2.1,bs="ts") + 
                    s(dens20,bs="ts") ,family = betar(),select=F,
                  data = alldata, method ="REML")
#^This one is the best
summary(sel_model)
AIC(sel_model)

alldata$fit=sel_model$fitted.values


r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.80
r_squared <- round(r2_cv_report, 2)

# edf Ref.df Chi.sq  p-value    
# s(CHELSA_bio15_1981.2010_V.2.1)   1.3136      9  53.97  < 2e-16 ***
#   s(CHELSA_vpd_max_1981.2010_V.2.1) 4.5813      9  33.55 4.12e-07 ***
#   s(dens20)                         0.9816      9  18.61 7.30e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.768   Deviance explained = 80.5%
# -REML = -35.018  Scale est. = 1         n = 63
# > AIC(sel_model)
# [1] -86.42245

# s(CHELSA_bio15_1981.2010_V.2.1)   2.2247      9  46.86  < 2e-16 ***
#   s(CHELSA_vpd_max_1981.2010_V.2.1) 5.5619      9  25.41 2.99e-05 ***
#   s(dens20)                         0.9841      9  17.21 1.16e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.774   Deviance explained = 82.1%
# NCV = -47.101  Scale est. = 1         n = 63
# > AIC(sel_model)
# [1] -84.86992


#Trying out the full model on this data
full_model=readRDS(file="GAMbetar_EcoPref_Model.RDS")
summary(full_model)
AIC(full_model)

alldata$fit <- predict(full_model, newdata = alldata, type = "response")



r2_cv_report <- cor(alldata$pref_scaled, alldata$fit)^2
r2_cv_report
#0.74
r_squared <- round(r2_cv_report, 2)

####Just senegal
SENplots=list()
SENalldata=subset(alldata, Country=="Senegal")
#full model
r2_cv_report <- cor(SENalldata$pref_scaled, SENalldata$fit)^2
r2_cv_report
#0.567
r_squared <- round(r2_cv_report, 2)
annotation_text1 <- expression(italic("p") ~ "< 0.001")
annotation_text2 <- bquote( R^2 * "=" * .(r_squared))
p1 <- ggplot(SENalldata, aes((pref_scaled*2)-1, (fit*2)-1, label=Pop)) +
  geom_segment(aes(x = -1, y = -1, xend = 1, yend = 1),
               linetype = "dashed", color = "red") +
geom_label_repel( 
    size = 2.5, 
    min.segment.length = 0, 
    box.padding = 0.0, 
    label.padding = 0.15, 
    point.padding = 0.0, 
    fill = "white", # box outline color 
    label.size = 0.01, # thickness of outline 
    segment.size = 0.35, 
    segment.color = "black", 
    max.overlaps = Inf)+
 # scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  labs(
    x = "Real or Geno Predicted Preference index",
    y = "Eco Predicted Preference index"
  ) +
  annotate("text", x = -.85, y = .9, label = "italic('p') < 0.001", parse = TRUE, size = 6) +
  annotate("text", x = -.85, y = .8, 
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
  file = paste0("./", "SenEcoComp_fullmodel.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p1)
dev.off()
SENplots[[1]]=p1
#
# vs new pred model
SENalldata$newfit=predict(sel_model, newdata = SENalldata, type = "response")
r2_cv_report <- cor(SENalldata$pref_scaled, SENalldata$newfit)^2
r2_cv_report
#0.78
r_squared <- round(r2_cv_report, 2)

p2 <- ggplot(SENalldata, aes((pref_scaled*2)-1, (newfit*2)-1, label=Pop)) +
  geom_segment(aes(x = -1, y = -1, xend = 1, yend = 1),
               linetype = "dashed", color = "red") +
  geom_label_repel( 
    size = 2.5, 
    min.segment.length = 0, 
    box.padding = 0.0, 
    label.padding = 0.15, 
    point.padding = 0.0, 
    fill = "white", # box outline color 
    label.size = 0.01, # thickness of outline 
    segment.size = 0.35, 
    segment.color = "black", 
    max.overlaps = Inf)+
  # scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  labs(
    x = "Real or Geno Predicted Preference index",
    y = "Eco Predicted Preference index"
  ) +
  annotate("text", x = -.85, y = .9, label = "italic('p') < 0.001", parse = TRUE, size = 6) +
  annotate("text", x = -.85, y = .8, 
           label = paste0("R^2 == ", signif(r_squared, 3)), parse = TRUE, size = 6) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  ) +
  ylim(-1, 1)
p2
cairo_pdf(
  file = paste0("./", "SenEcoComp_vpdmodel.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p2)
dev.off()
SENplots[[2]]=p2


#What about if we just used the actual pref data
SEN2alldata=subset(SENalldata, Dataset=="ii" |Dataset=="i")
r2_cv_report <- cor(SEN2alldata$pref_scaled, SEN2alldata$fit)^2
r2_cv_report
r_squared <- round(r2_cv_report, 2)

#0.8197476
p3 <- ggplot(SEN2alldata, aes((pref_scaled*2)-1, (fit*2)-1, label=Pop)) +
  geom_segment(aes(x = -1, y = -1, xend = 1, yend = 1),
               linetype = "dashed", color = "red") +
  geom_label_repel( 
    size = 2.5, 
    min.segment.length = 0, 
    box.padding = 0.0, 
    label.padding = 0.15, 
    point.padding = 0.0, 
    fill = "white", # box outline color 
    label.size = 0.01, # thickness of outline 
    segment.size = 0.35, 
    segment.color = "black", 
    max.overlaps = Inf)+
  # scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  labs(
    x = "Real Preference index",
    y = "Eco Predicted Preference index"
  ) +
  annotate("text", x = -.85, y = .9, label = "italic('p') < 0.001", parse = TRUE, size = 6) +
  annotate("text", x = -.85, y = .8, 
           label = paste0("R^2 == ", signif(r_squared, 3)), parse = TRUE, size = 6) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  ) +
  ylim(-1, 1)
p3
cairo_pdf(
  file = paste0("./", "SenEcoComp_fullmodelRealPref.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p3)
dev.off()
SENplots[[3]]=p3


SEN2alldata=subset(SENalldata, Dataset=="ii" |Dataset=="i")
r2_cv_report <- cor(SEN2alldata$pref_scaled, SEN2alldata$newfit)^2
r2_cv_report
#0.8249372
r_squared <- round(r2_cv_report, 2)

p4 <- ggplot(SEN2alldata, aes((pref_scaled*2)-1, (newfit*2)-1, label=Pop)) +
  geom_segment(aes(x = -1, y = -1, xend = 1, yend = 1),
               linetype = "dashed", color = "red") +
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
  # scale_shape_manual(values = c("0" = 21, "1" = 24))  +
  labs(
    x = "Real Preference index",
    y = "Eco Predicted Preference index"
  ) +
  annotate("text", x = -.85, y = .9, label = "italic('p') < 0.001", parse = TRUE, size = 6) +
  annotate("text", x = -.85, y = .8, 
           label = paste0("R^2 == ", signif(r_squared, 3)), parse = TRUE, size = 6) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  ) +
  ylim(-1, 1)

p4


cairo_pdf(
  file = paste0("./", "SenEcoComp_vpdmodelRealPref.pdf"),
  width = 6,
  height = 6,
  family = "Arial" )
print(p4)
dev.off()
SENplots[[4]]=p4

library(patchwork)


for (start in seq(1, length(SENplots), by = 4)) {
  end <- min(start + 3, length(SENplots))
  page_plots <- SENplots[start:end]
  pdf_name <- paste0("SuppPlot_SEN_EcoPreds", start, "_to_", end, ".pdf")
  pdf(pdf_name, width = 12, height = 12)  # adjust width/height as needed
  print(wrap_plots(page_plots, ncol = 2))
  dev.off()
}





library(ncf)

z <- resid(sel_model)

corres <- correlog(
  x = alldata$Long, 
  y = alldata$Lat, 
  z = z,
  increment = 100,   # km bin width
  resamp = 0,
  latlon = TRUE
)


plot(corres$mean.of.class, corres$cor, type='b',
     xlab='Distance (km)', ylab='Correlation', main='Correlogram')

# find first distance where correlation <= 0 (or <= 0.1 for conservative)
idx <- which(corres$cor <= 0)[1]
radius_km <- if(!is.na(idx)) corres$mean.of.class[idx] else max(corres$mean.of.class)


library(spdep)
library(geosphere)

# Residuals from your fitted GAM
res <- resid(selected_model2)

# Compute great-circle distance matrix (km)
coords <- cbind(alldata$Long, alldata$Lat)
D <- distm(coords) / 1000

radius_km <- radius_km
nb <- dnearneigh(coords, d1 = 0, d2 = radius_km, longlat = TRUE)

# Convert to weights list (row-standardized)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Moran's I test for spatial autocorrelation
mtest <- moran.test(res, lw, zero.policy = TRUE)
mtest


#145km radius neighborhoods for CV, 

radius_km <- 145
nei <- construct_nei(alldata, radius_km)

selected_model1 <- gam(pref_scaled ~  
                         s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                         s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                         s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                         s(CHELSA_pet_penman_mean_1981.2010_V.2.1, bs="ts")+
                         s(dens20, bs="ts"), family = betar(), select=F,
                       data= alldata, method="REML")
#Chec

full_model <- gam(pref_scaled ~  
                    s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                    s(CHELSA_bio19_1981.2010_V.2.1, bs="ts") +
                    s(CHELSA_npp_1981.2010_V.2.1, bs="ts")+
                    s(CHELSA_pet_penman_mean_1981.2010_V.2.1, bs="ts")+
                    s(dens20, bs="ts"), family = betar(), select=F,
                  data= alldata, method="REML")

summary(full_model)
AIC(full_model)


# s(CHELSA_bio15_1981.2010_V.2.1)           2.2233      9 22.235 1.68e-06 ***
#   s(CHELSA_bio19_1981.2010_V.2.1)           0.7396      9  2.809 0.045887 *  
#   s(CHELSA_npp_1981.2010_V.2.1)             0.8413      9  3.827 0.022201 *  
#   s(CHELSA_pet_penman_mean_1981.2010_V.2.1) 0.8088      9  3.395 0.033171 *  
#   s(dens20)                                 0.9394      9 10.756 0.000455 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.679   Deviance explained =   72%
# -REML = -23.995  Scale est. = 1         n = 63

#Make a raster
sel_model <- gam(pref_scaled ~ 
                   s(CHELSA_bio15_1981.2010_V.2.1,bs="ts") + 
                   s(CHELSA_vpd_max_1981.2010_V.2.1,bs="ts") + 
                   s(dens20,bs="ts") ,family = betar(),select=F,
                 data = alldata, method ="REML")
#^This one is the best
summary(sel_model)
AIC(sel_model)

# Load raster layers (replace with actual file paths)
bio15 <- rast(paste0(chelsa_path,"CHELSA_bio15_1981-2010_V.2.1.tif"))
vpdmax=rast(paste0(chelsa_path,"CHELSA_vpd_max_1981-2010_V.2.1.tif"))
dens20 <- rast("gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_2pt5_min.tif")
buffered_raster <- focal(dens20, w=9, fun=mean, na.policy="omit") #w of 9 means 6 cell radius which sould be about 20km with 2.5res

#I dont want to make the chelsa variables so coarse so going to extrapolate dens20 instead
#default is bilinear which is fine for continuous data
dens20_resampled <- resample(dens20, bio15)


predictors <-c(bio15, vpdmax, dens20_resampled)
names(predictors) <- c("CHELSA_bio15_1981.2010_V.2.1", 
                       "CHELSA_vpd_max_1981.2010_V.2.1", 
                       "dens20")
# Predict values over the raster
predicted_raster <- predict(predictors, sel_model, type = "response")

writeRaster(predicted_raster, "bio15_vpd_dens20_values.tif", overwrite=TRUE)
# 
# # Reverse the scaled transformation on the predicted raster
predicted_raster_pref <- (predicted_raster * 2) - 1
writeRaster(predicted_raster_pref, "bio15_vpd_dens20_values.tif", overwrite=TRUE)



#####################################################################################
####################### Make Maps ###################################################
#####################################################################################

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
merged_df$batch=as.numeric(ifelse(merged_df$Pop %in% dat$Pop, "1","0"))
merged_df$batch=as.numeric(ifelse(merged_df$Pop %in% predbehave1$Pop, "2", merged_df$batch))
merged_df$batch=as.numeric(ifelse(merged_df$Pop == "SKU" | merged_df$Pop == "LUA" | merged_df$Pop == "TWA" |
                                  merged_df$Pop == "YAO" | merged_df$Pop == "BUN" | merged_df$Pop == "KIC" |
                                  merged_df$Pop == "KAR", "2", merged_df$batch))


SEN=c("THI","STL","NGO","NGY","VLG", "KAN","FTK","SBR","RBT","KLK","DMS","TBA","KAF","DHR","MIN","LIG","KLD","BKJ","DBE","SCS","RNR","GNS","GLB","TCD","OSG","KAN","BTT","PKT","KED")
SEN[!(SEN %in% merged_df$Pop)]

SENmerged_df=subset(merged_df, Pop %in% SEN)

# Senegal subset

merged_df=subset(merged_df, !(Pop %in% SEN))
SEN[!(SEN %in% merged_df$Pop)]

  blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
  colors <- blue_to_red_palette(100)[as.numeric(cut(merged_df$pref, breaks = seq(-1, 1, length.out = 101)))]
  colors
  
  
  
  #####Map that shows predicted preference #####
  max_lat <- 40
  min_lat <- -45
  max_lon <- 50
  min_lon <- -30
  
  geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))
  
  world_map <- world(resolution = 3,
                     path = "data/")
  
  my_map <- crop(x = world_map, y = geographic_extent)
  
  bio15_bio19_dens20model_prefvalues=raster("bio15_bio19_bio9_bio17_dens20_npp_penman_model_prefvalues.tif")
  bio15_bio19_dens20model_prefvalues=raster("bio15_vpd_dens20_values.tif")
  
  bio15_bio19_dens20model_prefvalues_spat <- terra::rast(bio15_bio19_dens20model_prefvalues)
  
  terra::summary(bio15_bio19_dens20model_prefvalues_spat)

  
  crs(bio15_bio19_dens20model_prefvalues_spat) <- crs(my_map)
  
  my_map_raster <- terra::rasterize(my_map, bio15_bio19_dens20model_prefvalues_spat)
  
  
  breaks_raster <- seq(-1, 1, length.out = 101)
  
 

  #Downloaded land vector from 
  #https://www.naturalearthdata.com/downloads/110m-physical-vectors/
  #(physical labels)
  land <- terra::vect("C:/Users/james/OneDrive - UC San Diego/Documents/SAN_DIEGO/Rose/SahelProj/code/Environmental/ne_110m_geography_regions_polys/ne_110m_geography_regions_polys.shp")
  names(land)
  sahara <- land[land$NAME == "SAHARA", ]
  plot(land, col = "lightgrey")
  plot(sahara, col = "orange", add = TRUE)
  
  library(rworldmap)
  # Download African country boundaries
  africamap <- getMap(resolution = "low")  
  africa <- vect(subset(africamap, continent == "Africa"))
  niger=vect(subset(africamap, ADMIN ==
                      "Niger"))
  # Define extent: xmin, xmax = Niger longitude bounds, ymin = southern edge, ymax = 14
  niger_ext <- ext(xmin(niger), xmax(niger), ymin(niger), 15)
  
  # Crop Niger to this extent
  niger_crop <- crop(niger, niger_ext)
  unmask=c("Senegal","Burkina Faso","Nigeria","Cameroon")
  saharaafrica <- vect(subset(africamap, ADMIN %in% unmask))
  to_remove <- rbind(niger_crop, saharaafrica)
  sahara_masked <- erase(sahara, to_remove)
  
  # Plot to check
  plot(sahara, col="orange", add=TRUE)
  plot(to_remove, col="green", add=TRUE)
  plot(sahara_masked, col="red", add=TRUE)
  
  bio15_bio19_dens20model_prefvalues_cropped_africa_bbox <- mask(bio15_bio19_dens20model_prefvalues_spat, africa)
  
  
  
  
  fig3path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig3/"
  #pdf(paste0(fig3path,"PredictedPref_outset_Map_Plot.pdf"), width = 7, height = 7)
  pdf(paste0(fig3path,"PredictedPrefx2_outset_Map_Plot.pdf"), width = 7, height = 7)
  
  #plot base:
  plot(my_map, axes = T, col = "grey95",xlim = c(-31, 48.55),ylim=c(-33.5,30))
  
  
  plot(bio15_bio19_dens20model_prefvalues_cropped_africa_bbox, 
       axis=T, alpha=0.7, col = blue_to_red_palette(100),ylim=c(-35,20),breaks = breaks_raster, add = TRUE, legend=FALSE)
  
  plot(sahara_masked, col = "#E2CA76", add = TRUE)
  
  
  # plot(bio15_bio19_dens20model_prefvalues_cropped_africa_bbox, 
  #      axis=T, alpha=0.7, col = blue_to_red_palette(100),breaks = breaks_raster, add = TRUE)
  
  ##Adding points 
  
  merged_df$jitter_x=merged_df$Long
  merged_df$jitter_y=merged_df$Lat  
  
  library(dplyr)
  merged_df <- merged_df %>%
    mutate(jitter_x = ifelse(Pop =="STL" | Pop =="NGO"| Pop =="MIN", jitter_x - 4, jitter_x))%>%
    mutate(jitter_x = ifelse(Pop == "THI", jitter_x - 3, jitter_x))%>%
    mutate(jitter_x = ifelse(Pop == "KAN", jitter_x - 1, jitter_x))%>%
    mutate(jitter_x = ifelse(Pop =="OHI", jitter_x + 1, jitter_x))%>%
    mutate(jitter_y = ifelse(Pop =="SYL" , jitter_y - 1, jitter_y))%>%
    mutate(jitter_y = ifelse(Pop =="STL"|Pop =="KAN"|Pop =="ZDR"|Pop =="MAR"|Pop =="DRI"|Pop =="LPF"|Pop =="SEV" |Pop =="OHI" |Pop =="YAO" |Pop =="KAR", jitter_y + 2, jitter_y))%>%
    mutate(jitter_y = ifelse(Pop =="NGO"|Pop =="PKT"|Pop =="KED"|Pop =="BTT"|Pop =="KUM"|Pop =="KIN"|Pop =="BOA"|Pop =="KAD"|Pop =="LPV"|Pop =="RAB" | Pop =="KIC", jitter_y - 3, jitter_y))%>%
    mutate(jitter_y = ifelse(Pop =="MIN"|Pop =="ABK", jitter_y - 5, jitter_y))%>%
    mutate(jitter_x = ifelse(Pop =="KAN"|Pop =="MAR"|Pop =="LBV"|Pop =="ENT", jitter_x - 2, jitter_x))%>%
    mutate(jitter_x = ifelse(Pop =="KED"|Pop =="BOA"|Pop =="SHM"|Pop =="ABK", jitter_x - 3, jitter_x))%>%
    mutate(jitter_x = ifelse(Pop =="BTT"|Pop =="KIN"|Pop =="KBO"|Pop =="GND", jitter_x + 3, jitter_x))%>%
    mutate(jitter_y = ifelse(Pop =="BMK"|Pop =="SNM"|Pop =="KAK"|Pop =="KWA"|Pop =="GND"|Pop =="KAY" | Pop =="BUN", jitter_y + 3, jitter_y))%>%
    mutate(jitter_y = ifelse(Pop =="SEV", jitter_y + 3, jitter_y))%>%
    mutate(jitter_x = ifelse(Pop =="DRI"|Pop =="OGD"|Pop =="BNF"|Pop =="FCV"|Pop =="VMB", jitter_x + 2, jitter_x))
  
  
  
  
  jitter_x=merged_df$jitter_x
  jitter_y=merged_df$jitter_y
  
  
  
#pies  
  for (i in 1:nrow(merged_df)) {
    jitter_x=merged_df$jitter_x[i]
    jitter_y=merged_df$jitter_y[i]
    
   point_shape <- ifelse(merged_df$batch[i] == 0, 21, ifelse(merged_df$batch[i] == 1, 24,22))  # 22 = square, 21 = circle
    #point_shape=21
    segments(
      x0 = merged_df$Long[i], y0 = merged_df$Lat[i], 
      x1 = jitter_x, y1 = jitter_y, 
      col = "black",lwd=1
    )
     points(
       x = jitter_x, y = jitter_y,  #
       col = "black", 
       pch = point_shape,
       bg = colors[i],
       cex = 3.2
     )
   
    # Add the label
    text(
      x = jitter_x, 
      y = jitter_y,  
      labels = merged_df$Pop[i], 
      cex = 0.6, 
      col = "black"
    )
  }
  
  dev.off()
  

  
  #Senegal only
    #Download data with geodata's world function to use for our base map
  world_map <- world(resolution = 3,
                     path = "data/")
  
  #Define Senegal's bounding box
  senegal_extent <- ext(-17.5, -11.0, 12.0, 17.5)
  # Crop the world map to Senegal
  senegal_map <- crop(world_map, senegal_extent)
  # Plot the result to check
  #pdf("senegal.pdf", width = 3.5, height = 3.5)
  sen <- vect(subset(africamap, ADMIN == "Senegal"))
  
  bio15_bio19_dens20model_prefvalues_cropped_sen_bbox <- mask(bio15_bio19_dens20model_prefvalues_spat, sen)
  
  
  fig3path="C:/Users/james/OneDrive - UC San Diego/Rose Lab/2023_Fifer_Sahel_Genomics/Manuscript/Fig3/"
  pdf(paste0(fig3path,"PredictedPrefx2_inset_Map_Plot.pdf"), width = 3.5, height = 3.5)
  
  plot(senegal_map)
  plot(bio15_bio19_dens20model_prefvalues_cropped_sen_bbox, 
       axis=T, alpha=0.7, col = blue_to_red_palette(100),breaks = breaks_raster, add = TRUE, legend=FALSE)
  
  
  blue_to_red_palette <- colorRampPalette(c("#7698e5", "red"))
  colors <- blue_to_red_palette(100)[as.numeric(cut(SENmerged_df$pref, breaks = seq(-1, 1, length.out = 101)))]
  colors
  
 
  
  
  SENmerged_df$jitter_x=SENmerged_df$Long
  SENmerged_df$jitter_y=SENmerged_df$Lat  
  
  library(dplyr)
   SENmerged_df <- SENmerged_df %>%
   mutate(jitter_x = ifelse(Pop == "SBR"| Pop =="VLG"| Pop=="PKT" | Pop=="FTK" , jitter_x - .5, jitter_x))%>%
     mutate(jitter_x = ifelse(Pop == "OSG" |Pop =="KAF" | Pop =="VLG" | Pop =="THI", jitter_x - .2, jitter_x)) %>%
     mutate(jitter_x = ifelse(Pop == "MIN" | Pop =="GNS" | Pop =="KAF" | Pop =="KAN", jitter_x + .2, jitter_x)) %>%
     
     mutate(jitter_x = ifelse(Pop == "KED" |Pop =="GLB", jitter_x + .5, jitter_x)) %>%
     mutate(jitter_y = ifelse(Pop == "NGY" |Pop=="KLK" |Pop=="BTT" | Pop=="NGO" |Pop=="GLB" | Pop =="THI", jitter_y + .3, jitter_y)) %>%
     mutate(jitter_y = ifelse(Pop == "SBR" | Pop == "KLD" | Pop =="KAF", jitter_y - .5, jitter_y))%>%
     mutate(jitter_y = ifelse(Pop =="SCS", jitter_y - .75, jitter_y))%>%
     mutate(jitter_y = ifelse(Pop == "RBT" | Pop =="OSG" | Pop =="LIG" | Pop =="NGY" | Pop =="TCD", jitter_y + .5, jitter_y))
    
  # 
  
  jitter_x=SENmerged_df$jitter_x
  jitter_y=SENmerged_df$jitter_y
  
  
  
  
  for (i in 1:nrow(SENmerged_df)) {
    jitter_x=SENmerged_df$jitter_x[i]
    jitter_y=SENmerged_df$jitter_y[i]
    
    point_shape <- ifelse(SENmerged_df$batch[i] == 0, 21, ifelse(SENmerged_df$batch[i] == 1, 24,22))  # 22 = square, 21 = circle
    
    segments(
      x0 = SENmerged_df$Long[i], y0 = SENmerged_df$Lat[i], 
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
      labels = SENmerged_df$Pop[i], 
      cex = 0.6, 
      col = "black"
    )
  }
  
  dev.off()


  


