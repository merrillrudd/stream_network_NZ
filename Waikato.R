rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_network_NZ"

res_dir <- file.path(main_dir, "Waikato")
dir.create(res_dir, showWarnings=FALSE)

data_dir <- file.path(main_dir, "data_save")

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

R_dir <- file.path(res_dir, "R")
listR <- list.files(R_dir)
ignore <- sapply(1:length(listR), function(x) source(file.path(R_dir, listR[x])))

devtools::install_github("james-thorson/FishStatsUtils", ref="development")
library(VAST)

devtools::install_github("james-thorson/VAST", ref="development")
library(VAST)

library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)

#########################
## read in data
##########################

data <- load(file.path(data_dir, "nz_waikato_longfin_eel_downstream.rda"))
# data2 <- load(file.path(data_dir, "nz_waikato_longfin_eel_downstream_counts.rda"))

net_ex <- nz_waikato_longfin_eel_downstream[["network"]]
net_ex$width[which(is.na(net_ex$width))] <- median(net_ex$width, na.rm=TRUE)
net_ex2 <- net_ex %>% mutate("area_s" = width * length) %>% select('parent_s', 'child_s','dist_s', 'long','lat','area_s') %>% rename('Lat'=lat) %>% rename("Lon"=long)
net_ex2$Lat[which(net_ex2$parent_s == 0)] <- net_ex2$Lat[which(net_ex2$parent_s == 0)] + 0.00001
net_ex2$Long[which(net_ex2$parent_s == 0)] <- net_ex2$Long[which(net_ex2$parent_s == 0)] + 0.00001

obs_ex <- nz_waikato_longfin_eel_downstream[["observations"]] %>% 
          filter(data_type == "encounter") %>%
          select('child_i','fishmethod','data_value','width','year','easting','northing','agency','long','lat','dist_i') %>%
          rename('Lat'=lat) %>% rename("Lon"=long)
obs_ex$width[which(is.na(obs_ex$width))] <- median(obs_ex$width, na.rm=TRUE)
 
stream_network_eel_example <- list()
stream_network_eel_example[['network']] <- net_ex2
stream_network_eel_example[['observations']] <- obs_ex
save(stream_network_eel_example, file=file.path(data_dir, 'stream_network_eel_example.rda'))

# data <- data("stream_network_eel_example", package="FishStatsUtils")
# load(file.path("C:\\merrill\\FishStatsUtils\\examples\\nz_waikato_longfin_eel_downstream.rda"))
# stream_network_eel_example <- nz_waikato_longfin_eel_downstream
# save(stream_network_eel_example, file="C:\\merrill\\FishStatsUtils\\examples\\stream_network_eel_example.rda")
load(file.path("C:\\merrill\\FishStatsUtils\\examples\\stream_network_eel_example.rda"))

#########################
## Network
#########################

## downstream segments from all data available
network <- nz_waikato_longfin_eel_downstream[["network"]]

## add a small value to the latitude and longitude where we don't have updated location for the nodes emptying into the ocean
network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## remove NAs from width
network$width[which(is.na(network$width))] <- median(network$width, na.rm=TRUE)
network$area_s = network$width * network$length

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'area_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat) %>%
  rename('dist_s'=area_s)
Network_sz_LL$dist_s[which(Network_sz_LL$parent_s==0)] <- Inf

###########################
## Observations
###########################

## all data
obs <- nz_waikato_longfin_eel_downstream[["observations"]] %>%
	select(child_i, fishmethod, data_value, easting, northing, year, agency, data_type, pass, source, long, lat, width, dist_i)
obs$width[which(is.na(obs$width))] <- median(obs$width, na.rm=TRUE)
obs$area_i <- obs$width * obs$dist_i

## count only
obs_count <- obs %>% filter(data_type == "count")
obs_count$method_agency <- sapply(1:nrow(obs_count), function(x) paste0(obs_count[x,"fishmethod"], "_", obs_count[x,"agency"]))
obs_count$data_value <- as.integer(obs_count$data_value)

## count data turned into encounters
obs_count_enc <- obs_count
obs_count_enc$data_value <- sapply(1:nrow(obs_count_enc), function(x) ifelse(obs_count_enc$data_value[x]==0, 0, 1))
data_value <- obs_count_enc$data_value

set.seed(123)
devs <- rnorm(length(data_value), 0, 0.01)
data_value_new <- sapply(1:length(data_value), function(x) ifelse(data_value[x]==1, data_value[x]+devs[x], data_value[x]))
obs_count_enc$data_value <- data_value_new

## NZFFD data
obs_enc <- obs %>% filter(data_type == "encounter")
##### add small value to encounter observations
data_value <- obs_enc$data_value

set.seed(456)
devs <- rnorm(length(data_value), 0, 0.01)
data_value_new <- sapply(1:length(data_value), function(x) ifelse(data_value[x]==1, data_value[x]+devs[x], data_value[x]))
obs_enc$data_value <- data_value_new
obs_enc$method_agency <- sapply(1:nrow(obs_enc), function(x) paste0(obs_enc[x,"fishmethod"], "_", obs_enc[x,"agency"]))

# obs_all_count <- rbind.data.frame(obs_enc, obs_count)
## all encounter data - two datasets
obs_all_enc <- rbind.data.frame(obs_enc, obs_count_enc)

##### setup data frame
# Data_Geostat_all <- data.frame( "Catch_KG" = obs_all_count$data_value, 
#               "Year" = as.numeric(obs_all_count$year),
#                "Vessel" = obs_all_count$fishmethod,
#                "Vessel2" = obs_all_count$agency,
#                "Vessel3" = obs_all_count$method_agency, 
#                "AreaSwept_km2" = obs_all_count$dist_i, 
#                "Lat" = obs_all_count$lat, 
#                "Lon" = obs_all_count$long, 
#                "Pass" = obs_all_count$pass,
#                "Knot" = obs_all_count$child_i,
#                "Category" = obs_all_count$data_type,
#                "Vessel_Category" = paste0(obs_all_count$fishmethod,"_",obs_all_count$data_type))
# Data_Geostat_all$Data_type_num <- sapply(1:nrow(Data_Geostat_all), function(x) ifelse(Data_Geostat_all$Category[x]=="encounter",0,1))

## all encounter data- two datasets
Data_Geostat_all_enc <- data.frame( "Catch_KG" = obs_all_enc$data_value, 
              "Year" = as.numeric(obs_all_enc$year),
               "Vessel" = obs_all_enc$fishmethod,
               "Vessel2" = obs_all_enc$agency,
               "Vessel3" = obs_all_enc$method_agency, 
               "AreaSwept_km2" = obs_all_enc$area_i, 
               "Lat" = obs_all_enc$lat, 
               "Lon" = obs_all_enc$long, 
               "Pass" = obs_all_enc$pass,
               "Knot" = obs_all_enc$child_i,
               "Category" = obs_all_enc$data_type,
               "Vessel_Category" = paste0(obs_all_enc$fishmethod,"_",obs_all_enc$data_type))
Data_Geostat_all_enc$Data_type_num <- sapply(1:nrow(Data_Geostat_all_enc), function(x) ifelse(Data_Geostat_all_enc$Category[x]=="encounter",0,1))
Data_Geostat_all_enc$Vessel4 <- as.character(Data_Geostat_all_enc$Vessel3)
Data_Geostat_all_enc$Vessel4[which(grepl("fish&game",Data_Geostat_all_enc$Vessel3))] <- "Other"
Data_Geostat_all_enc$Vessel4[which(grepl("cawthron",Data_Geostat_all_enc$Vessel3))] <- "Other"
Data_Geostat_all_enc$Vessel4[which(grepl("Angling",Data_Geostat_all_enc$Vessel3))] <- "Other"
Data_Geostat_all_enc$Vessel4[which(grepl("boffa",Data_Geostat_all_enc$Vessel3))] <- "Other"
Data_Geostat_all_enc$Vessel4[which(grepl("Visual",Data_Geostat_all_enc$Vessel3))] <- "Other"


# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_all, byYear = TRUE, root = TRUE, FilePath = fig_dir, FileName = "Waikato_network")

## encounter data from NZFFD
Data_Geostat_enc <- data.frame( "Catch_KG" = obs_enc$data_value, 
              "Year" = as.numeric(obs_enc$year),
               "Vessel" = obs_enc$fishmethod,
               "Vessel2" = obs_enc$agency,
               "Vessel3" = obs_enc$method_agency, 
               "AreaSwept_km2" = obs_enc$area_i, 
               "Lat" = obs_enc$lat, 
               "Lon" = obs_enc$long, 
               "Pass" = obs_enc$pass,
               "Knot" = obs_enc$child_i,
               "Category" = "Longfin_eels",
               "Vessel_Category" = paste0(obs_enc$fishmethod, "_", obs_enc$data_type))
Data_Geostat_enc$Data_type_num <- 0
Data_Geostat_enc$Vessel4 <- as.character(Data_Geostat_enc$Vessel3)
Data_Geostat_enc$Vessel4[which(grepl("fish&game",Data_Geostat_enc$Vessel3))] <- "Other"
Data_Geostat_enc$Vessel4[which(grepl("cawthron",Data_Geostat_enc$Vessel3))] <- "Other"
Data_Geostat_enc$Vessel4[which(grepl("Angling",Data_Geostat_enc$Vessel3))] <- "Other"
Data_Geostat_enc$Vessel4[which(grepl("boffa",Data_Geostat_enc$Vessel3))] <- "Other"

# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_enc, byYear = TRUE, root = TRUE, FilePath = fig_dir, FileName = "Waikato_network_encounterdata")

## count data only
Data_Geostat_count <- data.frame( "Catch_KG" = obs_count$data_value, 
              "Year" = as.numeric(obs_count$year),
               "Vessel" = obs_count$fishmethod,
               "AreaSwept_km2" = obs_count$area_i, 
               "Lat" = obs_count$lat, 
               "Lon" = obs_count$long, 
               "Pass" = obs_count$pass,
               "Knot" = obs_count$child_i,
               "Category" = "Longfin_eels",
               "Data_type" = obs_count$data_type)
Data_Geostat_count$Data_type_num <- 0
# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_count, byYear = TRUE, root = TRUE, FilePath = fig_dir, FileName = "Waikato_network_countdata")

## encounters from count data
Data_Geostat_count_enc <- data.frame( "Catch_KG" = obs_count_enc$data_value, 
              "Year" = as.numeric(obs_count_enc$year),
               "Vessel" = obs_count_enc$fishmethod,
               "AreaSwept_km2" = obs_count_enc$area_i, 
               "Lat" = obs_count_enc$lat, 
               "Lon" = obs_count_enc$long, 
               "Pass" = obs_count_enc$pass,
               "Knot" = obs_count_enc$child_i,
               "Category" = "Longfin_eels",
               "Data_type" = obs_count_enc$data_type)
Data_Geostat_count_enc$Data_type_num <- 0

# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_all_enc, byYear = FALSE, byValue=FALSE, FilePath=fig_dir)
# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_all_enc, byYear = TRUE, byValue=TRUE, FilePath=fig_dir)
# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_all_enc, byYear = TRUE, byValue=FALSE, FilePath=fig_dir)


####################################
## habitat information
## treated as density covariates
###################################
## habitat data
hab <- nz_waikato_longfin_eel_downstream[['habitat']]
covar_toUse <- unique(hab$covariate) #c('MeanFlowCumecs','Dist2Coast_FromMid','loc_rnvar',"local_twarm",'DamAffected')
hab <- hab %>% filter(covariate %in% covar_toUse)

nodes <- network$child_s[order(network$child_s)]
years <- min(obs_all_enc$year):max(obs_all_enc$year)
years_count <- min(obs_count$year):max(obs_count$year)
years_enc <- min(obs_enc$year):max(obs_enc$year)
n_t <- length(years)
n_i_all <- nrow(obs_all_enc)
n_i_enc <- nrow(obs_enc)
n_i_count <- nrow(obs_count)

covar <- unique(hab$covariate)
n_x <- length(nodes)
n_p <- length(covar)

# for(i in 1:length(covar)){
#   p <- ggplot(hab %>% filter(covariate == covar[i])) +
#   geom_point(aes(x = easting, y = northing, color = value)) +
#   guides(color=guide_legend(title=covar[i])) +
#   scale_color_viridis_c() +
#   mytheme()
#   ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p, width=10, height=8)
# }

# hab2 <- hab %>% select(child_s, covariate, value) %>% spread(key = covariate, value = value)
# png(file.path(fig_dir, "Habitat_pairs.png"), width = 8, height=6, units = "in", res=200)
# pairs(hab2)
# dev.off()

X_gtp_input1 <- array(0, dim=c(n_x, n_t, n_p))
for(p in 1:n_p){
  psub <- hab %>% filter(covariate == covar[p])
  mat <- matrix(0, nrow=n_x, ncol = 1)
  mat[psub$child_s,1] <- psub$value
  if(covar[p]=="DamAffected"){
    X_gtp_input1[,,p] <- mat
  } else {
      mat_sd <- (mat - mean(mat, na.rm=TRUE))/sd(mat, na.rm=TRUE)
      X_gtp_input1[,,p] <- mat_sd
  }
}

## years since dam impact
X_choose <- X_gtp_input1[,,which(covar == "DamAffected")]
X_gtp1 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- years[x] - 1935
  return(sub)
})
X_gtp1_sd <- (X_gtp1 - mean(X_gtp1))/sd(X_gtp1)

## years since impact squared
X_gtp2 <- sapply(1:length(years), function(x){
  sub <- X_choose[,x]
  sub[which(sub == 1)] <- (years[x] - 1935)^2
  return(sub)
})
X_gtp2_sd <- (X_gtp2 - mean(X_gtp2))/sd(X_gtp2)

covar2 <- c(covar, "YearsSinceDam","YearsSinceDam2")[-which(covar=="DamAffected")]
n_p <- length(covar2)
X_gtp_input <- array(0, dim=c(n_x,n_t,n_p))
for(p in 1:(n_p)){
  ## skip dam affected
  if(p < length(covar)) X_gtp_input[,,p] <- X_gtp_input1[,,p]

  ## in place of dam affected, years since dam
  if(p == length(covar)) X_gtp_input[,,p] <- X_gtp1_sd

  ## additional covariate, years since dam squared
  if(p ==length(covar)+1) X_gtp_input[,,p] <- X_gtp2_sd
}

X_gtp_all <- X_gtp_input
X_gtp_count <- X_gtp_input[,which(years %in% years_count),]
X_gtp_enc <- X_gtp_input[,which(years %in% years_enc),]

## match habitat covariates to observations
## double check the indices will match up properly
X_itp_all <- array(0, dim=c(n_i_all,n_t,n_p))
for(i in 1:n_i_all){
  for(p in 1:n_p){
    child_i <- obs_all_enc$child_i[i]
    index <- which(nodes == child_i)
    X_itp_all[i,,p] <- X_gtp_input[index,,p]
  }
}

X_itp_enc <- array(0, dim=c(n_i_enc,length(years_enc),n_p))
for(i in 1:n_i_enc){
  for(p in 1:n_p){
    child_i <- obs_enc$child_i[i]
    index <- which(nodes == child_i)
    X_itp_enc[i,,p] <- X_gtp_input[index,which(years %in% years_enc),p]
  }
}

X_itp_count <- array(0, dim=c(n_i_count,length(years_count),n_p))
for(i in 1:n_i_count){
  for(p in 1:n_p){
    child_i <- obs_count$child_i[i]
    index <- which(nodes == child_i)
    X_itp_count[i,,p] <- X_gtp_input[index,which(years %in% years_count),p]
  }
}

## design matrix representing differences in catchability between data types
# Q_ik_all <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_all[,"Category"])[,-2,drop=FALSE]
Q_ik_study <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_all_enc[,"Category"])[,-1,drop=FALSE]

## design matrix representing differences in catchability between gears
vessels <- unique(Data_Geostat_enc$Vessel)
count_vessel <- sapply(1:length(vessels), function(x) nrow(Data_Geostat_enc %>% filter(Vessel == vessels[x])))
names(count_vessel) <- vessels
count_vessel[order(count_vessel)]
count_vessel[order(count_vessel)]/sum(count_vessel)

vessels <- unique(Data_Geostat_all_enc$Vessel)
count_vessel <- sapply(1:length(vessels), function(x) nrow(Data_Geostat_all_enc %>% filter(Vessel == vessels[x])))
names(count_vessel) <- vessels
count_vessel[order(count_vessel)]
count_vessel[order(count_vessel)]/sum(count_vessel)
Q_ik_all_vessel <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_all_enc[,"Vessel_Category"])[,-3,drop=FALSE]
Q_ik_enc_vessel <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_enc[,"Vessel_Category"])[,-2,drop=FALSE]


agencies <- unique(Data_Geostat_all_enc$Vessel2)
count_agency <- sapply(1:length(agencies), function(x) nrow(Data_Geostat_all_enc %>% filter(Vessel2 == agencies[x])))
names(count_agency) <- agencies
count_agency[order(count_agency)]
count_agency[order(count_agency)]/sum(count_agency)
Q_ik_all_agency <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_all_enc[,"Vessel2"])[,-7,drop=FALSE]

agencies <- unique(Data_Geostat_enc$Vessel2)
count_agency <- sapply(1:length(agencies), function(x) nrow(Data_Geostat_enc %>% filter(Vessel2 == agencies[x])))
names(count_agency) <- agencies
count_agency[order(count_agency)]
count_agency[order(count_agency)]/sum(count_agency)
Q_ik_enc_agency <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_enc[,"Vessel2"])[,-7,drop=FALSE]

va <- unique(Data_Geostat_enc$Vessel3)
count_va <- sapply(1:length(va), function(x) nrow(Data_Geostat_enc %>% filter(Vessel3 == va[x])))
names(count_va) <- va
count_va[order(count_va)]
count_va[order(count_va)]/sum(count_va)
Q_ik_enc_va <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_enc[,"Vessel4"])[,-6,drop=FALSE]

va <- unique(Data_Geostat_all_enc$Vessel3)
count_va <- sapply(1:length(va), function(x) nrow(Data_Geostat_all_enc %>% filter(Vessel3 == va[x])))
names(count_va) <- va
count_va[order(count_va)]
count_va[order(count_va)]/sum(count_va)
Q_ik_all_va <- ThorsonUtilities::vector_to_design_matrix( Data_Geostat_all_enc[,"Vessel4"])[,-6,drop=FALSE]

##################################
## save data used for model runs
##################################

# saveRDS(obs_all_count, file.path(res_dir, "observations_count.rds"))
saveRDS(obs_all_enc, file.path(res_dir, "observations_enc.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab, file.path(res_dir, "habitat.rds"))

##################################
## Models
##################################

###########################
## NZFFD_SST_AllGears_NoHab
###########################
path <- file.path(res_dir, "NZFFD_SST_AllGears_NoHab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc
save(Data, file=file.path(path, "Data_Geostat.RData"))

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL), 
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(1,6), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel")
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_SST_AllGears_Hab
###########################
path <- file.path(res_dir, "NZFFD_SST_AllGears_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)


# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_SST_Gear_Hab
###########################
path <- file.path(res_dir, "NZFFD_SST_Gear_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_ST_Gear_Hab
###########################
path <- file.path(res_dir, "NZFFD_ST_Gear_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_SST_Agency_Hab
###########################
path <- file.path(res_dir, "NZFFD_SST_Agency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_agency, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_agency, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_agency,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_SST_VesselAgency_Hab
###########################
path <- file.path(res_dir, "NZFFD_SST_VesselAgency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## NZFFD_ST_VesselAgency_Hab
###########################
path <- file.path(res_dir, "NZFFD_ST_VesselAgency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,dim(X_gtp_enc)[3]))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)


###########################
## NZFFD_SST_VesselAgency_NoHab
###########################
path <- file.path(res_dir, "NZFFD_SST_VesselAgency_NoHab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va, 
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_enc_va,
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_enc, X_itp = X_itp_enc,
                  model_args = list(Map = Map), 
                  optimize_args = list(startpar=fit1$parameter_estimates$par, newtonsteps = 3))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(1,6), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_SST_AllGears_NoHab
###########################
path <- file.path(res_dir, "AllEncounters_SST_AllGears_NoHab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_study, 
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(NA)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_study,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_study,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(1,6), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel")
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_SST_AllGears_Hab
###########################
path <- file.path(res_dir, "AllEncounters_SST_AllGears_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_study, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(NA)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_study,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_study,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_SST_Gear_Hab
###########################
path <- file.path(res_dir, "AllEncounters_SST_Gear_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_ST_Gear_Hab
###########################
path <- file.path(res_dir, "AllEncounters_ST_Gear_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_vessel, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_vessel,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)


###########################
## AllEncounters_SST_Agency_Hab
###########################
path <- file.path(res_dir, "AllEncounters_SST_Agency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_agency, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_agency,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_agency,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 


plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_SST_VesselAgency_Hab
###########################
path <- file.path(res_dir, "AllEncounters_SST_VesselAgency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)

###########################
## AllEncounters_ST_VesselAgency_Hab
###########################
path <- file.path(res_dir, "AllEncounters_ST_VesselAgency_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_va, 
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_maps(FileName=fig, plot_set=c(1,6,11,15), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)


###########################
## AllEncounters_SST_VesselAgency_NoHab
###########################
path <- file.path(res_dir, "AllEncounters_SST_VesselAgency_NoHab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_all_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Q_ik = Q_ik_all_va, 
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  run_model = FALSE)

Map <- input_list$tmb_list$Map
Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Q_ik = Q_ik_all_va,
                  # Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 


plot_maps(FileName=fig, plot_set=c(1,6), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot, category_names = "Longfin eel", covar_names = covar2)
plot_encounter_diagnostic(Report=fit$Report, Data_Geostat=Data, DirName=fig)
plot_range_index(Sdreport=fit$parameter_estimates$SD, Report=fit$Report, TmbData=fit$data_list, Year_Set=fit$year_labels, PlotDir=fig, category_names="Longfin eels")
Q <- plot_quantile_diagnostic(Report=fit$Report, TmbData=fit$data_list, DateFile=fig)
plot_residuals(Lat_i=Data[,"Lat"], Lon_i=Data[,"Lon"], TmbData=fit$data_list, spatial_list=fit$spatial_list, Report=fit$Report, Q=Q, savedir=fig, Year_Set=fit$year_labels, Years2Include=fit$years_to_plot, Cex=1.5)


df <- data.frame("Model" = c("NZFFD_SST_AllGears_Hab", "NZFFD_SST_AllGears_NoHab","NZFFD_SST_Gear_Hab","NZFFD_SST_Agency_Hab","NZFFD_SST_VesselAgency_Hab","NZFFD_SST_VesselAgency_NoHab","NZFFD_ST_VesselAgency_Hab","NZFFD_ST_Gear_Hab"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- df$AIC - min(df$AIC)
df[order(df$AIC),]

df <- data.frame("Model" = c("AllEncounters_SST_AllGears_Hab", "AllEncounters_SST_AllGears_NoHab", "AllEncounters_SST_Gear_Hab","AllEncounters_SST_Agency_Hab","AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_VesselAgency_NoHab","AllEncounters_ST_VesselAgency_Hab","AllEncounters_ST_Gear_Hab"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- df$AIC - min(df$AIC)
df[order(df$AIC),]

res1 <- readRDS(file.path(res_dir, "AllEncounters_SST_VESSELAGENCY_HAB", "Fit.rds"))
res2 <- readRDS(file.path(res_dir, "ALLEncounters_SST_AllGears_HAB", "Fit.rds"))
res3 <- readRDS(file.path(res_dir, "AllEncounters_SST_VESSELAGENCY_NOHAB", "Fit.rds"))
res4 <- readRDS(file.path(res_dir, "AllEncounters_ST_VesselAgency_HAB","Fit.rds"))
res5 <- readRDS(file.path(res_dir, "NZFFD_SST_VESSELAGENCY_HAB", "Fit.rds"))


logk <- c(res1$parameter_estimates$par[["logkappa1"]],
  res2$parameter_estimates$par[["logkappa1"]],
  res3$parameter_estimates$par[["logkappa1"]],
  res4$parameter_estimates$par[["logkappa1"]],
  res5$parameter_estimates$par[["logkappa1"]])

eps <- c(res1$parameter_estimates$par[["L_epsilon1_z"]],
  res2$parameter_estimates$par[["L_epsilon1_z"]],
  res3$parameter_estimates$par[["L_epsilon1_z"]],
  # res$parameter_estimates$par[["L_epsilon1_z"]])#,
  res5$parameter_estimates$par[["L_epsilon1_z"]])

omg <- c(res1$parameter_estimates$par[["L_omega1_z"]],
  res2$parameter_estimates$par[["L_omega1_z"]],
  res3$parameter_estimates$par[["L_omega1_z"]],
  res4$parameter_estimates$par[["L_omega1_z"]],
  res5$parameter_estimates$par[["L_omega1_z"]])

bet <- c(res1$parameter_estimates$par[["L_beta1_z"]],
  res2$parameter_estimates$par[["L_beta1_z"]],
  res3$parameter_estimates$par[["L_beta1_z"]],
  res4$parameter_estimates$par[["L_beta1_z"]],
  res5$parameter_estimates$par[["L_beta1_z"]])

mu <- c(res1$parameter_estimates$par[["Beta_mean1_c"]],
  res2$parameter_estimates$par[["Beta_mean1_c"]],
  res3$parameter_estimates$par[["Beta_mean1_c"]],
  res4$parameter_estimates$par[["Beta_mean1_c"]],
  res5$parameter_estimates$par[["Beta_mean1_c"]])


plot_compare(choose_param=1, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_ST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=1, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_VesselAgency_NoHab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=1, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_AllGears_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=1, model_vec=c("AllEncounters_SST_VesselAgency_Hab","NZFFD_SST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")


plot_compare(choose_param=2, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_ST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=2, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_VesselAgency_NoHab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=2, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_AllGears_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=2, model_vec=c("AllEncounters_SST_VesselAgency_Hab","NZFFD_SST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")


plot_compare(choose_param=3, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_ST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=3, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_VesselAgency_NoHab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=3, model_vec=c("AllEncounters_SST_VesselAgency_Hab","AllEncounters_SST_AllGears_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

plot_compare(choose_param=3, model_vec=c("AllEncounters_SST_VesselAgency_Hab","NZFFD_SST_VesselAgency_Hab"), resdir=res_dir, figdir=fig_dir, panel="category", category_names="Longfin eels")

###########################
## CountEnc_SST_AllGears_NoHab
###########################
path <- file.path(res_dir, "CountEnc_SST_AllGears_NoHab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_count_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

## L_beta1_z going to zero
## try random walk
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 
fit2$parameter_estimates$diagnostics

## L_omega1_z going to zero
FieldConfig = c("Omega1"=0, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit3 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit3$tmb_list$Obj) 
fit3$parameter_estimates$diagnostics

## L_epsilon1_z going to zero
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit4 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit4$tmb_list$Obj) 
fit4$parameter_estimates$diagnostics


fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(1,6), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot)

###########################
## CountEnc_SST_AllGears_Hab
###########################
path <- file.path(res_dir, "CountEnc_SST_AllGears_Hab")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_count_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

Xconfig_zcp_inp <- array(1, dim=c(2,1,length(covar2)))
Xconfig_zcp_inp[2,,] <- 0

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_count, X_itp = X_itp_count,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_count, X_itp = X_itp_count,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

newpar <- fit1$parameter_estimates$par
newpar[["logkappa1"]] <- 0.7

fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_count, X_itp = X_itp_count,
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar=newpar))
check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 
fit2$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "v_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  model_args = list(Map = Map),
                  Xconfig_zcp = Xconfig_zcp_inp, X_gtp = X_gtp_all, X_itp = X_itp_all,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_maps(FileName=fig, plot_set=c(1,6,11,13), Report=fit$Report, Sdreport=fit$parameter_estimates$SD, spatial_list=fit$spatial_list, TmbData=fit$data_list, Year_Set=fit$year_labels, Years2Include = fit$years_to_plot)






#################################
## CountsOnly_Explore_PosDist7
#################################
path <- file.path(res_dir, "CountsOnly_Explore_PosDist7")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_count

## first compponent is logit link zero inflation, second component is log link expectation for remainer
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=7,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

## first compponent is logit link zero inflation, second component is log link expectation for remainer
## often not enough information to estimate first component
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=3, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=7,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

## wrapper function to set up common settings
settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(startpar = fit1$parameter_estimates$par))

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz)

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_results(fit=fit, settings=settings, plot_set=c(1:3,6:7), working_dir=fig,
  year_labels=fit$year_labels, years_to_plot=fit$years_to_plot, use_biascorr=TRUE,
  category_names="Longfin eels")   

#################################
## CountsOnly_Explore_PosDist11
#################################
path <- file.path(res_dir, "CountsOnly_Explore_PosDist11")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_count

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=11,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

## L_omega1_z going very large
FieldConfig = c("Omega1"=0, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=7,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

## L_epsilon1_z going to very large
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=7,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz)

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_results(fit=fit, settings=settings, plot_set=c(1:3,6:7), working_dir=fig,
  year_labels=fit$year_labels, years_to_plot=fit$years_to_plot, use_biascorr=TRUE,
  category_names="Longfin eels")   

df <- data.frame("Model" = c("CountsOnly_Explore_PosDist7","CountsOnly_Explore_PosDist11"))
df$AIC <- NA
for(i in 1:nrow(df)){
	res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
	df$AIC[i] <- res$parameter_estimates$AIC
}

###########################
## EncounterCountsOnly_Explore
###########################
path <- file.path(res_dir, "EncounterCountsOnly_Explore")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig, showWarnings = FALSE)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_1_0.o"), to = path)

### choose dataset 
Data <- Data_Geostat_count_enc

## start with the best model that worked when PosDist = 7
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=0, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

## L_omega1_z going very small and L_epsilon1_z going very large
## try turning on temporal structure
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=1, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

## turn off L_beta1_z
FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=0, "Epsilon2"=0)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=0, "Epsilon2"=0)
ObsModel = cbind("PosDist"=2,"Link"=0)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
Options =  c("Calculate_Range"=1, 
            "Calculate_effective_area"=1)

settings <- make_settings(n_x = nrow(Network_sz), Region = "Stream_network", FieldConfig=FieldConfig, RhoConfig=RhoConfig, OverdispersionConfig=OverdispersionConfig, Options=Options, ObsModel=ObsModel, purpose = "index", fine_scale=FALSE )
settings$Method <- "Stream_network"
settings$grid_size_km <- 1

# check estimated parameters
input_list = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir=path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  run_model = FALSE)

# first model run
fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 
fit1$parameter_estimates$diagnostics

fit = fit_model( "settings"=settings, 
                  "Lat_i"=Data[,"Lat"], 
                  "Lon_i"=Data[,"Lon"], 
                  "t_iz"=Data[,'Year'], 
                  "c_i"=rep(0,nrow(Data)), 
                  "b_i"=Data[,'Catch_KG'], 
                  "a_i"=Data[,'AreaSwept_km2'], 
                  "e_i"=Data[,"Data_type_num"], 
                  working_dir = path,
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data[,"Lat"], "Lon"=Data[,"Lon"],"child_i"=Data[,"Knot"],"Area_km2"=Data[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  optimize_args = list(startpar=fit1$parameter_estimates$par))

saveRDS(fit, file.path(path, "Fit.rds"))    

fit <- readRDS(file.path(path, "Fit.rds")) 

plot_results(fit=fit, settings=settings, plot_set=c(1,6), working_dir=fig,
  year_labels=fit$year_labels, years_to_plot=fit$years_to_plot, use_biascorr=TRUE,
  category_names="Longfin eels")   
