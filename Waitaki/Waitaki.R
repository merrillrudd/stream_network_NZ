rm(list=ls())

################
## Directories
################

main_dir <- "/home/merrill/stream_network_NZ"
# R_dir <- file.path(main_dir, "R")
# R_files <- list.files(R_dir)
# readr <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

res_dir <- file.path(main_dir, "Waitaki")
dir.create(res_dir, showWarnings=FALSE)

data_dir <- file.path(res_dir, "data")

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
##remember to pull upstream development branches
devtools::install_github("james-thorson/VAST", INSTALL_opts="--no-staged-install")
devtools::install_github("merrillrudd/VASTPlotUtils")

library(VAST)
library(VASTPlotUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(foreach)
library(doParallel)

#########################
## read in data
##########################

load(file.path(data_dir, "nz_waitaki_longfin_eel.rda"))
network <- nz_waitaki_longfin_eel[["network"]]

## add a small value to the latitude and longitude where we don't have updated location for the nodes emptying into the ocean
# network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
# network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- nz_waitaki_longfin_eel[["observations"]] %>%
    dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value) %>%
    mutate('method_agency' = paste0(fishmethod, "_", agency)) %>%
    rename("Year"=year) %>%
    filter(fishmethod != "Angling")

## check that segment distance is more than 0.125
# obs$dist_i <- sapply(1:nrow(obs), function(x){
#   sub_net <- Network_sz %>% filter(child_s == obs$child_i[x])
#   if(sub_net$dist_s >= obs$dist_i[x]){
#     out <- obs$dist_i[x]
#   } else { out <- sub_net$dist_s }
#   return(out)
# })

##### add small value to encounter observations
set.seed(1234)
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = as.numeric(obs$Year),
               "Method" = obs$fishmethod,
               "Agency" = obs$agency,
               "Method_Agency" = obs$method_agency, 
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = "Longfin_eels")

## habitat data
hab <- nz_waitaki_longfin_eel[['habitat']]

covar_toUse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm",'DamAffected')
hab <- hab %>% filter(covariate %in% covar_toUse)

####################################
## habitat information
## treated as density covariates
###################################

nodes <- network$child_s[order(network$child_s)]
years <- min(obs$Year):max(obs$Year)
covar <- unique(hab$covariate)
n_x <- length(nodes)
n_t <- length(years)
n_p <- length(covar)
n_i <- nrow(obs)

hab_std <- lapply(1:length(covar), function(x){
  sub <- hab %>% 
    filter(covariate == covar[x]) %>%
    mutate(value_std = (value - mean(value))/sd(value))
  return(sub)
})
hab_std <- do.call(rbind, hab_std)

# p_habstd <- ggplot(hab_std) +
#   geom_point(aes(x = easting, y = northing, color = value_std)) +
#   facet_wrap(.~covariate, nrow = 2) +
#   scale_color_distiller(palette = "Spectral") +
#   scale_x_continuous(breaks = quantile(hab$easting, prob = c(0.1, 0.5, 0.95)), labels = round(quantile(hab$easting, prob = c(0.1,0.5,0.95)),0)) +
#   guides(color = guide_legend(title = "Standardised\nvalue")) +
#   theme_minimal() 
# ggsave(file.path(fig_dir, "Habitat_covariates_standardised.png"), p_habstd, width = 15, height = 8)

# for(i in 1:length(covar)){
#   p <- ggplot(hab_std %>% filter(covariate == covar[i])) +
#   geom_point(aes(x = easting, y = northing, color = value)) +
#   guides(color=guide_legend(title=covar[i])) +
#   scale_color_distiller(palette = "Spectral") +
#   theme_minimal()
#   ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p, width=10, height=8)
# }

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

## match habitat covariates to observations
## double check the indices will match up properly
X_itp_input <- array(0, dim=c(n_i,n_t,n_p))
for(i in 1:n_i){
  for(p in 1:n_p){
    child_i <- obs$child_i[i]
    index <- which(nodes == child_i)
    X_itp_input[i,,p] <- X_gtp_input[index,,p]
  }
}



####################################
## sampling information
## treated as catchability covariates
###################################

method_info <- Data_Geostat %>%
  group_by(Method) %>%
  summarise(Samples = length(Method), 
            Prop_samples = length(Method)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method))
method_info <- method_info[rev(order(method_info$Samples)),]
# Data_Geostat$Method2 <- as.character(Data_Geostat$Method)
# Data_Geostat$Method2[which(Data_Geostat$Method %in% c("Trap","Visual","Angling"))] <- "Trap_Visual_Angling"
Q_ik_method <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method"])[,-1,drop=FALSE]


agency_info <- Data_Geostat %>%
  group_by(Agency) %>%
  summarise(Samples = length(Agency), 
            Prop_samples = length(Agency)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Agency))
agency_info <- agency_info[rev(order(agency_info$Samples)),]
# Data_Geostat$Agency2 <- as.character(Data_Geostat$Agency)
# Data_Geostat$Agency2[which(Data_Geostat$Agency %in% c("fish&game","consultants"))] <- "fish&game_consultants"
Q_ik_agency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Agency"])[,-3,drop=FALSE]


method_agency_info <- Data_Geostat %>%
  group_by(Method_Agency) %>%
  summarise(Samples = length(Method_Agency), 
            Prop_samples = length(Method_Agency)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method_Agency))
method_agency_info <- method_agency_info[rev(order(method_agency_info$Samples)),]
# Data_Geostat$Method_Agency2 <- sapply(1:nrow(Data_Geostat), function(x){
#   if(grepl("Electric", Data_Geostat$Method_Agency[x]) == FALSE){
#     out <- as.character(Data_Geostat$Method[x])
#   } else{
#     out <- as.character(Data_Geostat$Method_Agency[x])
#   }
#   if(out %in% c("Trap","Visual","Angling")) out <- "Method_other"
#   if(grepl("fish&game",out) | grepl("consultants",out)) out <- "Electric fishing_agency_other"
#   return(out)
# })
# method_agency_info2 <- Data_Geostat %>%
#   group_by(Method_Agency2) %>%
#   summarise(Samples = length(Method_Agency2), 
#             Prop_samples = length(Method_Agency2)/nrow(Data_Geostat),
#             Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method_Agency2))
# method_agency_info2 <- method_agency_info2[rev(order(method_agency_info2$Samples)),]
Q_ik_methodagency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method_Agency"])[,-4,drop=FALSE]

##################################
## Models
##################################
##########################
## model1
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, fishing method covariate
###########################
path <- file.path(res_dir, "model1")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1)

plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")
plot_maps(plot_set = c(6), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model2
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, fishing method covariate
###########################
path <- file.path(res_dir, "model2")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))


plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model3
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## fishing method covariate
###########################
path <- file.path(res_dir, "model3")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model3b
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## fishing method covariate
## ST smoother -- random walk
###########################
path <- file.path(res_dir, "model3b")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_method

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6, 13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model4
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, agency covariate
###########################
path <- file.path(res_dir, "model4")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_agency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model5
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, method/agency covariate
###########################
path <- file.path(res_dir, "model5")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model6
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## method/agency covariate
###########################
path <- file.path(res_dir, "model6")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model7
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, method/agency covariate
###########################
path <- file.path(res_dir, "model7")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 5
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)

##########################
## model8
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## method/agency covariate
###########################
path <- file.path(res_dir, "model8")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- NULL
X_itp_inp <- NULL
Xconfig_zcp_inp <- NULL
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- 8
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = Data_inp[,"Method_Agency2"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.01, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model9
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate
###########################
path <- file.path(res_dir, "model9")
dir.create(path, showWarnings = FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.o"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_2_0.so"), to = path)

Data_inp <- Data_Geostat
X_gtp_inp <- X_gtp_input
X_itp_inp <- X_itp_input
Xconfig_zcp_inp <- array(1, dim = c(2,1,n_p))
Xconfig_zcp_inp[2,,] <- 0
Q_ik_inp <- NULL

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 1, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
ObsModel <- cbind("PosDist" = 1, "Link" = 0)
OverdispersionConfig <- c("Eta1" = 0, "Eta2" = 0)
Options <- c("Calculate_range" = 1,
            "Calculate_effective_area" = 1)

settings <- make_settings(Version = "VAST_v8_2_0", 
                          n_x = nrow(Network_sz), 
                          Region = "Stream_network",
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          ObsModel = ObsModel,
                          OverdispersionConfig = OverdispersionConfig,
                          Options = Options,
                          purpose = "index",
                          fine_scale = FALSE, 
                          bias.correct = FALSE)
settings$Method <- "Stream_network"
settings$grid_size_km <- 1


fit0 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = rep(0, nrow(Data_inp)),
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  run_model = FALSE)

Par <- fit0$tmb_list$Parameters
Map <- fit0$tmb_list$Map
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
Map[["lambda2_k"]] <- factor(rep(NA, length(Par[["lambda2_k"]])))

fit1 <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = rep(0, nrow(Data_inp)),
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(getsd = FALSE, newtonsteps = 0),
                  test_fit = FALSE)
check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj)

fit <- fit_model(settings = settings,
                  Lat_i = Data_inp[,"Lat"],
                  Lon_i = Data_inp[,"Lon"],
                  t_iz = Data_inp[,"Year"],
                  c_i = rep(0, nrow(Data_inp)),
                  b_i = Data_inp[,"Catch_KG"],
                  a_i = Data_inp[,"AreaSwept_km2"],
                  v_i = rep(0, nrow(Data_inp)),
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)



df <- data.frame("Model" = c("model1",
                            "model2",
                            "model3",
                            "model4",
                            "model5",
                            "model6",
                            "model7",
                            "model8",
                            "model9"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- df$AIC - min(df$AIC)
df[order(df$AIC),]


### MODEL A: SST, no habitat, gear/agency
modelA <- readRDS(file.path(res_dir, "model8", "Fit.rds"))

### MODEL B: ST, habitat, gear/agency
modelB <- readRDS(file.path(res_dir, "model5", "Fit.rds"))
