rm(list=ls())

################
## Directories
################

main_dir <- "/home/merrill/stream_network_NZ"
# R_dir <- file.path(main_dir, "R")
# R_files <- list.files(R_dir)
# readr <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

res_dir <- file.path(main_dir, "Waikato")
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

load(file.path(data_dir, "nz_waikato_longfin_eel.rda"))
network <- nz_waikato_longfin_eel[["network"]]

## add a small value to the latitude and longitude where we don't have updated location for the nodes emptying into the ocean
# network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
# network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- nz_waikato_longfin_eel[["observations"]]
obs$data_value <- sapply(1:nrow(obs), function(x) ifelse(obs$data_value[x] > 1, 1, obs$data_value[x]))
obs <- obs %>%
    # dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value) %>%
    mutate('method_agency' = paste0(fishmethod, "_", agency)) %>%
    rename("Year"=year)

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
hab <- nz_waikato_longfin_eel[['habitat']]

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

# plot_network(Network_sz_LL = Network_sz_LL, Data_Geostat = Data_Geostat_all, byYear = TRUE, root = TRUE, FilePath = fig_dir, FileName = "Waikato_network")
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
Q_ik_method <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method"])[,-2,drop=FALSE]

agency_info <- Data_Geostat %>%
  group_by(Agency) %>%
  summarise(Samples = length(Agency), 
            Prop_samples = length(Agency)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Agency))
agency_info <- agency_info[rev(order(agency_info$Samples)),]
Q_ik_agency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Agency"])[,-7,drop=FALSE]

method_agency_info <- Data_Geostat %>%
  group_by(Method_Agency) %>%
  summarise(Samples = length(Method_Agency), 
            Prop_samples = length(Method_Agency)/nrow(Data_Geostat),
            Prop_samples_w_encounters = length(which(Catch_KG > 0))/length(Method_Agency))
method_agency_info <- method_agency_info[rev(order(method_agency_info$Samples)),]
Q_ik_methodagency <- ThorsonUtilities::vector_to_design_matrix(Data_Geostat[,"Method_Agency"])[,-8,drop=FALSE]

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
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
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
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


# test1 <- readRDS(file.path(res_dir, "v1", "model1", "Fit.rds"))
##########################
## model2
## spatiotemporal, spatial, temporal####, WITH ST smoother (otherwise SD all NA)
## beta1 IID, lognormal dist, 
## habitat covariate, fishing method/agency covariate
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
Q_ik_inp <- Q_ik_methodagency

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 0, "Epsilon2" = 0)
RhoConfig <- c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 0, "Epsilon2" = 0)
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
                  v_i = Data_inp[,"Method_Agency"],
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
                  v_i = Data_inp[,"Method_Agency"],
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
                  v_i = Data_inp[,"Method_Agency"],
                  working_dir = path,
                  extrapolation_args_input = list(input_grid = cbind("Lat" = Data_inp[,"Lat"], "Lon" = Data_inp[,"Lon"], "child_i" = Data_inp[,"Knot"], "Area_km2" = Data_inp[,"AreaSwept_km2"])),
                  spatial_args = list(Network_sz_LL = Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, 
                  X_itp = X_itp_inp,
                  Q_ik = Q_ik_inp,
                  model_args = list(Map = Map, Par = Par),
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))
# check <- TMBhelper::Check_Identifiable(fit$tmb_list$Obj)

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)



##########################
## model3
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, agency covariate
###########################
path <- file.path(res_dir, "model3")
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
                  v_i = Data_inp[,"Agency"],
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
                  v_i = Data_inp[,"Agency"],
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
                  v_i = Data_inp[,"Agency"],
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

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model4
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, fishing method covariate
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
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
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
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model5
## spatiotemporal, spatial, temporal, 
## beta1 IID, lognormal dist, 
## fishing method covariate
###########################
path <- file.path(res_dir, "model5")
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
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
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
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model6
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## habitat covariate, fishing method * agency covariate
###########################
path <- file.path(res_dir, "model6")
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
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
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
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


##########################
## model7
## spatial, temporal, 
## beta1 IID, lognormal dist, 
## fishing method covariate
###########################
path <- file.path(res_dir, "model7")
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
Par[["logkappa1"]] <- log(1/median(Network_sz$dist_s))
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
                  # optimize_args = list(startpar = fit1$parameter_estimates$par),
                  test_fit = FALSE)
saveRDS(fit, file.path(path, "Fit.rds"))

fit <- readRDS(file.path(path, "Fit.rds"))

plot_maps(plot_set = c(1,6,13), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.5)
plot_maps(plot_set = c(1), fit = fit, Sdreport = fit$parameter_estimates$SD, TmbData = fit$data_list, spatial_list = fit$spatial_list, DirName = fig, category_names = "Longfin_eel", cex = 0.1, Panel = "Year")

map_list <- make_map_info(Region = settings$Region, spatial_list = fit$spatial_list, Extrapolation_List = fit$extrapolation_list)
Enc_prob <- plot_encounter_diagnostic(Report = fit$Report, Data_Geostat = Data_inp, DirName = fig)
plot_range_index(Report = fit$Report, TmbData = fit$data_list, Sdreport = fit$parameter_estimates$SD, Znames = colnames(fit$data_list$Z_xm), PlotDir = fig, Year_Set = fit$year_labels, use_biascorr = TRUE, category_names = "Longfin_eels")
plot_residuals(ObsModel = 1, fit = fit, Data = Data_inp, Network_sz_LL = Network_sz_LL, category_names = "Longfin_eels", FilePath = fig)


df <- data.frame("Model" = c("model1",
                            "model2",
                            "model3",
                            "model4",
                            "model5",
                            "model7"))
df$AIC <- NA
for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df$deltaAIC <- df$AIC - min(df$AIC)
df[order(df$AIC),]




### MODEL A: SST, habitat, gear*agency
modelA <- readRDS(file.path(res_dir, "model2", "Fit.rds"))

### MODEL B: SST, habitat, method
modelB <- readRDS(file.path(res_dir, "model1", "Fit.rds"))

### MODEL C: ST, habitat, method
modelC <- readRDS(file.path(res_dir, "model4", "Fit.rds"))


## compare maps
ep_byModel <- lapply(1:3, function(x){
  if(x == 1){
    Report <- modelA$Report
    year_labels = modelA$year_labels
    years_to_plot = modelA$years_to_plot
    spatial_list <- modelA$spatial_list
    name <- "Spatiotemporal variation,\n Habitat & Method*Agency covariates"
  }
  if(x == 2){
    Report <- modelB$Report
    year_labels = modelB$year_labels
    years_to_plot = modelB$years_to_plot
    spatial_list <- modelB$spatial_list
    name <- "Spatiotemporal variation,\n Habitat & Method covariates"
  }
  if(x == 3){
    Report <- modelC$Report
    year_labels = modelC$year_labels
    years_to_plot = modelC$years_to_plot
    spatial_list <- modelC$spatial_list
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  Array_xct = Report$R1_gcy
  dimnames(Array_xct) <- list(Node = 1:dim(Array_xct)[1], Category = "Longfin eels", Year = year_labels)
  xct <- reshape2::melt(Array_xct) %>% mutate(Model = name)
  xctll <- full_join(xct, cbind.data.frame("Node" = 1:spatial_list$n_g,spatial_list$latlon_g))
  return(xctll)
})
ep <- do.call(rbind, ep_byModel)

plot_ep <- ep %>% filter(Year %in% c(1965, 1995, 2018))

p <- ggplot(plot_ep) +
  geom_point(aes(x = Lon, y = Lat, color = value), cex = 0.5, alpha = 0.75) +
  scale_color_distiller(palette = "Spectral") +
  facet_grid(Year ~ Model) +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_encounter_prob_maps.png"), p, height = 10, width = 12)

## effective area occupied
eao_byModel <- lapply(2:3, function(x){
  if(x == 1){
    SD <- TMB::summary.sdreport(modelA$parameter_estimates$SD)
    TmbData <- modelA$data_list
    year_labels = modelA$year_labels
    years_to_plot = modelA$years_to_plot
    spatial_list <- modelA$spatial_list
    name <- "Spatiotemporal variation,\n Habitat & Method*Agency covariates"
  }
  if(x == 2){
    SD <- TMB::summary.sdreport(modelB$parameter_estimates$SD)
    TmbData <- modelB$data_list
    year_labels = modelB$year_labels
    years_to_plot = modelB$years_to_plot
    spatial_list <- modelB$spatial_list
     name <- "Spatiotemporal variation,\n Habitat & Method covariates"
  }
  if(x == 3){
    SD <- TMB::summary.sdreport(modelC$parameter_estimates$SD)
    TmbData <- modelC$data_list
    year_labels = modelC$year_labels
    years_to_plot = modelC$years_to_plot
    spatial_list <- modelC$spatial_list
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  EffectiveName = "effective_area_cyl"
  SD_effective_area_ctl = SD_log_effective_area_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
  use_biascorr = TRUE
    # Extract estimates
    SD_effective_area_ctl = SD_log_effective_area_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
    # Effective area
    if( use_biascorr==TRUE && "unbiased"%in%names(SD) ){
      SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(SD_effective_area_ctl)) ){
      message("Using bias-corrected estimates for effective area occupied (natural scale)...")
    }else{
      message("Not using bias-corrected estimates for effective area occupied (natural scale)...")
      SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Estimate','Std. Error')]
    }
    # Log-Effective area
    if( use_biascorr==TRUE && "unbiased"%in%names(SD) ){
      SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(SD_log_effective_area_ctl)) ){
      message("Using bias-corrected estimates for effective area occupied (log scale)...")
    }else{
      message("Not using bias-corrected estimates for effective area occupied (log scale)...")
      SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Estimate','Std. Error')]
    }

  Index_ctl=array(SD_log_effective_area_ctl[,,,'Estimate'],dim(SD_log_effective_area_ctl)[1:3])
  dimnames(Index_ctl) <- list(Category = "Longfin eel", Year = year_labels, Stratum = NA)

  sd_Index_ctl=array(SD_log_effective_area_ctl[,,,'Std. Error'],dim(SD_log_effective_area_ctl)[1:3])
  dimnames(sd_Index_ctl) <- list(Category = "Longfin eel", Year = year_labels, Stratum = NA)

  df1 <- reshape2::melt(Index_ctl) %>% rename("Estimate" = value)
  df2 <- reshape2::melt(sd_Index_ctl) %>% rename("SD" = value)
  df <- full_join(df1, df2) %>% mutate(Model = name)
  return(df)
})
eao <- do.call(rbind, eao_byModel)

p <- ggplot(eao) +
  geom_segment(aes(x = Year, xend = Year, y = Estimate - 1.96 * SD, yend = Estimate + 1.96 * SD), color = "red", lwd = 1.2) +
  geom_point(aes(x = Year, y = Estimate), color = "red", cex = 3) +
  geom_line(aes(x = Year, y = Estimate), color = "red") +
  coord_cartesian(ylim = c(0,max(eao$Estimate + 1.96 * eao$SD)*1.01)) +
  facet_grid(~Model) +
  ylab("Effective area occupied (km^2)") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_effective_area_occupied.png"), p, height = 6, width = 10)

### encounter diagnostic
ep_diag_byModel <- lapply(1:3, function(x){
  if(x == 1){
    Report <- modelA$Report
    name <- "Spatiotemporal variation,\n Habitat & Method*Agency covariates"
  }
  if(x == 2){
    Report <- modelB$Report
    name <- "Spatiotemporal variation,\n Habitat & Method covariates"
  }
  if(x == 3){
    Report <- modelC$Report
    name <- "No spatiotemporal variation,\n Habitat & Method covariates"
  }
  # Get bin for each datum
  cutpoints_z=seq(0,1,length=21)
  z_i = cut( Report$R1_i, breaks=cutpoints_z, include.lowest=TRUE )
  midpoints_z = rowMeans( cbind(cutpoints_z[-1],cutpoints_z[-length(cutpoints_z)]) )

  # Get encounter frequency for each bin
  freq_z = tapply( ifelse(Data_Geostat[,'Catch_KG']>0,1,0), INDEX=z_i, FUN=mean )

  # Get expectation given model
  num_z = tapply( Report$R1_i, INDEX=z_i, FUN=length )
  mean_z = tapply( Report$R1_i, INDEX=z_i, FUN=mean )
  var_z = tapply( Report$R1_i, INDEX=z_i, FUN=function(vec){sum(vec*(1-vec))} )
  sd_mean_z = sqrt(var_z / num_z^2)

  df_z <- data.frame('midpoint' = midpoints_z, 
                    'frequency' = freq_z, 
                    "num" = num_z,
                    'mean' = mean_z,
                    'var' = var_z,
                    'sd_mean' = sd_mean_z,
                    'Model' = name)
  return(df_z)
})
ep_diag <- do.call(rbind, ep_diag_byModel)

p <- ggplot(ep_diag %>% filter(is.na(frequency) == FALSE)) +
  geom_ribbon(aes(x = midpoint, ymin = mean - 1.96 * sd_mean, ymax = mean + 1.96 * sd_mean), fill = "red", color = NA, alpha = 0.25) +
  geom_line(aes(x = midpoint, y = mean), lwd = 1.5, color = "red") +
  geom_point(aes(x = midpoint, y = frequency), cex = 3) +
  geom_abline(aes(slope = 1, intercept = 0), lty = 2) +
  facet_grid(~Model) + 
  xlab("Observed encounter probability") + ylab("Predicted encounter probability") +
  theme_bw(base_size = 14)
ggsave(file.path(fig_dir, "Compare_encounter_diagnostic.png"), p, height = 6, width = 14)

















