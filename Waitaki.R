rm(list=ls())

################
## Directories
################

main_dir <- "C:\\merrill\\stream_networks"
R_dir <- file.path(main_dir, "R")
R_files <- list.files(R_dir)
readr <- sapply(1:length(R_files), function(x) source(file.path(R_dir, R_files[x])))

nz_dir <- file.path(main_dir, "NZ")

res_dir <- file.path(nz_dir, "Waitaki")
dir.create(res_dir, showWarnings=FALSE)

fig_dir <- file.path(res_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
##remember to pull upstream development branches

devtools::install_github("merrillrudd/FishStatsUtils", ref='stream')
devtools::install_github("james-thorson/VAST", ref="development")

library(VAST)
library(FishStatsUtils)
# library(StreamUtils)
library(TMB)
library(tidyverse)
library(RColorBrewer)
library(proj4)
library(RuddR)
library(foreach)
library(doParallel)

#########################
## read in data
##########################

data <- data("nz_waitaki_longfin_eel_downstream", package="FishStatsUtils")

network <- nz_waitaki_longfin_eel_downstream[["network"]]

## add a small value to the latitude and longitude where we don't have updated location for the nodes emptying into the ocean
network$lat[which(network$parent_s == 0)] <- network$lat[which(network$parent_s == 0)] + 0.00001
network$long[which(network$parent_s == 0)] <- network$long[which(network$parent_s == 0)] + 0.00001

## format network data
Network_sz = network %>% select(c('parent_s','child_s','dist_s'))
Network_sz_LL = network %>% select(c('parent_s', 'child_s', 'dist_s', 'lat', 'long')) %>%
  rename("Lon"=long, "Lat"=lat)

## make sure to use only encounter data
obs <- nz_waitaki_longfin_eel_downstream[["observations"]] %>%
    dplyr::filter(data_type=="encounter") %>%
    select(-data_type) %>%
    rename('present' = data_value) %>%
    mutate('vessel_agency' = paste0(fishmethod, "_", agency)) %>%
    mutate('fishmethod2' = ifelse(fishmethod != "Electric fishing", "Other", "Electric fishing")) %>%
    rename("Year"=year)

  # ### by year
  # years <- unique(obs$Year)[order(unique(obs$Year))]

  # Network_sz_LL_wYear <- lapply(1:length(years), function(x){
  #   out <- cbind.data.frame(Network_sz_LL, "Year"=years[x])
  #   return(out)
  # })
  # Network_sz_LL_wYear <- do.call(rbind, Network_sz_LL_wYear)

  # bb <- ggplot(Network_sz_LL_wYear) +
  #     geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", cex=0.5, alpha=0.6) +
  #     geom_point(data = obs, aes(x = long, y = lat, fill = factor(round(present,0))), cex=1.8, pch=22, alpha=0.6) +
  #     scale_fill_brewer(palette = "Set1") +
  #     scale_x_continuous(breaks = round(quantile(obs$long,prob=c(0.2,0.5,0.8)),0), labels=round(quantile(obs$long,prob=c(0.2,0.5,0.8)),0)) +
  #     xlab("Longtiude") + ylab("Latitude") +
  #     facet_wrap(~Year) +
  #     guides(fill = guide_legend(title = "Encounter")) +
  #     mytheme()
  # ggsave(file.path(fig_dir, "Network_byYear_byEncounter.png"), bb, width = 10, height = 8)

##### add small value to encounter observations
present <- obs$present
devs <- rnorm(length(present), 0, 0.01)
present_new <- sapply(1:length(present), function(x) ifelse(present[x]==1, present[x]+devs[x], present[x]))
obs$present <- present_new

##### setup data frame
Data_Geostat <- data.frame( "Catch_KG" = present_new, 
              "Year" = as.numeric(obs$Year),
               "Vessel" = obs$fishmethod,
               "Vessel2" = obs$vessel_agency, 
               "AreaSwept_km2" = obs$dist_i, 
               "Lat" = obs$lat, 
               "Lon" = obs$long, 
               "Pass" = 0,
               "Knot" = obs$child_i,
               "Category" = "Longfin_eels")

##### setup data frame
Data_Geostat_ef <- Data_Geostat %>% filter(Vessel == "Electric fishing")
obs_ef <- obs %>% filter(fishmethod == "Electric fishing")

###### remove angling
# Data_Geostat <- Data_Geostat %>% filter(Vessel != "Angling")
# obs <- obs %>% filter(fishmethod != "Angling")

## habitat data
hab <- nz_waitaki_longfin_eel_downstream[['habitat']]

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
n_i_ef <- nrow(obs_ef)

for(i in 1:length(covar)){
  p <- ggplot(hab %>% filter(covariate == covar[i])) +
  geom_point(aes(x = easting, y = northing, color = value)) +
  guides(color=guide_legend(title=covar[i])) +
  scale_color_viridis_c() +
  mytheme()
  ggsave(file.path(fig_dir, paste0("Habitat_covariate_", covar[i],".png")),p, width=10, height=8)
}

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

X_itp_ef <- array(0, dim=c(n_i_ef,n_t,n_p))
for(i in 1:n_i_ef){
  for(p in 1:n_p){
    child_i <- obs_ef$child_i[i]
    index <- which(nodes == child_i)
    X_itp_ef[i,,p] <- X_gtp_input[index,,p]
  }
}

# X_itp <- array(0, dim=c(n_i,n_t,n_p))
# for(i in 1:n_i){
#   for(p in 1:n_p){
#     child_i <- obs$child_i[i]
#     index <- which(nodes == child_i)
#     X_itp[i,,p] <- X_gtp_input[index,,p]
#   }
# }

####################################
## sampling information
## treated as catchability covariates
###################################

method <- unique(obs$fishmethod)
cmethod <- sapply(1:length(method), function(x) length(which(obs$fishmethod ==method[x])) / nrow(obs))
names(cmethod) <- method
cmethod[order(cmethod)]

Q_ik_method <- ThorsonUtilities::vector_to_design_matrix(obs $fishmethod )[,-2,drop=FALSE]
Q_ik_method2 <- ThorsonUtilities::vector_to_design_matrix( obs$fishmethod2 )[,-1,drop=FALSE]


agency <- unique(obs$agency)
cagency <- sapply(1:length(agency), function(x) length(which(obs$agency == agency[x])) / nrow(obs))
names(cagency) <- agency
cagency[order(cagency)]
## agency2
obs$agency <- as.character(obs$agency)
obs$agency2 <- sapply(1:nrow(obs), function(x) ifelse(obs$agency[x] %in% c("doc",'university'), obs$agency[x], "other"))
Q_ik_agency2 <- ThorsonUtilities::vector_to_design_matrix( obs$agency2 )[,-1,drop=FALSE]

metag <- unique(obs$vessel_agency)
cmetag <- sapply(1:length(metag), function(x) length(which(obs$vessel_agency == metag[x])) / nrow(obs))
names(cmetag) <- metag
cmetag[order(cmetag)]
obs$vessel_agency2 <- sapply(1:nrow(obs), function(x) ifelse(obs$vessel_agency[x] %in% c("Electric fishing_doc","Electric fishing_university"), obs$vessel_agency[x], "Other"))
Q_ik_vesselagency2 <- ThorsonUtilities::vector_to_design_matrix( obs$vessel_agency2 )[,-1,drop=FALSE]

# vessels <- unique(obs$fishmethod)
# count_vessels <- sapply(1:length(vessels), function(x) length(which(obs$fishmethod == vessels[x])))
# Q_ik_method <- ThorsonUtilities::vector_to_design_matrix( obs$fishmethod )[,-1,drop=FALSE]

# agency <- unique(obs_ef$agency)
# count_agency <- sapply(1:length(agency), function(x) length(which(obs_ef$agency == agency[x])))
# names(count_agency) <- agency

# obs_ef$agency2 <- sapply(1:nrow(obs_ef), function(x) ifelse(obs_ef$agency[x] %in% c('fish&game','council','consultants'), 'other', as.character(obs_ef$agency[x])))
# agency2 <- unique(obs_ef$agency2)
# count_agency2 <- sapply(1:length(agency2), function(x) length(which(obs_ef$agency2 == agency2[x])))
# names(count_agency2) <- agency2

# Q_ik_agency_ef <- ThorsonUtilities::vector_to_design_matrix( obs_ef$agency2 )[,-1,drop=FALSE]

##################################
## save data used for model runs
##################################

saveRDS(obs, file.path(res_dir, "observations.rds"))
saveRDS(network, file.path(res_dir, "network.rds"))
saveRDS(hab, file.path(res_dir, "habitat.rds"))


##################################
## Models
##################################
#########################
## SST_IID_ALL_HAB
#########################
path <- file.path(res_dir, "SST_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # save.image(file.path(nz_dir, "test", "NZ_eel_gear.Rdata"))

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map))#,
                # optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    



#########################
## ST_IID_ALL_HAB
#########################
path <- file.path(res_dir, "ST_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.3

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## ST_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "ST_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # save.image(file.path(nz_dir, "test", "NZ_eel_gear.Rdata"))

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  ## adjust starting value
  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- 0.3
  # newpar[["L_epsilon1_z"]] <- 0.01

  fit2 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0, startpar = newpar))
  check <- TMBhelper::Check_Identifiable(fit2$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## SST_RW_ALL_HAB
#########################
path <- file.path(res_dir, "SST_RW_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  newpar <- fit1$parameter_estimates$par
  newpar[["logkappa1"]] <- -1

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= newpar))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_RW_GEAR_HAB
#########################
path <- file.path(res_dir, "SST_RW_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar=fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## SST_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## ST_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "ST_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    


#########################
## SST_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "SST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## ST_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "ST_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8))    

#########################
## T_IID_ALL_HAB
#########################
path <- file.path(res_dir, "T_IID_ALL_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## T_IID_ALL_NOHAB
#########################
path <- file.path(res_dir, "T_IID_ALL_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }



  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  # model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                # model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## T_IID_GEAR_HAB
#########################
path <- file.path(res_dir, "T_IID_GEAR_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    

#########################
## T_IID_GEAR_NOHAB
#########################
path <- file.path(res_dir, "T_IID_GEAR_NOHAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(msubx[3] == "GEAR"){
      X_itp_inp <- X_itp_input
    }
    if(msubx[3] == "ALL"){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }

  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3)#,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_AGENCY_HAB
#########################
path <- file.path(res_dir, "SST_IID_AGENCY_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(grepl("EF", msubx[3])==FALSE){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


#########################
## SST_IID_GEARAGENCY_HAB
#########################
path <- file.path(res_dir, "SST_IID_GEARAGENCY_HAB")
dir.create(path, showWarnings=FALSE)
fig <- file.path(path, "figures")
dir.create(fig)

msub <- strsplit(path, "Waitaki/")[[1]][2]
msubx <- strsplit(msub, "_")[[1]]

ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.cpp"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.dll"), to = path)
ignore <- file.copy(from = file.path(res_dir, "VAST_v8_0_0.o"), to = path)

  if(msubx[3] == "GEAR"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_method2
  }
  if(msubx[3] == "ALL"){ 
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- NULL
  }
  if(msubx[3] == "GEARAGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_vesselagency2
  }
  if(msubx[3] == "AGENCY"){
    Data_Geostat_inp <- Data_Geostat
    Q_ik_inp <- Q_ik_agency2
  }

  if(msubx[4] == "HAB"){
    X_gtp_inp <- X_gtp_input
    Xconfig_zcp_inp <- array(1,dim=c(2,1,n_p))
    Xconfig_zcp_inp[2,,] <- 0
    if(grepl("EF", msubx[3])){
      X_itp_inp <- X_itp_ef
    }
    if(grepl("EF", msubx[3])==FALSE){
      X_itp_inp <- X_itp_input
    }
  } else{ 
    X_gtp_inp <- NULL
    X_itp_inp <- NULL
    Xconfig_zcp_inp <- NULL
  }


  settings <- model_setup(model_info = msubx, Method = "Stream_network")

  # check estimated parameters
  fit0 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  run_model = FALSE)

  Map <- fit0$tmb_list$Map
  Map[["lambda2_k"]] <- factor(rep(NA, length(Map[["lambda2_k"]])))

  # first model run
  fit1 = fit_model( "settings"=settings, 
                  "Lat_i"=Data_Geostat_inp[,"Lat"], 
                  "Lon_i"=Data_Geostat_inp[,"Lon"], 
                  "t_iz"=Data_Geostat_inp[,'Year'], 
                  "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                  "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                  "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                  "v_i"=Data_Geostat_inp[,'Vessel'], 
                  working_dir=path, 
                  extrapolation_args=list(
                    input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]), 
                    Network_sz_LL=Network_sz_LL),
                  Network_sz = Network_sz,
                  Xconfig_zcp = Xconfig_zcp_inp,
                  X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                  Q_ik = Q_ik_inp, 
                  model_args = list(Map = Map),
                  optimize_args = list(getsd=FALSE, newtonsteps=0))
  check <- TMBhelper::Check_Identifiable(fit1$tmb_list$Obj) 

  fit = fit_model( "settings"=settings, 
                "Lat_i"=Data_Geostat_inp[,"Lat"], 
                "Lon_i"=Data_Geostat_inp[,"Lon"], 
                "t_iz"=Data_Geostat_inp[,'Year'], 
                "c_i"=rep(0,nrow(Data_Geostat_inp)), 
                "b_i"=Data_Geostat_inp[,'Catch_KG'], 
                "a_i"=Data_Geostat_inp[,'AreaSwept_km2'], 
                "v_i"=Data_Geostat_inp[,'Vessel'], 
                working_dir=path, 
                extrapolation_args=list(
                  input_grid=cbind("Lat"=Data_Geostat_inp[,"Lat"], "Lon"=Data_Geostat_inp[,"Lon"],"child_i"=Data_Geostat_inp[,"Knot"],"Area_km2"=Data_Geostat_inp[,"AreaSwept_km2"]),
                  Network_sz_LL=Network_sz_LL),
                Network_sz = Network_sz,
                Xconfig_zcp = Xconfig_zcp_inp,
                X_gtp = X_gtp_inp, X_itp = X_itp_inp, 
                Q_ik = Q_ik_inp, 
                # newtonsteps=3,
                model_args = list(Map = Map),
                optimize_args = list(startpar= fit1$parameter_estimates$par))
  saveRDS(fit, file.path(path, "Fit.rds"))    

  fit <- readRDS(file.path(path, "Fit.rds"))    

  # Plot results
  plot_results( settings=settings, fit=fit, working_dir=fig, category_names="Longfin eels", strata_names="Longfin eels", plot_set=c(1,6,8,11,13,15), covar_names = covar2)    


df <- data.frame("Model"=c("SST_IID_AGENCY_HAB",
                          "SST_IID_ALL_HAB",
                          "SST_IID_ALL_NOHAB",
                          "SST_IID_GEAR_HAB",
                          "SST_IID_GEAR_NOHAB",
                          "SST_IID_GEARAGENCY_HAB",
                          "SST_RW_ALL_HAB",
                          "SST_RW_GEAR_HAB",
                          "ST_IID_ALL_HAB",
                          "ST_IID_ALL_NOHAB",
                          "ST_IID_GEAR_HAB",
                          "ST_IID_GEAR_NOHAB",
                          "T_IID_ALL_HAB",
                          "T_IID_ALL_NOHAB",
                          "T_IID_GEAR_HAB",
                          "T_IID_GEAR_NOHAB"))
df$AIC <- NA

for(i in 1:nrow(df)){
  res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
  df$AIC[i] <- res$parameter_estimates$AIC
}
df <- df %>% mutate(dAIC = AIC - min(AIC, na.rm=TRUE))
df[order(df$dAIC),]


gamma <- list()
for(i in 1:nrow(df)){
    res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
    sum <- summary(res$parameter_estimates$SD)
    sum2 <- sum[grepl('gamma1',rownames(sum)),]
    if(nrow(sum2)>0){
      gamma[[i]] <- data.frame("Model"=df[i,"Model"], "gamma"=sum2[,"Estimate"])
    } else{ gamma[[i]] <- NULL}
}


par <- list()
for(i in 1:nrow(df)){
    res <- readRDS(file.path(res_dir, df[i,"Model"], "Fit.rds"))
    dfx <- data.frame(Model = df[i,"Model"], 'par' = res$parameter_estimates$par, "names"=names(res$parameter_estimates$par))
    par[[i]] <- dfx
}




