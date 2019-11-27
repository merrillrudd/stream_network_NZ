rm(list=ls())

################
## Directories
################

nz_dir <- "/home/merrill/stream_network_NZ"

data_dir <- file.path(nz_dir, "data")

fig_dir <- file.path(nz_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################
library(tidyverse)
library(proj4)

# ################
## Load data
################

## all network
network_raw <- readRDS(file.path(data_dir, "REC2.4fromGDB.rds"))
# network_raw <- REC2.4fromGDB

load(file.path(data_dir, "FINAL_REC2_FOR_PREDICTIONS.Rdata"))
network_old <- REC2

## look at old network for width info
choose <- network_old %>% select("nzsegment", "FWENZ_DSDamAffected", "REC1_WidthHUC.MALF_cumecs") %>%
	rename('DamAffected'=FWENZ_DSDamAffected, "width"=REC1_WidthHUC.MALF_cumecs) %>%
	mutate('width' = width / 1000)

## list of covariates to consider initially at least
covariates <- read.csv(file.path(data_dir, "longfin_covariates_network.csv"), header=TRUE, stringsAsFactors=FALSE)
covar_toUse <- covariates[which(covariates$toUse==1),"x"]

## filter important things from network and adjust distance to km
network <- network_raw %>%
	select(c('CatName','nzsegment','fnode','tnode','Shape_Leng', 'upcoordX', 'upcoordY', 'downcoordX', 'downcoordY','NextDownID','Headwater','REC2_TerminalSegment', covar_toUse)) %>%
	rename('parent_s' = tnode, 'child_s' = fnode, 'length'=Shape_Leng, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX,'NextDownSeg'=NextDownID, 'Headwater'=Headwater) %>%
	mutate('length' = length / 1000)

## remove NAs from network
network <- network[-which(is.na(network$child_s)),]

network_combine <- left_join(network, choose)
sapply(1:ncol(network_combine), function(x) length(which(is.na(network_combine[,x])))/nrow(network_combine))

## add dam affected to covariates to look for
covar_toUse <- c(covar_toUse, "width", "DamAffected")

## network
network_reformat <- network_combine %>% 
	select('CatName', 'nzsegment','parent_s', 'child_s', 'length', 'width','easting_child','northing_child', "NextDownSeg", covar_toUse)  %>%
	rename('easting'=easting_child, 'northing'=northing_child) %>%
	filter(easting != 0) %>%
	mutate('dist_s' = length )#* width)

## identify root nodes
root_nodes <- which(network_reformat$parent_s %in% network_reformat$child_s == FALSE)
true_root_node <- which(network_reformat$NextDownSeg==-1)

root_list <- lapply(1:length(root_nodes), function(x){
	sub <- network_reformat[root_nodes[x],]
	df <- sub %>% mutate('child_s'=parent_s) %>% 
		mutate('parent_s'=0) %>% 
		mutate('dist_s'=Inf) 
	return(df)
})
roots <- do.call(rbind, root_list)

child_roots <- unique(roots$child_s)
root_byChild <- lapply(1:length(child_roots), function(x){
	sub <- roots %>% filter(child_s == child_roots[x])
	return(sub)
})
ii <- sapply(1:length(root_byChild), function(x) nrow(root_byChild[[x]]))
root_single <- root_byChild[which(ii==1)]
root_multi <- root_byChild[which(ii > 1)]
multi_to_single <- lapply(1:length(root_multi), function(x){
	sub <- root_multi[[x]]
	if(all(sub$parent_s == sub$parent_s[1]) & all(sub$child_s == sub$child_s[1])){
		out <- sub[1,]
	} else {
		out <- NULL
	}
	return(out)
})
any(is.null(multi_to_single))
root_single2 <- do.call(rbind, root_single)
root_single3 <- do.call(rbind, multi_to_single)
root_toUse <- unique(rbind.data.frame(root_single2, root_single3))

network_all <- rbind.data.frame(network_reformat, unique(root_toUse))
nrow(network_all)
length(unique(network_all$child_s))
nrow(unique(network_all %>% select(easting,northing)))

e_mult <- names(table(network_all$easting))[which(table(network_all$easting)>1)]
e_uni <- network_all %>% filter(easting %in% e_mult == FALSE)
set.seed(123)
e_rep <- network_all %>% filter(easting %in% e_mult) %>% mutate(easting = easting + runif(length(easting),-1,1)) 

network_all2 <- rbind.data.frame(e_uni, e_rep)
max(table(network_all2$easting))

n_mult <- names(table(network_all2$northing))[which(table(network_all2$northing)>1)]
n_uni <- network_all2 %>% filter(northing %in% n_mult == FALSE)
set.seed(456)
n_rep <- network_all2 %>% filter(northing %in% n_mult) %>% mutate(northing = northing + runif(length(northing),-1,1))

network_all3 <- rbind.data.frame(n_uni, n_rep)
max(table(network_all3$northing))

#############################
## latitude and longitude
#############################

## function to calculate latitude and longitude from eastings and northings
calc_NZ_latlon <- function(northing, easting){
	proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
	p <- project(matrix(c(easting, northing),nrow=1), proj=proj4string, inv=T)
	colnames(p) <- c('long', 'lat')
	return(p)
}

## latitude and longitude for child nodes in network
network_ll_child <- lapply(1:nrow(network_all3), function(x){
	p <- calc_NZ_latlon(northing = network_all3$northing[x], easting = network_all3$easting[x])
	return(p)
})
network_ll_child <- do.call(rbind, network_ll_child)
# network_ll_child <- data.frame(network_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

## attach latitude and longtiude to network
network_full <- cbind.data.frame(network_all3, network_ll_child)
nrow(network_full)
nrow(unique(network_full))
nrow(network_full %>% select('lat','long'))
nrow(network_full %>% select('easting','northing'))

## all observations
# load(file.path(data_dir, "Diadromous fish dataset.Rdata"))
# obs_raw <- NZFFD.REC2.Diad.EF
###############
## NZFFD ######
###############
obs_raw <- read.csv(file.path(data_dir, "NZFFD_Joined_REC2_fulldataset.csv"))

## filter information we need from observations, do some renaming, and label encounter data
obs_enc <- obs_raw %>% 
	select(c('catchname','nzsegment', 'fishmeth','angdie', 'upcoordX','downcoordX','upcoordY','downcoordY','y', 'org.groups')) %>%
	rename('catchment'=catchname,'fishmethod'=fishmeth, 'present'=angdie, 'northing_child'=upcoordY, 'easting_child'=upcoordX, 'northing_parent'=downcoordY, 'easting_parent'=downcoordX, 'year'=y, 'agency'=org.groups) %>%
	mutate('year' = as.numeric(as.character(year))) %>%
	na.omit() %>%
	# mutate('fishmethod' = 'ef') %>%
	mutate('data_type'='encounter') %>%
	rename('data_value'='present') %>%
	mutate('pass'=0) %>%
	mutate('source'='NZFFD') %>%
	select(-c(easting_parent, northing_parent)) %>%
	rename('easting'=easting_child, 'northing'=northing_child)

obs_ll_child <- lapply(1:nrow(obs_enc), function(x){
	p <- calc_NZ_latlon(northing = obs_enc$northing[x], easting = obs_enc$easting[x])
	return(p)
})
obs_ll_child <- do.call(rbind, obs_ll_child)
# obs_ll_child <- data.frame(obs_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

obs_enc <- cbind.data.frame(obs_enc, obs_ll_child)
all(obs_enc$nzsegment %in% network_full$nzsegment)
# obs <- obs %>% filter(nzsegment %in% network_full$nzsegment == TRUE)

## observations
network_sz <- network_full %>% select('CatName','nzsegment','parent_s','child_s','width')
obs_reformat <- inner_join(network_sz, obs_enc, by='nzsegment') %>% filter(parent_s !=0) # %>% select(-c('catchment','northing', 'easting'))
obs_reformat$data_value <- sapply(1:nrow(obs_reformat), function(x){
	if(obs_reformat$data_type[x]!="encounter") out <- obs_reformat$data_value[x]
	if(obs_reformat$data_type[x]=="encounter") out <- ifelse(obs_reformat$data_value[x]==FALSE, 0, 1)
	return(out)
})
obs_reformat <- obs_reformat %>%
	mutate('length' = 125 / 1000) %>%
	mutate('dist_i' = length) %>%
	select(-catchment)

# ## waikato densities 
# waikato_dens_raw <- read.csv(file.path(data_dir, "Waikato Abundance.REC.csv"), stringsAsFactors=FALSE)

# waikato_dens <- waikato_dens_raw %>%
# 	select("Sample.Date_fish",'nzsegment', "E.NZTM", "N.NZTM", "average_measured_stream_channel_width_m", "Angdie.ALL") 
# waikato_dens$Sample.Date_fish <- as.character(waikato_dens$Sample.Date_fish)
# waikato_dens$year <- sapply(1:nrow(waikato_dens), function(x) strsplit(waikato_dens$Sample.Date_fish[x], "/")[[1]][3])
# waikato_dens <- waikato_dens %>%
# 	select(-Sample.Date_fish) %>%
# 	rename('easting'=E.NZTM, "northing"=N.NZTM, 'width' = average_measured_stream_channel_width_m, 'count'=Angdie.ALL) %>%
# 	mutate(width = width/1000) %>%
# 	mutate(length = 150/1000) %>%
# 	mutate(dist_i = length) %>% #* length) %>%
# 	mutate('data_type' = "count") %>%
# 	rename('data_value' = count)

# obs_ll_child <- lapply(1:nrow(waikato_dens), function(x){
# 	p <- calc_NZ_latlon(northing = waikato_dens$northing[x], easting = waikato_dens$easting[x])
# 	return(p)
# })
# obs_ll_child <- do.call(rbind, obs_ll_child)
# # obs_ll_child <- data.frame(obs_ll_child) #%>% dplyr::rename('long_child'=long, 'lat_child'=lat)

# waikato_dens_ll <- cbind.data.frame(waikato_dens, obs_ll_child)
# all(waikato_dens_ll$nzsegment %in% network_full$nzsegment)

# waikato_reformat <- inner_join(network_sz, waikato_dens_ll, by="nzsegment") %>% filter(parent_s != 0)
# waikato_reformat <- waikato_reformat %>% 
# 	select(-width.x) %>%
# 	rename('width' = width.y) %>%
# 	mutate('fishmethod'=unique(obs_reformat$fishmethod)[grepl("Electric",unique(obs_reformat$fishmethod))]) %>%
# 	mutate('agency'=unique(obs_reformat$agency)[grepl('council',unique(obs_reformat$agency))]) %>%
# 	mutate(pass = 0) %>%
# 	mutate(source = "Waikato_region")

obs_all <- rbind.data.frame(obs_reformat)#, waikato_reformat)

obs_full <- obs_all %>% 
			rename('parent_i' = parent_s, 'child_i' = child_s)


## select habitat data from network separately
hab_full <- network_full %>% 
		dplyr::select('nzsegment', 'parent_s','child_s', covar_toUse, easting, northing) %>%
		tidyr::gather(key = covariate, value = value, covar_toUse[1]:covar_toUse[length(covar_toUse)])


saveRDS(obs_full, file.path(data_dir, "NZ_observations.rds"))
saveRDS(network_full, file.path(data_dir, "NZ_network.rds"))
saveRDS(hab_full, file.path(data_dir, "NZ_habitat.rds"))

nzmap <- ggplot(network_full) +
		geom_point(aes(x = long, y = lat), cex=0.2) +
		xlab("Longitude") + ylab("Latitude") +
		theme_minimal()
ggsave(file.path(fig_dir, "NZmap.png"), nzmap)

obsmap <- ggplot() +
		geom_point(data=network_full, aes(x = long, y = lat), col = "black", cex=0.2) +
		geom_point(data=obs_full %>% filter(data_type=="encounter"), aes(x = long, y = lat, color = data_type)) +
		xlab("Longitude") + ylab("Latitude") +
		theme_minimal()
ggsave(file.path(fig_dir, "NZmap_obs_encounter.png"), obsmap)

## save rda
nz_longfin_eel <- list()
nz_longfin_eel$observations <- obs_full
nz_longfin_eel$network <- network_full
nz_longfin_eel$habitat <- hab_full
save(nz_longfin_eel, file=file.path(data_dir, "nz_longfin_eel.rda"))

