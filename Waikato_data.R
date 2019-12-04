rm(list=ls())

################
## Directories
################

nz_dir <- "/home/merrill/stream_network_NZ"
sub_dir <- file.path(nz_dir, "Waikato")
dir.create(sub_dir, showWarnings = FALSE)

data_dir <- file.path(nz_dir, "data")
data_dir2 <- file.path(sub_dir, "data")
dir.create(data_dir2, showWarnings = FALSE)

fig_dir <- file.path(sub_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

library(tidyverse)
library(proj4)
library(akima)

# ################
## Load data
################
obsfull <- readRDS(file.path(data_dir, "NZ_observations.rds"))
netfull <- readRDS(file.path(data_dir, "NZ_network.rds"))
habfull <- readRDS(file.path(data_dir, "NZ_habitat.rds"))

#############################
## subset Waikato catchment
#############################

network_sub <- netfull %>% filter(grepl("aikato", CatName))

obs_sub1 <- obsfull %>% filter(grepl("aikato", CatName))
obs_sub2 <- obsfull %>% filter(grepl("aikato", source))
obs_sub <- rbind.data.frame(obs_sub1, obs_sub2)

all(obs_sub2$nzsegment %in% netfull$nzsegment)
all(obs_sub$nzsegment %in% network_sub$nzsegment)

## filter Waikato observations
obs_sub <- obs_sub %>% filter(nzsegment %in% network_sub$nzsegment)
all(obs_sub$nzsegment %in% network_sub$nzsegment)

hab_sub <- habfull %>% filter(child_s %in% network_sub$child_s == TRUE)
covar <- unique(hab_sub$covariate)
covar_toUse <- c('MeanFlowCumecs','Dist2Coast_FromMid','loc_elev','loc_slope','loc_rnvar',"local_twarm",'DamAffected')
all(covar_toUse %in% covar)

#############################
## format
#############################

## rename nodes
nodes <- unique(c(network_sub$child_s, network_sub$parent_s))
inodes <- seq_along(nodes)

net_parents <- sapply(1:nrow(network_sub), function(x){
  if(network_sub$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub$parent_s[x])]
  if(network_sub$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
net_children <- sapply(1:nrow(network_sub), function(x) inodes[which(nodes == network_sub$child_s[x])])

network_sub$parent_s <- net_parents
network_sub$child_s <- net_children

obs_parents <- sapply(1:nrow(obs_sub), function(x){
  if(obs_sub$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub$parent_i[x])]
  if(obs_sub$parent_i[x] == 0) new_node <- 0
  return(new_node)  
})
obs_children <- sapply(1:nrow(obs_sub), function(x) inodes[which(nodes == obs_sub$child_i[x])])

obs_sub$parent_i <- obs_parents
obs_sub$child_i <- obs_children

hab_parents <- sapply(1:nrow(hab_sub), function(x){
  if(hab_sub$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub$parent_s[x])]
  if(hab_sub$parent_s[x] == 0) new_node <- 0
  return(new_node)
})
hab_children <- sapply(1:nrow(hab_sub), function(x) inodes[which(nodes == hab_sub$child_s[x])])

hab_sub$parent_s <- hab_parents
hab_sub$child_s <- hab_children

sapply(1:length(covar_toUse), function(x){
	sub <- hab_sub %>% filter(covariate == covar_toUse[x])
	find <- any(is.na(sub$value))
	names(find) <- covar_toUse[x]
	return(find)
})

hab_sub2 <- lapply(1:length(covar_toUse), function(x){
	sub <- hab_sub %>% filter(covariate == covar_toUse[x])
	# any(is.na(sub$value))
	if(any(is.na(sub$value))){

		if(covar_toUse[x]!="DamAffected"){
			interp_east <- sub$easting[which(is.na(sub$value)==FALSE)]
			interp_north <- sub$northing[which(is.na(sub$value)==FALSE)]
			interp_z <- sub$value[which(is.na(sub$value)==FALSE)]	

			find_df <- data.frame('east' = sub$easting[which(is.na(sub$value))], 'north' = sub$northing[which(is.na(sub$value))])	

			east <- sub$easting[order(sub$easting)]
			north <- sub$northing[order(sub$northing)]
			# mat2 <- zoo::na.approx(object = mat)
			compute <- akima::interp(x = interp_east, y = interp_north, z = interp_z, xo=east, yo=north, extrap=TRUE)
			mat2 <- compute$z	

			vals <- sapply(1:nrow(find_df), function(y){
				mat2[which(compute$x == find_df$east[y]), which(compute$y == find_df$north[y])]
			})	

			inp_vals <- sub$value
			inp_vals[which(is.na(inp_vals))] <- vals	

			sub$value <- inp_vals	

			if(length(which(is.na(sub$value)))==1){
				xx <- sub[(which(is.na(sub$value))-5):(which(is.na(sub$value))+5),]
				xx2 <- xx[order(xx$easting),]
				val_inp <- median(xx$value, na.rm=TRUE)
				sub$value[which(is.na(sub$value))] <- val_inp
			}
			if(length(which(is.na(sub$value)))>1){
				val_inp <- median(sub$value, na.rm=TRUE)
				sub$value[which(is.na(sub$value))] <- val_inp
			}
		}
		if(covar_toUse[x]=="DamAffected"){
			inp <- sub$value
			inp[which(is.na(inp))] <- 2
			ggplot(sub) + geom_point(aes(x = easting, y = northing, color = factor(inp)))

			ina <- which(is.na(sub$value))
			new <- rep(NA, length(ina))
			for(i in 1:length(ina)){
				sub2 <- sub[ina[i],]
				up <- sub %>% filter(parent_s == sub2$child_s)
				down <- sub %>% filter(child_s == sub2$parent_s)
				check <- rbind.data.frame(up, down)
				if(any(is.na(check$value))) check <- check %>% filter(is.na(value)==FALSE)
				if(all(check$value == 0)) new[i] <- 0
				if(all(check$value == 1)) new[i] <- 1
			}

			sub$value[ina] <- new
		}
	}
	return(sub)
})
check <- sapply(1:length(hab_sub2), function(x) any(is.na(hab_sub2[[x]]$value)))
all(check == FALSE)
hab_sub2 <- do.call(rbind, hab_sub2)

find0 <- sapply(1:length(covar_toUse), function(x){
	sub <- hab_sub2 %>% filter(covariate == covar_toUse[x])
	any(sub$value==0)
})
names(find0) <- covar_toUse

hab_sub3 <- lapply(1:length(covar_toUse), function(x){
	sub <- hab_sub2 %>% filter(covariate == covar_toUse[x])
	# any(is.na(sub$value))
	if(any(sub$value == 0)){

		if(covar_toUse[x]!="DamAffected"){
			interp_east <- sub$easting[which(sub$value != 0)]
			interp_north <- sub$northing[which(sub$value != 0)]
			interp_z <- sub$value[which(sub$value != 0)]	

			find_df <- data.frame('east' = sub$easting[which(sub$value == 0)], 'north' = sub$northing[which(sub$value == 0)])	

			east <- sub$easting[order(sub$easting)]
			north <- sub$northing[order(sub$northing)]
			# mat2 <- zoo::na.approx(object = mat)
			compute <- akima::interp(x = interp_east, y = interp_north, z = interp_z, xo=east, yo=north, extrap=TRUE)
			mat2 <- compute$z	

			vals <- sapply(1:nrow(find_df), function(y){
				mat2[which(compute$x == find_df$east[y]), which(compute$y == find_df$north[y])]
			})	

			inp_vals <- sub$value
			inp_vals[which(inp_vals == 0)] <- vals	

			sub$value <- inp_vals	

			if(length(which(sub$value == 0))==1){
				xx <- sub[(which(sub$value == 0)-5):(which(sub$value == 0)+5),]
				xx2 <- xx[order(xx$easting),]
				val_inp <- median(xx$value, na.rm=TRUE)
				sub$value[which(sub$value == 0)] <- val_inp
			}
		}
	}
	print(any(sub$value == 0))
	return(sub)
})
check <- sapply(1:length(hab_sub3), function(x) any(hab_sub3[[x]]$value == 0))
all(check == FALSE)
hab_sub3 <- do.call(rbind, hab_sub3)

saveRDS(obs_sub, file.path(data_dir2, "Waikato_observations.rds"))
saveRDS(network_sub, file.path(data_dir2, "Waikato_network.rds"))
saveRDS(hab_sub3, file.path(data_dir2, "Waikato_habitat.rds"))


## save rda
nz_waikato_longfin_eel <- list()
nz_waikato_longfin_eel$observations <- obs_sub
nz_waikato_longfin_eel$network <- network_sub
nz_waikato_longfin_eel$habitat <- hab_sub3
save(nz_waikato_longfin_eel, file=file.path(data_dir2, "nz_waikato_longfin_eel.rda"))



obs_sub <- readRDS(file.path(data_dir2, "Waikato_observations.rds"))
network_sub <- readRDS(file.path(data_dir2, "Waikato_network.rds"))
hab_sub <- readRDS(file.path(data_dir2, "Waikato_habitat.rds"))

# catchmap <- ggplot() +
# 		geom_point(data=network_sub, aes(x = long, y = lat), col="gray") +
# 		geom_point(data=obs_sub, aes(x = long, y = lat, fill = data_type), pch=22, alpha=0.6) +
# 		scale_fill_brewer(palette = "Set1") +
# 		xlab("Longitude") + ylab("Latitude") +
# 		# guides(fill = FALSE) +
# 		mytheme()
# ggsave(file.path(fig_dir, "Waikato_map.png"), catchmap)

catchmap2 <- ggplot() +
		geom_point(data=netfull, aes(x = long, y = lat), col = "black", cex=0.2) +
		geom_point(data = network_sub, aes(x = long, y = lat), col = "gray") +
		# geom_point(data=obsfull %>% filter(data_type=="encounter"), aes(x = long, y = lat, fill=data_type), pch=22, alpha=0.6) +
		xlab("Longitude") + ylab("Latitude") +
		# scale_fill_brewer(palette = "Set1") +
		theme_minimal()
ggsave(file.path(fig_dir, "Waikato_on_NZ.png"), catchmap2)

# #########################################
# ## Waikato catchment downstream segments
# #########################################

# # loc_df <- hab_sub %>% select('easting', 'northing')
# # loc_mat <- as.matrix(loc_df)
# # loc <- st_linestring(loc_mat)
# # loc_simp <- st_simplify(loc, dTolerance = 5)

# obs_child <- unique(obs_sub$child_i)

# obs_enc <- obs_sub %>% filter(data_type == "encounter")
# obs_count <- obs_sub %>% filter(data_type == "count")

# #################################
# ## network with all observations
# net_obs <- network_sub %>% filter(child_s %in% obs_child)
# nextdown <- network_sub %>% filter(child_s %in% net_obs$parent_s)
# save <- rbind.data.frame(net_obs,nextdown)
# for(i in 1:100){
#   nextdown <- network_sub %>% filter(child_s %in% nextdown$parent_s)
#   save <- unique(rbind.data.frame(save, nextdown))
#   print(nrow(save))
# }
# network_sub2 <- save

# hab_sub2 <- hab_sub %>% filter(child_s %in% network_sub2$child_s)

# ## rename nodes
# nodes <- unique(c(network_sub2$child_s, network_sub2$parent_s))
# inodes <- seq_along(nodes)

# net_parents <- sapply(1:nrow(network_sub2), function(x){
#   if(network_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub2$parent_s[x])]
#   if(network_sub2$parent_s[x] == 0) new_node <- 0
#   return(new_node)
# })
# net_children <- sapply(1:nrow(network_sub2), function(x) inodes[which(nodes == network_sub2$child_s[x])])

# network_sub2$parent_s <- net_parents
# network_sub2$child_s <- net_children

# obs_parents <- sapply(1:nrow(obs_sub), function(x){
#   if(obs_sub$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub$parent_i[x])]
#   if(obs_sub$parent_i[x] == 0) new_node <- 0
#   return(new_node)  
# })
# obs_children <- sapply(1:nrow(obs_sub), function(x) inodes[which(nodes == obs_sub$child_i[x])])

# obs_sub2 <- obs_sub
# obs_sub2$parent_i <- obs_parents
# obs_sub2$child_i <- obs_children

# hab_parents <- sapply(1:nrow(hab_sub2), function(x){
#   if(hab_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub2$parent_s[x])]
#   if(hab_sub2$parent_s[x] == 0) new_node <- 0
#   return(new_node)
# })
# hab_children <- sapply(1:nrow(hab_sub2), function(x) inodes[which(nodes == hab_sub2$child_s[x])])

# hab_sub2$parent_s <- hab_parents
# hab_sub2$child_s <- hab_children

# network_all <- network_sub2
# obs_all <- obs_sub2
# hab_all <- hab_sub2

# ###############################
# ## network with count data 
# net_obs <- network_sub %>% filter(child_s %in% obs_count$child_i)
# nextdown <- network_sub %>% filter(child_s %in% net_obs$parent_s)
# save <- rbind.data.frame(net_obs,nextdown)
# for(i in 1:500){
#   nextdown <- network_sub %>% filter(child_s %in% nextdown$parent_s)
#   save <- unique(rbind.data.frame(save, nextdown))
#   print(nrow(save))
# }
# network_sub2 <- save

# hab_sub2 <- hab_sub %>% filter(child_s %in% network_sub2$child_s)

# ## rename nodes
# nodes <- unique(c(network_sub2$child_s, network_sub2$parent_s))
# inodes <- seq_along(nodes)

# net_parents <- sapply(1:nrow(network_sub2), function(x){
#   if(network_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == network_sub2$parent_s[x])]
#   if(network_sub2$parent_s[x] == 0) new_node <- 0
#   return(new_node)
# })
# net_children <- sapply(1:nrow(network_sub2), function(x) inodes[which(nodes == network_sub2$child_s[x])])

# network_sub2$parent_s <- net_parents
# network_sub2$child_s <- net_children

# obs_sub <- obs_count
# obs_parents <- sapply(1:nrow(obs_sub), function(x){
#   if(obs_sub$parent_i[x] != 0) new_node <- inodes[which(nodes == obs_sub$parent_i[x])]
#   if(obs_sub$parent_i[x] == 0) new_node <- 0
#   return(new_node)  
# })
# obs_children <- sapply(1:nrow(obs_sub), function(x) inodes[which(nodes == obs_sub$child_i[x])])

# obs_sub2 <- obs_sub
# obs_sub2$parent_i <- obs_parents
# obs_sub2$child_i <- obs_children

# hab_parents <- sapply(1:nrow(hab_sub2), function(x){
#   if(hab_sub2$parent_s[x] != 0) new_node <- inodes[which(nodes == hab_sub2$parent_s[x])]
#   if(hab_sub2$parent_s[x] == 0) new_node <- 0
#   return(new_node)
# })
# hab_children <- sapply(1:nrow(hab_sub2), function(x) inodes[which(nodes == hab_sub2$child_s[x])])

# hab_sub2$parent_s <- hab_parents
# hab_sub2$child_s <- hab_children

# network_count <- network_sub2
# obs_count <- obs_sub2
# hab_count <- hab_sub2

# network_all$network <- "all"
# network_count$network <- "counts"

# net_compare <- rbind.data.frame(network_all, network_count)
# p <- ggplot(net_compare) +
# 	geom_point(aes(x = long, y = lat, color = network), alpha=0.5) +
# 	scale_color_brewer(palette = "Set1") +
# 	mytheme()
# ggsave(file.path(fig_dir, "Waikato_network_all_vs_counts.png"), p)
# all(network_count$child_s %in% network_all$child_s)


# saveRDS(obs_all, file.path(data_dir, "Waikato_observations_downstreamOnly.rds"))
# saveRDS(network_all, file.path(data_dir, "Waikato_network_downstreamOnly.rds"))
# saveRDS(hab_all, file.path(data_dir, "Waikato_habitat_downstreamOnly.rds"))
# saveRDS(obs_count, file.path(data_dir, "Waikato_observations_downstreamOnly_count.rds"))
# saveRDS(network_count, file.path(data_dir, "Waikato_network_downstreamOnly_count.rds"))
# saveRDS(hab_count, file.path(data_dir, "Waikato_habitat_downstreamOnly_count.rds"))

# obs_sub2 <- readRDS(file.path(data_dir, "Waikato_observations_downstreamOnly.rds"))
# network_sub2 <- readRDS(file.path(data_dir, "Waikato_network_downstreamOnly.rds"))

# 		l2 <- lapply(1:nrow(network_sub2), function(x){
# 			parent <- network_sub2$parent_s[x]
# 			find <- network_sub2 %>% filter(child_s == parent)
# 			if(nrow(find)>0) out <- cbind.data.frame(network_sub2[x,], 'long2'=find$long, 'lat2'=find$lat)
# 			if(nrow(find)==0) out <- cbind.data.frame(network_sub2[x,], 'long2'=NA, 'lat2'=NA)
# 			# if(nrow(find)>0) out <- cbind.data.frame(network_sub2[x,], 'long2'=find$long, 'lat2'=find$lat)
# 			# if(nrow(find)==0) out <- cbind.data.frame(network_sub2[x,], 'long2'=NA, 'lat2'=NA)
# 			return(out)
# 		})
# 		l2 <- do.call(rbind, l2)

# catchmap <- ggplot() +
# 		geom_point(data=network_sub2, aes(x = long, y = lat), col="gray") +
# 		geom_segment(data=l2, aes(x = long2,y = lat2, xend = long, yend = lat), col="gray") +
# 		geom_point(data=obs_sub2, aes(x = long, y = lat, fill = data_type), pch=22, alpha=0.6) +
# 		scale_fill_brewer(palette = "Set1") +
# 		xlab("Longitude") + ylab("Latitude") +
# 		# guides(fill = FALSE) +
# 		mytheme()
# ggsave(file.path(fig_dir, "Waikato_map_downstream.png"), catchmap)


# data_dir2 <- file.path(nz_dir, "data_save")
# ## save rda
# nz_waikato_longfin_eel_downstream <- list()
# nz_waikato_longfin_eel_downstream$observations <- obs_all
# nz_waikato_longfin_eel_downstream$network <- network_all
# nz_waikato_longfin_eel_downstream$habitat <- hab_all
# save(nz_waikato_longfin_eel_downstream, file=file.path(data_dir2, "nz_waikato_longfin_eel_downstream.rda"))

# nz_waikato_longfin_eel_downstream_counts <- list()
# nz_waikato_longfin_eel_downstream_counts$observations <- obs_count
# nz_waikato_longfin_eel_downstream_counts$network <- network_count
# nz_waikato_longfin_eel_downstream_counts$habitat <- hab_count
# save(nz_waikato_longfin_eel_downstream_counts, file=file.path(data_dir2, "nz_waikato_longfin_eel_downstream_counts.rda"))





