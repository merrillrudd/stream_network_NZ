rm(list=ls())

################
## Directories
################

nz_dir <- "/home/merrill/stream_network_NZ"
sub_dir1 <- file.path(nz_dir, "Waitaki")
sub_dir2 <- file.path(nz_dir, "Waikato")

data_dir <- file.path(nz_dir, "data")
data_dir1 <- file.path(sub_dir1, "data")
data_dir2 <- file.path(sub_dir2, "data")

fig_dir <- file.path(nz_dir, "figures")
dir.create(fig_dir, showWarnings=FALSE)

#################
## Packages
#################

library(tidyverse)
library(RColorBrewer)

# ################
## Load data
################
load(file.path(data_dir, 'nz_longfin_eel.rda'))
load(file.path(data_dir1, 'nz_waitaki_longfin_eel.rda'))
load(file.path(data_dir2, 'nz_waikato_longfin_eel.rda'))


taki_net <- nz_waitaki_longfin_eel[["network"]] %>% mutate("Catchment" = "Waitaki")
kato_net <- nz_waikato_longfin_eel[["network"]] %>% mutate("Catchment" = "Waikato")
reg_net <- rbind.data.frame(taki_net, kato_net)
nz_net <- nz_longfin_eel[["network"]] %>% mutate("Catchment" = "Other") %>% filter(lat %in% reg_net$lat == FALSE) %>% filter(long %in% reg_net$long == FALSE)
nz_reg_net <- rbind.data.frame(reg_net, nz_net)

taki_obs <- nz_waitaki_longfin_eel[["observations"]] %>% 
	mutate(Catchment = "Waitaki") %>%
	rename(present = data_value)
kato_obs <- nz_waikato_longfin_eel[["observations"]] %>% 
	mutate(Catchment = "Waikato") %>%
	rename(present = data_value)
reg_obs <- rbind.data.frame(taki_obs, kato_obs)

map <- ggplot() +
	geom_point(data = nz_net, aes(x = long, y = lat), pch = ".") +
	geom_point(data = reg_net, aes(x = long, y = lat, color = Catchment), cex=0.5) +
	scale_color_brewer(palette = "Set1") +
	xlab("Longitude") + ylab("Latitude") +
	theme_bw()
ggsave(file.path(fig_dir, "Region_map.png"), map, height = 10, width = 11)

## map with regions
regmap <- ggplot() +
	geom_point(data = reg_net, aes(x = long, y = lat), color = "gray") +
	geom_point(data = reg_obs, aes(x = long, y = lat, fill = year), pch=21, cex=5) +
	scale_fill_distiller(palette = "RdBu") +
	facet_wrap(Catchment~., scales = "free") +
	guides(fill=guide_legend(title="Year")) +
	xlab("Longitude") + ylab("Latitude") +
	theme_bw()
ggsave(file.path(fig_dir, "Region_map_observations.png"), regmap, height = 10, width = 20)

takimap <- ggplot() + 
	geom_point(data = taki_net, aes(x = long, y = lat), color = "gray", cex = 0.5) +
	geom_point(data = taki_obs, aes(x = long, y = lat, fill = factor(present)), pch = 21, cex = 4, alpha = 0.75) +
	scale_fill_viridis_d() +
	facet_wrap(year~.) +
	guides(fill=guide_legend(title="Present")) +
	theme_bw()
ggsave(file.path(fig_dir, "Waitaki_map_byYear.png"), takimap, height = 15, width = 18)

katomap <- ggplot() + 
	geom_point(data = kato_net, aes(x = long, y = lat), color = "gray", cex = 0.5) +
	geom_point(data = kato_obs, aes(x = long, y = lat, fill = factor(present)), pch = 21, cex = 4, alpha = 0.75) +
	scale_fill_viridis_d() +
	facet_wrap(year~.) +
	guides(fill=guide_legend(title="Present")) +
	theme_bw()
ggsave(file.path(fig_dir, "Waikato_map_byYear.png"), katomap, height = 15, width = 18)
