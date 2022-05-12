# Teresita M. Porter, April 9, 2022

# Figure 1 - sampling map and mapped gradients

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('Fig1_mapped_gradients.R')

# try to improve reproducibility using groundhog
library(groundhog)

groundhog.library("stringr", '2022-05-03', tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", '2022-05-03', tolerate.R.version='4.1.1') # dcast
groundhog.library("ggplot2", '2022-05-03', tolerate.R.version='4.1.1') # ggplot
groundhog.library("gridExtra", '2022-05-03', tolerate.R.version='4.1.1') # grid.arrange
groundhog.library("cowplot", '2022-05-03', tolerate.R.version='4.1.1') # get_legend
groundhog.library("data.table", '2022-05-03', tolerate.R.version='4.1.1') # setDT
groundhog.library("stringr", '2022-05-03', tolerate.R.version='4.1.1') # str_split
groundhog.library("ggmap", '2022-05-03', tolerate.R.version='4.1.1')
groundhog.library("ggrepel", '2022-05-03', tolerate.R.version='4.1.1') # geom_label_repel

# devtools::install_github('oswaldosantos/ggsn')
groundhog.library("ggsn", '2022-05-03', tolerate.R.version='4.1.1')

############################
# get sewer outflow locations, stations, and don river in one df

sites <- read.csv('Infiles/SitesCSOs.csv', header=T)
names(sites)[1] <- "Station"
sites$Location <- factor(sites$Location)

############################
# main map

bbox <- c(bottom = min(sites$Lat)-0.01, top = max(sites$Lat)+0.01 , right = max(sites$Long)+0.02, left = min(sites$Long)-0.02)

base <- get_stamenmap(bbox = bbox, zoom = 15, maptype = 'terrain') 

map <- ggmap(base, extent = "device") +
  # Don River 43.65857926616557, -79.354480730911
  geom_point(data=sites[35,], aes(x=Long, y=Lat), fill="blue", shape=25, color="blue", cex=2.5) + # plot the points
  # stations
  geom_point(data=sites[1:25,], aes(x=Long, y=Lat, color=Location), cex=2.5) + # plot the points
  geom_label_repel(data=sites[1:25,], aes(x=Long, y=Lat, label=Station), 
                   box.padding = 0.75, min.segment.length = 0, max.overlaps=Inf,
                   size=2.5) +
  # CSOs
  geom_point(data=sites[26:34,], aes(x=Long, y=Lat), fill="black", colour="black", shape=25, cex=2.5) + # plot the points
  labs(x="", y="", title="") + # label the axes
  theme_bw() + 
  theme(
    text = element_text(size = 9),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.text = element_text(size=9),
    legend.position = "right",
    legend.key = element_rect(colour = "white"),
    legend.title = element_blank(),
    legend.spacing.y = unit(0.01, 'cm')
  ) 

map

map <- map + 
  scalebar(x.min = min(sites$Long)-0.02, x.max = max(sites$Long)+0.012,
           y.min = min(sites$Lat)-0.01,  y.max = max(sites$Lat)+0.005, location = "topright",
           dist = 1, transform = TRUE, model = 'WGS84', dist_unit = "km",
           height = 0.02, st.bottom = TRUE, st.dist = 0.02, st.size = 3.5,
           box.fill = c("black", "white"), box.color = "black", border.size = 0.5,
           st.color = "black")
map






############################
# inset map

bbox2 <- c(bottom = min(sites$Lat)-2.8, top = max(sites$Lat)+2.8 , right = max(sites$Long)+4.7, left = min(sites$Long)-3.9)

base2 <- get_stamenmap(bbox = bbox2, zoom = 9, maptype = 'terrain') 

map2 <- 
  ggmap(base2, extent = "device") +
  # labs(title="Toronto") +
  # Toronto harbour  43.63270406524798, -79.37190629436049
   geom_point(data=sites[1:25,], x=-79.37190629436049, y=43.63270406524798, 
              color="black", fill="transparent",shape=22, cex=3) + # plot the points
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.line = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.title= element_blank()
  )
map2

# put inset map in the right place
map3 <- map +
  inset(ggplotGrob(map2), xmin = -79.421, xmax = -79.3825, ymin = 43.65, ymax = 43.665)
map3


#######################################
# Clean up water physical-chemical measurements

W <-read.csv(file='Infiles/metadata.csv',head=TRUE)
names(W)[1] <- "SampleName"

# Remove redundant part of names
W$SampleName <- gsub("TorontoHarbour_", "", W$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
W.1<-data.frame(W, do.call(rbind, str_split(W$SampleName,"_")), stringsAsFactors = FALSE)
names(W.1)[48:50] <- c("Station", "Replicate","Location")

# drop invariant vars, ex. "b_Fluoride"
W.1 <- W.1[,-26]

# no b_temp available for 6699 (use same as 1370, 16.0C)
W.1$b_temp[73:75] <- 16.0

# Average values across 3 reps per station, but only for the not significantly correlated
# Max 13 vars for hier.part, but one paper suggest do not use > 9 vars
# only keep vars not strongly correlated
keep <- c("b_Ammonia_N", "b_Sulfate", "b_Nitrogen_Total",
          "b_Magnesium", "b_pH", "b_SpCond", 
          "b_temp")
water.env <- data.frame(setDT(W.1)[, lapply(.SD, mean),.SDcols=keep, by=Station])

# add water vars to sites
water.map <- merge(sites, water.env, by="Station", all.x=TRUE)

# map for predictors chosen by HP NH3, temp, TN, pH, Mg

# plot b_Ammonia_N
p1 <- ggplot(water.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=water.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=water.map[11:35,], aes(x=Long, y=Lat, size=b_Ammonia_N, color=b_Ammonia_N), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Blues", direction=1) +
  ggtitle("Ammonia (W)") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p1

# plot b_temp
# no b_temp available for 6699 (use same as 1370, 16.0C)
p2 <- ggplot(water.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=water.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=water.map[11:35,], aes(x=Long, y=Lat, size=b_temp, color=b_temp), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Blues", direction=1) +
  ggtitle("Temperature") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p2

# plot b_Nitrogen_Total
p3 <- ggplot(water.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=water.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=water.map[11:35,], aes(x=Long, y=Lat, size=b_Nitrogen_Total, color=b_Nitrogen_Total), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Blues", direction=1) +
  ggtitle("Total Nitrogen") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p3

# plot b_pH
p4 <- ggplot(water.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=water.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=water.map[11:35,], aes(x=Long, y=Lat, size=b_pH, color=b_pH), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Blues", direction=1) +
  ggtitle("pH") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p4

# plot b_Magnesium
p5 <- ggplot(water.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=water.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=water.map[11:35,], aes(x=Long, y=Lat, size=b_Magnesium, color=b_Magnesium), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Blues", direction=1) +
  ggtitle("Magnesium") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p5




############################
# Read in metadata (metals, aromatics, pesticides)
S <-read.table(file='Infiles/TH_extra_data.txt', head=TRUE, sep="\t")
names(S)[1] <- "Station"

# replace NA with zero
S[is.na(S)] <- 0

# change 1765 to 1765A so vals match up across water & sediment samples
S$Station <- gsub("1765", "1765A", S$Station)

# Sum the values acccross each variable type
S.2 <- reshape2::dcast(S, Station ~ VariableType, value.var = "Value", fun.aggregate = sum)

# remove spaces from names
names(S.2) <- c("Station", "Aromatic", "Halogenated_Substituted_Aliphatic",
                "Halogenated_Substituted_Aromatic", "Halogenated_monoaromatic",
                "Herbicides", "Insecticides","Metal",
                "Metalloids", "Organochlorine", "PCB", "Pesticides",
                "Rare_Earth")
# drop invariant vars #Metalloids
S.2 <- S.2[,-9]

# drop predictors that are strongly correlated
# remove Rare Earth, Herbicides, Pesticides, PCB, Halogenated_substituted_aromatic
S.2 <- S.2[,-c(4,6,10,11,12)]

# [1] "Aromatic"                         
# [2] "Halogenated_Substituted_Aliphatic"
# [3] "Halogenated_monoaromatic"         
# [4] "Insecticides"                     
# [5] "Metal"                            
# [6] "Organochlorine" 

# add water vars to sites
sed.map <- merge(sites, S.2, by="Station", all.x=TRUE)

# map for sig vars from HP only

# plot Aromatic hydrocarbons
p6 <- ggplot(sed.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Aromatic, color=Aromatic), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("Aromatic hydrocarbons (S)") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p6

# plot Halogenated_Substituted_Aliphatic
p7 <- ggplot(sed.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Halogenated_Substituted_Aliphatic, color=Halogenated_Substituted_Aliphatic), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("HalSubAli") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p7

# plot Halogenated_monoaromatic
p8 <- ggplot(sed.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Halogenated_monoaromatic, color=Halogenated_monoaromatic), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("HalMon") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p8

# plot Insecticides
p9 <- ggplot(sed.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Insecticides, color=Insecticides), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("Insecticides") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p9

# plot Metals
p10 <- ggplot(sed.map, aes(x=Long, y=Lat,)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Metal, color=Metal), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("Metals (S)") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p10

# plot Organochlorines
p11 <- ggplot(sed.map, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=sed.map[11:35,], aes(x=Long, y=Lat, size=Organochlorine, color=Organochlorine), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Oranges", direction=1) +
  ggtitle("Organochlorines") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p11


###############################
# read in diversity metrics
d <- read.csv("Outfiles/supportingFiles/diversity.metrics.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1)
setDT(d, keep.rownames = TRUE)[]
names(d)[1] <- "Station"

# add lat long
diversity <- merge(d, sites, by="Station", all.x=TRUE)


# plot mi.richness.dna
p12 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=mi.richness.dna, color=mi.richness.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("mi.richness.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p12

# plot Thoracosphaeraceae.dna
p13 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Thoracosphaeraceae.dna, color=Thoracosphaeraceae.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Thoracosphaeraceae.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p13

# plot Haptoria.dna
p14 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Haptoria.dna, color=Haptoria.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Haptoria.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p14

# plot Candonidae.dna
p15 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Candonidae.dna, color=Candonidae.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Candonidae.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p15

# plot Cyprididae.dna
p16 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Cyprididae.dna, color=Cyprididae.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Cyprididae.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p16

# plot Limnesiidae.morph
p17 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Limnesiidae.morph, color=Limnesiidae.morph), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Limnesiidae.morph") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p17

# plot Bosminidae.dna
p18 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Bosminidae.dna, color=Bosminidae.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Bosminidae.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p18

# plot Naididae.dna
p19 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Naididae.dna, color=Naididae.dna), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Naididae.dna") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p19

# plot Hydridae.morph
p20 <- ggplot(diversity, aes(x=Long, y=Lat)) +
  # sewer outflows
  geom_point(data=sed.map[1:9,], aes(x=Long, y=Lat), pch = 25, fill="black", color = "black", size=4) +
  # stations
  geom_point(data=diversity, aes(x=Long, y=Lat, size=Hydridae.morph, color=Hydridae.morph), show.legend=FALSE) +
  # scale_color_viridis() +
  scale_color_distiller(palette="Greys", direction=1) +
  ggtitle("Hydridae.morph") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    legend.position="bottom",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title=element_text(size=10),
    axis.text=element_blank(),
    axis.ticks = element_blank(),
    axis.title=element_blank())
p20

metals.row <- plot_grid(p10, p12, p13, p14, nrow=1)
metals.row

aromatic.row <- plot_grid(p6, p15, p16, p17, nrow=1)
aromatic.row

ammonia.row <- plot_grid(p1, p18, p19, p20, nrow=1)
ammonia.row

g <- plot_grid(map3, metals.row, aromatic.row, ammonia.row, ncol=1, rel_heights=c(2,1,1,1))
ggsave("Outfiles/Fig1_mapped_gradients.pdf", g, width = 8, height = 10)











