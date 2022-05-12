# Teresita M. Porter, May 5, 2022
# Fig 3 PCA & NMDS

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('Fig2_richness.R')

library(groundhog)

groundhog.library("icesTAF", "2022-05-04", tolerate.R.version='4.1.1') # mkdir
groundhog.library("scales", "2022-05-04", tolerate.R.version='4.1.1')
groundhog.library("cowplot", "2022-05-04", tolerate.R.version='4.1.1')
groundhog.library("maps", "2022-05-04", tolerate.R.version='4.1.1') #needed for goeveg
groundhog.library("goeveg", "2022-05-04", tolerate.R.version='4.1.1') # scree
groundhog.library("stringr", "2022-05-04", tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", "2022-05-04", tolerate.R.version='4.1.1') # dcast
groundhog.library("vegan", "2022-05-04", tolerate.R.version='4.1.1') # rrarefy
groundhog.library("data.table", "2022-05-04", tolerate.R.version='4.1.1') # setDT
groundhog.library("plyr", "2022-05-04", tolerate.R.version='4.1.1') # ddply
groundhog.library("gridExtra", "2022-05-04", tolerate.R.version='4.1.1') #grid.arrange
groundhog.library("moments", "2022-05-04", tolerate.R.version='4.1.1') #skewness
groundhog.library("dplyr", "2022-05-04", tolerate.R.version='4.1.1') # group_by

#############################################
# Process DNA, macroinvert phyla, family rank, pooled reps

# Read in COI.csv from MetaWorks v1.0.0
F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE) # limited to Arthropoda
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE) # Arthropoda, Nematoda(remove because not found from Morph), Cnidaria, Annelida, Platyhelminthes, Mollusca

# edit column headers so matrices can be merged
# GlobalESV
names(F)[1] <- "GlobalESV"
names(M)[1] <- "GlobalESV"

# Seqfield
names(F)[4] <- "ORF/ESVSeq"
names(M)[4] <- "ORF/ESVSeq"

# combine into 1 matrix
COI <- rbind(F, M)

# limit to macroinvertebrate phyla detected using morphology for a fair comparison
morphPhyla <- c("Arthropoda", "Annelida", "Mollusca", "Cnidaria", "Platyhelminthes")
COI2 <- COI[COI$Phylum %in% morphPhyla,]

# filter for 95% correct species, 99% correct genus & up
COI2$Species <- ifelse(COI2$sBP >= 0.7, COI2$Species, "")
COI2$Genus <- ifelse(COI2$gBP >= 0.3, COI2$Genus, "")
COI2$Family <- ifelse(COI2$fBP >= 0.2, COI2$Family, "")

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
COI2$SampleName <- gsub("TorontoHarbour_", "", COI2$SampleName)

# Reformat negative control (blanks) filenames
COI2$SampleName <- gsub("BCBl_COI18Sp","BCBl__BoxCorerBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct1_BCBl__BoxCorerBlank_S77","BCBl__BoxCorerBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct3_BCBl__BoxCorerBlank_S80","BCBl__BoxCorerBlank", COI2$SampleName)
COI2$SampleName <- gsub("HoBl_COI18Sp","HoBl__HoseBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct1_HoBl__HoseBlank_S76","HoBl__HoseBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct3_HoBl__HoseBlank_S79","HoBl__HoseBlank", COI2$SampleName)
COI2$SampleName <- gsub("StBl_COI18Sp","StBl__StrawBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct1_StBl__StrawBlank_S78","StBl__StrawBlank", COI2$SampleName)
COI2$SampleName <- gsub("Oct3_StBl__StrawBlank_S81","StBl__StrawBlank", COI2$SampleName)

# Track negative control samples (ESVs 68/2517 = 0.027 or 2.7%; reads 69480/1565582 or 0.04 or 4%)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- COI2[(COI2$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from neg controls and neg control samples from main dataset
COI3 <- COI2[!COI2$GlobalESV %in% remove,]
COI4 <- COI3[!(COI3$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_COI_S78", "Opt_2_COI_S79", "Opt_3_COI_S80", 
                         "Opt_4_COI_S81","Opt_5_COI_S82","Opt_1_18S_S83",
                         "Opt_2_18S_S84", "Opt_3_18S_S85", "Opt_4_18S_S86",
                         "Opt_5_18S_S87")
COI5 <- COI4[!(COI4$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get station & replicate fields
COI5 <- data.frame(COI5, do.call(rbind, str_split(COI5$SampleName,"_")), stringsAsFactors = FALSE)
names(COI5)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
COI5$Station <- gsub("1271","1371", COI5$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
COI6 <- data.frame(COI5, do.call(rbind, str_split(COI5$GlobalESV,"_")), stringsAsFactors = FALSE)
names(COI6)[37:38] <- c("Amplicon", "Zotu")

# Species

COI9 <- COI6[!COI6$Species=="",]

# remove extra columns 
species.dna <- COI9[,c("Station", "Species", "ESVsize")]

# pivot to make sample x ESV matrix 
species.dna.df <- reshape2::dcast(species.dna, Station ~ Species, value.var = "ESVsize", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each species, per sample
# denominator should be all reads per station even if unassigned
denominator <- data.frame(COI6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotalReads"
species.dna.prop <- species.dna.df[,-1]/denominator$TotalReads
species.dna.prop$Station <- species.dna.df$Station

# edit location to specify DNA samples
species.dna.prop$Station <- paste(species.dna.prop$Station, "DNA", sep=".")

# now bring in morph, convert to prop, and combine with DNA for NMDS & PCA
#############################################
# Morphology, family rank, pooled replicates
# Parse morphology results

trad <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# drop entries without a Family level identification
trad4 <- trad[trad$Species!="",]

# replace space with hyphen
trad4$Species <- gsub(" ", "_", trad4$Species)

# rename 1765 from morph as 1765A to match DNA and siteTable
trad4$station <- gsub("1765", "1765A", trad4$station)

# pivot to make ESV matrix 
species.morph.df <- reshape2::dcast(trad4, station ~ Species, value.var = "individuals", fun.aggregate = sum)

# convert to proportions, for each station, show proportion of individuals recovered per species
# denominator should be all indv per station, even if unassigned
denominator <- data.frame(trad %>% group_by(station) %>% dplyr::summarize(sum(individuals)))
names(denominator)[2] <- "TotIndv"
species.morph.prop <- species.morph.df[,-1]/denominator$TotIndv
species.morph.prop$Station <- species.morph.df$station

# edit show show stations from morph
species.morph.prop$Station <- paste(species.morph.prop$Station, "Morph", sep=".")



###################
# combine, automatically add new cols and rows
dna.morph <- plyr::rbind.fill(species.dna.prop, species.morph.prop)

# move Station to rownames
rownames(dna.morph) <- dna.morph$Station
dna.morph$Station <- NULL

# change NA to zeros
dna.morph[is.na(dna.morph)] <- 0



###################
# directly compare DNA amd morph & species rank using PCA

my.pca <- rda(dna.morph)

pca1var <- summary(my.pca)[[6]]$importance[2]
# 0.2874014

pca2var <- summary(my.pca)[[6]]$importance[5]
# [1] 0.1919886

# total variation explained by first 2 axes
pca1var + pca2var
# [1] 0.4793899


# calculate scores, type 2 scaling, correlation plots
# angles between vectors represent correlations
# choices 1 (species)
species_scor = data.frame(scores(my.pca, display=c("sp"), scaling=2))
setDT(species_scor, keep.rownames = TRUE)[]
names(species_scor)[1] <- "Species"
# filter for clarity
species_scor <- species_scor[abs(species_scor$PC1)>= 0.2 |
                               abs(species_scor$PC2)>=0.2,]

#sites
sites_scor = data.frame(scores(my.pca, display=c("si"), scaling=2))
sites_scor$Station.Type <- rownames(sites_scor)
sites_scor <- data.frame(sites_scor, do.call(rbind, str_split(sites_scor$Station.Type,"\\.")), stringsAsFactors = FALSE)
names(sites_scor)[4:5] <- c("Station", "Type")


# map locations to stations
sites <- read.csv('Infiles/SitesCSOs.csv', header=T)
names(sites)[1] <- "Station"

# remove CSOs
sites <- sites[!sites$Key=="CSO",]
sites <- sites[!sites$Key=="River",]
# remove Key field
sites <- sites[,-3]

# left join to add locations to stations (always check for problems here)
site_scores.1 <- merge(sites_scor, sites, by="Station", all.x = TRUE)

# create factors
site_scores.1$Location <- factor(site_scores.1$Location, 
                                 levels=unique(site_scores.1$Location),
                                 labels=c("Inner Harbour West", "Bathurst Quay", "Inner Harbour Centre",
                                          "Queen Elizabeth Docks", "Inner Harbour East", "Keating Channel",
                                          "Ship Channel", "Ship Turning Basin", "Toronto Island"))

chulls.location <- ddply(site_scores.1, .(Location), function(site_scores.1) site_scores.1[chull(site_scores.1$PC1, site_scores.1$PC2), ])


pca <- ggplot(site_scores.1, aes(x = PC1, y = PC2))+
  ggtitle("") +
  labs(x= paste("PC1 (", round(pca1var*100,1), "%)", sep=""), y=paste("PC2 (", round(pca2var*100,1),"%)", sep="")) +
  # coord_cartesian(x = c(-1, 1.25), y = c(-1, 1))+
  scale_fill_manual(values = c("darkgrey", "white")) +
  scale_shape_manual(values = c(21, 22)) +
  # geom_polygon(data=chulls.location, aes(x=PC1, y=PC2, fill=Type), alpha=0.5, show.legend = FALSE) +
  geom_point(data=site_scores.1, aes(x = PC1, y = PC2, fill=Type, shape=Type), size = 3) +
  geom_segment(data = species_scor,
               aes(x = 0, xend = PC1,
                   y = 0, yend = PC2),
               arrow = NULL, 
               colour = "black") + 
  geom_text(data = species_scor,
            aes(x= PC1*1.05, y = PC2*1.05, label = Species),
            vjust = 0.5,
            hjust = 1,
            size = 3, 
            col="black") +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=5),
    legend.position="none") +
  guides(fill=guide_legend(ncol=1),
         color=NULL)
pca





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
keep <- c("b_Ammonia_N", "b_Sulfate", "b_Nitrogen_Total",
          "b_Magnesium", "b_pH", "b_SpCond", 
          "b_temp")
W.2 <- data.frame(setDT(W.1)[, lapply(.SD, mean),.SDcols=keep, by=Station])
rownames(W.2) <- W.2$Station
W.2$Station <- NULL

# test each env var for normality
hist(W.2$b_Ammonia_N)
shapiro.test(W.2$b_Ammonia_N)
skewness(W.2$b_Ammonia_N)
skewness(sqrt(W.2$b_Ammonia_N))
skewness(log10(W.2$b_Ammonia_N))
skewness(1/(W.2$b_Ammonia_N)) # this one

hist(W.2$b_Sulfate)
shapiro.test(W.2$b_Sulfate)
skewness(W.2$b_Sulfate)
skewness(sqrt(W.2$b_Sulfate))
skewness(log10(W.2$b_Sulfate))
skewness(1/(W.2$b_Sulfate)) # this one

hist(W.2$b_Nitrogen_Total)
shapiro.test(W.2$b_Nitrogen_Total)
skewness(W.2$b_Nitrogen_Total)
skewness(sqrt(W.2$b_Nitrogen_Total))
skewness(log10(W.2$b_Nitrogen_Total))
skewness(1/(W.2$b_Nitrogen_Total)) # this one

hist(W.2$b_Magnesium) # this one
shapiro.test(W.2$b_Magnesium)
skewness(W.2$b_Magnesium)
skewness(sqrt(max(W.2$b_Magnesium+1) - W.2$b_Magnesium))
skewness(log10(max(W.2$b_Magnesium+1) - W.2$b_Magnesium))
skewness(1/(max(W.2$b_Magnesium+1) - W.2$b_Magnesium)) # this one

hist(W.2$b_pH) # this one
shapiro.test(W.2$b_pH)
skewness(W.2$b_pH)

hist(W.2$b_SpCond)
shapiro.test(W.2$b_SpCond)
skewness(W.2$b_SpCond)
skewness(sqrt(W.2$b_SpCond))
skewness(log10(W.2$b_SpCond))
skewness(1/(W.2$b_SpCond)) # this one

hist(W.2$b_temp)
shapiro.test(W.2$b_temp)
skewness(W.2$b_temp)
skewness(sqrt(max(W.2$b_temp+1) - W.2$b_temp)) # this one

# create df.envir with transformed vars and NAs filled from nearest station
water.env <- cbind(1/(W.2$b_Ammonia_N), 
                   1/(W.2$b_Sulfate),
                   1/(W.2$b_Nitrogen_Total),
                   1/(max(W.2$b_Magnesium+1) - W.2$b_Magnesium),
                   W.2$b_pH,
                   1/(W.2$b_SpCond),
                   sqrt(max(W.2$b_temp+1) - W.2$b_temp)
)

water.envir <- data.frame(apply(water.env[,], 2, scale))
names(water.envir) <- c("NH3_N", "SO4", "TN", "Mg",
                        "pH", "SpCond", "Temp")
rownames(water.envir) <- rownames(W.2)
water.envir <- setDT(water.envir, keep.rownames = TRUE)[]
names(water.envir)[1] <- "Station"

#######################################
# Clean up sediment metal/aromatics/pesticides measurements

S <-read.table(file='Infiles/TH_extra_data.txt', head=TRUE, sep="\t")
names(S)[1] <- "Station"

# replace NA with zero
S[is.na(S)] <- 0

# change 1765 to 1765A so vals match up across water & sediment samples
S$Station <- gsub("1765", "1765A", S$Station)

# Sum the values acccross each variable type
S.2 <- reshape2::dcast(S, Station ~ VariableType, value.var = "Value", fun.aggregate = sum)
row.names(S.2) <- S.2$Station
S.2$Station <- NULL

# remove spaces from names
names(S.2) <- c("Aromatic", "Halogenated_Substituted_Aliphatic",
                "Halogenated_Substituted_Aromatic", "Halogenated_monoaromatic",
                "Herbicides", "Insecticides","Metal",
                "Metalloids", "Organochlorine", "PCB", "Pesticides",
                "Rare_Earth")
# drop invariant vars #Metalloids
S.2 <- S.2[,-8]

# drop predictors that are strongly correlated
# remove Rare Earth, Herbicides, Pesticides, PCB, Halogenated_substituted_aromatic
S.2 <- S.2[,-c(3,5,9,10,11)]

# [1] "Aromatic"                         
# [2] "Halogenated_Substituted_Aliphatic"
# [3] "Halogenated_monoaromatic"         
# [4] "Insecticides"                     
# [5] "Metal"                            
# [6] "Organochlorine" 

# test for normality

hist(S.2$Aromatic)
shapiro.test(S.2$Aromatic)
skewness(S.2$Aromatic)
skewness(sqrt(S.2$Aromatic))
skewness(log10(S.2$Aromatic))
skewness(1/(S.2$Aromatic)) # this one

hist(S.2$Halogenated_Substituted_Aliphatic)
shapiro.test(S.2$Halogenated_Substituted_Aliphatic)
skewness(S.2$Halogenated_Substituted_Aliphatic)
skewness(sqrt(S.2$Halogenated_Substituted_Aliphatic))
skewness(log10(S.2$Halogenated_Substituted_Aliphatic)) # this one

hist(S.2$Halogenated_monoaromatic)
shapiro.test(S.2$Halogenated_monoaromatic)
skewness(S.2$Halogenated_monoaromatic)
skewness(sqrt(S.2$Halogenated_monoaromatic)) # this one

hist(S.2$Insecticides)
shapiro.test(S.2$Insecticides)
skewness(S.2$Insecticides)
skewness(sqrt(S.2$Insecticides))
skewness(log10(S.2$Insecticides))
skewness(1/(S.2$Insecticides)) # this one

hist(S.2$Metal) # this one
shapiro.test(S.2$Metal)
skewness(S.2$Metal)

hist(S.2$Organochlorine)
shapiro.test(S.2$Organochlorine)
skewness(S.2$Organochlorine) 
skewness(sqrt(S.2$Organochlorine)) 
skewness(log10(S.2$Organochlorine))
skewness(1/(S.2$Organochlorine)) # this one

# create df.envir with transformed vars and NAs filled from nearest station
sed.env <- cbind(
  1/(S.2$Aromatic),
  log10(S.2$Halogenated_Substituted_Aliphatic),
  sqrt(S.2$Halogenated_monoaromatic),
  1/(S.2$Insecticides),
  S.2$Metal,
  1/(S.2$Organochlorine)
)

# column standardization (accounts for different units accross vars)
# standardize using z-scores
# 'scale' centers=TRUE, scales=TRUE columns of matrix
# 'apply' 2-columns, 1-rows
sed.envir <- data.frame(apply(sed.env[,], 2, scale))
names(sed.envir) <- c(
  "Aro",
  "HalSubAli", 
  "HalMon",
  "Insect",
  "Metal",
  "Orgchl")
rownames(sed.envir) <- rownames(S.2)
sed.envir <- setDT(sed.envir, keep.rownames = TRUE)[]
names(sed.envir)[1] <- "Station"


# combine water and sed
envs <- cbind(water.envir, sed.envir)


#############################
# trad taxa

# Station 1772 is a severe outlier due to high abundance of "Dreissena polymorpha"
# drop to visualize the rest of the plot

trad.tmp <- trad[!trad$station=="1772",]

# rename 1765 from morph as 1765A to match DNA and siteTable
trad.tmp$station <- gsub("1765", "1765A", trad.tmp$station)

# pivot to make sample x ESV matrix 
taxa.df <- reshape2::dcast(trad.tmp, station ~ taxa, value.var = "individuals", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each sv, per sample
# denominator should be all indv per station
taxa.prop <- taxa.df[,-1]/rowSums(taxa.df[,-1])
rownames(taxa.prop) <- taxa.df$station

# Scree plots to determine number of dimensions to use for NMDS, use k=2!
mkdir("Outfiles")
mkdir("Outfiles/supportingFiles/")
pdf("Outfiles/supportingFiles/Scree_taxa.pdf")
# check dims
dimcheckMDS(taxa.prop)
dev.off()

# Do 3 dimensional NMDS  
set.seed(12345)
nmds3 <- metaMDS(taxa.prop, k=3, trymax=100)
# stress = 0.08736722 

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("Outfiles/supportingFiles/stressplot_taxa.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n", main="SSU")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.975, non-metric fit R2 = 0.992

# Create grouping matrix for samples by grabbing row names from above matrix
names <- data.frame(row.names(taxa.prop), stringsAsFactors = FALSE)
names(names) <- c("Station")

# # use Station instead of station for easier comparisons with morph info
# siteTable <- read.csv('Infiles/SitesCSOs.csv', header=T)
# names(siteTable)[1] <- "Station"
# # remove CSOs
# siteTable <- siteTable[!siteTable$Key=="CSO",]
# siteTable <- siteTable[!siteTable$Key=="River",]
# # remove Key field
# siteTable <- siteTable[,-3]

# Add location
names.2 <- merge(names, sites, by="Station", all.x=TRUE)

# Grab sites/species scores from NMDS output
df.sites <- data.frame(vegan::scores(nmds3, display = "sites"))
df.sites <- setDT(df.sites, keep.rownames = TRUE)[]
names(df.sites)[1] <- "Station"


# look at species too?
df.species <- data.frame(vegan::scores(nmds3, display = "species"))

# centroids for locations
set.seed(12345)
centroids <- envfit(nmds3, names.2[,c(1:2)], perm = 999, na.rm=TRUE)
factors.centroids <- data.frame(centroids[2]$factors[1])
factors.centroids <- setDT(factors.centroids, keep.rownames = TRUE)[]
names(factors.centroids)[1] <- "Location"
factors.centroids <- factors.centroids[c(25:33),]
names(factors.centroids)[2:3] <- c("NMDS1", "NMDS2")
factors.centroids$Location <- gsub("Location", "", factors.centroids$Location)

# create factors
factors.centroids$Location <- factor(factors.centroids$Location,
                       levels = c("Inner Harbour West","Bathurst Quay","Inner Harbour Centre","Queen Elizabeth Docks","Inner Harbour East","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"),
                       labels = c("Inner Harbour (West)","Bathurst Quay","Inner Harbour (Centre)","Queen Elizabeth Docks","Inner Harbour (East)","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"))


#######################################
# fit environmental variables
#######################################

# also remove 1772 from envs
envs.tmp <- envs[!envs$Station=="1772",]
set.seed(12345)
fit3 <- envfit(nmds3, envs.tmp[,-1], perm = 999, na.rm=TRUE)

# vectors
fit.vectors3 <- fit3[[1]]
fit.pvals3 <- fit3[[1]]$pvals
fit.r3 <- fit3[[1]]$r
vectors.df3 <- as.data.frame(fit.vectors3[[1]])
vectors.df3$pvals <- fit.pvals3
vectors.df3$r <- fit.r3
vectors.df.sig3 <- vectors.df3[vectors.df3$pvals <0.05 &
                                vectors.df3$r >= 0.30,]
vectors.df.sig3$x <- 0
vectors.df.sig3$y <- 0

# plot to visualize how close DNA and morph are
#############################################
# NMDS color by location and fit sediment vars
# chulls.location <- ddply(factors., .(Location), function(gg6) gg6[chull(gg6$NMDS1, gg6$NMDS2), ])

# NMDS plot, color by site
tax <- ggplot(data=factors.centroids, aes(x=NMDS1, y=NMDS2)) + 
  geom_text(aes(label=Location), size=3) +
  # geom_point(data=factors.df3, aes(x=NMDS1, y=NMDS2, col=Location, fill=Location), size = 3, show.legend = TRUE) +
  # geom_polygon(data=chulls.location, aes(x=NMDS1, y=NMDS2, fill=Location), alpha=0.5, show.legend = FALSE) +
  # geom_text(data=gg6, aes(x=NMDS1, y=NMDS2, col=Location, label=Station), show.legend = FALSE) +
  ggtitle("") +
  coord_cartesian(x = c(-0.5, 0.5), y = c(-0.5, 1.25))+
  geom_segment(data=vectors.df.sig3, 
               mapping=aes(x=x,y=y,xend=NMDS1, yend=NMDS2), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.5,
               color="black") +
  geom_text(vectors.df.sig3,
            mapping=aes(x=NMDS1+0.03, y=NMDS2+0.03, label=rownames(vectors.df.sig3)),
            #  position=position_jitter(),
            hjust="outward", vjust="outward", 
            size=2.5, color="black") +
  # annotate("text", x=Inf, y=Inf, label="Stress = 0.14, R2 = 0.96", hjust = 1.1, vjust = 1.3, size=2.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=5),
    legend.position="none") +
  guides(fill=guide_legend(ncol=1),
         color=NULL)
tax




# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
# sort dna.morph so that the order is the same as in gg where the metadata is
sor <- vegdist(taxa.prop, "bray", binary=TRUE)


# Put it all in one df for ggplot
gg6 <- merge(df.sites, names.2, by="Station")

# create factors
gg6$Station <- factor(gg6$Station,
                      levels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","6699"),
                      labels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","6699"))
gg6$Location <- factor(gg6$Location,
                      levels = c("Inner Harbour West","Bathurst Quay","Inner Harbour Centre","Queen Elizabeth Docks","Inner Harbour East","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"),
                      labels = c("Inner Harbour (West)","Bathurst Quay","Inner Harbour (Centre)","Queen Elizabeth Docks","Inner Harbour (East)","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"))


# # Calculate beta dispersion (homogeneity needed for adonis)
target <- rownames(taxa.prop)
gg6 <- gg6[match(target, gg6$Station),]
#ensure gg$Type is in the same order as sor
bd.location <- betadisper(sor, as.factor(gg6$Location))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
anova(bd.location) # 1.892e-06 ***

pdf("Outfiles/supportingFiles/Taxa_betadisp.pdf")
boxplot(bd.location, main="Type")
dev.off()


# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis(sor ~ Location, data=gg6, permutations=999)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# Location   8    1.5892 0.19865   1.258 0.40154  0.079 .
# Residuals 15    2.3686 0.15791         0.59846         
# Total     23    3.9578                 1.00000 










# COI ESVs

# remove extra columns 
sv.dna <- COI6[,c("Station", "GlobalESV", "ESVsize")]

# pivot to make sample x ESV matrix 
sv.dna.df <- reshape2::dcast(sv.dna, Station ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each sv, per sample
# denominator should be total reads per station, even if not assigned
denominator <- data.frame(COI6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotReads"
sv.dna.prop <- sv.dna.df[,-1]/denominator$TotReads
sv.dna.prop$Station <- sv.dna.df$Station

# edit location to specify DNA samples
sv.dna.prop$Station <- paste(sv.dna.prop$Station, "DNA", sep=".")

# move Station to rownames
rownames(sv.dna.prop) <- sv.dna.prop$Station
sv.dna.prop$Station <- NULL

# change NA to zeros
sv.dna.prop[is.na(sv.dna.prop)] <- 0

# save for RDA
write.csv(sv.dna.prop, "Outfiles/supportingFiles/COI.ESV.prop.csv", quote=FALSE, row.names = TRUE)

# Scree plots to determine number of dimensions to use for NMDS, use k=3!
pdf("Outfiles/supportingFiles/Scree_COI_sv_pooledreps.pdf")
# check dims
dimcheckMDS(sv.dna.prop)
dev.off()

# Do 3 dimensional NMDS this time to keep stress < 2
set.seed(12345)
nmds3 <- metaMDS(sv.dna.prop, k=3, trymax=100)
# stress = 0.1407833 

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("Outfiles/supportingFiles/stressplot_COI_sv_pooledreps.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n", main="SSU")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.858, non-metric fit R2 = 0.98

# Create grouping matrix for samples by grabbing row names from above matrix
names <- data.frame(row.names(sv.dna.prop), stringsAsFactors = FALSE)
names.1 <- data.frame(names, do.call(rbind, strsplit(names$row.names.sv.dna.prop.,'\\.')), stringsAsFactors = FALSE)
names(names.1) <- c("Station.Type", "Station", "Type")

# # use Station instead of station for easier comparisons with morph info
# siteTable <- read.csv('SitesCSOs.csv', header=T)
# names(siteTable)[1] <- "Station"

# Add location
names.2 <- merge(names.1, sites, by="Station", all.x=TRUE)
rownames(names.2) <- names.2$Station.Type
names.2$Station.Type <- NULL

# Grab sites/species scores from NMDS output
df <- data.frame(vegan::scores(nmds3, display = "sites"))

# look at species too?
df.2 <- data.frame(vegan::scores(nmds3, display = "species"))
df.2 <- setDT(df.2, keep.rownames = TRUE)[]
names(df.2)[1] <- "GlobalESV"

# centroids for locations
set.seed(12345)
centroids <- envfit(nmds3, names.2[,c(1,3)], perm = 999, na.rm=TRUE)
factors.centroids <- data.frame(centroids[2]$factors[1])
factors.centroids <- setDT(factors.centroids, keep.rownames = TRUE)[]
names(factors.centroids)[1] <- "Location"
factors.centroids <- factors.centroids[c(26:34),]
names(factors.centroids)[2:3] <- c("NMDS1", "NMDS2")
factors.centroids$Location <- gsub("Location", "", factors.centroids$Location)

# centroids
factors.centroids$Location <- factor(factors.centroids$Location,
                       levels = c("Inner Harbour West","Bathurst Quay","Inner Harbour Centre","Queen Elizabeth Docks","Inner Harbour East","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"),
                       labels = c("Inner Harbour (West)","Bathurst Quay","Inner Harbour (Centre)","Queen Elizabeth Docks","Inner Harbour (East)","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"))


#######################################
# fit environmental variables
#######################################
set.seed(12345)
fit <- envfit(nmds3, envs[,-1], perm = 999, na.rm=TRUE)

# vectors
fit.vectors <- fit[[1]]
fit.pvals <- fit[[1]]$pvals
fit.r <- fit[[1]]$r
vectors.df <- as.data.frame(fit.vectors[[1]])
vectors.df$pvals <- fit.pvals
vectors.df$r <- fit.r
vectors.df.sig <- vectors.df[vectors.df$pvals <0.05 &
                                             vectors.df$r >= 0.30,]
vectors.df.sig$x <- 0
vectors.df.sig$y <- 0

# plot to visualize how close DNA and morph are
#############################################
# NMDS color by location and fit sediment vars
# chulls.location <- ddply(gg4, .(Location), function(gg4) gg4[chull(gg4$NMDS1, gg4$NMDS2), ])

# NMDS plot, color by site
sv.tmp <- ggplot(data=factors.centroids, aes(x=NMDS1, y=NMDS2)) + 
  geom_text(aes(label=Location)) +
  # geom_point(data=gg4, aes(x=NMDS1, y=NMDS2, col=Location, fill=Location), size = 3, show.legend = TRUE) +
  #  geom_text(data=gg, aes(x=NMDS1, y=NMDS2, col=Location, label=Station), show.legend = FALSE) +
  # geom_polygon(data=chulls.location, aes(x=NMDS1, y=NMDS2, fill=Location), alpha=0.5, show.legend = FALSE) +
  ggtitle("Macroinvertebrate ESVs") +
  geom_segment(data=vectors.df.sig, 
               mapping=aes(x=x,y=y,xend=NMDS1, yend=NMDS2), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.5,
               color="black") +
  geom_text(vectors.df.sig,
            mapping=aes(x=NMDS1+0.03, y=NMDS2+0.03, label=rownames(vectors.df.sig)),
            #  position=position_jitter(),
            hjust="outward", vjust="outward", 
            size=2.5, color="black") +
  annotate("text", x=Inf, y=Inf, label="Stress = 0.03, R2 = 0.99", hjust = 1.1, vjust = 1.3, size=2.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=5),
    legend.position="bottom") +
  guides(fill=guide_legend(ncol=1),
         color=NULL)
sv.tmp

l <- get_legend(sv.tmp)

sv <- ggplot(data=factors.centroids, aes(x=NMDS1, y=NMDS2)) + 
  geom_text(aes(label=Location), size=3) +
  # geom_point(data=gg4, aes(x=NMDS1, y=NMDS2, col=Location, fill=Location), size = 3, show.legend = TRUE) +
  #  geom_text(data=gg, aes(x=NMDS1, y=NMDS2, col=Location, label=Station), show.legend = FALSE) +
  # geom_polygon(data=chulls.location, aes(x=NMDS1, y=NMDS2, fill=Location), alpha=0.5, show.legend = FALSE) +
  ggtitle("") +
  coord_cartesian(x = c(-1.5, 2.25), y = c(-1, 1))+
  geom_segment(data=vectors.df.sig, 
               mapping=aes(x=x,y=y,xend=NMDS1, yend=NMDS2), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.5,
               color="black") +
  geom_text(vectors.df.sig,
            mapping=aes(x=NMDS1+0.03, y=NMDS2+0.03, label=rownames(vectors.df.sig)),
            #  position=position_jitter(),
            hjust="outward", vjust="outward", 
            size=2.5, color="black") +
  # annotate("text", x=Inf, y=Inf, label="Stress = 0.03, R2 = 0.99", hjust = 1.1, vjust = 1.3, size=2.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=5),
    legend.position="none") +
  guides(fill=guide_legend(ncol=1),
         color=NULL)
sv




# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
# sort dna.morph so that the order is the same as in gg where the metadata is
sor <- vegdist(sv.dna.prop, "bray", binary=TRUE)

# # Put it all in one df for ggplot
gg4 <- merge(df, names.2, by="row.names")

# create factors
gg4$Station <- factor(gg4$Station,
                       levels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","1772","6699"),
                       labels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","1772","6699"))


# Calculate beta dispersion (homogeneity needed for adonis)
target <- rownames(sv.dna.prop)
gg4 <- gg4[match(target, gg4$Row.names),]
# ensure gg$Type is in the same order as sor
bd.location <- betadisper(sor, as.factor(gg4$Location))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
anova(bd.location) # 2.069e-07 ***

pdf("Outfiles/supportingFiles/dnamorph_COI_sv_betadisp.pdf")
boxplot(bd.location, main="Type")
dev.off()


# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis(sor ~ Location, data=gg4, permutations=999)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   8    3.4614 0.43268  1.3894 0.40992  0.001 ***
# Residuals 16    4.9827 0.31142         0.59008           
# Total     24    8.4441                 1.00000  












# 18S ESVs

S <- read.csv(file="Infiles/results_18S_v4.1.csv", head=TRUE) # includes non-arthropods too
names(S)[1] <- "GlobalESV"
names(S)[4] <- "ORF/ESVSeq"
# SuperKingdom/Domain
names(S)[9] <- "SuperKingdom"
names(S)[10] <- "SuperKingdomRank"
names(S)[11] <- "skBP"

# add Species, SpeciesRank, sBP to S
S$Species <- ""
S$SpeciesRank <- ""
S$sBP <- 0

# Curate SILVA taxa higher level taxa are equivalent to NCBI taxa at this high level
S$Kingdom <- gsub("Animalia", "Metazoa", S$Kingdom)
S$Kingdom <- gsub("Chloroplastida", "Viridiplantae", S$Kingdom)

# exclude macro inverts (only a couple found anyways)
morphPhyla <- c("Arthropoda", "Annelida", "Mollusca", "Cnidaria", "Platyhelminthes")
SSU <- S[!S$Phylum %in% morphPhyla,]

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
S$SampleName <- gsub("TorontoHarbour_", "", S$SampleName)

# Reformat negative control (blanks) filenames
S$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", S$SampleName)
S$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", S$SampleName)
S$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", S$SampleName)
S$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", S$SampleName)
S$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", S$SampleName)
S$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", S$SampleName)


# Track negative control samples (ESVs 336/6602 = 0.05 or 5%; reads 69480/1285245 = 0.05 or 5%)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- S[(S$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from neg controls and neg control samples from main dataset
S3 <- S[!S$GlobalESV %in% remove,]
S4 <- S3[!(S3$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_18S_S83", "Opt_2_18S_S84", "Opt_3_18S_S85", 
                         "Opt_5_18S_S87","Opt_4_18S_S86")
S5 <- S4[!(S4$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get station & replicate fields
S5 <- data.frame(S5, do.call(rbind, str_split(S5$SampleName,"_")), stringsAsFactors = FALSE)
names(S5)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
S5$Station <- gsub("1271","1371", S5$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
S6 <- data.frame(S5, do.call(rbind, str_split(S5$GlobalESV,"_")), stringsAsFactors = FALSE)
names(S6)[37:38] <- c("Amplicon", "Zotu")

# remove extra columns 
sv.dna2 <- S6[,c("Station", "GlobalESV", "ESVsize")]

# pivot to make sample x ESV matrix 
sv.dna.df2 <- reshape2::dcast(sv.dna2, Station ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each sv, per sample
# denominator should be total reads per station
denominator <- data.frame(S6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotReads"
sv.dna.prop2 <- sv.dna.df2[,-1]/denominator$TotReads
sv.dna.prop2$Station <- sv.dna.df2$Station

# edit location to specify DNA samples
sv.dna.prop2$Station <- paste(sv.dna.prop2$Station, "DNA", sep=".")

# morph can't produce data at this resolution usually

# move Station to rownames
rownames(sv.dna.prop2) <- sv.dna.prop2$Station
sv.dna.prop2$Station <- NULL

# change NA to zeros
sv.dna.prop2[is.na(sv.dna.prop2)] <- 0

# Scree plots to determine number of dimensions to use for NMDS, use k=3!
pdf("Outfiles/supportingFiles/Scree_18S.pdf")
# check dims
dimcheckMDS(sv.dna.prop2)
dev.off()

# Do 3 dimensional NMDS this time to keep stress < 2
set.seed(12345)
nmds3 <- metaMDS(sv.dna.prop2, k=3, trymax=100)
# stress = 0.06361157 

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("Outfiles/supportingFiles/stressplot_18S.pdf")
stressplot(nmds3)
gof <-goodness(nmds3)
gof
plot(nmds3, display = "sites", type="n", main="SSU")
points(nmds3, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.989, non-metric fit R2 = 0.996

# Create grouping matrix for samples by grabbing row names from above matrix
names <- data.frame(row.names(sv.dna.prop2), stringsAsFactors = FALSE)
names.1 <- data.frame(names, do.call(rbind, strsplit(names$row.names.sv.dna.prop2.,'\\.')), stringsAsFactors = FALSE)
names(names.1) <- c("Station.Type", "Station", "Type")

# # use Station instead of station for easier comparisons with morph info
# siteTable <- read.csv('Sites.csv', header=T)
# names(siteTable)[1] <- "Station"

# Add location
names.2 <- merge(names.1, sites, by="Station", all.x=TRUE)
rownames(names.2) <- names.2$Station.Type
names.2$Station.Type <- NULL

# Grab sites/species scores from NMDS output
df <- data.frame(vegan::scores(nmds3, display = "sites"))

# look at species too?
df.2 <- data.frame(vegan::scores(nmds3, display = "species"))
df.2 <- setDT(df.2, keep.rownames = TRUE)[]
names(df.2)[1] <- "GlobalESV"

# centroids for locations
set.seed(12345)
centroids <- envfit(nmds3, names.2[,c(1,3)], perm = 999, na.rm=TRUE)
factors.centroids <- data.frame(centroids[2]$factors[1])
factors.centroids <- setDT(factors.centroids, keep.rownames = TRUE)[]
names(factors.centroids)[1] <- "Location"
factors.centroids <- factors.centroids[c(26:34),]
names(factors.centroids)[2:3] <- c("NMDS1", "NMDS2")
factors.centroids$Location <- gsub("Location", "", factors.centroids$Location)

# centroids
factors.centroids$Location <- factor(factors.centroids$Location,
                       levels = c("Inner Harbour West","Bathurst Quay","Inner Harbour Centre","Queen Elizabeth Docks","Inner Harbour East","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"),
                       labels = c("Inner Harbour (West)","Bathurst Quay","Inner Harbour (Centre)","Queen Elizabeth Docks","Inner Harbour (East)","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"))


#######################################
# fit environmental variables
#######################################
set.seed(12345)
fit2 <- envfit(nmds3, envs[,-1], perm = 999, na.rm=TRUE)

# vectors
fit.vectors2 <- fit2[[1]]
fit.pvals2 <- fit2[[1]]$pvals
fit.r2 <- fit2[[1]]$r
vectors.df2 <- as.data.frame(fit.vectors2[[1]])
vectors.df2$pvals <- fit.pvals2
vectors.df2$r <- fit.r2
vectors.df.sig2 <- vectors.df2[vectors.df2$pvals <0.05 &
                               vectors.df2$r >= 0.30,]
vectors.df.sig2$x <- 0
vectors.df.sig2$y <- 0

# plot to visualize how close DNA and morph are
#############################################
# NMDS color by location and fit sediment vars
# chulls.location <- ddply(gg5, .(Location), function(gg5) gg5[chull(gg5$NMDS1, gg5$NMDS2), ])

# NMDS plot, color by site
sv2 <- ggplot(data=factors.centroids, aes(x=NMDS1, y=NMDS2)) + 
  geom_text(aes(label=Location), size=3) +
  # geom_point(data=gg5, aes(x=NMDS1, y=NMDS2, col=Location, fill=Location), size = 3, show.legend = TRUE) +
  # geom_text(data=gg, aes(x=NMDS1, y=NMDS2, col=Location, label=Station), show.legend = FALSE) +
  # geom_polygon(data=chulls.location, aes(x=NMDS1, y=NMDS2, fill=Location), alpha=0.5, show.legend = FALSE) +
  ggtitle("") +
  coord_cartesian(x = c(-2.25, 1), y = c(-1.0, 1.5))+
  geom_segment(data=vectors.df.sig2, 
               mapping=aes(x=x,y=y,xend=NMDS1, yend=NMDS2), 
               arrow=arrow(angle=25, length=unit(0.10, "inches")), 
               size=0.5,
               color="black") +
  geom_text(vectors.df.sig2,
            mapping=aes(x=NMDS1+0.03, y=NMDS2+0.03, label=rownames(vectors.df.sig2)),
             position=position_jitter(width = 0.15, height = 0.15),
            # hjust="outward", vjust="outward", 
            size=2.5, color="black") +
  # annotate("text", x=Inf, y=Inf, label="Stress = 0.06, R2 = 0.99", hjust = 1.1, vjust = 1.3, size=2.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size=9),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=5),
    legend.position="none") +
  guides(fill=guide_legend(ncol=1),
         color=NULL)
sv2




# Create distance matrix based on P-A data using binary Bray Curtis (Sorensen) dissimilarity
# sort dna.morph so that the order is the same as in gg where the metadata is
sor2 <- vegdist(sv.dna.prop2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
target2 <- rownames(sv.dna.prop2)

 # Put it all in one df for ggplot
 gg5 <- merge(df, names.2, by="row.names")

# create factors
 gg5$Station <- factor(gg5$Station,
                       levels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","1772","6699"),
                       labels = c("1012","1346","1349","1355","1358","1359","1362","1363","1364","1365","1366","1370","1371","1372","1375","1379A","1380","1381","1384","1401","1403","1761","1765A","1772","6699"))
 gg5$Location <- factor(gg5$Location,
                                     levels = c("Inner Harbour West","Bathurst Quay","Inner Harbour Centre","Queen Elizabeth Docks","Inner Harbour East","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"),
                                     labels = c("Inner Harbour (West)","Bathurst Quay","Inner Harbour (Centre)","Queen Elizabeth Docks","Inner Harbour (East)","Keating Channel","Ship Channel","Ship Turning Basin","Toronto Island"))


gg5 <- gg5[match(target, gg5$Row.names),]
# ensure gg$Type is in the same order as sor
bd.location <- betadisper(sor2, as.factor(gg5$Location))

# check for heterogeneity of beta dispersions within groups 
set.seed(1234)
anova(bd.location) # 3.664e-08 ***

pdf("Outfiles/supportingFiles/18S_sv_betadisp.pdf")
boxplot(bd.location, main="Type")
dev.off()


# Use ADONIS to check for significant interactions between disturbance and development stage (within layers)
adonis(sor2 ~ Location, data=gg5, permutations=999)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location   8    3.0050 0.37563  1.6027 0.44487  0.001 ***
# Residuals 16    3.7499 0.23437         0.55513           
# Total     24    6.7549                 1.00000   








top <- plot_grid(pca, tax, sv, sv2, labels="auto")

ggsave("Outfiles/Fig3_PCA_NMDS.pdf", top, width = 8, height = 8)


