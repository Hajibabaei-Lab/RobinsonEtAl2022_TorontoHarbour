# Teresita Porter, May 6, 2022
# Fig 5 hierarchical partitioning

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('R_scripts/Fig5_hierarchical_partitioning.R')

library(groundhog)
groundhog.library("ggpubr", "2022-05-01", tolerate.R.version='4.1.1') # stat_compare_means

groundhog.library("hier.part", "2022-05-04", tolerate.R.version='4.1.1')
groundhog.library("vegan", "2022-05-04", tolerate.R.version='4.1.1')  #package written for vegetation analysis
groundhog.library("stats", "2022-05-04", tolerate.R.version='4.1.1')    # e.g for hclust() function
groundhog.library("stringr", "2022-05-04", tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", "2022-05-04", tolerate.R.version='4.1.1') # dcast
groundhog.library("data.table", "2022-05-04", tolerate.R.version='4.1.1') # setDT
groundhog.library("moments", "2022-05-04", tolerate.R.version='4.1.1') # skewness
groundhog.library("dplyr", "2022-05-04", tolerate.R.version='4.1.1') # group_by, summarize
groundhog.library("RColorBrewer", "2022-05-04", tolerate.R.version='4.1.1')

#############################################
# First clean up the COI macroinvertebrate data

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

# Track negative control samples
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





################################################
# Now get species richness for COI metabarcoding data

COI7 <- COI6[!COI6$Species=="",]

# remove extra columns 
species.dna <- COI7[,c("Station", "Species", "ESVsize")]

# pivot to make sample x ESV matrix 
species.dna.df <- reshape2::dcast(species.dna, Station ~ Species, value.var = "ESVsize", fun.aggregate = sum)

# calculate richness
mi.richness.dna <- data.frame(mi.richness.dna=specnumber(species.dna.df[,-1]))
rownames(mi.richness.dna) <- species.dna.df$Station


##############################################
# Now get family level diversity indicators

COI8 <- COI6[!COI6$Family=="",]

# remove extra columns 
family.dna <- COI8[,c("Station", "Family", "ESVsize")]

# pivot to make sample x ESV matrix 
family.dna.df <- reshape2::dcast(family.dna, Station ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# convert to percent
# denominator should be all reads per station (not just the ones assigned)
denominator <- data.frame(COI6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "sum"
family.dna.pct <- family.dna.df[,-1]/denominator$sum * 100

# sort by relabund and keep top 5
sort.coi <- data.frame(relabund=colSums(family.dna.pct), family=names(family.dna.pct))
sort.coi <- sort.coi[order(-sort.coi$relabund),]
target <- head(sort.coi$family,5)

mi.family.dna.pct <- family.dna.pct[,target]
rownames(mi.family.dna.pct) <- family.dna.df$Station
names(mi.family.dna.pct) <- paste(names(mi.family.dna.pct), "dna", sep=".")


##############################################
# Now get order level diversity indicators to calc EPT (not possible anymore after removing ESVs & samples from blanks)

# remove extra columns 
order.dna <- COI6[,c("Station", "Order", "ESVsize")]

# pivot to make sample x ESV matrix 
order.dna.df <- reshape2::dcast(order.dna, Station ~ Order, value.var = "ESVsize", fun.aggregate = sum)

# convert to percent
mi.order.dna.pct <- order.dna.df[,-1]/rowSums(order.dna.df[,-1]) * 100

# just keep Ephemeroptera, Plecoptera, Trichoptera 
# no longer present after removing blanks and optimization samples)
# not detected much using morph either so just leave this out





###########################################
# now clean up morphology data

trad <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# rename 1765 from morph as 1765A to match DNA and siteTable
trad$station <- gsub("1765", "1765A", trad$station)



################################################
# Now get species richness for morphology data

# drop entries without a Species level identification
trad.species <- trad[trad$Species!="",]

# pivot to make ESV matrix 
species.morph.df <- reshape2::dcast(trad.species, station ~ Species, value.var = "individuals", fun.aggregate = sum)

# get richness
mi.richness.morph <- data.frame(mi.richness.morph=specnumber(species.morph.df[,-1]))
rownames(mi.richness.morph) <- species.morph.df$station




#############################################
# Now get family level diversity indicators

# drop entries without a Family level identification
trad.family <- trad[trad$Family!="",]

# pivot to make ESV matrix 
family.morph.df <- reshape2::dcast(trad.family, station ~ Family, value.var = "individuals", fun.aggregate = sum)

# convert to percent
# denominator should be total number of reads from each station before filtering out taxa that couldn't be id'd to family
denominator <- data.frame(trad %>% group_by(station) %>% dplyr::summarize(sum(individuals)))
names(denominator) <- c("Station", "Individuals")
family.morph.pct <- family.morph.df[,-1]/denominator[,2] * 100

# sort by relabund and keep top 5
sort.morph <- data.frame(relabund=colSums(family.morph.pct), family=names(family.morph.pct))
sort.morph <- sort.morph[order(-sort.morph$relabund),]
target2 <- head(sort.morph$family,5)

mi.family.morph.pct <- family.morph.pct[,target2]
rownames(mi.family.morph.pct) <- family.morph.df$Station
names(mi.family.morph.pct) <- paste(names(mi.family.morph.pct), "morph", sep=".")

###########################################
# Calculate order level EPT (virtually nothing, only found from 2 stations,  0.09 - 0.2 % )

# drop entries without a Order level identification
trad.order <- trad[trad$Order!="",]

# pivot to make ESV matrix 
order.morph.df <- reshape2::dcast(trad.order, station ~ Order, value.var = "individuals", fun.aggregate = sum)

# convert to percent
order.morph.pct <- order.morph.df[,-1]/rowSums(order.morph.df[,-1]) * 100
order.morph.pct$Station <- order.morph.df$station

# Just keep Ephemeroptera & Trichoptera (no Plecoptera)
EPT.morph.pct <- order.morph.pct[,c("Ephemeroptera", "Trichoptera")]

# only detected at low level from two samples, don't bother using





#############################################
# Next clean up the 18S eukaryote (non-macroinvertebrate) data

# Read in COI.csv from MetaWorks v1.0.0
S <- read.csv(file="Infiles/results_18S_v4.1.csv", head=TRUE, stringsAsFactors = FALSE) # includes non-arthropods too

# edit column headers so matrices can be merged
# GlobalESV
names(S)[1] <- "GlobalESV"

# Seqfield
names(S)[4] <- "ORF/ESVSeq"

# SuperKingdom/Domain
names(S)[9] <- "SuperKingdom"
names(S)[10] <- "SuperKingdomRank"
names(S)[11] <- "skBP"

# add Species, SpeciesRank, sBP to S
S$Species <- ""
S$SpeciesRank <- ""
S$sBP <- 0

# keep all taxa for COI, but filter 18S
# sBP > 0.70, gBP > 0.20
S$Species <- ifelse(S$sBP >=0.70, S$Species, "")
S$Genus <- ifelse(S$gBP >=0.20, S$Genus, "")

# # combine into 1 matrix
# COI <- rbind(F, M)

# focus on eukaryotes (only a couple macroinverts found using this marker, so just exclude those)
morphPhyla <- c("Arthropoda", "Annelida", "Mollusca", "Cnidaria", "Platyhelminthes")
S2 <- S[!S$Phylum %in% morphPhyla,]

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
S2$SampleName <- gsub("TorontoHarbour_", "", S2$SampleName)

# Reformat negative control (blanks) filenames
S2$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", S2$SampleName)
S2$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", S2$SampleName)
S2$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", S2$SampleName)
S2$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", S2$SampleName)
S2$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", S2$SampleName)
S2$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", S2$SampleName)


# Track negative control samples
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- S2[(S2$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from neg controls and neg control samples from main dataset
S3 <- S2[!S2$GlobalESV %in% remove,]
S4 <- S3[!(S3$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_4_18S_S86", "Opt_1_18S_S83", "Opt_2_18S_S84", 
                         "Opt_5_18S_S87","Opt_3_18S_S85")
S5 <- S4[!(S4$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get station & replicate fields
S5 <- data.frame(S5, do.call(rbind, str_split(S5$SampleName,"_")), stringsAsFactors = FALSE)
names(S5)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
S5$Station <- gsub("1271","1371", S5$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
S6 <- data.frame(S5, do.call(rbind, str_split(S5$GlobalESV,"_")), stringsAsFactors = FALSE)
names(S6)[37:38] <- c("Amplicon", "Zotu")





################################################
# Now get genus richness for 18S metabarcoding data

# filter to ensure at least 80% correct, genus 0.70 (family 0.20)
S7 <- S6[S6$gBP >= 0.70,]

# remove extra columns 
genus.dna <- S7[,c("Station", "Genus", "ESVsize")]

# pivot to make sample x ESV matrix 
genus.dna.df <- reshape2::dcast(genus.dna, Station ~ Genus, value.var = "ESVsize", fun.aggregate = sum)

# calculate richness
euk.richness.dna <- data.frame(euk.richness.dna=specnumber(genus.dna.df[,-1]))
rownames(euk.richness.dna) <- genus.dna.df$Station


##############################################
# Now get family level diversity indicators

# filter to ensure at least 80% correct, (genus 0.70) family 0.20
S8 <- S6[S6$fBP >= 0.20,]

# remove if contains undef or \\d+ 
S9 <- S8[!grepl("undef|\\d+", S8$Family),]

# Handle Cryptomycota
S9$Family <- gsub("Cryptomycota_Incertae_Sedis_Incertae_Sedis_Incertae_Sedis", "Cryptomycota_Incertae_Sedis", S9$Family)

# remove extra columns 
family.dna2 <- S9[,c("Station", "Family", "ESVsize")]

# pivot to make sample x ESV matrix 
family.dna.df2 <- reshape2::dcast(family.dna2, Station ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# convert to percent
# denominator should be from matrix before filtering out taxa we couldn't identify to a good family S6 for each station
denominator <- data.frame(S6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotalReads" 
euk.family.dna.pct <- family.dna.df2[,-1]/denominator[,2] * 100

# just keep top 5 most abundant
tail(sort(colSums(euk.family.dna.pct)),5)
# Hypotrichia  Oligohymenophorea Thoracosphaeraceae 
# 67.73614           80.74114           97.34158 
# Haptoria        Prostomatea 
# 240.77569          407.18394 

# top 1 or maybe 3?
euk.family.dna.pct <- euk.family.dna.pct[,c("Prostomatea","Haptoria","Thoracosphaeraceae","Oligohymenophorea","Hypotrichia")]
rownames(euk.family.dna.pct) <- family.dna.df2$Station
names(euk.family.dna.pct) <- paste(names(euk.family.dna.pct) , "dna", sep=".")





########################################
# compile a matrix of diverersity metrics

diversity.metrics <- cbind(mi.richness.dna, mi.richness.morph, euk.richness.dna,
                           mi.family.dna.pct, mi.family.morph.pct, euk.family.dna.pct)

# print out for mapping
write.csv(diversity.metrics, "Outfiles/supportingFiles/diversity.metrics.csv", quote=FALSE)

# Standardize (entries are transformed relative to other entries) response variables (rows = Stations, cols = community metrics such as richness & indicator taxon occurrence )
diversity.metrics.std <- decostand (diversity.metrics, method = "hellinger")
diversity.metrics.std





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

hist(S.2$PCB)
shapiro.test(S.2$PCB)
skewness(S.2$PCB)
skewness(sqrt(S.2$PCB)) 
skewness(log10(S.2$PCB)) # this one

hist(S.2$Pesticides)
shapiro.test(S.2$Pesticides)
skewness(S.2$Pesticides)
skewness(sqrt(S.2$Pesticides)) # this one

hist(S.2$Rare_Earth)
shapiro.test(S.2$Rare_Earth)
skewness(S.2$Rare_Earth) 
skewness(sqrt(max(S.2$Rare_Earth+1)-S.2$Rare_Earth)) # this one

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




########################################
# Hierarchical partitioning of the variance in the diversity metrics
# from Caroline Emilson Aug. 26/20

# assess water vars first
response.v=names(diversity.metrics.std)
# initialize list to collect data
set.seed(1234)
results <- list()
for (i in response.v){
  hp<-hier.part(diversity.metrics.std[,i], water.envir[,-1], barplot=TRUE)
  # 10 reps for debugging, 100 reps for testing, 1000 for real
  rand<-rand.hp(diversity.metrics.std[,i], water.envir[,-1], num.reps=1000, gof="logLik")
  hp.results<-cbind(hp$IJ, hp$I.perc, rand$Iprobs)
  hp.results$response <- i
  setDT(hp.results, keep.rownames = TRUE)[]
  names(hp.results)[1] <- "predictors"
  results[[i]] <- hp.results
}

# turn list into df for ggplot
results.df <- rbindlist(results)

# rename col
names(results.df)[5] <- "I.perc"

# create I/J col
results.df$IdivJ <- results.df$I/results.df$J

# only keep if significant
results.df.sig <- results.df[results.df$sig95=="*",]

# only keep key cols
water.hp <- results.df.sig[,c("predictors", "I.perc", "IdivJ", "response")]
# order, descending
water.hp <- water.hp[order(-water.hp$I.perc),]

# save this under new name so doesn't get overwritten
WaterHP <- water.hp

# Print for later analysis, keeping all sig results even if J > I
write.csv(WaterHP, file="Outfiles/TableS5_water_diversity.csv", quote=FALSE, row.names=FALSE) # keeps the rownames

# only keep if independent contributions exceed joint, do this filtering before plotting
WaterHP <- WaterHP[abs(WaterHP$IdivJ)>=1.0,]





# assess sediment vars next
response.v=names(diversity.metrics.std)
# initialize list to collect data
set.seed(1234)
results <- list()
for (i in response.v){
  hp<-hier.part(diversity.metrics.std[,i], sed.envir[,-1], barplot=TRUE)
  # 10 reps for debugging, 100 reps for testing, 1000 for real
  rand<-rand.hp(diversity.metrics.std[,i], sed.envir[,-1], num.reps=1000, gof="logLik")
  hp.results<-cbind(hp$IJ, hp$I.perc, rand$Iprobs)
  hp.results$response <- i
  setDT(hp.results, keep.rownames = TRUE)[]
  names(hp.results)[1] <- "predictors"
  results[[i]] <- hp.results
}

# turn list into df for ggplot
results.df <- rbindlist(results)

# rename col
names(results.df)[5] <- "I.perc"

# create I/J col
results.df$IdivJ <- results.df$I/results.df$J

# only keep if significant
results.df.sig <- results.df[results.df$sig95=="*",]

# only keep key cols
sed.hp <- results.df.sig[,c("predictors", "I.perc", "IdivJ", "response")]
# order, descending
sed.hp <- sed.hp[order(-sed.hp$I.perc),]

# save this under new name so doesn't get overwritten
SedimentHP <- sed.hp

# Print for later analysis, keep all even if J > I
write.csv(SedimentHP, file="Outfiles/TableS6_sediment_diversity.csv", quote=FALSE, row.names=FALSE) # keeps the rownames

# only keep if independent contributions exceed joint
# SedimentHP <- SedimentHP[abs(SedimentHP$IdivJ)>=1.0,]



################
# # Visualize
# 
# # ###### start here if program crashes
# WaterHP <- read.csv("Outfiles/TableS5_water_diversity.csv", header=TRUE)
# SedimentHP <- read.csv("Outfiles/TableS6_sediment_diversity.csv", header=TRUE)

# prep to combine matrices
WaterHP$Predictor <- "Water"
SedimentHP$Predictor <- "Sediment"
# combine
combined <- rbind(WaterHP, SedimentHP)

# tweak for easier splitting
combined$response <- gsub("mi.richness", "mi_richness", combined$response)
combined$response <- gsub("euk.richness", "euk_richness", combined$response)

# Split to get separate type field
combined2 <- data.frame(combined, do.call(rbind, str_split(combined$response,"\\.")), stringsAsFactors = FALSE)
names(combined2)[6:7] <- c("DiversityIndicator", "Type")

# "Prostomatea","Haptoria","Thoracosphaeraceae","Oligohymenophorea","Hypotrichia")]

# edit dna to be COI or 18S and morph to be Morphology
combined2$Type <- ifelse(grepl("dna", combined2$Type) & grepl("mi_richness|Naididae|Candonidae|Cicadidae|Chironomidae|Dreissenidae|Pennaridae|Hydridae|Limnesiidae", combined2$DiversityIndicator), "COI", 
                    ifelse(grepl("dna", combined2$Type) & grepl("euk_richness|Prostomatea|Haptoria|Thoracosphaeraceae|Oligohymenophorea|Hypotrichia", combined2$DiversityIndicator), 
                         "18S", "Morphology"))

# sort predictors by descending median I.perc
summary <- data.frame(combined2 %>% group_by(predictors) %>% dplyr::summarize(median(I.perc)))
summary <- summary[order(-summary$median.I.perc.),]

# create factors
combined2$Predictor <- factor(combined2$Predictor, levels=c("Water", "Sediment"))
combined2$predictors <- factor(combined2$predictors, levels=summary$predictors)
combined2$Type <- factor(combined2$Type, levels=c("COI", "Morphology", "18S"))

# pick 3 colors
blues <- brewer.pal(9,"Blues")
# [1] "#F7FBFF" "#DEEBF7" "#C6DBEF" "#9ECAE1" "#6BAED6" "#4292C6" "#2171B5"
# [8] "#08519C" "#08306B"
oranges <- brewer.pal(9,"Oranges")
# [1] "#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913" "#D94801"
# [8] "#A63603" "#7F2704"

p1 <- ggplot(combined2, aes(x=response, y=I.perc, label=response, group=predictors)) + 
  geom_bar(stat="identity", aes(fill=predictors)) +
  # geom_text(aes(color=Type, hjust=-0.05, vjust=-0.2), position=position_jitter(width=0.3)) +
  # ylim(c(0,100)) +
  coord_flip() +
  labs(x="Diversity metrics", y="Predictor independent contribution (%)") +
  scale_fill_manual(values=c(blues[6:9], oranges[1:9])) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        # legend.position = c(0.5, 0.95),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
p1

# geom_bar()# plot it
# check normality
hist(combined2$I.perc)
shapiro.test(combined2$I.perc) # use wilcox

my_comparisons <- list(c("COI", "Morphology"), c("Morphology", "18S"), c("COI", "18S"))

p2 <- ggplot(combined2, aes(x=Type, y=I.perc)) + 
  geom_boxplot() +
  ylim(c(0,100)) +
  labs(x="Diversity metric responses", y="Variance explained (%)") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons = my_comparisons,
                     method="wilcox.test")
p2
# The independent contribution of predictors tended to explain
# higher levels of variance in COI diversity metrics












#################################
# Now get FFG's for COI metabarcodes

# Read in metabarcoding data from MetaWorks v1.0.0

F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE) # already limited to Arthropoda?
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE) # includes non-arthropods too
# B <- read.csv(file="results_BR5_OPTonly.csv", head=TRUE) # optimization samples
# S <- read.csv(file="results_18S_v4.1.csv", head=TRUE) # includes non-arthropods too

# edit column headers so matrices can be merged
# GlobalESV
names(F)[1] <- "GlobalESV"
names(M)[1] <- "GlobalESV"
# names(S)[1] <- "GlobalESV"

# Seqfield
names(F)[4] <- "ORF/ESVSeq"
names(M)[4] <- "ORF/ESVSeq"
# names(S)[4] <- "ORF/ESVSeq"

# combine into 1 matrix
A <- rbind(F, M)

# limit to macroinvertebrates phyla because need to calc %EPT and %Chironomid
# i.e. the ones picked up by morphology in 2018
morphPhyla <- c("Arthropoda", "Annelida", "Mollusca", "Cnidaria", "Platyhelminthes")
A.1 <- A[A$Phylum %in% morphPhyla,]

# filter for 95% correct species, 99% correct genus & up
A.1$Species <- ifelse(A.1$sBP >= 0.7, A.1$Species, "")
A.1$Genus <- ifelse(A.1$gBP >= 0.3, A.1$Genus, "")
A.1$Family <- ifelse(A.1$fBP >= 0.2, A.1$Family, "")

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
A.1$SampleName <- gsub("TorontoHarbour_", "", A.1$SampleName)

# Reformat negative control (blanks) filenames
A.1$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", A.1$SampleName)
A.1$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", A.1$SampleName)
A.1$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", A.1$SampleName)
A.1$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", A.1$SampleName)
A.1$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", A.1$SampleName)
A.1$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", A.1$SampleName)

# Track negative control samples
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- A.1[(A.1$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)

# Remove ESVs from blanks from the rest of the dataset
A.2 <- A.1[!A.1$GlobalESV %in% remove,]
# now remove the blanks
A.2 <- A.2[!A.2$SampleName %in% negativeControls,]

# Remove optimization samples
optimizationSamples <- c("Opt_3_COI_S80", "Opt_1_COI_S78", "Opt_2_COI_S79", 
                         "Opt_5_COI_S82","Opt_4_COI_S81","Opt_1_18S_S83",
                         "Opt_2_18S_S84", "Opt_3_18S_S85", "Opt_1_18S_S83")
A.3 <- A.2[!(A.2$SampleName %in% optimizationSamples),]

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.3 <- data.frame(A.3, do.call(rbind, str_split(A.3$SampleName,"_")), stringsAsFactors = FALSE)
names(A.3)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
A.3$Station <- gsub("1271","1371", A.3$Station)

# Split up GlobalESV with pkg 'stringr'
A.4<-data.frame(A.3, do.call(rbind, str_split(A.3$GlobalESV,"_")), stringsAsFactors = FALSE)
names(A.4)[37:38] <- c("Amplicon", "Zotu")

# filter
A.5 <- A.4[!A.4$Family=="",]

# pivot to make Family matrix 
pivot <- reshape2::dcast(A.5, Station ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(pivot) <- pivot$Station
pivot$Station <- NULL

# 61 families



############################
# map taxonomy, need family names
df.t <- as.data.frame(t(pivot))
# setDT(df.t, keep.rownames = TRUE)[]
# names(df.t)[1] <- "GlobalESV"






########
# add FFG for families

# read in FFG info for each family (for now, only based on EPA freshwater dbase)
# wait to see if EU dbase account is approved then fill in more blanks
EPA <- read.csv('Infiles/EPA_families.csv', header=T)

# get target list of taxa
keep_for_EPA_check <- names(pivot)

# only annotated 10 families using EPA
FFG.1 <- EPA[EPA$Family %in% keep_for_EPA_check,]

# Import Moog FFGs from Freshwaterecology.info
Moog <- read.csv('Infiles/Freshwaterecol_Moog.csv', header=T)

# get target list of taxa
# remove families that I already have EPA FFG for
keep_for_Moog_check <- keep_for_EPA_check[!(keep_for_EPA_check %in% FFG.1$Family)]

# annotated 12 more families using Moog
FFG.2 <- Moog[Moog$Family %in% keep_for_Moog_check,]

# Import Tachet FFGs from Freshwaterecology.info
Tachet <- read.csv('Infiles/Freshwaterecol_Tachet.csv', header=T)

# get target list of taxa
# remove families that I already have EPA + Moog FFG for
keep_for_Tachet_check <- keep_for_Moog_check[!(keep_for_Moog_check %in% FFG.2$Family)]

# annotated 1 more families using Tachet
FFG.3 <- Tachet[Tachet$Family %in% keep_for_Tachet_check,]

# combine these FFG annotations and ignore what we can't annotate in an automated manner
# add 'Other' FFG category (was only present with Moog)
FFG.1$Other <- NA
FFG.3$Other <- NA
DNA.FFG <- rbind(FFG.1, FFG.2, FFG.3)

# set NA to 0
DNA.FFG[is.na(DNA.FFG)] <- 0

# 58 families annotated out of 382 (!) 0.15 or 15% annotated, poor
# only proceed with the 58 families we could functionally assign

# add FFG info to metabarcoding data

# matrix multiplication
# df.tax2 x FFG (do each FFG type at a time)
# get total number of reads associated with an FFG type per station
DNA.ffg.CF <- df.t * DNA.FFG[,2]
sum.ffg.CF <- colSums(DNA.ffg.CF)
DNA.ffg.CG <- df.t * DNA.FFG[,3]
sum.ffg.CG <- colSums(DNA.ffg.CG)
DNA.ffg.HB <- df.t * DNA.FFG[,4]
sum.ffg.HB <- colSums(DNA.ffg.HB)
DNA.ffg.PA <- df.t * DNA.FFG[,5]
sum.ffg.PA <- colSums(DNA.ffg.PA)
DNA.ffg.PR <- df.t * DNA.FFG[,6]
sum.ffg.PR <- colSums(DNA.ffg.PR)
DNA.ffg.SH <- df.t * DNA.FFG[,7]
sum.ffg.SH <- colSums(DNA.ffg.SH)
DNA.ffg.Other <- df.t * DNA.FFG[,8]
sum.ffg.Other <- colSums(DNA.ffg.Other)

# for each station record sum of reads belonging to each FFG type
DNA.FFG.df <- data.frame(CF=sum.ffg.CF, CG=sum.ffg.CG, HB=sum.ffg.HB, PA=sum.ffg.PA, PR=sum.ffg.PR, SH=sum.ffg.SH, Other=sum.ffg.Other)

# convert FFG reads to pct for each station
# denominator should be total number of reads even if we couldn't assign to family
denominator <- data.frame(A.4 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <-  "sum"
DNA.FFG.pct <- DNA.FFG.df/denominator$sum * 100

# sanity check, each row should sum to < 1 because 
# we were only able to assign function to a fraction of families, 
# and ESVs to families
rowSums(DNA.FFG.pct)
# good




#############################################
# Parse morphology results
# Read in Morphology counts
B <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# rename 1765 in morph to match 1765A in DNA
B$station <- gsub("1765", "1765A", B$station)

# drop entries without a Family level identification
B.2 <- B[B$Family!="",]


# pivot to make ESV matrix 
morph <- reshape2::dcast(B.2, station ~ Family, value.var = "individuals", fun.aggregate = sum)
rownames(morph) <- morph$station
morph$station <- NULL



# family FFG

# get target list of tax 29
keep_for_EPA_check <- names(morph[,-1])

# only annotated 13 families using EPA
FFG.4 <- EPA[EPA$Family %in% keep_for_EPA_check,]

# get target list of taxa
# remove families that I already have EPA FFG for
keep_for_Moog_check <- keep_for_EPA_check[!(keep_for_EPA_check %in% FFG.4$Family)]

# annotated 11 more families using Moog
FFG.5 <- Moog[Moog$Family %in% keep_for_Moog_check,]

# get target list of taxa
# remove families that I already have EPA + Moog FFG for
keep_for_Tachet_check <- keep_for_Moog_check[!(keep_for_Moog_check %in% FFG.5$Family)]

# annotated 0 more families using Tachet
FFG.6 <- Tachet[Tachet$Family %in% keep_for_Tachet_check,]

# combine these FFG annotations and ignore what we can't annotate in an automated manner
FFG.4$Other <- NA
morph.FFG <- rbind(FFG.4, FFG.5)
# set NA to 0
morph.FFG[is.na(morph.FFG)] <- 0



# annotated 24 out of 29 morph families

# add FFG info to morph data
morph.t <- data.frame(t(morph))
# data.frame(setDT(morph.t, keep.rownames = TRUE)[])
# names(morph.t)[1] <- "Family"
# rownames(morph.t) <- morph.t$Family


# add FFG info to metabarcoding data

# matrix multiplication
# df.tax2 x FFG (do each FFG type at a time)
# get total number of reads associated with an FFG type per station
morph.ffg.CF <- morph.t * morph.FFG[,2]
sum.ffg.CF <- colSums(morph.ffg.CF)
morph.ffg.CG <- morph.t * morph.FFG[,3]
sum.ffg.CG <- colSums(morph.ffg.CG)
morph.ffg.HB <- morph.t * morph.FFG[,4]
sum.ffg.HB <- colSums(morph.ffg.HB)
morph.ffg.PA <- morph.t * morph.FFG[,5]
sum.ffg.PA <- colSums(morph.ffg.PA)
morph.ffg.PR <- morph.t * morph.FFG[,6]
sum.ffg.PR <- colSums(morph.ffg.PR)
morph.ffg.SH <- morph.t * morph.FFG[,7]
sum.ffg.SH <- colSums(morph.ffg.SH)
morph.ffg.Other <- morph.t * morph.FFG[,8]
sum.ffg.Other <- colSums(morph.ffg.Other)

# for each station record sum of reads belonging to each FFG type
morph.FFG.df <- data.frame(CF=sum.ffg.CF, CG=sum.ffg.CG, HB=sum.ffg.HB, PA=sum.ffg.PA, PR=sum.ffg.PR, SH=sum.ffg.SH, Other=sum.ffg.Other)

# convert FFG reads to pct for each station
morph.FFG.pct <- morph.FFG.df/rowSums(morph.FFG.df) * 100

# sanity check, each row should sum to 1
rowSums(morph.FFG.pct)
# good

########################################
# compile a matrix of functional metrics

# fix names
names(DNA.FFG.pct) <- c("CF.dna", "CG.dna", "HB.dna", "PA.dna", "PR.dna", "SH.dna", "Other.dna")
names(morph.FFG.pct) <- c("CF.morph", "CG.morph", "HB.morph", "PA.morph", "PR.morph", "SH.morph", "Other.morph")

functional.metrics <- cbind(DNA.FFG.pct, morph.FFG.pct)

# print out for mapping
write.csv(functional.metrics, "Outfiles/supportingFiles/functional.metrics.csv", quote=FALSE)

# Standardize (entries are transformed relative to other entries) response variables (rows = Stations, cols = community metrics such as richness & indicator taxon occurrence )
functional.metrics.std <- decostand (functional.metrics, method = "hellinger")
functional.metrics.std


# test for differences between means
cf <- data.frame(cbind(functional.metrics$CF.dna, functional.metrics$CF.morph))
names(cf) <- c("DNA", "morph")
cf.melt <- reshape2::melt(cf)
compare_means(value ~ variable, cf.melt, method="wilcox.test")
# n/s

cg <- data.frame(cbind(functional.metrics$CG.dna, functional.metrics$CG.morph))
names(cg) <- c("DNA", "morph")
cg.melt <- reshape2::melt(cg)
compare_means(value ~ variable, cg.melt, method="wilcox.test")
# n/s

hb <- data.frame(cbind(functional.metrics$HB.dna, functional.metrics$HB.morph))
names(hb) <- c("DNA", "morph")
hb.melt <- reshape2::melt(hb)
compare_means(value ~ variable, hb.melt, method="wilcox.test")
# n/s

pa <- data.frame(cbind(functional.metrics$PA.dna, functional.metrics$PA.morph))
names(pa) <- c("DNA", "morph")
pa.melt <- reshape2::melt(pa)
compare_means(value ~ variable, pa.melt, method="wilcox.test")
# *** p.adj 0.000034

pr <- data.frame(cbind(functional.metrics$PR.dna, functional.metrics$PR.morph))
names(pr) <- c("DNA", "morph")
pr.melt <- reshape2::melt(pr)
compare_means(value ~ variable, pr.melt, method="wilcox.test")
# * p.adj 0.016

sh <- data.frame(cbind(functional.metrics$SH.dna, functional.metrics$SH.morph))
names(sh) <- c("DNA", "morph")
sh.melt <- reshape2::melt(sh)
compare_means(value ~ variable, sh.melt, method="wilcox.test")
# n/s

other <- data.frame(cbind(functional.metrics$Other.dna, functional.metrics$Other.morph))
names(other) <- c("DNA", "morph")
other.melt <- reshape2::melt(other)
compare_means(value ~ variable, other.melt, method="wilcox.test")
# * p.adj 0.035


########################################
# Hierarchical partitioning of the variance in the functional metrics
# from Caroline Emilson Aug. 26/20

# assess water vars first
response.v=names(functional.metrics.std)
# initialize list to collect data
set.seed(1234)
results <- list()
for (i in response.v){
  hp<-hier.part(functional.metrics.std[,i], water.envir[,-1], barplot=TRUE)
  # 10 reps for debugging, 100 reps for testing, 1000 for real
  rand<-rand.hp(functional.metrics.std[,i], water.envir[,-1], num.reps=1000, gof="logLik")
  hp.results<-cbind(hp$IJ, hp$I.perc, rand$Iprobs)
  hp.results$response <- i
  setDT(hp.results, keep.rownames = TRUE)[]
  names(hp.results)[1] <- "predictors"
  results[[i]] <- hp.results
}

# turn list into df for ggplot
results.df <- rbindlist(results)

# rename col
names(results.df)[5] <- "I.perc"

# create I/J col
results.df$IdivJ <- results.df$I/results.df$J

# only keep if significant
results.df.sig <- results.df[results.df$sig95=="*",]

# only keep key cols
water.hp2 <- results.df.sig[,c("predictors", "I.perc", "IdivJ", "response")]
# order, descending
water.hp2 <- water.hp2[order(-water.hp2$I.perc),]

# save this under new name so doesn't get overwritten
WaterHP2 <- water.hp2

# Print for later analysis, keep all sig even if J > I
write.csv(WaterHP2, file="Outfiles/TableS5_water_functional.csv", quote=FALSE, row.names=FALSE) # keeps the rownames

# only plot if independent contributions exceed joint
WaterHP2 <- WaterHP2[abs(WaterHP2$IdivJ)>=1.0,]





# assess sediment vars next
response.v=names(functional.metrics.std)
# initialize list to collect data
set.seed(1234)
results <- list()
for (i in response.v){
  hp<-hier.part(functional.metrics.std[,i], sed.envir[,-1], barplot=TRUE)
  # 10 reps for debugging, 100 reps for testing, 1000 for real
  rand<-rand.hp(functional.metrics.std[,i], sed.envir[,-1], num.reps=1000, gof="logLik")
  hp.results<-cbind(hp$IJ, hp$I.perc, rand$Iprobs)
  hp.results$response <- i
  setDT(hp.results, keep.rownames = TRUE)[]
  names(hp.results)[1] <- "predictors"
  results[[i]] <- hp.results
}

# turn list into df for ggplot
results.df <- rbindlist(results)

# rename col
names(results.df)[5] <- "I.perc"

# create I/J col
results.df$IdivJ <- results.df$I/results.df$J

# only keep if significant
results.df.sig <- results.df[results.df$sig95=="*",]

# only keep key cols
sed.hp2 <- results.df.sig[,c("predictors", "I.perc", "IdivJ", "response")]
# order, descending
sed.hp2 <- sed.hp2[order(-sed.hp2$I.perc),]

# save this under new name so doesn't get overwritten
SedimentHP2 <- sed.hp2
# Print for later analysis, keep all sig even if J > I
write.csv(SedimentHP2, file="Outfiles/TableS6_sediment_functional.csv", quote=FALSE, row.names=FALSE) # keeps the rownames

# only plot if independent contributions exceed joint
# SedimentHP2 <- SedimentHP2[abs(SedimentHP2$IdivJ)>=1.0,]







################
# Visualize

# # ###### start here if program crashes
# WaterHP2 <- read.csv("Outfiles/TableS5_water_functional.csv", header=TRUE)
# SedimentHP2 <- read.csv("Outfiles/TableS6_sediment_functional.csv", header=TRUE)

# prep to combine matrices
WaterHP2$Predictor <- "Water"
SedimentHP2$Predictor <- "Sediment"
# combine
combined3 <- rbind(WaterHP2, SedimentHP2)

# Split to get separate type field
combined4 <- data.frame(combined3, do.call(rbind, str_split(combined3$response,"\\.")), stringsAsFactors = FALSE)
names(combined4)[6:7] <- c("FunctionalIndicator", "Type")

# edit dna to be COI or 18S and morph to be Morphology
combined4$Type <- ifelse(grepl("dna", combined4$Type), "DNA", "Morphology")

# sort predictors by descending median I.perc
summary2 <- data.frame(combined4 %>% group_by(predictors) %>% dplyr::summarize(median(I.perc)))
summary2 <- summary2[order(-summary2$median.I.perc.),]

# create factors
combined4$Predictor <- factor(combined4$Predictor, levels=c("Water", "Sediment"))
combined4$predictors <- factor(combined4$predictors,
                               levels=summary2$predictors)
combined4$Type <- factor(combined4$Type, levels=c("DNA", "Morphology"),
                         labels=c("COI", "Morphology"))



p3 <- ggplot(combined4, aes(x=response, y=I.perc, label=response, group=predictors)) + 
  geom_bar(stat="identity", aes(fill=predictors)) +
  # geom_text(aes(color=Type, hjust=-0.05, vjust=-0.2), position=position_jitter(width=0.3)) +
  # ylim(c(0,100)) +
  coord_flip() +
  labs(x="Functional metrics", y="Predictor independent contribution (%)") +
  # scale_color_manual(values=c(mypalette[c(2,3,2)])) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        # legend.position = c(0.5, 0.95),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2))
p3

# plot it
# check normality
hist(combined4$I.perc)
shapiro.test(combined2$I.perc) # use wilco      x

my_comparisons <- list(c("COI", "Morphology"))


# plot it
p4 <- ggplot(combined4, aes(x=Type, y=I.perc)) + 
  geom_boxplot() +
  ylim(c(1,100)) +
  labs(x="Functional metric responses", y="Variance explained (%)") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  stat_compare_means(comparisons=my_comparisons, label.y = 95,
    method="wilcox.test")
p4
# The independent contribution of predictors tended to explain
# higher levels of variance in COI diversity metrics





# put combined2 and combined4 together to get a harmonized legend
names(combined2)[6] <- "Indicator" 
names(combined4)[6] <- "Indicator" 
combined2$ResponseType <- "Diversity"
combined4$ResponseType <- "Functional"
combined5 <- rbind(combined2, combined4)


combined5$predictors <- factor(combined5$predictors,
                               levels=c(unique(combined5$predictors[combined5$Predictor=="Water"]), 
                                        unique(combined5$predictors[combined5$Predictor=="Sediment"])))

names(combined5)[1] <- "Predictors"

p5 <- ggplot(combined5, aes(x=response, y=I.perc, label=response, group=Predictors)) + 
  geom_bar(stat="identity", aes(fill=Predictors)) +
  # geom_text(aes(color=Type, hjust=-0.05, vjust=-0.2), position=position_jitter(width=0.3)) +
  # ylim(c(0,100)) +
  facet_wrap(~ResponseType, scales="free_y") +
  coord_flip() +
  labs(x="Response metrics", y="Predictor independent contribution (%)") +
  scale_fill_manual(values=c(blues[5:9], oranges[2:9])) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        strip.text= element_blank(),
        # legend.title = element_blank()
        ) +
  guides(fill = guide_legend(ncol = 1))
p5


ggsave("Outfiles/Fig5_hierarchical_partitioning.pdf", p5, width = 8, height = 8)


