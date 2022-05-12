# Teresita M. Porter, May 6, 2022
# Fig 4 relative abundance

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('R_scripts/Fig4_relative_abundance.R')

library(groundhog)

groundhog.library("ggpubr", '2022-05-01', tolerate.R.version='4.1.1') # stats_compare_means
# comes with ggplot
groundhog.library("stringr", '2022-05-04', tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", '2022-05-04', tolerate.R.version='4.1.1') # dcast
groundhog.library("data.table", '2022-05-04', tolerate.R.version='4.1.1') # setDT
groundhog.library("cowplot", '2022-05-04', tolerate.R.version='4.1.1')

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

# Family
COI7 <- COI6[!COI6$Family=="",]

# remove extra columns 
family.dna <- COI7[,c("Station", "Family", "ESVsize")]

# pivot to make sample x ESV matrix 
family.dna.df <- reshape2::dcast(family.dna, Station ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each family, per sample
# denominator should be total number of reads per station even if not tax assigned
denominator <- data.frame(COI6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotalReads"
family.dna.pct <- family.dna.df[,-1]/denominator$TotalReads * 100
family.dna.pct$Station <- family.dna.df$Station

# edit location to specify DNA samples
family.dna.pct$Station <- paste(family.dna.pct$Station, "COI", sep=".")

df <- reshape2::melt(family.dna.pct, id.vars=c("Station"))
names(df)[2:3] <- c("Family", "RelAbund")

# figure out top 5 families
df2 <- data.frame(df %>% group_by(Family) %>% dplyr::summarize(sum(RelAbund)))
df2 <- df2[order(-df2$sum.RelAbund.),]
keep <- head(df2$Family, 5)

COI.family <- df[df$Family %in% keep,]
COI.family$Type <- "COI"


# now bring in morph, convert to prop, and combine with dna for nmds
#############################################
# Morphology, family rank, pooled replicates
# Parse morphology results

trad <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# drop entries without a Family level identification
trad2 <- trad[trad$Family!="",]

# rename 1765 from morph as 1765A to match DNA and siteTable
trad2$station <- gsub("1765", "1765A", trad2$station)

# pivot to make ESV matrix 
family.morph.df <- reshape2::dcast(trad2, station ~ Family, value.var = "individuals", fun.aggregate = sum)

# convert to proportions, for each station, show proportion of individuals recovered per family
# denominator should be all individuals even if not assigned to family
denominator <- data.frame(trad %>% group_by(station) %>% dplyr::summarize(sum(individuals)))
names(denominator)[2] <- "TotalIndv"
family.morph.pct <- family.morph.df[,-1]/denominator$TotalIndv * 100
family.morph.pct$Station <- family.morph.df$station

# edit show show stations from morph
family.morph.pct$Station <- paste(family.morph.pct$Station, "Morph", sep=".")

df3 <- reshape2::melt(family.morph.pct, id.vars=c("Station"))
names(df3)[2:3] <- c("Family", "RelAbund")

# figure out top 5 families
df4 <- data.frame(df3 %>% group_by(Family) %>% dplyr::summarize(sum(RelAbund)))
df4 <- df4[order(-df4$sum.RelAbund.),]
keep <- head(df4$Family,5)

trad.family <- df3[df3$Family %in% keep,]
trad.family$Type <- "Morphology"

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

# genus 0.7, family 0.2 (80% correct)
sv.dna3 <- S6[S6$fBP >= 0.2, ]
# remove if undef
sv.dna3 <- sv.dna3[!grepl("undef", sv.dna3$Family),]

# pivot to make sample x ESV matrix 
sv.dna.df2 <- reshape2::dcast(sv.dna3, Station ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# Convert to proportions, proportion of reads assigned to each sv, per sample
# denominator should be total number of reads even if can't id to family
denominator <- data.frame(S6 %>% group_by(Station) %>% dplyr::summarize(sum(ESVsize)))
names(denominator)[2] <- "TotalReads"
sv.dna.pct <- sv.dna.df2[,-1]/denominator$TotalReads * 100
sv.dna.pct$Station <- sv.dna.df2$Station

# melt for ggplot
ssu.df <- reshape2::melt(sv.dna.pct, id.vars=c("Station"))
names(ssu.df)[2:3] <- c("Family", "RelAbund")

# filter to jsut keep the top 5 most abundant
ssu.df.summary <- data.frame(ssu.df %>% group_by(Family) %>% dplyr::summarize(sum(RelAbund)))
ssu.df.summary <- ssu.df.summary[order(-ssu.df.summary$sum.RelAbund.),]
keep <- head(ssu.df.summary$Family,5)

ssu.df2 <- ssu.df[ssu.df$Family %in% keep,]

# figure out what the top 5 are
ssu.df3 <- data.frame(ssu.df %>% group_by(Family) %>% dplyr::summarize(sum(RelAbund)))
ssu.df3 <- ssu.df3[order(-ssu.df3$sum.RelAbund.),]
keep <- head(ssu.df3$Family,5)

SSU.family <- ssu.df2[ssu.df2$Family %in% keep,]
SSU.family$Type="18S"












# put them together
all <- rbind(COI.family, trad.family, SSU.family)

# sort families by relative abundance
all2 <- data.frame(all %>% group_by(Family) %>% dplyr::summarize(sum(RelAbund)))
all2 <- all2[order(all2$sum.RelAbund.), ]

# create factors
all$Family <- factor(all$Family, levels=all2$Family)

# check for normality
hist(all$RelAbund)
shapiro.test(all$RelAbund)
# W = 0.50461, p-value < 2.2e-16
# use wilcox test

# compare means between just the macroinverts
macroinverts <- all[all$Type=="COI" |
                      all$Type=="Morphology",]
comps <- data.frame(compare_means(RelAbund ~ Type, macroinverts, group.by = "Family", method="wilcox.test"))
# Family      .y.     group2 group1            p   p.adj p.format p.signif   method
# 1 Chironomidae RelAbund Morphology    COI 3.894625e-02 3.9e-02    0.039        * Wilcoxon
# 2     Naididae RelAbund Morphology    COI 3.638461e-10 7.3e-10  3.6e-10     **** Wilcoxon


# check chironomids
chi <-all[all$Family=="Chironomidae",]
chi %>% group_by(Type) %>% dplyr::summarize(median(RelAbund))
# Type       `median(RelAbund)`
# <chr>                   <dbl>
#   1 COI                      9.19
# 2 Morphology               4.64

# check naididae
nai <-all[all$Family=="Naididae",]
nai %>% group_by(Type) %>% dplyr::summarize(median(RelAbund))
# Type       `median(RelAbund)`
# <chr>                   <dbl>
#   1 COI                      31.8
# 2 Morphology               90.6

# assess Pearson corrleations across stations
# rownames still match
dna <- all[all$Type=="COI",]
morph <- all[all$Type=="Morphology",]

# Calculating Pearson's product-moment correlation
p <- cor.test(dna$RelAbund, morph$RelAbund, method = "pearson", conf.level = 0.95)
p.corr <- round(as.numeric(p[4]), 3)
# 0.707
p.pval <- round(as.numeric(p[3]), 10)
# 0.0
# 95 percent confidence interval:
# 0.6062252 0.7848300

# plot 


p <- ggplot(all, aes(x=Family, y=RelAbund)) + 
  geom_boxplot() +
  coord_flip() +
  facet_wrap(~Type, scales="free_y") +
  # geom_point(position=position_jitterdodge(jitter.width = 0.15), size = 0.5) +
  labs(y="", x="") +
  scale_fill_manual(values=c("dark grey", "white")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank()) 
p














############################
# Functional relative abundances

F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE) # already limited to Arthropoda?
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE) # includes non-arthropods too

# edit column headers so matrices can be merged
# GlobalESV
names(F)[1] <- "GlobalESV"
names(M)[1] <- "GlobalESV"

# Seqfield
names(F)[4] <- "ORF/ESVSeq"
names(M)[4] <- "ORF/ESVSeq"

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

# Throw out optimization samples
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

# work with individuals for now, convert to proportions before combining with morph

# pivot
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

# normally, each would be < 100, 
# but because some taxa can be assigned to more than one FFG,
# may be > 100
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
# setDT(morph.t, keep.rownames = TRUE)[]
# names(morph.t)[1] <- "Family"


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

# convert FFG reads to proportions for each station
# denominator should be total number of indivdiuals even if not assigned
denominator <- data.frame(B %>% group_by(station) %>% dplyr::summarize(sum(individuals)))
names(denominator)[2] <- "RelAbund"

morph.FFG.pct <- morph.FFG.df/denominator$RelAbund * 100

# normally would sum to < 1
# but because some families can be assigned to more than one FFG
# can be ? one
rowSums(morph.FFG.pct)
# good

# put both matrices together
setDT(DNA.FFG.pct, keep.rownames = TRUE)[]
names(DNA.FFG.pct)[1] <- "Station"
DNA.FFG.pct$Type <- "DNA"

setDT(morph.FFG.pct, keep.rownames = TRUE)[]
names(morph.FFG.pct)[1] <- "Station"
morph.FFG.pct$Type <- "Morphology"

all2 <- rbind(DNA.FFG.pct, morph.FFG.pct)

# melt for ggplot
all.long <- reshape2::melt(all2, id.vars=c("Station", "Type"))

# fix Station
all.long$Station <- gsub("X", "", all.long$Station)

# check for normality
hist(all.long$value)
shapiro.test(all.long$value)
# W = 0.66105, p-value < 2.2e-16
# use Wilcox test

# sort by relative abundance
all.long2 <- data.frame(all.long %>% group_by(variable) %>% dplyr::summarize(sum(value)))
all.long2 <- all.long2[order(all.long2$sum.value.),]

# create factor
all.long$variable <- factor(all.long$variable, levels=all.long2$variable)

# compare means between just the macroinverts
comps2 <- data.frame(compare_means(value ~ Type, all.long, group.by = "variable", method="wilcox.test"))
# variable   .y. group1     group2            p   p.adj p.format p.signif   method
# 1       CF value    DNA Morphology 9.065050e-01 1.00000    0.907       ns Wilcoxon
# 2       CG value    DNA Morphology 8.626083e-01 1.00000    0.863       ns Wilcoxon
# 3       HB value    DNA Morphology 1.260844e-01 0.50000    0.126       ns Wilcoxon
# 4       PA value    DNA Morphology 3.757504e-05 0.00023  3.8e-05     **** Wilcoxon
# 5       PR value    DNA Morphology 1.828004e-05 0.00013  1.8e-05     **** Wilcoxon
# 6       SH value    DNA Morphology 6.975870e-01 1.00000    0.698       ns Wilcoxon
# 7    Other value    DNA Morphology 1.394648e-02 0.07000    0.014        * Wilcoxon


# check PA
pa <-all.long[all.long$variable=="PA",]
pa %>% group_by(Type) %>% dplyr::summarize(median(value))
# Type       `median(value)`
# <chr>                <dbl>
#   1 DNA                 0.0311
# 2 Morphology          0     

# check PR
pr <-all.long[all.long$variable=="PR",]
pr %>% group_by(Type) %>% dplyr::summarize(median(value))
# Type       `median(value)`
# <chr>                <dbl>
#   1 DNA                   10.8
# 2 Morphology            97.3

# check Other
oth <-all.long[all.long$variable=="Other",]
oth %>% group_by(Type) %>% dplyr::summarize(median(value))
# Type       `median(value)`
# <chr>                <dbl>
#   1 DNA                  0.182
# 2 Morphology           2.17 

# assess Pearson corrleations across stations
# rownames still match
dna2 <- all.long[all.long$Type=="DNA",]
morph2 <- all.long[all.long$Type=="Morphology",]

# Calculating Pearson's product-moment correlation
p2 <- cor.test(dna2$value, morph2$value, method = "pearson", conf.level = 0.95)
p.corr2 <- round(as.numeric(p2[4]), 3)
# 0.078
p.pval2 <- round(as.numeric(p2[3]), 10)
# 0.3063534
# 95 percent confidence interval:
#  -0.07140243  0.22352892


# plot
p2 <- ggplot(all.long, aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(y="Relative Abundance (%)", x="") +
  coord_flip() +
  facet_wrap(~Type, scales="free_y") +
  scale_fill_manual(values=c("darkgrey", "white")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        # strip.background = element_blank(),
        # strip.text.x = element_text(angle = 0, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank()) 
p2


# bottom <- plot_grid(p2, NULL, rel_widths = c(2,1))

g <- plot_grid(p , p2, nrow=2, labels="auto", rel_heights = c(5,6))
ggsave("Outfiles/Fig4_relative_abundance.pdf", g, width = 8, height = 8)

