# Teresita M. Porter, May 5, 2022
# Fig 2 alpha and gamma diversity

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('Fig2_richness.R')

library(groundhog)

groundhog.library("stringr", '2022-05-03', tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", '2022-05-03', tolerate.R.version='4.1.1') # dcast
groundhog.library("vegan", '2022-05-03', tolerate.R.version='4.1.1') # rrarefy
groundhog.library("ggplot2", '2022-05-03', tolerate.R.version='4.1.1') # ggplot
groundhog.library("cowplot", '2022-05-03', tolerate.R.version='4.1.1')
groundhog.library("ggpubr", '2022-05-03', tolerate.R.version='4.1.1') # stat_compare_means

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

# filter genus 0.70bp 80% correct
S7 <- S6[S6$gBP >= 0.70,]

# remove extra columns 
genus.dna <- S7[,c("Station", "Genus", "ESVsize")]

# pivot to make sample x ESV matrix 
genus.dna.df <- reshape2::dcast(genus.dna, Station ~ Genus, value.var = "ESVsize", fun.aggregate = sum)







# create richness df
# median richness
richness.species.dna <- specnumber(species.dna.df)
median.species.dna <- median(richness.species.dna)
# sd.species.dna <- sd(richness.species.dna)

richness.species.morph <- specnumber(species.morph.df)
median.species.morph <- median(richness.species.morph)
# sd.species.morph <- sd(richness.species.morph)

richness.genus.dna <- specnumber(genus.dna.df)
median.genus.dna <- median(richness.genus.dna)
# sd.genus.dna <- sd(richness.genus.dna)

gg <- data.frame(COI=richness.species.dna, 
                 Morphology=richness.species.morph, 
                 SSU=richness.genus.dna)
gg.melt <- reshape2::melt(gg)
names(gg.melt)<- c("Method", "Richness")

# check normality
shapiro.test(gg.melt$Richness)
# W = 0.77812, p-value = 2.773e-09

res <- data.frame(compare_means(Richness ~ Method, gg.melt, method="wilcox.test", 
              p.adjust.method = "holm"))
res

# .y.     group1     group2            p   p.adj p.format p.signif   method
# 1 Richness        COI Morphology 2.667090e-04 2.7e-04  0.00027      *** Wilcoxon
# 2 Richness        COI        SSU 1.396603e-09 4.1e-09  1.4e-09     **** Wilcoxon
# 3 Richness Morphology        SSU 1.381525e-09 4.1e-09  1.4e-09     **** Wilcoxon

# create factors
gg.melt$Method <- factor(gg.melt$Method, levels=c("SSU", "COI", "Morphology"),
                         labels=c("18S", "COI", "Morphology"))

# plot it 
p1 <- ggplot(gg.melt, aes(x=Method, y=Richness)) +
  geom_boxplot() +
  labs(x="", y="Alpha diversity") +
theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        # strip.background = element_blank(),
        # strip.text.x = element_text(angle = 0, hjust = 0),
        legend.position = "bottom",
        legend.title = element_blank()) 
p1
  


# total richness
total.species.dna <- length(unique(species.dna$Species))
total.species.morph <- length(unique(trad4$Species))
total.genus.dna <- length(unique(genus.dna$Genus))

gg2 <- data.frame(Total=c(total.species.dna, total.species.morph, total.genus.dna),
                  Type=c("COI", "Morphology", "18S"))

# create factors
gg2$Type <- factor(gg2$Type, levels=c("18S", "COI", "Morphology"))

# plot it 
p2 <-ggplot(gg2, aes(x=Type, y=Total)) +
  geom_bar(stat="identity") +
  labs(x="", y="Gamma diversity") +
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

g <- plot_grid(p1, p2, nrow=1, labels="auto")
ggsave("Outfiles/Fig2_richness.pdf", g, width = 8, height = 4)
