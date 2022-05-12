# Teresita M. Porter, May 6, 2022

# Table S3 Taxon checklist

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('R_scripts/TableS3_checklist.R')

library(groundhog)

groundhog.library("stringr", "2022-05-04", tolerate.R.version='4.1.1')
groundhog.library("dplyr", "2022-05-04", tolerate.R.version='4.1.1')

#############################################
# Parse morphology results
# Read in Morphology counts
B <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Get unique taxa
Morph <- data.frame(unique(B[c(5:6,8:11)]))

# add col to indicate morph presence
Morph$Morph <- 1




#########################################
# Process metabarcoding data

# Read in COI.csv from MetaWorks v1.0.0
F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE, stringsAsFactors = FALSE) # already limited to Arthropoda
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE, stringsAsFactors = FALSE) # includes non-arthropods too
S <- read.csv(file="Infiles/results_18S_v4.1.csv", head=TRUE, stringsAsFactors = FALSE) # includes non-arthropods too

# edit column headers so matrices can be merged
# GlobalESV
names(F)[1] <- "GlobalESV"
names(M)[1] <- "GlobalESV"
names(S)[1] <- "GlobalESV"

# Seqfield
names(F)[4] <- "ORF/ESVSeq"
names(M)[4] <- "ORF/ESVSeq"
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
# gBP > 0.70, fBP > 0.20
S$Genus <- ifelse(S$gBP >=0.70, S$Genus, "")
S$Family <- ifelse(S$fBP >=0.20, S$Family, "")

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

# remove GlobalESVs from rest of data
A.2 <- A.1[!(A.1$GlobalESV %in% remove),]
# now remove the negative control samples
A.2 <- A.2[!(A.2$SampleName %in% negativeControls),]

# Throw out optimization samples
optimizationSamples <- c("Opt_3_COI_S80", "Opt_1_COI_S78", "Opt_2_COI_S79", 
                         "Opt_5_COI_S82","Opt_4_COI_S81","Opt_1_18S_S83")
A.3 <- A.2[!(A.2$SampleName %in% optimizationSamples),]

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.4 <- data.frame(A.3, do.call(rbind, str_split(A.3$SampleName,"_")), stringsAsFactors = FALSE)
names(A.4)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
A.4$Station <- gsub("1271","1371", A.4$Station)

# Split up GlobalESV with pkg 'stringr'
A.5 <- data.frame(A.4, do.call(rbind, str_split(A.4$GlobalESV,"_")), stringsAsFactors = FALSE)
names(A.5)[37:38] <- c("Amplicon", "Zotu")

# get unique taxa
DNA <- data.frame(unique(A.5[,c(15,18,21,24,27,30)]))

DNA <- data.frame(lapply(DNA, function(x) { gsub("_", " ", x) }))

# Add a col to indicate DNA presence
DNA$DNA <- 1


###########
# put dna and morph together

# merge morph with metabarcoding
merged <- merge(Morph, DNA, by=c("Phylum", "Class", "Order", "Family", "Genus", "Species"), all=TRUE)

# change NAs to 0
merged[is.na(merged)] <- 0
# merged$Genus <- gsub("0","", merged$Genus)
# merged$Species <- gsub("0", "", merged$Species)

# sort for easier reading
merged.sorted <- merged[order(merged$Phylum, merged$Class, merged$Order, merged$Family, merged$Genus, merged$Species),]

# # test
# merged.sorted$sum <- merged.sorted$DNA+merged.sorted$Morph
# 
# merged.sorted[order(-merged.sorted$sum),]

# print it
write.csv(merged.sorted, "Outfiles/COI_checklist_2022-04-11.csv", quote = FALSE, row.names = FALSE)


# some stats for the text
# minus one because sometimes comes up as ""
length(unique(Morph$Species))
# 48-1
length(unique(DNA$Species))
# 89-1


length(unique(Morph$Genus))
# 78-1
length(unique(DNA$Genus))
# 80-1

length(unique(Morph$Family))
# 31-1
length(unique(DNA$Family))
# 62-1






# now just the 18S (non-macroinverts)

# limit to macroinvertebrates phyla because need to calc %EPT and %Chironomid
# i.e. the ones picked up by morphology in 2018
morphPhyla <- c("Arthropoda", "Annelida", "Mollusca", "Cnidaria", "Platyhelminthes")
S.1 <- S[!S$Phylum %in% morphPhyla,]

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
S.1$SampleName <- gsub("TorontoHarbour_", "", S.1$SampleName)

# Reformat negative control (blanks) filenames
S.1$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", S.1$SampleName)
S.1$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", S.1$SampleName)
S.1$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", S.1$SampleName)
S.1$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", S.1$SampleName)
S.1$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", S.1$SampleName)
S.1$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", S.1$SampleName)

# Track negative control samples
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- S.1[(S.1$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)

# remove GlobalESVs from rest of data
S.2 <- S.1[!(S.1$GlobalESV %in% remove),]
# now remove the negative control samples
S.2 <- S.2[!(S.2$SampleName %in% negativeControls),]

# Throw out optimization samples
optimizationSamples <- c("Opt_1_18S_S83", "Opt_2_18S_S84", "Opt_3_18S_S85", 
                         "Opt_5_18S_S87", "Opt_4_18S_S86")
S.3 <- S.2[!(S.2$SampleName %in% optimizationSamples),]

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
S.4 <- data.frame(S.3, do.call(rbind, str_split(S.3$SampleName,"_")), stringsAsFactors = FALSE)
names(S.4)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
S.4$Station <- gsub("1271","1371", S.4$Station)

# Split up GlobalESV with pkg 'stringr'
S.5 <- data.frame(S.4, do.call(rbind, str_split(S.4$GlobalESV,"_")), stringsAsFactors = FALSE)
names(S.5)[37:38] <- c("Amplicon", "Zotu")

# get unique taxa
DNA <- data.frame(unique(S.5[,c(15,18,21,24,27)]))

DNA <- data.frame(lapply(DNA, function(x) { gsub("_", " ", x) }))

# Add a col to indicate DNA presence
DNA$DNA <- 1


###########
# put dna and morph together


# sort for easier reading
DNA.sorted <- DNA[order(DNA$Phylum, DNA$Class, DNA$Order, DNA$Family, DNA$Genus),]

# # test
# merged.sorted$sum <- merged.sorted$DNA+merged.sorted$Morph
# 
# merged.sorted[order(-merged.sorted$sum),]

# print it
write.csv(DNA.sorted, "Outfiles/18S_checklist_2022-04-11.csv", quote = FALSE, row.names = FALSE)






# some stats for the text

genus.df <- DNA.sorted[!grepl("undef|\\d+|uncultured", DNA.sorted$Genus),]
genus.df <- genus.df[!genus.df$Genus=="",]
length(unique(genus.df$Genus))
# [1] 271

family.df <- DNA.sorted[!grepl("undef|\\d+|uncultured", DNA.sorted$Family),]
family.df <- family.df[!family.df$Family=="",]
length(unique(family.df$Family))
# [1] 111

order.df <- DNA.sorted[!grepl("undef|\\d+|uncultured", DNA.sorted$Order),]
order.df <- order.df[!order.df$Order=="",]
length(unique(order.df$Order))
# [1] 144

class.df <- DNA.sorted[!grepl("undef|\\d+|uncultured", DNA.sorted$Class),]
class.df <- class.df[!class.df$Class=="",]
length(unique(class.df$Class))
# [1] 77

phylum.df <- DNA.sorted[!grepl("undef|\\d+|uncultured", DNA.sorted$Phylum),]
phylum.df <- phylum.df[!phylum.df$Phylum=="",]
length(unique(phylum.df$Phylum))
# 36

# count up number of taxa per phylum, show in descending order
DNA.sorted %>% group_by(Phylum) %>% dplyr::tally(sort=TRUE)

# What's up with Archaeorhizomyces
length(unique(S.5$GlobalESV[S.5$Genus=="Archaeorhizomyces"]))
# 12 ESVs
sum(S.5$ESVsize[S.5$Genus=="Archaeorhizomyces"])
# 117

