# Teresita M. Porter, May 6, 2022

# FigS1 resolution rarefaction

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('R_scripts/FigS1_resolution_rarefaction.R')

library(groundhog)

groundhog.library("vegan", "2022-05-04", tolerate.R.version='4.1.1') # rarecurve
groundhog.library("stringr", "2022-05-04", tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", "2022-05-04", tolerate.R.version='4.1.1') # dcast
groundhog.library("purrr", "2022-05-04", tolerate.R.version='4.1.1') # for map_dfr
groundhog.library("ggplot2", "2022-05-04", tolerate.R.version='4.1.1') # ggplot
groundhog.library("scales", "2022-05-04", tolerate.R.version='4.1.1') # comma
groundhog.library("cowplot", "2022-05-04", tolerate.R.version='4.1.1') # get_legend
groundhog.library("dplyr", "2022-05-04", tolerate.R.version='4.1.1') # map_dfr

#####################################################################

# Read in COI.csv from MetaWorks v1.0.0
F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE) # already limited to Arthropoda?
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE) # includes non-arthropods too
S <- read.csv(file="Infiles/results_18S_v4.1.csv", head=TRUE) # includes non-arthropods too

# edit column headers so matrices from different amplicons (using different taxonomies) can be merged
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

# create lineage
# keep unique lineages only
# recode species = 1, genus = 2, family = 3, order = 4, class = 5, phylum = 6, kingdom = 7, superkingdom = 8
# add to df

# F230R

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
F$SampleName <- gsub("TorontoHarbour_", "", F$SampleName)

# Reformat negative control (blanks) filenames
F$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", F$SampleName)
F$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", F$SampleName)
F$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", F$SampleName)
F$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", F$SampleName)
F$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", F$SampleName)
F$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", F$SampleName)

# Track negative control sample ESVs (0.02)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- F[(F$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
F.1 <- F[! F$GlobalESV %in% remove,]
# remove samples from df
F.1 <- F.1[!(F.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_3_COI_S80", "Opt_1_COI_S78", "Opt_2_COI_S79", 
                         "Opt_5_COI_S82","Opt_4_COI_S81")
F.2 <- F.1[!(F.1$SampleName %in% optimizationSamples),]

# just keep unique GlobalESV & Species
F.3 <- unique(F.2[,"Species", drop=FALSE])
F.3$resolution <- 1
F.3$type <- "DNA COI F230R"

names(F.3) <- c("Taxon", "Resolution", "Type")



# ml-jg

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
M$SampleName <- gsub("TorontoHarbour_", "", M$SampleName)

# Reformat negative control (blanks) filenames
M$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", M$SampleName)
M$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", M$SampleName)
M$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", M$SampleName)
M$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", M$SampleName)
M$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", M$SampleName)
M$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", M$SampleName)

# Track negative control sample ESVs (0.005)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- M[(M$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
M.1 <- M[! M$GlobalESV %in% remove,]
# remove samples from df
M.1 <- M.1[!(M.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_COI_S78", "Opt_4_COI_S81", "Opt_2_COI_S79", 
                         "Opt_3_COI_S80","Opt_5_COI_S82", "Opt_1_18S_S83")
M.2 <- M.1[!(M.1$SampleName %in% optimizationSamples),]

# just keep unique GlobalESV & Species
M.3 <- unique(M.2[,"Species",drop=FALSE])
M.3$resolution <- 1
M.3$Type <- "DNA COI ml-jg"

names(M.3) <- c("Taxon", "Resolution", "Type")





# 18S

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

# Track negative control sample ESVs (0.01)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- S[(S$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
S.1 <- S[!S$GlobalESV %in% remove,]
# remove samples from df
S.1 <- S.1[!(S.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_18S_S83", "Opt_2_18S_S84", "Opt_3_18S_S85", 
                         "Opt_5_18S_S87","Opt_4_18S_S86")
S.2 <- S.1[!(S.1$SampleName %in% optimizationSamples),]

# genus cutoff 0.70 for 80% correct
# family cutoff 0.20 for 80% correct
# c - o no cutoff needed for 90% correct
# k no cutoff needed for 95% correct
# sk no cutoff needed for 99% correct
S.2$taxon <- ifelse(!grepl("undef|\\d+|uncultured", S.2$Genus) & !S.2$gBP<0.70, paste("g__",S.2$Genus,sep=""),
                    ifelse(!grepl("undef|\\d+|uncultured", S.2$Family) & !S.2$fBP<0.20, paste("f__",S.2$Family,sep=""),
                           ifelse(!grepl("undef|\\d+|uncultured", S.2$Order) , paste("o__",S.2$Order,sep=""),
                                  ifelse(!grepl("undef|\\d+|uncultured", S.2$Class) , paste("c__",S.2$Class,sep=""),
                                         ifelse(!grepl("undef|\\d+|uncultured", S.2$Phylum) , paste("p__",S.2$Phylum,sep=""),
                                                ifelse(!grepl("undef|\\d+|uncultured", S.2$Kingdom) , paste("k__",S.2$Kingdom,sep=""), 
                                                       paste("sk__",S.2$SuperKingdom,sep="")))))))

# add a column for resolution
S.2$resolution <- ifelse(grepl("g__", S.2$taxon), 2,
                         ifelse(grepl("f__", S.2$taxon), 3,
                                ifelse(grepl("o__", S.2$taxon), 4,
                                       ifelse(grepl("c__", S.2$taxon), 5,
                                              ifelse(grepl("p__", S.2$taxon), 6,
                                                     ifelse(grepl("k__", S.2$taxon), 7, 8))))))

# just keep GlobalESV, taxon, resolution              
S.3 <- unique(S.2[,c("taxon", "resolution")])  
S.3$Type <- "DNA 18S"

names(S.3)<- c("Taxon", "Resolution", "Type")


##########################
# Read in morphology, 2018 only
# Raw data file: Invertebrate_2018_morphology.xlsx

B <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Correct a typo Mollusca;Gastropoda;Littoridinomorpha  - morphology
# should probably be Mollusca;Gastropoda;Littorinimorpha – DNA
B$Order <- gsub("Littoridinomorpha", "Littorinimorpha", B$Order)

# replace blank with NA
B[B==""] <- NA

# Create taxon column
B$taxon <- ifelse(!is.na(B$Species), paste("s__",B$Species,sep=""),
                  ifelse(!is.na(B$Genus), paste("g__",B$Genus,sep=""),
                         ifelse(!is.na(B$Family), paste("f__",B$Family,sep=""),
                                ifelse(!is.na(B$Order), paste("o__",B$Order,sep=""),
                                       ifelse(!is.na(B$Class), paste("c__",B$Class,sep=""), 
                                              paste("p__",B$Phylum,sep=""))))))

# add resolution
B$resolution <- ifelse(grepl("s__", B$taxon), 1,
                       ifelse(grepl("g__", B$taxon), 2,
                              ifelse(grepl("f__", B$taxon), 3,
                                     ifelse(grepl("o__", B$taxon), 4,
                                            ifelse(grepl("c__", B$taxon), 5, 6)))))


# just keep taxon and resolution
B.3 <- unique(B[,c("taxon" ,"resolution")])
B.3$Type <- "Morphology"
names(B.3) <- c("Taxon", "Resolution", "Type")

# put it all together
all <- rbind(F.3, M.3, S.3, B.3)

# create factors
all$Type <- factor(all$Type, levels=c("DNA COI F230R", "DNA COI ml-jg", "DNA 18S", "Morphology"))

# ggplot
p1 <- ggplot(all, aes(x=Type, y=Resolution)) + 
  geom_boxplot() +
  labs(y="", x="") +
  scale_y_continuous(name="",
                     breaks=c(1:7), 
                     labels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 
p1

resolution_figure <- p1











###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in COI.csv from MetaWorks v1.0.0
F <- read.csv(file="Infiles/results_F230R_10Aug20.csv", head=TRUE) # already limited to Arthropoda?
M <- read.csv(file="Infiles/results_ml-jg_10Aug20.csv", head=TRUE) # includes non-arthropods too
S <- read.csv(file="Infiles/results_18S_v4.1.csv", head=TRUE) # includes non-arthropods too

# edit column headers so matrices from different amplicons (using different taxonomies) can be merged
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

# Do rarefaction separately for COI-F230R, COI-mljg, 18S and morph


# F230R

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
F$SampleName <- gsub("TorontoHarbour_", "", F$SampleName)

# Reformat negative control (blanks) filenames
F$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", F$SampleName)
F$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", F$SampleName)
F$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", F$SampleName)
F$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", F$SampleName)
F$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", F$SampleName)
F$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", F$SampleName)

# Track negative control sample ESVs (0.02)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- F[(F$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
F.1 <- F[! F$GlobalESV %in% remove,]
# remove samples from df
F.1 <- F.1[!(F.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_3_COI_S80", "Opt_1_COI_S78", "Opt_2_COI_S79", 
                         "Opt_5_COI_S82","Opt_4_COI_S81")
F.2 <- F.1[!(F.1$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get fields for station and replicate
F.3<-data.frame(F.2, do.call(rbind, str_split(F.2$SampleName,"_")), stringsAsFactors = FALSE)
names(F.3)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
F.3$Station <- gsub("1271","1371", F.3$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
F.4<-data.frame(F.3, do.call(rbind, str_split(F.3$GlobalESV,"_")), stringsAsFactors = FALSE)
names(F.4)[37:38] <- c("Amplicon", "Zotu")

# # map locations to stations
# S <- read.csv('Sites.csv', header=T)
# names(S)[1] <- "Station"
# # left join to add locations to stations (always check for problems here)
# F.5 <- merge(x = A.4, y = S, by="Station", all.x = TRUE)

# # Create new field StationReplicateLocation
# A.5$StationReplicateLocation <- paste(A.5$Station, A.5$Replicate, A.5$Location, sep="_")

# pivot to make sample x ESV matrix
dna <- reshape2::dcast(F.4, Station ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(dna[,-1], 
                           sample=500, 
                           step=50, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
F.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# # Add sample names to vegan output (df) (rownames)
# sample_names <- rownames(dna)
# names(rare.df) <- sample_names

# Map rownames to vegan output (df)
F.df <- map_dfr(F.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# add col for marker
F.df$Type <- "DNA COI-F230R"


# color by station
p1 <- ggplot(data = F.df) +
  ggtitle("DNA COI F230R") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU), size=0.1) +
  # geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=3, override.aes = list(size = 2))) 
p1







# ml-jg

# Tweak sample names to have a consistent number of fields
# Remove TorontoHarbour prefix
M$SampleName <- gsub("TorontoHarbour_", "", M$SampleName)

# Reformat negative control (blanks) filenames
M$SampleName <- gsub("Oct1_BCBl_COI18Sp_S77","BCBl__BoxCorerBlank", M$SampleName)
M$SampleName <- gsub("Oct3_HoBl_COI18Sp_S79","HoBl__HoseBlank", M$SampleName)
M$SampleName <- gsub("Oct3_BCBl_COI18Sp_S80","BCBl__BoxCorerBlank", M$SampleName)
M$SampleName <- gsub("Oct1_HoBl_COI18Sp_S76","HoBl__HoseBlank", M$SampleName)
M$SampleName <- gsub("Oct1_StBl_COI18Sp_S78","StBl__StrawBlank", M$SampleName)
M$SampleName <- gsub("Oct3_StBl_COI18Sp_S81","StBl__StrawBlank", M$SampleName)

# Track negative control sample ESVs (0.005)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- M[(M$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
M.1 <- M[! M$GlobalESV %in% remove,]
# remove samples from df
M.1 <- M.1[!(M.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_COI_S78", "Opt_4_COI_S81", "Opt_2_COI_S79", 
                         "Opt_3_COI_S80","Opt_5_COI_S82", "Opt_1_18S_S83")
M.2 <- M.1[!(M.1$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get fields for station and replicate
M.3<-data.frame(M.2, do.call(rbind, str_split(M.2$SampleName,"_")), stringsAsFactors = FALSE)
names(M.3)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
M.3$Station <- gsub("1271","1371", M.3$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
M.4<-data.frame(M.3, do.call(rbind, str_split(M.3$GlobalESV,"_")), stringsAsFactors = FALSE)
names(M.4)[37:38] <- c("Amplicon", "Zotu")

# pivot to make sample x ESV matrix
dna <- reshape2::dcast(M.4, Station ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(dna[,-1], 
                           sample=500, 
                           step=50, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
M.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Map rownames to vegan output (df)
M.df <- map_dfr(M.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# add col for marker
M.df$Type <- "DNA COI-mljg"

# list -> df


# color by station
p2 <- ggplot(data = M.df) +
  ggtitle("DNA COI ml-jg") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU), size=0.1) +
  # geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=3, override.aes = list(size = 2))) 
p2







# 18S

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

# Track negative control sample ESVs (0.01)
negativeControls <- c("BCBl__BoxCorerBlank","HoBl__HoseBlank","StBl__StrawBlank")
neg <- S[(S$SampleName %in% negativeControls),]
remove <- unique(neg$GlobalESV)
# remove ESVs from df
S.1 <- S[!S$GlobalESV %in% remove,]
# remove samples from df
S.1 <- S.1[!(S.1$SampleName %in% negativeControls),]

# Drop optimization samples
optimizationSamples <- c("Opt_1_18S_S83", "Opt_2_18S_S84", "Opt_3_18S_S85", 
                         "Opt_5_18S_S87","Opt_4_18S_S86")
S.2 <- S.1[!(S.1$SampleName %in% optimizationSamples),]

# Split up SampleName with pkg 'stringr' to get fields for station and replicate
S.3<-data.frame(S.2, do.call(rbind, str_split(S.2$SampleName,"_")), stringsAsFactors = FALSE)
names(S.3)[33:36] <- c("Station", "Replicate", "COI18Sp","IlluminaSample")

# correct typo, change Station 1271 to 1371
S.3$Station <- gsub("1271","1371", S.3$Station)

# Split up GlobalESV with pkg 'stringr' to get amplicon field
S.4<-data.frame(S.3, do.call(rbind, str_split(S.3$GlobalESV,"_")), stringsAsFactors = FALSE)
names(S.4)[37:38] <- c("Amplicon", "Zotu")

# pivot to make sample x ESV matrix
dna <- reshape2::dcast(S.4, Station ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(dna[,-1], 
                           sample=500, 
                           step=50, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
S.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Map rownames to vegan output (df)
S.df <- map_dfr(S.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# add col for marker
S.df$Type <- "DNA 18S"


# color by station
p3 <- ggplot(data = S.df) +
  ggtitle("DNA 18S") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU), size=0.1) +
  # geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma,
                     breaks=c(30000, 60000)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=3, override.aes = list(size = 2))) 
p3



##########################
# Read in morphology, 2018 only
# Raw data file: Invertebrate_2018_morphology.xlsx

B <- read.table("Infiles/morph_new.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)

# Correct a typo Mollusca;Gastropoda;Littoridinomorpha  - morphology
# should probably be Mollusca;Gastropoda;Littorinimorpha – DNA
B$Order <- gsub("Littoridinomorpha", "Littorinimorpha", B$Order)

# replace blank with NA
B[B==""] <- NA

# Create lineage column
B$lineage <- paste(B$Phylum, B$Class, B$Order, B$Family, B$Genus, B$Species, sep="|")

# pivot to create station x lineage matrix
morph <- reshape2::dcast(B, station+core ~ lineage, value.var = "individuals", fun.aggregate = sum)
rownames(morph) <- paste(morph$station, morph$core, sep="_")
morph$station <- NULL
morph$core <- NULL



# set random seed
set.seed(1234)

# calculcate rarecurve with vegan
morph.df <- rarecurve2(morph, 
                       sample=15, 
                       step=5, 
                       label=T)

# Reformat vegan list as df (cols OTU, raw.read)
morph.df <- lapply(morph.df, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], indvs = rownames(b))
  b$indvs <- as.numeric(gsub("N", "",  b$indvs))
  return(b)
})


# Map rownames to vegan output (df)
morph.df <- map_dfr(morph.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Rename cols
names(morph.df) <- c("Station", "Taxa", "Individuals")

# color by station
p4 <- ggplot(data = morph.df) +
  ggtitle("Morphology") +
  labs(x="Individuals", y="Taxa") +
  geom_point(aes(x = Individuals, y = Taxa), size=0.1) +
  # geom_vline(xintercept = morph.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma,
                     breaks=c(1000, 2000)) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")
p4


g <- plot_grid(p1, p2, p3, p4, nrow = 1)
g

rarefaction_figure <- g




pg <- plot_grid(resolution_figure, rarefaction_figure, ncol=1, labels="auto")
ggsave("Outfiles/FigS1_resolution_rarefaction.pdf", pg, width = 8, height = 8)



