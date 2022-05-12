# Teresita Porter, May 6, 2022

# FigS2 predictor correlations

# Run at command line to identify libraries not used
# Run this before using groundhog
# funchir::stale_package_check('R_scripts/FigS2_predictor_correlations.R')

library(groundhog)

groundhog.library("moments", "2022-05-04", tolerate.R.version='4.1.1') #skewness
groundhog.library("stats", "2022-05-04", tolerate.R.version='4.1.1')    # e.g for hclust() function
groundhog.library("stringr", "2022-05-04", tolerate.R.version='4.1.1') # str_split
groundhog.library("reshape2", "2022-05-04", tolerate.R.version='4.1.1') # dcast
groundhog.library("Hmisc", "2022-05-04", tolerate.R.version='4.1.1') # rcorr
groundhog.library("corrplot", "2022-05-04", tolerate.R.version='4.1.1') # corrplot

#######################################
# Read in metadata (water and sediment predictors)
C <-read.csv(file='Infiles/metadata.csv',head=TRUE)
names(C)[1] <- "SampleName"

# Tweak sample names to have a consistent number of fields
C$SampleName <- gsub("TorontoHarbour_", "", C$SampleName)
C$SampleName <- gsub("BCBl_COI18Sp","BCBl__BoxCorerBlank", C$SampleName)
C$SampleName <- gsub("HoBl_COI18Sp","HoBl__HoseBlank", C$SampleName)
C$SampleName <- gsub("StBl_COI18Sp","StBl__StrawBlank", C$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
C.1<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")), stringsAsFactors = FALSE)
names(C.1)[48:50] <- c("Station", "Replicate","Location")

# only keep bottom water and sediment vars
keep <- names(C.1)[c(9:15,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46)]

# remove SampleName
C.2 <- unique(C.1[,-c(1,49:50)])

# move Station to rownames
rownames(C.2) <- C.2$Station
C.2$Station <- NULL

# only keep bottom water and sediment
C.3 <- C.2[,-c(1:7,15:16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46)]


###################################
# Read in sediment predictors (metals, aromatics, pesticides)
D <-read.table(file='Infiles/TH_extra_data.txt', head=TRUE, sep="\t", stringsAsFactors = FALSE)
names(D)[1] <- "Station"

# replace NA with zero
D[is.na(D)] <- 0

# change 1765 to 1765A so vals match up across water & sediment samples
D$Station <- gsub("1765", "1765A", D$Station)

# Pool the VariableTypes by summing the individual chemicals
D.1 <- reshape2::dcast(D, Station ~ VariableType, value.var = "Value", fun.aggregate = sum)
row.names(D.1) <- D.1$Station
D.1$Station <- NULL

# remove spaces from names
names(D.1) <- c("Aromatic", 
                  "Halogenated_Substituted_Aliphatic",
                  "Halogenated_Substituted_Aromatic", 
                  "Halogenated_monoaromatic",
                  "Herbicides",
                  "Insecticides",
                  "Metal", 
                  "Metalloids", 
                  "Organochlorine", 
                  "PCB", 
                  "Pesticides",
                  "Rare_Earth")
# drop invariant vars, i.e. Metalloids
D.1 <- D.1[,-8]

##########################################
# combine vars for one analysis
merged <- cbind(C.3, D.1)

# test each env var for normality
# bottom water vars
hist(merged$b_temp)
shapiro.test(merged$b_temp)
skewness(merged$b_temp)
skewness(sqrt(max(merged$b_temp+1) - merged$b_temp))
skewness(log10(max(merged$b_temp+1) - merged$b_temp))
skewness(1/(max(merged$b_temp+1) - merged$b_temp)) # this one

hist(merged$b_DO_p) # ok

hist(merged$b_DO) # ok

hist(merged$b_SpCond)
shapiro.test(merged$b_SpCond)
skewness(merged$b_SpCond)
skewness(sqrt(merged$b_SpCond))
skewness(log10(merged$b_SpCond))
skewness(1/(merged$b_SpCond)) # this one

hist(merged$b_TDS)
shapiro.test(merged$b_TDS)
skewness(merged$b_TDS)
skewness(sqrt(merged$b_TDS))
skewness(log10(merged$b_TDS))
skewness(1/(merged$b_TDS)) # this one

hist(merged$b_SAL)
shapiro.test(merged$b_SAL)
skewness(merged$b_SAL)
skewness(sqrt(merged$b_SAL))
skewness(log10(merged$b_SAL))
skewness(1/(merged$b_SAL)) # this one

hist(merged$b_pH) # ok

hist(merged$b_Ammonia_N)
shapiro.test(merged$b_Ammonia_N)
skewness(merged$b_Ammonia_N)
skewness(sqrt(merged$b_Ammonia_N))
skewness(log10(merged$b_Ammonia_N))
skewness(1/(merged$b_Ammonia_N)) # this one

hist(merged$b_Calcium)
shapiro.test(merged$b_Calcium)
skewness(merged$b_Calcium)
skewness(sqrt(merged$b_Calcium))
skewness(log10(merged$b_Calcium))
skewness(1/(merged$b_Calcium)) # this one

hist(merged$b_Chloride)
shapiro.test(merged$b_Chloride)
skewness(merged$b_Chloride)
skewness(sqrt(merged$b_Chloride))
skewness(log10(merged$b_Chloride))
skewness(1/(merged$b_Chloride)) # this one

hist(merged$b_Magnesium) 
shapiro.test(merged$b_Magnesium)
skewness(merged$b_Magnesium)
skewness(sqrt(max(merged$b_Magnesium+1) - merged$b_Magnesium))
skewness(log10(max(merged$b_Magnesium+1) - merged$b_Magnesium))
skewness(1/(max(merged$b_Magnesium+1) - merged$b_Magnesium)) # this one

hist(merged$b_Nitrate_Nitrite_N) 
shapiro.test(merged$b_Nitrate_Nitrite_N)
skewness(merged$b_Nitrate_Nitrite_N)
skewness(sqrt(merged$b_Nitrate_Nitrite_N)) # this one
skewness(log10(merged$b_Nitrate_Nitrite_N))
skewness(1/(merged$b_Nitrate_Nitrite_N)) 

hist(merged$b_Nitrogen_Total)
shapiro.test(merged$b_Nitrogen_Total)
skewness(merged$b_Nitrogen_Total)
skewness(sqrt(merged$b_Nitrogen_Total))
skewness(log10(merged$b_Nitrogen_Total))
skewness(1/(merged$b_Nitrogen_Total)) # this one

hist(merged$b_Nitrogen_Total_Dissolved)
shapiro.test(merged$b_Nitrogen_Total_Dissolved)
skewness(merged$b_Nitrogen_Total_Dissolved)
skewness(sqrt(merged$b_Nitrogen_Total_Dissolved))
skewness(log10(merged$b_Nitrogen_Total_Dissolved))
skewness(1/(merged$b_Nitrogen_Total_Dissolved)) # this one

hist(merged$b_Phosphorus_SRP)
shapiro.test(merged$b_Phosphorus_SRP)
skewness(merged$b_Phosphorus_SRP)
skewness(sqrt(merged$b_Phosphorus_SRP))
skewness(log10(merged$b_Phosphorus_SRP)) # this one
skewness(1/(merged$b_Phosphorus_SRP)) 

hist(merged$b_Phosphorus_Total)
shapiro.test(merged$b_Phosphorus_Total)
skewness(merged$b_Phosphorus_Total)
skewness(sqrt(merged$b_Phosphorus_Total))
skewness(log10(merged$b_Phosphorus_Total))
skewness(1/(merged$b_Phosphorus_Total)) # this one

hist(merged$b_Potassium)
shapiro.test(merged$b_Potassium)
skewness(merged$b_Potassium)
skewness(sqrt(merged$b_Potassium))
skewness(log10(merged$b_Potassium))
skewness(1/(merged$b_Potassium)) # this one

hist(merged$b_Silica)
shapiro.test(merged$b_Silica)
skewness(merged$b_Silica)
skewness(sqrt(merged$b_Silica))
skewness(log10(merged$b_Silica))
skewness(1/(merged$b_Silica)) # this one

hist(merged$b_Sodium)
shapiro.test(merged$b_Sodium)
skewness(merged$b_Sodium)
skewness(sqrt(merged$b_Sodium))
skewness(log10(merged$b_Sodium))
skewness(1/(merged$b_Sodium)) # this one

hist(merged$b_Sulfate)
shapiro.test(merged$b_Sulfate)
skewness(merged$b_Sulfate)
skewness(sqrt(merged$b_Sulfate))
skewness(log10(merged$b_Sulfate))
skewness(1/(merged$b_Sulfate)) # this one


# sediment vars
hist(merged$Aromatic)
shapiro.test(merged$Aromatic)
skewness(merged$Aromatic)
skewness(sqrt(merged$Aromatic))
skewness(log10(merged$Aromatic))
skewness(1/(merged$Aromatic)) # this one

hist(merged$Halogenated_Substituted_Aliphatic)
shapiro.test(merged$Halogenated_Substituted_Aliphatic)
skewness(merged$Halogenated_Substituted_Aliphatic)
skewness(sqrt(merged$Halogenated_Substituted_Aliphatic))
skewness(log10(merged$Halogenated_Substituted_Aliphatic)) # this one

hist(merged$Halogenated_Substituted_Aromatic)
shapiro.test(merged$Halogenated_Substituted_Aromatic)
skewness(merged$Halogenated_Substituted_Aromatic)
skewness(sqrt(merged$Halogenated_Substituted_Aromatic)) # this one

hist(merged$Halogenated_monoaromatic)
shapiro.test(merged$Halogenated_monoaromatic)
skewness(merged$Halogenated_monoaromatic)
skewness(sqrt(merged$Halogenated_monoaromatic)) # this one

hist(merged$Herbicides)
shapiro.test(merged$Herbicides)
skewness(merged$Herbicides)
skewness(sqrt(merged$Herbicides)) # this one

hist(merged$Insecticides)
shapiro.test(merged$Insecticides) # this one

hist(merged$Metal)
shapiro.test(merged$Metal) # ok

hist(merged$Organochlorine)
shapiro.test(merged$Organochlorine)
skewness(merged$Organochlorine)
skewness(sqrt(merged$Organochlorine))
skewness(log10(merged$Organochlorine))
skewness(1/(merged$Organochlorine)) # this one

hist(merged$PCB)
shapiro.test(merged$PCB)
skewness(merged$PCB)
skewness(sqrt(merged$PCB)) 
skewness(log10(merged$PCB)) # this one

hist(merged$Pesticides)
shapiro.test(merged$Pesticides)
skewness(merged$Pesticides) 
skewness(sqrt(merged$Pesticides)) # this one

hist(merged$Rare_Earth)
shapiro.test(merged$Rare_Earth)
skewness(merged$Rare_Earth) 
skewness(sqrt(max(merged$Rare_Earth+1)-merged$Rare_Earth)) # this one

# create df.envir with transformed vars and NAs filled from nearest station
df.env <- cbind(
  1/(max(C.2$b_temp+1) - C.2$b_temp), 
  C.2$b_DO_p,
  C.2$b_DO, 
  1/(C.2$b_SpCond),
  1/(C.2$b_TDS),
  1/(C.2$b_SAL),
  C.2$b_pH, 
  1/(C.2$b_Ammonia_N),
  1/(C.2$b_Calcium),
  1/(C.2$b_Chloride),
  1/(max(C.2$b_Magnesium+1) - C.2$b_Magnesium),
  sqrt(C.2$b_Nitrate_Nitrite_N),
  1/(C.2$b_Nitrogen_Total),
  1/(C.2$b_Nitrogen_Total_Dissolved),
  log10(C.2$b_Phosphorus_SRP),
  1/(C.2$b_Phosphorus_Total), 
  1/(C.2$b_Potassium), 
  1/(C.2$b_Silica),
  1/(C.2$b_Sodium),
  1/(C.2$b_Sulfate),
                 
  1/(merged$Aromatic), 
  log10(merged$Halogenated_Substituted_Aliphatic), 
  sqrt(merged$Halogenated_Substituted_Aromatic), 
  sqrt(merged$Halogenated_monoaromatic),
  sqrt(merged$Herbicides),
  merged$Insecticides,
  merged$Metal, 
  1/(merged$Organochlorine), 
  log10(merged$PCB), 
  sqrt(merged$Pesticides), 
  sqrt(max(merged$Rare_Earth+1)-merged$Rare_Earth))

# column standardization (accounts for different units accross vars)
# standardize using z-scores
# 'scale' centers=TRUE, scales=TRUE columns of matrix
# 'apply' 2-columns, 1-rows
df.envir <- data.frame(apply(df.env[,], 2, scale))
names(df.envir) <- c(
"b_Temp", "b_DO_p", "b_DO", "b_SpCond", "b_TDS", "b_SAL",
"b_pH", "b_NH3_N", "b_Ca","b_Cl", "b_Mg", "b_NO3_NO2_N",
"b_TN", "b_TN_Diss", "b_SRP_P", "b_TP", "b_K", "b_Silica",
"b_Na","b_SO4",
                     
"Aromatic", "Halogenated_Substituted_Aliphatic", 
"Halogenated_Substituted_Aromatic", "Halogenated_monoaromatic",
"Herbicides", "Insecticides",
"Metal", "Organochlorine", "PCB",
"Pesticides", "Rare_Earth")

rownames(df.envir) <- rownames(C.3)

# plot

pdf("Outfiles/FigS2_predictor_correlations.pdf") 

vars <- df.envir

vars.rcorr <- rcorr(as.matrix(vars), type = "pearson")

# correlations > 0.7 may indicate multicollinearity, 
# so only plot those ones if they are significant
is.na(vars.rcorr$r) <- abs(vars.rcorr$r) < 0.7
vars.rcorr$r[is.na(vars.rcorr$r)] <- 0

# Insignificant correlations are leaved blank
corrplot(vars.rcorr$r, 
         #      type="lower", 
         #       addrect = 5,
         number.cex=0.5,
         order="hclust", 
         tl.col="black",
         tl.cex=0.5,
         cl.cex=0.5,
         p.mat = vars.rcorr$P, sig.level = 0.05, insig = "blank")

dev.off()

# Size of circles reflects size of Pearson correlation, 
# only plotted if > 0.7
# Only significant results shown else blank (p-value < 0.05)



