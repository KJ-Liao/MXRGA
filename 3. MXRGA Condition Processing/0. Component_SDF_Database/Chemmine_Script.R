library("ChemmineR")
sdf_database <- read.SDFset("MERGE_SDF.txt")
valid <- validSDF(sdf_database)
sdfset <- sdf_database[valid]

# MCS Tanimoto
library(fmcsR)

# Empty Matrix
DISM <- c() 

for(i in c(1:length(sdfset))){
  Dists <- c()
  for(j in c(1:length(sdfset))){
    
    MCS <-fmcs(sdfset[i], sdfset[j], au=0, bu=0)
    result <- unname(MCS@stats[4])
    Dists <- c(Dists, result) 
  }
  DISM <- rbind(DISM, Dists) 
}

write.csv(DISM, "MCS_DISM.csv")

#
#
#

# AP Tanimoto
apset <- sdf2ap(sdfset)

# Empty Matrix
DISM <- c() 

for(i in c(1:length(apset))){
  Dists <- c()
  for(j in c(1:length(apset))){
    result <- cmp.similarity(apset[i], apset[j])
    Dists <- c(Dists, result) 
  }
  DISM <- rbind(DISM, Dists) 
}

write.csv(DISM, "AP_DISM.csv")