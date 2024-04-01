library("ChemmineR")
sdf_database <- read.SDFset("MERGE_SDF.txt")
valid <- validSDF(sdf_database)
sdfset <- sdf_database[valid]

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
