library(tmvtnorm)
library(kinship2)
args = commandArgs(trailingOnly=TRUE)

#trait <- "AD"
#h2 <- 0.75
trait <- as.character(args[1])
h2 <- as.numeric(args[2])
#batch <- as.numeric(args[3])

prev <- read.table(paste("prev_files/",trait,".age_prev",sep=""),header = T)

ped <- read.table(paste("ped_files/",trait,".ped.gz",sep=""),header = T,colClasses = c(rep("character",5),rep("numeric",2)))
fams <- unique(ped$FID)

# 1 == female, 0 == male
ped$FATHER[which(is.na(ped$FATHER))] <- ""
ped$MOTHER[which(is.na(ped$MOTHER))] <- ""

                                        # change sex to: 1 == male, 2 == female, 3 == unknown
ped$SEX[which(is.na(ped$SEX))] <- "unknown"
ped$SEX[which(ped$SEX == "1")] <- "female"
ped$SEX[which(ped$SEX == "0")] <- "male"
big_ped <- pedigree(id = ped$IID,
                    dadid = ped$FATHER,
                    momid = ped$MOTHER,
                    sex = ped$SEX,
                    famid = ped$FID,
                    affected = ped$PHENO)


liab_pheno <- data.frame()
for(fam in fams){
    fam_info <- ped[ped$FID == fam,]
    fam_ped <- pedigree(id = fam_info$IID,
                        dadid = fam_info$FATHER,
                        momid = fam_info$MOTHER,
                        sex = fam_info$SEX,
                        famid = fam_info$FID,
                        affected = fam_info$PHENO)
    #fam_ped <- big_ped[fam]
    sigma <- 2*kinship(fam_ped)
    sigma <- h2 * sigma
    diag(sigma) <- rep(1,nrow(sigma))
    ids <- fam_ped$id

    #if(NA %in% fam_info$PHENO){
                                        #next
    
                                        #}
    prevs <- rep(0,length(ids))

    #upper
    l_upper <- rep(0,nrow(fam_info))
    l_lower <- rep(0,nrow(fam_info))
    for(i in 1:nrow(fam_info)){
        age <- fam_info$AGE[i]
        # ugh. missing age. fudge for now
        if(is.na(age)){ age <- 50}
        sex <- fam_info$SEX[i]
        pheno <- fam_info$PHENO[i]
        # if missing phenotypes in family, exclude. fix later
        if(age < min(prev$AGE)){ age <- min(prev$AGE)}
        if(age > max(prev$AGE)){ age <- max(prev$AGE)}
        k <- sum(prev[prev$AGE == age,c(2,3)])/2
        if(sex == "male"){ k <- prev[prev$AGE == age,4]}
        if(sex == "female"){ k <- prev[prev$AGE == age,5]}
        t <- qnorm(1-k)
        if(is.na(pheno)){
            l_upper[i] <- Inf
            l_lower[i] <- -Inf
        }
        else{
            if(pheno == 1){
                l_upper[i] <- Inf
                l_lower[i] <- t
                if(k == 0){
                    all_prev <- c(prev[,2],prev[,3])
                    min_k <- min(all_prev[all_prev!=0])
                    t <- qnorm(1-min_k)
                    l_lower[i] <- t
                }
            }
            if(pheno == 0){
                l_upper[i] <- t
                l_lower[i] <- -Inf
            }
        }
    }
    liab <- rtmvnorm(1000,
                     mean = rep(0,length(ids)),
                     sigma = sigma,
                     lower = l_lower,
                     upper = l_upper,
                     algorithm = "gibbs",
                     burn.in.samples = 100,
                     thinning = 5)
    liab <- colSums(liab)/1000
    liab_pheno <- rbind(liab_pheno,data.frame(ids,liab))
    if(nrow(liab_pheno) %% 500 == 0){ print(nrow(liab_pheno))}
}

write.table(liab_pheno,paste("liab_files/",trait,".liab",sep=""),col.names=F,row.names=F,quote=F)
