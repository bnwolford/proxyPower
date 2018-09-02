#!/usr/bin/Rscript

Revised from make_liabs.R script by Jimmy Liu (https://github.com/jimmyzliu/pmliability)

#===============================================================================
# Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#==============================================================================

################### Read Command Line Parameters #############################
suppressPackageStartupMessages(library(optparse))
optionList <- list(
    make_option(c("-h", "--heritability"), type="numeric", help="Estimated eritability of the trait"),
    make_option(c("-v","--prevalence"), type="character", help="File with prevalence by age and sex"),
    make_option(c("-p","--ped"),type="character",help=".ped file"),
    make_option(c("-o","--output"),type="character",help="Output file name")
	        )

parser <- OptionParser(
    usage="%prog -f file -a alpha",
        option_list=optionList
	    )

# a hack to fix a bug in optparse that won't let you use positional args
# if you also have non-boolean optional args:
getOptionStrings <- function(parserrObj){
optionStrings <- character()
        for(item in   parserObj@options) {
	        optionStrings <- append(optionStrings, c(item@short_flag, item@long_flag))
		} 
		optionStrings
		}

optStrings <- getOptionStrings(parser)
arguments <- parse_args(parser, positional_arguments=TRUE)


########################### Load functions and libraries ###################################################

library(tmvtnorm)
library(kinship2)
library(data.table)

############################ Main Commands ##################################################



h2 <- arguments$options$heritability
pedigree <-arguments$options$ped
prev <-arguments$options$prevalence
out <- arguments$options$output


prev <- fread(prev,header=T)
ped <- fread(pedigree,header=T)
fams <- unique(ped$FID)

# 1 == female, 0 == male
ped$FATHER[which(is.na(ped$FATHER))] <- ""
ped$MOTHER[which(is.na(ped$MOTHER))] <- ""

# change sex to: 1 == male, 2 == female, 3 == unknown
ped$SEX[which(is.na(ped$SEX))] <- "unknown"
ped$SEX[which(ped$SEX == "2")] <- "female"
ped$SEX[which(ped$SEX == "1")] <- "male"

#must have both mom and dad in ped
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
